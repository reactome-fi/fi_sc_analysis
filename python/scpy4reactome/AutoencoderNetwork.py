import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras.layers import Input, Dense
from keras.models import Model
from keras.callbacks import EarlyStopping
from keras import regularizers

import os
import numpy as np
import pandas as pd
from resVAE.resvae import resVAE
import resVAE.utils as cutils
from resVAE.config import config
import scanpy as sc


def gene_autoencoder(adata, n_dim=200, epochs=30):
    """
    resources:
    https://bit.ly/2Rce0YE
    https://bit.ly/3if980t
    https://www.tensorflow.org/api_docs/python/
    https://github.com/theislab/dca
    :param adata: preprocessed 10X genomics counts
    """
    # autoencoder architecture & model
    encoding_dim = n_dim  # ~ 200 gene sets via https://bit.ly/2Rce0YE
    input_layer = Input(shape=(adata.n_vars,), name='count')  # input layer is all preprocessed 10X gene counts

    # TODO:
    # - add another dense layer for clustered genes based on chemical and physical properties learned from MSigDB-CGP per https://bit.ly/2Rce0YE
    #   section: Incorporate gene sets into the encoder layer
    #   these are human annotations. https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C2
    #   homology options  https://www.biostars.org/p/209855/
    #   autoenoder architecture should be: input -> MSigDB-CGP (genesets) -> ~200 geneset codings ("superset") -> MSigDB-CGP (genesets) -> output
    # - How to track metadata or "labeled data" pg 576 oreilly to extract and reconstruct anndata via numpy array?
    encoded = Dense(encoding_dim, activation='relu')(input_layer)
    decoded = Dense(adata.n_vars, activation='sigmoid')(encoded)

    autoencoder = Model(input_layer, decoded)
    encoder = Model(input_layer, encoded)

    encoded_input = Input(shape=(n_dim,))
    decoder_layer = autoencoder.layers[-1]
    decoder = Model(encoded_input, decoder_layer(encoded_input))

    # compile model param vals via https://bit.ly/2Rce0YE
    # Section: Establish and train the gene superset autoencoder
    optimizer_instance = tf.keras.optimizers.SGD(learning_rate=0.05,
                                                 momentum=0.9,
                                                 nesterov=True,
                                                 name='SGD',
                                                 decay=10-6)

    autoencoder.compile(optimizer=optimizer_instance, loss='mean_squared_error')

    # fit model
    # ~ %5-%25 set for validation. choose %5 via https://bit.ly/2Rce0YE
    # anndata where anndata.var holds gene and anndata.obs holds cell information
    # subset_cells = adata.obs.sample(frac=0.5, replace=False, random_state=1234).index.values
    # x_val = adata[subset_cells, :]
    # x_train = adata[~adata.obs.index.isin(x_val.obs.index.values), :]

    batch_size = n_dim
    validation_split = 0.2
    verbose = True
    callbacks = []
    early_stop = 3
    early_stop_cb = EarlyStopping(monitor='val_loss', patience=early_stop, verbose=verbose)
    callbacks.append(early_stop_cb)

    # cluster label example
    # sc.pp.neighbors(adata)
    # sc.tl.leiden(adata)

    # RE-VISIT: haults due to layer dims
    labels = np.array(adata.obs.leiden.values)
    autoencoder_history = autoencoder.fit([adata.X.toarray(), labels],
                                       epochs=epochs,
                                       batch_size=batch_size,
                                       shuffle=True,
                                       validation_split=validation_split,
                                       callbacks=callbacks,
                                       # validation_data=(x_val.X.toarray(), x_val.X.toarray()),
                                       verbose=verbose)

    # compressed representation layer - codings - superset
    encoded_genes = encoder.predict(adata.X.toarray())
    # output layer
    decoded_genes = decoder.predict(encoded_genes)

    # TODO: return dictionary - or gene supersets depending on usecase
    # return {'encoded_genes': encoded_genes, 'loss': autoencoder_history, 'decoded_genes': decoded_genes}
    adata.uns['loss'] = autoencoder_history
    adata.obsm['X_codings'] = encoded_genes

    return adata


def resVAE_autoencoder(adata=adata, model_path='data/models'){
    """
    Example of resVAE initial workable workflow with AnnData 
    
    # python 3.7 env setup ------------------------
    conda env list
    conda create -n scfipy python=3.7 
    conda activate scfipy
    pip install pkginfo>=1.4.2 bleach>=2.1.0 docutils>=0.13.1 tqdm>=4.14 ipython argparse scanpy scvelo --upgrade 
    # for resVAE ----------------------------------
    pip install tensorflow==1.15 keras==2.3.1
    git clone https://github.com/lab-conrad/resVAE.git
    cd resVAE
    vim setup.py # comment out #'tensorflow-gpu==1.15.4', wq 
    python setup.py install 
    cd data 
    mkdir models 
    modify config file for your model 
    pip list 
    """

    # -----------------------------------
    # resVAE autoencoder
    # -----------------------------------
    adata = adata[:, adata.var.highly_variable]
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata) # example clusters to work with

    categorical_clusters, l_encoder = cutils.one_hot_encoder(np.array(adata.obs.leiden.values))

    # config would already be loaded or we could possibly fine tune the params currently the loss does not converge to 0
    config['INPUT_SHAPE'] = (adata.X.toarray().shape[1], categorical_clusters.shape[1])
    config['ENCODER_SHAPE'] = [600, 100]
    config['DECODER_SHAPE'] = [100, 600]
    config['DECODER_REGULARIZER'] = 'var_l1_l2'
    config['DECODER_REGULARIZER_INITIAL'] = 1e-6
    config['LATENT_SCALE'] = 3
    config['DECODER_RELU_THRESH'] = 0
    config['OPTIMIZER'] = 'SGD'

    reactome_resvae = resVAE(model_dir=model_path, config=config)
    reactome_resvae.genes = adata.var.index

    reactome_resvae.describe_model('encoder')

    reactome_resvae.add_latent_labels(l_encoder.classes_)

    reactome_resvae.compile()

    reactome_resvae.fit(exprs=adata.X.toarray(), classes=categorical_clusters, model_dir=model_path)

    # -----------------------------------
    # save model
    # -----------------------------------
    reactome_resvae.save_model_new(model_path)
    # reactome_resvae.load_model_new(model_dir=model_path, model_name='my_model')

    # -----------------------------------
    # map back to genes and NN nuerons
    # -----------------------------------
    weights_clusters = reactome_resvae.get_latent_to_gene(normalized=True)
    weights_neurons_1 = reactome_resvae.get_neuron_to_gene(normalized=True, initial_layer=1)
    weights_neurons_2 = reactome_resvae.get_neuron_to_gene(normalized=False, initial_layer=2)
    weights_latent_neurons_1 = reactome_resvae.get_latent_to_neuron(normalized=True, target_layer=1)
    weights_latent_neurons_2 = reactome_resvae.get_latent_to_neuron(normalized=True, target_layer=2)
    biases = reactome_resvae.get_gene_biases(relative=True)

    return {reactome_resvae=reactome_resvae, weights_clusters=weights_clusters}
}

