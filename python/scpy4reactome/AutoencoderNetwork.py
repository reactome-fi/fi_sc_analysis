import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras.layers import Input, Dense
from keras.models import Model
from keras.callbacks import EarlyStopping
from keras import regularizers


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
    subset_cells = adata.obs.sample(frac=0.5, replace=False, random_state=1234).index.values
    x_val = adata[subset_cells, :]
    x_train = adata[~adata.obs.index.isin(x_val.obs.index.values), :]

    batch_size = n_dim
    validation_split = 0.1
    verbose = True
    callbacks = []
    early_stop = 3
    early_stop_cb = EarlyStopping(monitor='val_loss', patience=early_stop, verbose=verbose)
    callbacks.append(early_stop_cb)

    autoencoder_history = autoencoder.fit(x_train.X.toarray(), x_train.X.toarray(),
                                          epochs=epochs,
                                          batch_size=batch_size,
                                          shuffle=True,
                                          validation_split=validation_split,
                                          callbacks=callbacks,
                                          validation_data=(x_val.X.toarray(), x_val.X.toarray()),
                                          verbose=verbose)
    # compressed representation layer - codings - superset
    encoded_genes = encoder.predict(adata.X.toarray())
    # output layer
    decoded_genes = decoder.predict(encoded_genes)

    # TODO: return dictionary - or gene supersets depending on usecase
    # return {'encoded_genes': encoded_genes, 'loss': autoencoder_history, 'decoded_genes': decoded_genes}
    adata.uns['loss'] = autoencoder_history
    adata.obsm['codings'] = encoded_genes

    return adata
