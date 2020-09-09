import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from keras.layers import Input, Dense
from keras.models import Model
from keras.callbacks import EarlyStopping
from keras import regularizers


def gene_autoencoder(adata):
    """
    resources:
    https://bit.ly/2Rce0YE
    https://bit.ly/3if980t
    https://www.tensorflow.org/api_docs/python/
    https://github.com/theislab/dca

    :param adata: preprocessed 10X genomics counts 
    """
    # autoencoder architecture & model
    encoding_dim = 32  # adata.X class 'numpy.float32
    input_layer = Input(shape=(adata.n_vars,), name='count') # TODO: clustered gene set vs counts as input layer
    # TODO: add activity_regularizer=regularizers.l1(10) # default for l1 or l2 for smaller val_loss convergence?
    # encoded = Dense(encoding_dim, activation='relu', activity_regularizer=regularizers.l1(10 - 6))(input_layer)
    encoded = Dense(encoding_dim, activation='relu')(input_layer)
    decoded = Dense(adata.n_vars, activation='sigmoid')(encoded)

    autoencoder = Model(input_layer, decoded)
    encoder = Model(input_layer, encoded)

    encoded_input = Input(shape=(32,))
    decoder_layer = autoencoder.layers[-1]
    decoder = Model(encoded_input, decoder_layer(encoded_input))

    # compile model
    # TODO:  # add default param vals for smaller val_loss convergence? ex. https://bit.ly/2Rce0YE
    optimizer_instance = tf.keras.optimizers.SGD(learning_rate=0.05,
                                                 momentum=0.9,
                                                 nesterov=True,
                                                 name='SGD',
                                                 decay=10-6)
    autoencoder.compile(optimizer='SGD', loss='mean_squared_error')

    # fit model
    # TODO: set aside x_test to avoid over fitting in autoencoder.fit()
    # ~ %5-%25
    # int((anndata.n_vars * 0.25))

    inputs = adata.X.toarray()
    output = adata.X.toarray()

    batch_size = 32
    validation_split = 0.1
    verbose = True
    callbacks = []
    early_stop = 3
    early_stop_cb = EarlyStopping(monitor='val_loss', patience=early_stop, verbose=verbose)
    callbacks.append(early_stop_cb)

    autoencoder_loss = autoencoder.fit(inputs, output,
                                       epochs=200,
                                       batch_size=batch_size,
                                       shuffle=True,
                                       validation_split=validation_split,
                                       callbacks=callbacks,
                                       # validation_data=(x_test, x_test),
                                       verbose=verbose)
    # compressed representation layer
    encoded_genes = encoder.predict(adata.X.toarray())
    # output layer
    decoded_genes = decoder.predict(encoded_genes)

    # TODO: return dictionary - or gene supersets depending on usecase
    return(encoded_genes)