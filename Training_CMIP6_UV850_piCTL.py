import tensorflow as tf
from tensorflow.keras import layers, models, optimizers
import numpy as np
import xarray as xr
from tqdm import tqdm
import os
from tensorflow.keras.callbacks import ReduceLROnPlateau, EarlyStopping, ModelCheckpoint, History
import matplotlib.pyplot as plt

#GPU
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    tf.config.experimental.set_visible_devices(gpus[1], 'GPU')
    tf.config.experimental.set_memory_growth(gpus[1], True)


def highestPowerof2(n):
    res = 0
    for i in range(n, 0, -1):
        if (i & (i - 1)) == 0:  # Vérifie si i est une puissance de 2
            res = i
            break
    return res

shape_inputs = [40, 42]
size = min(highestPowerof2(shape_inputs[0]), highestPowerof2(shape_inputs[1]))
diff_lat = 40 - size + 1
diff_lon = 42 - size + 1

def build_unet(input_shape=(40, 42, 2)):
    inputs = layers.Input(shape=input_shape)

    # Initial Conv
    x = layers.Conv2D(32, (diff_lat, diff_lon), padding='valid')(inputs)
    x = layers.BatchNormalization()(x)
    x = tf.keras.activations.tanh(x)

    # Encoder path
    enc1 = block_conv(x, 32)
    enc1_1 = layers.MaxPooling2D(pool_size=2)(enc1)
    enc2 = block_conv(enc1_1, 64)
    enc2_1 = layers.MaxPooling2D(pool_size=2)(enc2)
    enc3 = block_conv(enc2_1, 128)
    enc3_1 = layers.MaxPooling2D(pool_size=2)(enc3)

    # Bottleneck
    bottleneck = layers.Conv2D(128, 3, padding='same')(enc3_1)

    # Decoder path
    dec3 = layers.Conv2DTranspose(128, 3, strides=2, padding='same')(bottleneck)
    dec3 = layers.Concatenate()([dec3, enc3])
    dec3 = block_conv(dec3, 128)

    dec2 = layers.Conv2DTranspose(64, 3, strides=2, padding='same')(dec3)
    dec2 = layers.Concatenate()([dec2, enc2])
    dec2 = block_conv(dec2, 64)

    dec1 = layers.Conv2DTranspose(32, 3, strides=2, padding='same')(dec2)
    dec1 = layers.Concatenate()([dec1, enc1])
    dec1 = block_conv(dec1, 32)

    # Output layer
    outputs = layers.Conv2D(1, 3, padding='same')(dec1)

    return models.Model(inputs, outputs)



def block_conv(x, filters):
    x = layers.Conv2D(filters, 3, padding='same')(x)
    x = layers.BatchNormalization()(x)
    x = layers.ReLU()(x)
    x = layers.Conv2D(filters, 3, padding='same')(x)
    x = layers.BatchNormalization()(x)
    x = layers.ReLU()(x)
    return x

# ---------------------------- DataLoader ---------------------
class DataLoader(tf.keras.utils.Sequence):
    def __init__(self, ua_file, va_file, tas_file, batch_size=32):
        self.ua_data = xr.open_dataset(ua_file)['ua'][:,0,:,:]
        self.va_data = xr.open_dataset(va_file)['va'][:,0,:,:]
        self.tas_data = xr.open_dataset(tas_file)['tas']
        self.batch_size = batch_size

    def __len__(self):
        return int(np.ceil(len(self.tas_data) / self.batch_size))

    def __getitem__(self, idx):
        ua = self.ua_data[idx * self.batch_size:(idx + 1) * self.batch_size].values
        va = self.va_data[idx * self.batch_size:(idx + 1) * self.batch_size].values
        tas = self.tas_data[idx * self.batch_size:(idx + 1) * self.batch_size].values

        # Combine slp and second_var along the channel dimension
        input_stack = np.stack((ua, va), axis=-1)
        input_tensor = tf.convert_to_tensor(input_stack,dtype=tf.float32)
        tas_tensor = tf.convert_to_tensor(np.expand_dims(tas, axis=-1),dtype=tf.float32)

        return input_tensor, tas_tensor

# Configuration
model = build_unet()
model.compile(optimizer=optimizers.Adam(learning_rate=1e-3), loss='mse',metrics=['mse'])

ua_file_train = 'ua_850_day_CMIP6_all_piControl_1.3x1.3_fillmiss_CR_TRAIN.nc'
va_file_train = 'va_850_day_CMIP6_all_piControl_1.3x1.3_fillmiss_CR_TRAIN.nc'

ua_file_val = 'ua_850_day_CMIP6_all_piControl_1.3x1.3_fillmiss_CR_VAL.nc'
va_file_val = 'va_850_day_CMIP6_all_piControl_1.3x1.3_fillmiss_CR_VAL.nc'

tas_file_train = 'Tano_day_CMIP6_all_piControl_1.3x1.3_TRAIN.nc'
tas_file_val = 'Tano_day_CMIP6_all_piControl_1.3x1.3_VAL.nc'


# Datasets
dataset_train = DataLoader(ua_file_train, va_file_train, tas_file_train, batch_size=2500)
dataset_val = DataLoader(ua_file_val,va_file_val,tas_file_val, batch_size=2500)


callback = [ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=2, verbose=1), EarlyStopping(monitor='val_loss', patience=10, verbose=1),
                 ModelCheckpoint('UNET_CMIP6_piCTL_UV850_FineTuning_ERA519401969', monitor='val_loss', verbose=1, save_best_only=True),History()]

tf.config.run_functions_eagerly(True)

history = model.fit(dataset_train,epochs=100,callbacks = callback, validation_data=dataset_val)







