import tensorflow as tf
from tensorflow.keras import layers, models, optimizers
import numpy as np
import xarray as xr
from tqdm import tqdm
import os
from tensorflow.keras.models import load_model, Model
from tensorflow.keras.layers import Conv2D, Conv3D, MaxPooling2D, MaxPooling3D, Flatten, Dense, UpSampling2D, Conv2DTranspose,Conv3DTranspose, Reshape, concatenate, BatchNormalization, Activation
from tensorflow.keras.layers import Input,  LeakyReLU, Concatenate, Dropout, Conv1D
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.callbacks import ReduceLROnPlateau, EarlyStopping, ModelCheckpoint, History
from tensorflow.keras.models import Sequential
from tensorflow.keras.regularizers import l2
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
from tensorflow.keras.models import clone_model


#GPU
gpus = tf.config.experimental.list_physical_devices('GPU')
if gpus:
    tf.config.experimental.set_visible_devices(gpus[0], 'GPU')
    tf.config.experimental.set_memory_growth(gpus[0], True)

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


era5_tas = xr.open_dataset('Tano_era5_day_tas_1940_1969_1.4x1.4_s30.nc').tas
era5_U = xr.open_dataset('era5_day_ua_850_1940_1969_1.4x1.4_CR19401969.nc').ua[:,0,:,:]
era5_V = xr.open_dataset('era5_day_va_850_1940_1969_1.4x1.4_CR19401969.nc').va[:,0,:,:]

input_stack = np.stack((era5_U, era5_V), axis=-1)

base_model = load_model('Training_CMIP6_piCTL_all_UV850_base')
model = clone_model(base_model)
model.set_weights(base_model.get_weights())
model.compile(loss='mse',optimizer=Adam(learning_rate=1e-3),metrics=['mse'])

callback = [ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=4, verbose=1), EarlyStopping(monitor='val_loss', patience=10, verbose=1),
		 ModelCheckpoint('Training_CMIP6_piCTL_all_UV850_FineTuning_ERA519401969_V2', monitor='val_loss', verbose=1, save_best_only=True),History()]

X_train_era,X_test_era,Y_train_era,Y_test_era=train_test_split(input_stack,era5_tas,test_size=0.1)

Y_train_era=np.expand_dims(Y_train_era,axis=-1)
Y_test_era=np.expand_dims(Y_test_era,axis=-1) 

tf.config.run_functions_eagerly(True)
history =  model.fit(x=X_train_era, y=Y_train_era, batch_size=128,epochs=100,shuffle=True,callbacks=callback, validation_data=(X_test_era,Y_test_era))
