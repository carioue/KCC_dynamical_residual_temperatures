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
    
#Models=['BCC-CSM2-MR','CanESM5','CESM2','CESM2-WACCM','CMCC-CM2-SR5','CMCC-ESM2','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','EC-Earth3-CC','EC-Earth3-Veg','EC-Earth3-Veg-LR',
#	'FGOALS-g3','IITM-ESM','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','MIROC-ES2L','MIROC6','MPI-ESM1-2-LR','TaiESM1','UKESM1-0-LL']

Models=['CNRM-CM6-1']
for m in Models :
	print(m)

	tas = xr.open_dataset('Tano_day_'+str(m)+'_1.3x1.3_noleapday_1850_1899_allmembers.nc').tas
	ua = xr.open_dataset('ua_day_'+str(m)+'_1.3x1.3_noleapday_1850_1899_allmembers_fillmiss.nc').ua[:,0,:,:]
	va = xr.open_dataset('va_day_'+str(m)+'_1.3x1.3_noleapday_1850_1899_allmembers_fillmiss.nc').va[:,0,:,:]
    
	wgt=np.cos(np.deg2rad(ua.lat))
	fldmean_ua = ua.weighted(wgt).mean(dim=['lon','lat'])
	timmean_fldmean_ua = fldmean_ua.mean(dim='time')
	fldstd_ua = ua.weighted(wgt).std(dim=['lon','lat'])
	timstd_fldstd_ua = fldstd_ua.std(dim='time')

	fldmean_va = va.weighted(wgt).mean(dim=['lon','lat'])
	timmean_fldmean_va = fldmean_va.mean(dim='time')
	fldstd_va = va.weighted(wgt).std(dim=['lon','lat'])
	timstd_fldstd_va = fldstd_va.std(dim='time')

	ua_norm = (ua-timmean_fldmean_ua)/timstd_fldstd_ua
	va_norm = (va-timmean_fldmean_va)/timstd_fldstd_va

	input_stack = np.stack((ua_norm, va_norm), axis=-1)

	base_model = load_model('Training_CMIP6_piCTL_all_UV850_base')
	model = clone_model(base_model)
	model.set_weights(base_model.get_weights())
	model.compile(loss='mse',optimizer=Adam(learning_rate=1e-3),metrics=['mse'])

	callback = [ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=4, verbose=1), EarlyStopping(monitor='val_loss', patience=8, verbose=1),
			 ModelCheckpoint('UNET_'+str(m)+'_1850_1899', monitor='val_loss', verbose=1, save_best_only=True),History()]

	X_train,X_test,Y_train,Y_test=train_test_split(input_stack,tas,test_size=0.2)

	Y_train=np.expand_dims(Y_train,axis=-1)
	Y_test=np.expand_dims(Y_test,axis=-1) 
	    
	taille_data= len(tas)
	if taille_data<20000:
		batch=128
	if taille_data>=20000 and taille_data<100000:
		batch=512
	if taille_data>=100000  and taille_data<200000:
		batch=1024
	if taille_data>=200000:
		batch=2048

	tas.close()
	ua.close()
	va.close()

	tf.config.run_functions_eagerly(True)
	history =  model.fit(x=X_train, y=Y_train, batch_size=batch,epochs=100,shuffle=True,callbacks=callback, validation_data=(X_test,Y_test))
	
	train_loss = history.history['loss']
	val_loss = history.history['val_loss']

	try:
	# Tracer les courbes de la loss d'entraînement et de validation
	    plt.figure(figsize=(12, 8))
	    plt.plot(train_loss, label='Training Loss')
	    plt.plot(val_loss, label='Validation Loss')
	    plt.xlabel('Epochs')
	    plt.ylabel('Loss')
	    plt.title('Training and Validation Loss Over Epochs')
	    plt.legend()
	    plt.grid(True)

	# Enregistrer l'image dans un fichier
	    plt.savefig('UNET_'+str(m)+'_1850_1899/training_validation_loss.pdf')
	except:
	    print('pb fig')

	
	

