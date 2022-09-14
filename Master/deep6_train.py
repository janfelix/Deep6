###Jan Felix Finke, 2022
###deep6 version 1
###License: GNU AGPLv3
#!/usr/bin/env python3

import os
import sys
import optparse
import numpy as np
import tensorflow as tf 
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
from tensorflow.keras import layers
from tensorflow.keras.layers import Average
from tensorflow.keras.callbacks import ModelCheckpoint, EarlyStopping
from tensorflow.keras.optimizers import Adam
import h5py
import sklearn
from sklearn.metrics import roc_auc_score 
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from deep6_functions import layerx

###input arguments
inargs = optparse.OptionParser()
inargs.add_option("-l", "--length", action = "store", type = int, dest = "mlength", help = "fragment length")
inargs.add_option("-t", "--train", action = "store", type = "string", dest = "tdir", help = "directory for training data")
inargs.add_option("-v", "--val", action = "store", type = "string", dest = "vdir", help = "directory for validation data")
inargs.add_option("-o", "--out", action = "store", type = "string", dest = "odir", help = "output directory")
inargs.add_option("-n", "--nlayer", action = "store", type = int, dest = "nlayer", default=525, help = "neurons dense layer")
inargs.add_option("-k", "--ksize", action = "store", type = int, dest = "ksize", default=10, help = "kernel size")
inargs.add_option("-e", "--epochs", action = "store", type = int, dest = "epochs", default=40, help = "max number of epochs")
(options, args) = inargs.parse_args()
if (options.mlength is None or options.tdir is None or options.vdir is None or options.odir is None ) :
	sys.stderr.write("Input arg missing:")
	inargs.print_help()
	sys.exit(0)

#set variables
mlength = options.mlength
tdir = options.tdir
vdir = options.vdir
odir = options.odir
if not os.path.exists(odir):
    os.makedirs(odir)
nlayer = options.nlayer
ksize = options.ksize
epochs = options.epochs

###load and prep data
print("loading training data from: "+str(tdir))
train_cfw = [ x for x in os.listdir(tdir) if 'eukaryote' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
eukRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
eukRef_trcbw = np.load(os.path.join(tdir, train_cbw))
train_cfw = [ x for x in os.listdir(tdir) if 'prokaryote' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
proRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
proRef_trcbw = np.load(os.path.join(tdir, train_cbw))
train_cfw = [ x for x in os.listdir(tdir) if 'duplodnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
duploRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
duploRef_trcbw = np.load(os.path.join(tdir, train_cbw))
train_cfw = [ x for x in os.listdir(tdir) if 'varidnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
variRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
variRef_trcbw = np.load(os.path.join(tdir, train_cbw))
train_cfw = [ x for x in os.listdir(tdir) if 'monodnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
monoRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
monoRef_trcbw = np.load(os.path.join(tdir, train_cbw))
train_cfw = [ x for x in os.listdir(tdir) if 'riboviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
riboRef_trcfw = np.load(os.path.join(tdir, train_cfw))
train_cbw = train_cfw.replace('cfw', 'cbw')
riboRef_trcbw = np.load(os.path.join(tdir, train_cbw))

print("concatenate training data")
con_trcfw = np.concatenate((eukRef_trcfw, proRef_trcfw, duploRef_trcfw, variRef_trcfw, monoRef_trcfw, riboRef_trcfw), axis=0)
con_trcbw = np.concatenate((eukRef_trcbw, proRef_trcbw, duploRef_trcbw, variRef_trcbw, monoRef_trcbw, riboRef_trcbw), axis=0)
labels_trcfw = ["euk"] * eukRef_trcfw.shape[0] + ["pro"] * proRef_trcfw.shape[0] + ["duplo"] * duploRef_trcfw.shape[0] + ["vari"] * variRef_trcfw.shape[0] + ["mono"] * monoRef_trcfw.shape[0] + ["ribo"] * riboRef_trcfw.shape[0]
del eukRef_trcfw, proRef_trcfw, duploRef_trcfw, variRef_trcfw, monoRef_trcfw, riboRef_trcfw
del eukRef_trcbw, proRef_trcbw, duploRef_trcbw, variRef_trcbw, monoRef_trcbw, riboRef_trcbw

#build and print to output binary classification labels
lb = preprocessing.LabelBinarizer()
lb_trcfw = lb.fit_transform(labels_trcfw)
print('training labels class order: '+str(lb.classes_)+'\n')

print("shuffling training data")
index_trcfw = list(range(0, con_trcfw.shape[0]))
np.random.seed(0)
np.random.shuffle(index_trcfw)
con_trcfw_shuf = con_trcfw[np.ix_(index_trcfw, range(con_trcfw.shape[1]), range(con_trcfw.shape[2]))]
del con_trcfw
con_trcbw_shuf = con_trcbw[np.ix_(index_trcfw, range(con_trcbw.shape[1]), range(con_trcbw.shape[2]))]
del con_trcbw
lb_trcfw_shuf = lb_trcfw[index_trcfw]
del lb_trcfw

print("loading validation data from: "+str(vdir))
val_cfw = [ x for x in os.listdir(vdir) if 'eukaryote' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
eukRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
eukRef_valcbw = np.load(os.path.join(vdir, val_cbw))
val_cfw = [ x for x in os.listdir(vdir) if 'prokaryote' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
proRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
proRef_valcbw = np.load(os.path.join(vdir, val_cbw))
val_cfw = [ x for x in os.listdir(vdir) if 'duplodnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
duploRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
duploRef_valcbw = np.load(os.path.join(vdir, val_cbw))
val_cfw = [ x for x in os.listdir(vdir) if 'varidnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
variRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
variRef_valcbw = np.load(os.path.join(vdir, val_cbw))
val_cfw = [ x for x in os.listdir(vdir) if 'monodnaviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
monoRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
monoRef_valcbw = np.load(os.path.join(vdir, val_cbw))
val_cfw = [ x for x in os.listdir(vdir) if 'riboviria' in x and '_length:'+str(mlength)  in x and 'cfw.npy' in x][0]
riboRef_valcfw = np.load(os.path.join(vdir, val_cfw))
val_cbw = val_cfw.replace('cfw', 'cbw')
riboRef_valcbw = np.load(os.path.join(vdir, val_cbw))

print("concatenate validation data")
con_valcfw = np.concatenate((eukRef_valcfw, proRef_valcfw, duploRef_valcfw, variRef_valcfw, monoRef_valcfw, riboRef_valcfw), axis=0)
con_valcbw= np.concatenate((eukRef_valcbw, proRef_valcbw, duploRef_valcbw, variRef_valcbw, monoRef_valcbw, riboRef_valcbw), axis=0)
labels_valcfw = ["euk"] * eukRef_valcfw.shape[0] + ["pro"] * proRef_valcfw.shape[0] + ["duplo"] * duploRef_valcfw.shape[0] + ["vari"] * variRef_valcfw.shape[0] +["mono"] * monoRef_valcfw.shape[0] + ["ribo"] * riboRef_valcfw.shape[0]
del eukRef_valcfw, proRef_valcfw, duploRef_valcfw, variRef_valcfw, monoRef_valcfw, riboRef_valcfw
del eukRef_valcbw, proRef_valcbw, duploRef_valcbw, variRef_valcbw, monoRef_valcbw, riboRef_valcbw

#build and print to output binary classification labels
lb = preprocessing.LabelBinarizer()
lb_valcfw = lb.fit_transform(labels_valcfw)
print('validation label class order: '+str(lb.classes_)+'\n')

#build model
mname = os.path.join( odir, 'model_'+str(mlength)+'_fn'+str(500)+'_dn'+str(nlayer)+'_ks'+str(ksize)+'.h5')
if os.path.isfile(mname):
    sys.exit("model name already exists")
else :
    input_fw = layers.Input(shape=(None, 4))
    input_rv = layers.Input(shape=(None, 4))
    hlayers = [
        layers.Conv1D(filters = 500, kernel_size=ksize, activation='relu'),
        layers.GlobalMaxPooling1D(),
        layers.Dropout(0.1),
        layers.Dense(nlayer, activation='relu'),
        layers.Dropout(0.1),
        layers.Dense(6, activation='softmax')
    ]
    output_fw = layerx(input_fw, hlayers)     
    output_rv = layerx(input_rv, hlayers)
    output = Average()([output_fw, output_rv])
    model = Model(inputs=[input_fw, input_rv], outputs=output)
    model.compile(Adam(lr=0.001), loss='categorical_crossentropy', metrics=['accuracy'])
    model.summary()

batchs = int(con_trcfw_shuf.shape[0]/(1000000/mlength)) #based on data set size and mlength (min.= 1)
if batchs ==0:
    batchs=1
checkpoint = keras.callbacks.ModelCheckpoint(filepath=mname, verbose=1,save_best_only=True)
earlystopping = keras.callbacks.EarlyStopping(monitor='val_accuracy', min_delta=0.0001, patience=6, verbose=1)

#fit model
print("fitting model: "+str(mlength)+'_fn'+str(500)+'_dn'+str(nlayer)+'_ks'+str(ksize)+'_ep'+str(epochs))
model.fit(x = [con_trcfw_shuf, con_trcbw_shuf], y = lb_trcfw_shuf, \
            batch_size=batchs, epochs=epochs, verbose=2, \
            validation_data=([con_valcfw, con_valcbw], lb_valcfw), \
            callbacks=[checkpoint, earlystopping])
          
###evaluation of AUC
print('model evaluation, class order: '+str(lb.classes_)+'and averaged auroc scores:'+'\n')

# train data
subset = 'train'
pred_tr = model.predict([con_trcfw_shuf, con_trcbw_shuf], batch_size=1)
auc = sklearn.metrics.roc_auc_score(lb_trcfw_shuf, pred_tr, multi_class='ovo')
print('auroc_'+subset+'='+str(auc)+'\n')

# val data
subset = 'val'
pred_val = model.predict([con_valcfw, con_valcbw], batch_size=1)
auc = sklearn.metrics.roc_auc_score(lb_valcfw, pred_val, multi_class='ovo')
print('auroc_'+subset+'='+str(auc)+'\n')
np.savetxt(os.path.join(mname.replace('.h5', '_val_pred.txt')), np.transpose(pred_val))
np.savetxt(os.path.join(mname.replace('.h5', '_val_true.txt')), np.transpose(lb_valcfw))

#get confusion matrix, and precision and recall, must convert labels to numeric
class_true = np.argmax(lb_valcfw, axis=1)
class_predicted = np.argmax(pred_val, axis=1)
print(sklearn.metrics.confusion_matrix(class_true, class_predicted))
print(sklearn.metrics.classification_report(class_true, class_predicted, digits=3))
print("deep down!")