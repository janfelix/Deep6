###Jan Felix Finke, 2022
###deep6 version 1
##License: GNU AGPLv3

#encode function 
def encodex(seq) : 
    code = list()
    dict = {'A':[1,0,0,0], 'a':[1,0,0,0],'C':[0,1,0,0], 'c':[0,1,0,0],'G':[0,0,1,0], 'g':[0,0,1,0],'T':[0,0,0,1], 't':[0,0,0,1]}
    for n in seq:
         if n in dict.keys():
            code.append(dict[n])
         else:
            code.append([1/4, 1/4, 1/4, 1/4])
    return code

#function to create model output layers
def layerx(input, hidden):
    for h in hidden:
        input = h(input)
    return input

#function for predicting class scores, mdict and outfile is added during prediction
import numpy as np

def predx(cont) :
    print('predicting: '+str(cont[0]))
    if len(cont[1]) < 500 :
        models = mdict['250']
    elif len(cont[1]) < 1000 and len(cont[1]) >= 500 :
        models = mdict['500']
    elif len(cont[1]) < 1500 and len(cont[1]) >= 1000:
        models = mdict['1000']
    else:
        models = mdict['1500']
    scores = models.predict([np.array([cont[1]]), np.array([cont[2]])], batch_size=1)[0]
        
    results = open(outfile, 'a')
    writef = results.write('\t'.join([str(cont[0]), str(len(cont[1])), str(float(scores[0])), str(float(scores[1])), str(float(scores[2])), str(float(scores[3])), str(float(scores[4])), str(float(scores[5]))])+'\n')
    flushf = results.flush()
    results.close()
    
    return [str(cont[0]), str(len(cont[1])), float(scores[0]), float(scores[1]), float(scores[2]), float(scores[3]), float(scores[4]), float(scores[5])] 