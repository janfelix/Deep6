###Jan Felix Finke, 2022
###deep6 version 1
###License: GNU AGPLv3
#!/usr/bin/env python3

import os, sys, optparse, warnings
import h5py
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import tensorflow as tf 
from tensorflow import keras
from tensorflow.keras.models import load_model
os.environ['KERAS_BACKEND'] = 'theano'
import deep6_functions as d6

###input arguments
inargs = optparse.OptionParser()
inargs.add_option("-i", "--infile", action = "store", type = "string", dest = "ifile", help = "input file")
inargs.add_option("-l", "--minlength", action = "store", type = int, dest = "mlength", help = "min length")
inargs.add_option("-m", "--mod", action = "store", type = "string", dest = "mdir", help = "model directory")
inargs.add_option("-o", "--out", action = "store", type = "string", dest = "odir", help = "output directory")
(options, args) = inargs.parse_args()
if (options.ifile is None or options.mlength is None or options.mdir is None or options.odir is None) :
	sys.stderr.write("Input arg missing:")
	inargs.print_help()
	sys.exit(0)

#set variables
ifile = options.ifile
mlength = options.mlength
mdir = options.mdir
odir = options.odir
if not os.path.exists(odir):
    os.makedirs(odir)

#load models and add to function module
print("load models")
mdict = {}
warnings.filterwarnings('ignore', 'Error in loading the saved optimizer ')
for clength in ['250', '500', '1000', '1500'] :
    mname = [ x for x in os.listdir(mdir) if 'model_' in x and clength in x and x.endswith(".h5") ][0]
    mdict[clength] = load_model(os.path.join(mdir, mname))
d6.mdict = mdict

#create output file and add to function module
print("create output file")
outfile = os.path.join(odir, os.path.basename(ifile)+'_predict_'+str(mlength)+'bp_deep6.txt')
results = open(outfile, 'w')
writef = results.write('\t'.join(['name', 'length', 'duplo', 'euk', 'mono', 'pro', 'ribo', 'vari'])+'\n')
results.close()
d6.outfile = outfile

#load sequences
print("load sequences")
record = list(SeqIO.parse(ifile, "fasta"))

#encode and predict contigs in batches of 100
print("start predicting")
cfw = []
cbw = []
sname = []
for seqx in range(len(record)):
    if len(str(record[seqx].seq)) < mlength:
        continue
    cfw.append(d6.encodex(str(record[seqx].seq)))
    cbw.append(d6.encodex(str(record[seqx].seq.reverse_complement())))
    sname.append(str(record[seqx].id)) 
    if len(cfw) > 0 and len(cfw) % 100 == 0:
        z=list(zip(sname,cfw,cbw))
        resx = list(map(d6.predx, z))
        cfw = []
        cbw = []
        sname = []
        z = None

#encode and predict the last batch of contigs
if len(cfw) > 0:
    z=list(zip(sname,cfw,cbw))
    resx = list(map(d6.predx, z))
    cfw = []
    cbw = []
    sname = []
    z = None
print("done predicting")
