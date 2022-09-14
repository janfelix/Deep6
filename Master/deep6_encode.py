###Jan Felix Finke, 2022
###deep6 version 1
##License: GNU AGPLv3
#!/usr/bin/env python3

import os, sys
import numpy as np
import optparse
import random
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from deep6_functions import encodex

###input arguments
inargs = optparse.OptionParser()
inargs.add_option("-i", "--infile", action = "store", type = "string", dest = "ifile", help = "input file")
inargs.add_option("-l", "--length", action = "store", type = int, dest = "mlength", help = "fragment length")
inargs.add_option("-c", "--contigclass", action = "store", type = "string", dest = "cclass", help = "class: prokaryote, eukaryote, duplodnaviria, varidnaviria, monodnaviria, riboviria")
(options, args) = inargs.parse_args()
if (options.ifile is None or options.mlength is None or options.cclass is None ) :
	sys.stderr.write("Input arg missing:")
	inargs.print_help()
	sys.exit(0)

#set variables
ifile = options.ifile
mlength = options.mlength
cclass = options.cclass
fname = os.path.splitext((os.path.basename(ifile)))[0]
idir = os.path.dirname(ifile)

#create train and val data directories
tdir = os.path.join(idir, "train_encode")
if not os.path.exists(tdir):
    os.makedirs(tdir)
vdir = os.path.join(idir, "val_encode")
if not os.path.exists(vdir):
    os.makedirs(vdir)

#open file and encode by blocks of mlength, forward, reverse complement and for a fasta file of sequence blocks and seq labels
record = list(SeqIO.parse(ifile, "fasta"))

#index and split 90/10 into training and validation data
x = list(range(0, len(record)))
random.shuffle(x)
y = int(len(record)*0.9)
subset1 = x[0:y]
subset2 = x[y:]

#Training data
subset = "train"
cfw = []
cbw = []
sname = []
filenum = 0
#encode and contigs in exact fragments of contiglength, trim end for fw and beginning for revcomp sqeuences, save in batches of 2000000
for seqx in subset1:
    if len(str(record[seqx].seq)) < mlength:
        continue
    seqfw=str(record[seqx].seq)[0:(mlength*(len(str(record[seqx].seq))//mlength))]
    seqbw=str(record[seqx].seq.reverse_complement())[0:(mlength*(len(str(record[seqx].seq))//mlength))]
    cfw.extend(list(map(encodex, [seqfw[i:i+mlength] for i in range(0,len(seqfw),mlength)])))
    cbw.extend(list(map(encodex, [seqbw[i:i+mlength] for i in range(0,len(seqbw),mlength)])))
    for s in range(len([str(record[seqx].seq)[i:i+mlength] for i in range(0,len(str(record[seqx].seq)),mlength)])):
        sname.append(SeqRecord(Seq([str(record[seqx].seq)[i:i+mlength] for i in range(0,len(str(record[seqx].seq)),mlength)][s]),id= record[seqx].id)) 
    if len(cfw) > 0 and len(cfw) % 2000000 == 0:
        filenum += 1
        cfwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cfw.npy"
        np.save(os.path.join(tdir, cfwname), np.array(cfw))
        cbwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cbw.npy"
        np.save(os.path.join(tdir, cbwname), np.array(cbw))
        fastaname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+".fasta"
        SeqIO.write(sname, os.path.join(tdir, fastaname), "fasta")
        cfw = []
        cbw = []
        sname = []
#encode and save the last batch of contigs
if len(cfw) > 0:    
    filenum += 1
    cfwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cfw.npy"
    np.save(os.path.join(tdir, cfwname), np.array(cfw))
    cbwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cbw.npy"
    np.save(os.path.join(tdir, cbwname), np.array(cbw))
    fastaname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+".fasta"
    SeqIO.write(sname, os.path.join(tdir, fastaname), "fasta")
    cfw = []
    cbw = []
    sname = []



#Validation data
subset = "val"
cfw = []
cbw = []
sname = []
filenum = 0
#encode and contigs in exact fragments of contiglength, trim end for fw and beginning for revcomp sqeuences, save in batches of 2000000
for seqx in subset2:
    if len(str(record[seqx].seq)) < mlength:
        continue
    seqfw=str(record[seqx].seq)[0:(mlength*(len(str(record[seqx].seq))//mlength))]
    seqbw=str(record[seqx].seq.reverse_complement())[0:(mlength*(len(str(record[seqx].seq))//mlength))]
    cfw.extend(list(map(encodex, [seqfw[i:i+mlength] for i in range(0,len(seqfw),mlength)])))
    cbw.extend(list(map(encodex, [seqbw[i:i+mlength] for i in range(0,len(seqbw),mlength)])))
    for s in range(len([str(record[seqx].seq)[i:i+mlength] for i in range(0,len(str(record[seqx].seq)),mlength)])):
        sname.append(SeqRecord(Seq([str(record[seqx].seq)[i:i+mlength] for i in range(0,len(str(record[seqx].seq)),mlength)][s]),id= record[seqx].id)) 
    if len(cfw) > 0 and len(cfw) % 2000000 == 0:
        filenum += 1
        cfwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cfw.npy"
        np.save(os.path.join(vdir, cfwname), np.array(cfw))
        cbwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cbw.npy"
        np.save(os.path.join(vdir, cbwname), np.array(cbw))
        fastaname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+".fasta"
        SeqIO.write(sname, os.path.join(vdir, fastaname), "fasta")
        cfw = []
        cbw = []
        sname = []
#encode and save the last batch of contigs
if len(cfw) > 0:    
    filenum += 1
    cfwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cfw.npy"
    np.save(os.path.join(vdir, cfwname), np.array(cfw))
    cbwname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+"_cbw.npy"
    np.save(os.path.join(vdir, cbwname), np.array(cbw))
    fastaname = cclass+"_"+subset+"_"+fname+"#"+str(filenum)+"_length:"+str(mlength)+"_seqs:"+str(len(cfw))+".fasta"
    SeqIO.write(sname, os.path.join(vdir, fastaname), "fasta")
    cfw = []
    cbw = []
    sname = []