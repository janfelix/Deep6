###Author: Jan Felix Finke, 2022
###See license.txt for licensing



###create environment and install packages
#using virtual environment
module load python/3.6
virtualenv --no-download ~/dsix #create the empty environment
source ~/dsix/bin/activate
pip install --no-index --upgrade pip
pip install keras==2.2.4 numpy scipy pandas sklearn biopython
pip install --no-index tensorflow_gpu
pip install 'h5py==2.10.0' --force-reinstall
pip install 'scipy==1.4.1' --force-reinstall
#using conda
conda create --name dsix --file dsix_env.txt
#or
conda env create -f dsix.yml



###Useage:
#predictsequences with default model, e.g. minimum legnth set to 250nt
python Master/deep6.py -i inputfile.fasta -l 250 -m Models -o outdir
#Options:
  -h, --help            show this help message and exit
  -i IFILE, --infile=IFILE
                        input file
  -l MLENGTH, --minlength=MLENGTH
                        min length
  -m MDIR, --mod=MDIR   model directory
  -o ODIR, --out=ODIR   output directory



#batch encode sequences for custom training datasets, e.g. duplodnaviria with length set to 250nt
python Master/deep6_encode.py -i .../inputfile -l 250 -c duplodnaviria
#Options:
  -h, --help            show this help message and exit
  -i IFILE, --infile=IFILE
                        input file
  -l MLENGTH, --length=MLENGTH
                        fragment length
  -c CCLASS, --contigclass=CCLASS
                        class: prokaryote, eukaryote, duplodnaviria,
                        varidnaviria, monodnaviria, riboviria

#train custom model, e.g. fragent length set to 250nt, 525 nerons, kernel size 10 and 40 epochs
python Master/deep6_train.py -l 250 -t .../train_encode -v .../val_encode -o outdir -n 525 -k 10 -e 40

#Options:
  -h, --help            show this help message and exit
  -l MLENGTH, --length=MLENGTH
                        fragment length
  -t TDIR, --train=TDIR
                        directory for training data
  -v VDIR, --val=VDIR   directory for validation data
  -o ODIR, --out=ODIR   output directory
  -n NLAYER, --nlayer=NLAYER
                        neurons dense layer
  -k KSIZE, --ksize=KSIZE
                        kernel size
  -e EPOCHS, --epochs=EPOCHS
                        max number of epochs



###See reference_taxa.txt for overview of taxa for training datasets
