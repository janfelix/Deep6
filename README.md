Author: Jan Felix Finke, 2022

See [LICENSE](https://github.com/janfelix/Deep6/blob/main/LICENSE) for licensing

### Setup Environment and Install Packages
These instructions are focused on Linux based servers, especially the Tensorflow-gpu package can be difficult to implement for other systems. The sequence prediction is fully functional, custom model training when using different versions of tensorflow-gpu might be impaired.

#### Using virtual environment
```
module load python/3.6
virtualenv --no-download ~/dsix #create the empty environment
source ~/dsix/bin/activate

pip install --no-index --upgrade pip
pip install keras==2.2.4 numpy scipy pandas sklearn biopython
pip install --no-index tensorflow_gpu
pip install 'h5py==2.10.0' --force-reinstall
pip install 'scipy==1.4.1' --force-reinstall
```
#### Using `conda`

For Linux based servers:

Create a conda environment and install required packages. This is the recommended method, but can take a few minutes:

`conda create -n dsix python=3.6 numpy pandas h5py biopython scipy keras scikit-learn tensorflow-gpu`

Or using the yml setup file for the specific package version. This can create sytem based conflicts:

`conda create -n dsix -f dsix.yml` 

For Mac OSX systems:

There is no conda installer for tensorflow-gpu 2.6.0 for Mac OSX, the following installs tesnorflow-gpu 1.1.0 through pypi:
```
conda create -n dsix python=3.6 numpy pandas h5py biopython scipy keras scikit-learn
conda activate dsix
pip3 install tensorflow tensorflow_gpu
```

### Usage:

Predict sequences with default model, e.g. minimum legnth set to 250nt

`python Master/deep6.py -i inputfile.fasta -l 250 -m Models -o outdir`
##### Options:
```
  -h, --help            show this help message and exit
  -i IFILE, --infile=IFILE
                        input file
  -l MLENGTH, --minlength=MLENGTH
                        min length
  -m MDIR, --mod=MDIR   model directory
  -o ODIR, --out=ODIR   output directory
```

Batch encode sequences for custom training datasets, e.g. duplodnaviria with length set to 250nt

`python Master/deep6_encode.py -i ../inputfile -l 250 -c duplodnaviria`
##### Options:
```
  -h, --help            show this help message and exit
  -i IFILE, --infile=IFILE
                        input file
  -l MLENGTH, --length=MLENGTH
                        fragment length
  -c CCLASS, --contigclass=CCLASS
                        class: prokaryote, eukaryote, duplodnaviria,
                        varidnaviria, monodnaviria, riboviria
```
Train custom model, e.g. fragment length set to 250nt, 525 neurons, kernel size 10 for 40 epochs

`python Master/deep6_train.py -l 250 -t ../train_encode -v ../val_encode -o outdir -n 525 -k 10 -e 40`

##### Options:
```
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
```

#### See [refseq_taxa.txt](https://github.com/janfelix/Deep6/blob/main/Master/refseq_taxa.txt) for overview of taxa for training datasets
