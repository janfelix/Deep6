#### Author: Jan Felix Finke, 2022

#### See [LICENSE](https://github.com/janfelix/Deep6/blob/main/LICENSE) for licensing.

### Description:
Deep6 is a reference-independent and alignment-free tool to predict viral sequences in short metatrascriptome data, it can process input contigs as short as 250 nt. Deep6 is a Convolutional Neural Networ (CNN) with 500 convolutions and 525 dense layers. Four models for different length ranges of query sequences are provided, the repository also includes scripts to custom train alternative models. The primary prediction script selects the appropriate model for each input sequence length and automatically encodes and predicts the sequence. For each sequence group scores are calculated, based on those the final sequence prediction is derived in the downstream analysis. 

To train custom models, the batch-encoding script splits sets of CDS into 90% training data and 10% validation data, chunks of CDS of appropriate length for the model in overlapping forward and reverse order are then one-hot encoded. In a next step the training script feeds the encoded training and validation data into the CNN and saves the best model. During training models are optimized for accuracy; model performance is assessed by training and validation area under receiver operating characteristic curve (AUROC), average accuracy, and group precision, recall and derived F1-scores. Additionally, a confusion matrix of sequence predictions is produced for the final model.


### Setup Environment and Install Packages:
Clone the "Deep6" repository using the code button on top or in your terminal using: `git clone https://github.com/janfelix/Deep6.git`. The Master directory includes the license file, this readme file, a conda yml file and the python scripts for sequence prediction, batch encoding and custom model training. The provided R-script processes group scores per sequence into a final prediction for easy down stream processing. The Model directory contains the pre-trained default models for marine samples. 

The installation instructions are focused on Linux systems and use Python 3.6, numpy, pandas, h5py, biopython, scipy, keras, tensorflow-gpu and scikit-learn. The tensorflow-gpu package can be difficult to implement for non Linux systems, custom model training when using different versions of tensorflow-gpu might be impaired. Especially for custom model training it is advised to use gpu equipped servers, any reasonably scaled model training will exceed cpu systems.

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

For Linux systems:

Create a conda environment and install required packages. This is the recommended method, but can take a few minutes. Conda will need the following channels: default, anaconda, conda-forge. If necessary add channels e.g. `conda config --add channels conda-forge`. This also works for Windows 64-bit systems.

`conda create -n dsix python=3.6 numpy pandas h5py biopython scipy keras scikit-learn tensorflow-gpu`

Or using the yml setup file for the specific package versions.

`conda env create -n dsix -f dsix_linux.yml` 

For Mac OSX systems:

There is no conda installer for tensorflow-gpu 2.6.0 for Mac OSX, the following installs tensorflow-gpu 1.1.0 through pypi. Conda will need the following channels: default, anaconda, conda-forge, pypi. If necessary add channels e.g. `conda config --add channels conda-forge`.
```
conda create -n dsix python=3.6 numpy pandas h5py biopython scipy keras scikit-learn
conda activate dsix
pip install tensorflow tensorflow_gpu
```
Or using the yml setup file for the specific package versions, using pypi.

`conda env create -n dsix -f dsix_mac.yml`

### Usage:

#### Predict sequences with default model, e.g. with minimum length set to 250nt

`python Deep6/Master/deep6.py -i ./inputfile.fasta -l 250 -m Deep6/Models -o outdir`

Options: -i defines the input file with contigs in fasta format; -l is the minimum contig length to analyze; -m defines the directory with the models; -o defines the directory to save the outputfile.
```
  -h, --help            show this help message and exit
  -i IFILE, --infile=IFILE
                        input file
  -l MLENGTH, --minlength=MLENGTH
                        min length
  -m MDIR, --mod=MDIR   model directory
  -o ODIR, --out=ODIR   output directory
```

#### Batch encode sequences for custom training datasets, e.g. duplodnaviria with fragment length set to 250nt

`python Deep6/Master/deep6_encode.py -i Deep6/Master/test_data.fasta -l 250 -c duplodnaviria`

Options: -i defines the input file with sequences to encode in fasta format, one file per training class; -l is the fragment length to encode, -c denotes the corresponding class; encoded files for training and validation are saved in the input file directories.
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

#### Train custom model, e.g. fragment length set to 250nt, 525 neurons, kernel size 10 for 40 epochs

`python Deep6/Master/deep6_train.py -l 250 -t ./train_encode -v ./val_encode -o outdir -n 525 -k 10 -e 40`

Options: -l defines the fragment length for the specific model, corrsponding to the fragment lengths defined during batch encoding; -t defines the directory for training data; -v defines the directory for validation data; -n defines the number of dense layers (default 525); -e defines the maximum number of epochs to run model optimization.
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

### Test run:

To test the installation and demonstrate the function of Deep6 analyze the provided [test_data.fasta](https://github.com/janfelix/Deep6/blob/main/Master/test_data.fasta) file. This file contains one reference sequences per prokaryote, eukaryote, duplodnaviria, varidnaviria, monodnaviria, riboviria.

`python Deep6/Master/deep6.py -i Deep6/Master/test_data.fasta -l 250 -m Deep6/Models -o ./`

### Interpreting Results:

The provided R-script [deep6_interpetration.R](https://github.com/janfelix/Deep6/blob/main/Master/deep6_interpetration.R) assigns group prediction for each sequence based on the group with the highest score, if the score also is 1.25x the median score. It also produces an overview of the contig prediction for the analysed input file. It should be noted that this is only a suggested approach to interpret the predictions. Depending on the research question it also feasible to only differentiate between e.g. any cellular vs. any viral or to consider e.g. second highest scores.

### Supplementary files:
 
See [refseq_taxa.txt](https://github.com/janfelix/Deep6/blob/main/Master/refseq_taxa.txt) for an overview of taxa for training datasets.
