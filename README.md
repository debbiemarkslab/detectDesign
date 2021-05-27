# About
You can find single Cas9 (CRISPR associated protein) guides by just looking for potential binding sites that meet minimum requirements for binding to DNA. The problem detectDesign solves is finding Cas9 guides for use in a paired dCas9-kinase sensor system. The software can find paired sgRNA protospacers within linker range on a given target and rank the pairs based on matches to non-intended protospacers.

Full user manual with example: https://docs.google.com/document/d/1syjWFL8aQuPjSYoZhMGAh0bqcl6YhPHQR2l0A6avz5M/edit?usp=sharing

# Installation
## Required software
detectDesign was developed using Python 3.6.8, Python 2 is not supported.
- Python 3.6+


## Required Python packages
matplotlib==3.0.3
numpy==1.17.3
pandas==1.0.1 
seaborn==0.8.1 
regex==2018.8.29 


## Cloud computing setup (SLURM/O2)
To move files into O2
Get filezilla: https://filezilla-project.org/
Host info: https://wiki.rc.hms.harvard.edu/display/O2/File+Transfer

##### Open interactive node (8G was good for setting up the env)
`srun --pty -p interactive --mem 8G -t 0-06:00 /bin/bash`

##### Load conda package to use package manager commands
`module load conda2`

##### Creates a conda environment that has the packages necessary (skip if done before)
`conda create -n detectDesign-env-04_07_2021 python=3.7 matplotlib==3.0.3 seaborn==0.8.1 pandas==1.0.1 regex==2018.8.29 numpy==1.17.3`

##### Activate the conda environment
`conda activate detectDesign-env-04_07_2021`

##### Add permissions to execute, read, and write.
`chmod 777 ./run_detectDesign.py`

##### Run the module 
`sbatch -p short -t 0-04:00:00 --mem 32000MB -o detectDesign_%j.out -e detectDesign_err_%j.err --wrap="./run_detectDesign.py"`