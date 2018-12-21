# CheckM

[![version status](https://img.shields.io/pypi/v/checkm-genome.svg)](https://pypi.python.org/pypi/checkm-genome)

### Quick start with `conda` (recommended)
The easiest way to install CheckM is with Conda.  


```
# Create a new environment with Python 2.7  
conda create -n checkm python=2.7 -y

# Install CheckM and dependencies
conda install -c bioconda pplacer hmmer prodigal checkm-genome

# Test the install
checkm -h   # Confirm CheckM is installed.
hmmalign -h # Confirm HMMER is installed.
prodigal -h # Confirm Prodigal is installed.
pplacer -h  # Confirm pplacer is installed.
```


Please see the project home page for usage details and installation instructions:
https://github.com/Ecogenomics/CheckM/wiki

Note: we do not recommend installing CheckM from the master branch. This may be unstable. Please install an official release of CheckM or use `pip` or `conda`.

Copyright © 2014 Donovan Parks, Connor Skennerton, Michael Imelfort. See LICENSE for further details.
