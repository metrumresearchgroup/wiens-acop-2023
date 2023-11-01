################################################################################
## Python Setup Script
## 
## Only run this the first time that you set up this workflow. This will install
## all the python and TensorFlow/Keras components that you need to run the rest
## of the scripts in this repo.
################################################################################
library(reticulate)

## First step: install the latest version of python
## this call also returns the path that you will use to for reticulate
path_to_python <- install_python(version='3.7.7')
## on my machine the resulting path is: 
# /data/home/samc/.pyenv/versions/3.7.7/bin/python3.7

## Next, we create the virtual environment of python used to run Keras
virtualenv_create("r-reticulate", python = path_to_python)

## this seems to install/update the virtual environment
## instead of any other versions of python on the disk (of which there are many)
use_virtualenv("r-reticulate")

## Install both the core Keras library as well as the TensorFlow backend to this
## environment. This allows you to use the install_keras() function
tensorflow::install_tensorflow(envname = "r-reticulate",
                               version="2.11",
                               extra_packages = "tensorflow_probability")
keras::install_keras(envname = "r-reticulate")

## use the following command line call to upgrade to
## the latest version of tensorflow_probability
system("pip install --upgrade tensorflow")
system("pip install --upgrade tensorflow_probability")

## Optional for plotting
system("pip install pydot")
#system("sudo apt install graphviz")