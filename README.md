# PhasingPackage
Simple tool for phasing human PB-oocyte trios, as a programming assignment for the candidate informatician at Hoffmann Lab. It is the pilot stage of the phasing algorithm for phasing human oocytes from the SNP array data. This is the phasing by a common reference. The algorithm is described in Ottolini et al. 2016 and Ottolini et al. 2015. It reaches the step 45 of the Protocol described in Ottolini et al. 2016, previous to remove the common crossover. The program follows all the rules described in the Protocolo of Ottolini et al. 2016, until the step 45, except for the ones regarding the choosing of a reference.


# Example of usage

## Download and install the package
```
library(devtools)
install_github("ffalfred/PhasingPackage")
library(PhasingPackage)
```
## Usage
```
PhasingPackage::function_phasing(python_path,input,min,output)
```
### Options
* **python_path**: The path to the executable of python required to use the python code that uses the package. Python has to have installed the packages: math, os, numpy, pandas, matplotlib, operator, argparse.
* **input**: file containing the input data in format .txt tablimited.
* **min**: amount of SNPs required to define an haploblock.
* **output**: directory where to save the results. The directory has to exist previously.

### Example
```
PhasingPackage::function_phasing("/usr/bin/python3","~/Documents/meiomap/input_data_trios.txt",20,"~/Documents/meiomap")
```

## Author
Alfred Ferrer Florensa ; f.f.alfred@gmail.com 
