# Quantum Thermochemistry Calculator (QTC) 
[![Build Status](https://travis-ci.org/keceli/QTC.svg?branch=master)](https://travis-ci.org/keceli/QTC)

QTC includes modules that integrates open babel with quantum chemistry calculations and generates NASA polynomials in different formats.

It depends on:
  * Open Babel for cheminformatics (read/write chemical identifiers)
  * MOPAC, NWChem, Gaussian, Molpro for quantum chemistry calculations
  * MESS for calculating partititon function
  * PAC99, thermp for format conversions

## Installation
`conda` simplifies the installation process
If you don't have it you can install it without root privilages
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh
```
After the installation is completed you can install QTC with:
```
conda env create -f environment.yml
```
Activate `qtc-env` environment.
```
source activate qtc-env
```

## Run QTC
To run a series of quantum chemistry calculations for a water molecule (SMILES for H2O is O) with different packages
```
python src/qtc.py -i O -k "opt/mp2/dz/gaussian,freq/ccsd/dz/molpro,energy/ccsd(t)/adz/nwchem" -Q
```
To run thermochemistry calculations at the end, add `-T`
```
python src/qtc.py -i O -k 'opt/mp2/dz/gaussian,freq/ccsd/dz/molpro,energy/ccsd(t)/adz/nwchem' -Q -T
```
If you add -J, QTC generates a JSON file containing all the results
```
python src/qtc.py -i O -k 'opt/mp2/dz/gaussian,freq/ccsd/dz/molpro,energy/ccsd(t)/adz/nwchem' -Q -T -J
```
