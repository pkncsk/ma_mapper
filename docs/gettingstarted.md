# Getting started
## Quick links
- [About MA Mapper](#about-ma-mapper)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Prerequiste](#prerequiste)
## About MA Mapper

A python package for mapping genome-wide data onto multiple alignment of transposable elements. The main function of this package is to work as a wrapper for other biological data file handlers such as biopython to steamline multiple alignment data overlay. The output of this package is a numerical numpy matrix, which can be easily used for visualization by matplotlib or seaborn. 

## Dependencies

### general system
- python ^3.10
- compress-pickle ^2.1.0
- h5py ^3.7.0
- numpy ^1.21.5
- pandas ^1.3.5
- scipy ^1.7.3

### biological data file handler
- biopython ^1.83
- pybedtools ^0.10.0
- pybigwig ^0.3.23
- pysam ^0.22.0
- cyvcf2 ^0.30.28

### visualization 
- logomaker ^0.8
- matplotlib ^3.8.4

## Installation

There is a plan for a proper release but as of now, this package can be installed with pip:
```bash
pip install git+https://github.com/pkncsk/ma_mapper
```
## Prerequiste
The general workflow of MA Mapper 