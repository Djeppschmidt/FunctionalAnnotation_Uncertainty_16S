# FunctionalAnnotation_Uncertainty_16S
This repository provides a customizable pipeline to access and summarize all whole genome annotations for prokaryotic genome assemblies hosted on GenBank and RefSeq. It also contains some functions for visualizing important characteristics of your dataset prior to annotation. It has been formatted as an R package, for easy installation. This pipeline does not include functionality to re-annotate the genomes, rather to summarize existing annotations. This approach requires far less memory and computational power. To install:

# Pipeline process:

A tutorial can be found here: https://djeppschmidt.github.io/ProkaryoteGeneAnnotation_tutorial

# Additional information

Coming soon: Tools for predicting functional gene abundance perform poorly in general. Chapter 2 of my dissertation covers why this is, and what we can do about it. The core analysis with figures will be presented here. I use this package as a custom annotation pipeline that allows me to "peek under the hood" to understand the confidence behind gene predictions; then I compare it to other gene prediction software that exists (PICRUSt, Tax4Fun2), and demonstrate clear shortcomings of all of these methods for predicting gene abundances for a number of genes that are important for biogeochemical cycling in soils.
