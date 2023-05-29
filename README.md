README.md
# killifish-aging-aggregates

This repository contains data files and scripts used to analyze the quantitative TMT-based multiplexing mass spectrometry data and quantitative microscopy data for the manuscripts "Tissue-specific landscape of protein aggregation and quality control in an aging vertebrate” and "Identification of protein aggregates with prion-like and phase separation properties in the aging vertebrate brain". 

All the data files and scripts are organized by their origins into 3 folders: Microscopy, MSAnalysis, and MSResults. 

The Microscopy folder contains CellProfiler pipeline used to identify and quantify cells and fluorescent foci in S. cerevisiae, example quantification results, and scripts for visualization. 

The MSResults folder contains summarized data files from our TMT-based multiplexing mass spectrometry study where we profiled tissue lysate (TL) and aggregate (AGG) proteome from 7 tissues (brain, gut, heart, liver, muscle, skin, and testis) from young (3.5 months), old (7 months), and old telomerase deficient (7 months) male African killifish (N. furzeri). All raw mass spectrometry reads as well as processed datasets can be found in the MassIVE database (https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp) under ID MSV000086315. We also included the killifish protein fasta file and NCBI annotations including human homologs used in our analysis. 

The MSAnalyasis folder contains scripts used to analyze and visualize our mass spectrometry results as well as result files. All python scripts are compatible with Python 3.8.5 and Python 2.7.15. All r scripts are compatible with R 3.6.3.

Please refer to README file within each folder for list of scripts and files and their brief summary. Each script also includes a header that provide more description. 

[![DOI](https://zenodo.org/badge/328281013.svg)](https://zenodo.org/badge/latestdoi/328281013)
