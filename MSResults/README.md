README.md
# killifish-aging-aggregates/MSResults

This MSResults folder contains the following files and scripts:

1. BGLPM folder contains all the Byonic(v2.6.49) output files (14 files) organized by tissue origin (brain, gut, heart, liver, muscle, skin, and testis) and sample types (tissue lysate (TL) or aggregate (AGG)) from our TMT-based multiplexing mass spectrometry study. 

2. BGLPMOut folder contains the result files obtained by running the TissueTMT.py file (make sure the utility file TMTUtl.py is included). The input files are in the BGLPM folder. The output kf_combined_prop.csv is the main result table that serves as input for all scripts in the MSAnalysis folder. 

3. longest-Protein-seqs_Nfu-ncbi100_andMito.fasta is the killifish protein fasta file used in our study.

4. NCBI annotations of killifish proteins and homologs information are tabulated in nfur_updated_symbols_20161013.csv. We used human homolog information to infer the function and cellular localization of killifish proteins.

All raw mass spectrometry reads as well as processed datasets can be found in the MassIVE database (https://massive.ucsd.edu/ProteoSAFe/static/massive.jsp) under ID MSV000086315.
