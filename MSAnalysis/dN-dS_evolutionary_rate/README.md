### Protein families for dN, dS analysis

> This pipeline is based on codes and daatsets that have been previously verified for Valenzano et al 2015 and Wagner et al. 2019. 

#### Steps:

1. get 1-to-1 ortholog clusters using 6 fish genomes. Identity default is 25%, default coverage is 50%.
2. get coding sequences for these clusters
3. align codons using PRANK
4. convert alignment to PAML format
5. run YN00 to get dN, DS
6. Parse YN00 output to get the data

The output of YN00 anaysis is Guidance_YN00_Out.txt
Use parse_values.pl and Guidance_YN00_Out.txt to parse dNdS values and generate the final output.
