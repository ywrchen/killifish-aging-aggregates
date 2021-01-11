README.md
# killifish-aging-aggregates/MSAnalysis

This MSAnalysis folder contains the following files and scripts grouped by the main objective.

QCs-Descriptive-Plots

1. Reproducibility.py
--Description: Assess reproducibility by plotting one sample against another.
--Input files: kf_combined_prop.csv

2. PCA.py
--Description: PCA analysis
--Input files: kf_combined_prop.csv

3. Overlap7Venn.R
--Description: Make Venn diagram to visualize shared and tissue-specific proteins and those with age-associated increase for tissue lysate (TL), aggregate (AGG), and aggregation propensity (Prop).
--Input files: kf_combined_prop.csv


Differential-Expression-OldvsYoung-TertvsOld

4. SigHits.py
--Description: Create table to filter the proteins based on different criteria (z-score, logFC, p value etc). These tables are input for correlomics, biophysical analysis, and hit-specific heatmap. Section 4, 5, and 6 are enrichment analysis for proteins with putative prion-domains (predicted from PLAAC, http://plaac.wi.mit.edu/), intrinsically disordered regions (predicted from DISOPRED2, http://bioinf.cs.ucl.ac.uk/psipred/?disopred=1), and neurodegenerative disease association respectively. The list of proteins with neurodegenerative disease association (20201115_CuratedND.csv) is manually curated from databases including ALzGene (http://www.alzgene.org/default.asp), FTDGene (http://www.molgen.vib-ua.be/FTDmutations),  and Global Variome shared LOVD (https://databases.lovd.nl/shared/genes?search_diseases_=PARK#id=0&order=id_%2CASC&search_diseases_=PARK&page_size=1000&page=1). Neurovegetative diseases that are considered include AD (Alzheimer’s disease), PD (Parkinson’s disease), ALS (Amyotrophic lateral sclerosis, Lou Gehrig’s disease), SCA (Spinocerebellar ataxia), and SMA (Spinal muscular atrophy).
--Input files: kf_combined_prop.csv and BrainFig/20201115_CuratedND.csv

5. ProteinBarChat.py
--Description: Box plot on protein abundance and/or aggregation propensity for different age/disease groups. Section 2, part 4,5,6, and 7 are the relevant code.
--Input files: kf_combined_prop.csv

6. ManhattanPlot.py
--Description: Manhattan plot on age-associated changes. Run section 4.1 part 1. To generate the plot based on z-score (z-score of log2_FC on the y-axis). This is used for the TERT section. Run section 7 to generate plot by sample type (TL, AGG and Prop). 
--Input files: kf_combined_prop.csv

7. TertAnalysis.py
--Description:  Create a heatmap of fraction of terms that met specific z-cutoff across different tissues and sample types. This is used to correlate with proliferation index of different tissues.
--Input files: kf_combined_prop.csv

Tissue-Specificity

8. SharedVsTissueSpecific.py
--Description: Section 1: Diagonal heatmap that show overlap among proteins identified in AGG or WCL, or proteins with significant age-associated changes in each category. (“Pyramid plot”). Section 2: Sort the significant TL/AGG/Prop hits (either increase or decrease in Old compared to Young) by shared vs. tissue specific. The output table is used as input for TLAGGPropHeatmap.py to generate shared vs. tissue-specific heatmap. Section 3: run this to analyze the expression of a protein across multiple tissues.
--Input files: kf_combined_prop.csv

9. TLAGGPropHeatmap.py
--Description: Based on the sorted significant hits (sort by tissue-specific and shared entries), overlay their logFC to generate heatmap.
--Input files: output csv files from SharedVsTissueSpecific.py

Cross-Datasets-Comparisons

10. ProteomeComp.py
--Description: Compare our killifish dataset with two papers that profiled aggregates in old and young worms: David et al., 2010 and Walther et al., 2015 as well as the proteomic study that examed aged mice and killifish brain aggregates by Kelmer Sacramento et al., 2020. Read the script for more details on specifically which result tables were used for different comparison. 
--Input files: in ProteomeComp folder where each sub-folder contain relevant comparison with a specific paper (David2010_kenyon, Ori2020_killifish, and Walther2015_worm). kf_combined_prop.csv

Gene-Set-Enrichment-Analysis

11. getGeneLists.py
--Description: compute rankstats of all protein detected in each tissue sample. Rank-stats is defined as -log10(p-value)*log2FC on TL/AGG/Prop for each age pair comparison. The result is used as input for the GSEA analysis. 
--Input files: kf_combined_prop.csv

12. GSEA_v1.R
--Description: script to do enrichment analysis using GSEA (GO, MSigDb, KEGG, and DO genet sets based on human homlogs). It generates all combinations of input files and does enrichment and sorts results in the
--Input files: all input data files are in GSEA/GeneSets and are output from getGeneLists.py. The collection of MSigDB input are in GSEA/MSigDb_Collections_20190318 downloaded from their database on March 18, 2019.

13. EnrichmentCleanup.py
--Description: Combine the GSEA results. Part 1 is the main section I used to generate the final figures to preserve as many significant terms as possible. Part 2 is written to filter out gene set if it is part of another significant gene set (complete subset of another gene set or exact duplicate are effectively removed, in case of complete duplicate the term with the most significant enrichment is preserved). Note that part 2 takes a while to run especially on MSigDB (budget ~30 min per age comparison).
--Input files: GSEA output files in GSEA/Results section.

14. BidirectionalHeatmap.py
--Description: Blanket/check plot on enrichment analysis. The tissue-specific enrichment terms are placed on the top followed by shared terms (terms shared among the most tissues are placed the closet to the bottom). The rank-stats defined as -log10(p-value)*log2FC on TL/AGG/Prop for each age pair comparison (computed in getGeneLists.py) for each terms were visualized.  
--Input files: output files from EnrichmentCleanup.py

Killifish-GO-Parsers

15. GoOboParser.py
--Description: This script was used to parse uniport and Reactome cellular compartment information.
--Input files: CellularLocalization/Uniprot/….

16. CellularCompartment_Final.py
--Description: Cellular compartment enrichment analysis. Count number of proteins that belong to certain compartments. The proteins can either be detected proteins in TL or AGG (section 3), or AGG (or Prop) hits. For the hits, z-score cutoff consistent with the rest of the hits analysis is used (either none or 0.75, section 2).
--Input files: CellularLocalization/Uniprot/ killifish_human_CM.csv and kf_combined_prop.csv 

Biophysical-Properties

Properties folders contain input files from other proteome-wide computational analysis (PLACC, DISOPRED, localCIDER, and evolutionary analysis).

17. Cider.py
--Description: Script used to obtain all the localCIDER output. Note that this is done on the entire killifish proteome. Use section 2 to run localCider on a few sequences. localCIDER is developed by the Pappu lab: https://github.com/Pappulab/localCIDER
--Input files: Killifish_DISOPRED_PLAAC.csv

18. Correlomics.py
--Description: 
----Section 1: merge the sequence feature results from various analysis into one table.
----Section 2: use this to generate histogram on various parameters without or with weights (use individual part of the code as needed).
----Section 3: compute the correlation between evolutionary rates, expression data from this study, protein disorder metrics etc. 
----Section 4: use this section to run simulation in which proteins are randomly sampled (equal weights or weights by their counts) and their physical properties were retrieved to approximate a true population ensemble. The significance (p-value) of how the mean of a query property differs between hit and population is computed based on one-tail test. The null is to assume the hit and population are similar. The output here is used to generate heatmap. Weights based on young protein expression counts were also considered when comparing AGG hits. Budget ~30 min – 1hr for this section.
----Section 5: test if the distribution (equal weights or weights by their counts) of metric of interest for a given tissue sample (TL/AGG/Prop) follows normal distribution.
----Section 6: similar to section 4, but instead the “hits” here are the highly abundant proteins in young samples. This is to infer features of aggregation-prone proteins in young killifish. Budget ~30 min – 1hr for this section.
----Section 7: test if certain biophysical features are predictive of query parameter (either aggregation propensity or simply aggregate abundance) for a particular age group. 
--Input files: kf_combined_prop.csv or alternatively Properties/ kf_combined_prop_metrics.csv

19. CorrelomicsHeatmap.py
--Description: heatmap visualization of biophysical analysis in Correlomics_Final.py. 
--Input files: Make sure to run Correlomics.py to get the desired input files before running this.

20. AA_Heatmap.py
--Description: Visualization of protein sequence properties (disorder, putative prion region, hydropathy, charged regions, and aromatic residues).
--Input files: Various inputs are in Properties folder. The important one is Properties/ kf_combined_prop_metrics.csv

21. delta.py
--Description: Protein charge patterning calculator. The exact rationale is available in https://doi.org/10.1073/pnas.1304749110. The color coding is a simple way to visualize the positive and negative residues. 
--Input files: The input example peptide sequences were in the script.

22. SysSeedingCorrelomics.py
--Description: 
----Section 1: append various biophysical properties to proteins tested in the yeast single color overexpression assay (test for autonomous aggregating proteins). 
----Section 2: output the protein sequences on protein of interests.
----Section 3: visualization of seeding results as scatterplot (use aggregation status as hue).
----Section 4: visualization of the quantified yeast single color overexpression assay (fraction of cells with puncta during induction is generated).
--Input files: SysSeedingAll_dotplot.csv (this contains the quantified fraction of foci in one color microscopy images). SysSeedingAll.csv is the file used for generated example figures that highlight charge distribution. Note these two input files are almost identical except that dotplot only include those with age-associated increase in AGG or Prop whereas SysSeedingAll.csv contain a few controls a well. Run section 1 to initialize the result table. 

23. SysSeedingSVM.py
--Description: Support vector machine classifier to classify aggregate and diffuse proteins based on yeast in vivo aggregation status. 
----Section 1: build 2 parameter model (go through all possible 2 parameter combinations) and define hyperplane (support vector machine). 
----Section 2: visualize some of the best SVM classifiers by specifying x and y parameters (scaled so the x and y axis are standardized).
----Section 3: visualize some of the best SVM classifiers by specifying x and y parameters (unscaled so the x and y axis are in the original untransformed coordinates/range).
----Section 4: visualize some of the best SVM classifiers by specifying x and y parameters (scaled so the x and y axis are standardized), define the hyperplane explicitly.
----Section 5: fine-tuning a model to find the best parameter that yield the best hyperplane with the highest accuracy.
----Section 6: overlay age-associated aggregates on the best linear SVM classifier using the two predictive parameters.
----Section 7: visualize the performance of all two-parameter SVM classifiers generated in section 1. The accuracy for each model is visualized as heatmap.
----Section 8: compute pairwise correlation coefficients on all metrics tested.
--Input files: SysSeeding/SysSeedingAll_metrics.csv

Protein-Quality-Control-Sets

24. GeneSetFilter.py
--Description: Python script used to extract the MS data on gene/protein set of interest. The output were used in downstream visualization in bubble plot and protein complex cartoon (the fold change and p-value as used to extrapolate color for these plots). 
----Section 1: run this on a geneset_category(i.e. chaperones, ribosomes, proteasomes) of choice on  the sample type of choice (i.e. AGG, TL, and Prop). 
----Section 2: run this to combine all genesets with killifish identifiers. The output is used for Cytoscape visualization.
--Input files: kf_combined_prop.csv and files in GeneSets/Killifish/. All these tables were curated. Proteasome and ribosomes were listed. Proteins involved in degradation pathways (ie. ubiquitin proteasomal pathway and lysosomal/autophagy pathway) were downloaded from KEGG. Chaperones list were carefully curated from multiple sources.
--Output files: Various table organized by sample type (i.e. AGG, TL, and Prop).

25. BubblePlot.py
--Description: Bubble plot on age-associated changes with a curated list of genes. The size of each circle is reflective of their p-value and color is indicative of log2FC. Run different sections on different genesets.
--Input files: output files from GeneSetFilter.py in the format of GeneSets/Killifish/Chaperones_AGG/…flat.csv and GeneSets/Killifish/Proteasome_AGG/…flat.csv

26. ColorExtrapolation.py
--Description: extrapolate color (in hex value) from a linear colormap based on log2FC or rankstats (-np.log(pval)*log2FC) values. The output values are used to fill the cartoon for protein complexes. The colorscale here is consistent with that in BubblePlot.
--Input files: various geneset results that are output from GeneSetFilter.py. 


OMIM-Analysis

27. OMIM.py 
--Description: Identify age-associated aggregates that were also associated with Mendelian diseases (obtained from OMIM). Those that were significant across multiple tissues were visualized using ProteinBarChart.py
--Input files: OMIM/genemap2_melt.csv is re-formatted from genemap2.txt downloaded from the OMIM database (retrieved on April 10, 2019). 

Evolutionary Analysis

28. parse_values.pl 
--Description: script from Param Singh for evolutionary rate analysis. This pipeline is based on codes and datasets that have been previously verified for Valenzano et al 2015 and Wagner et al. 2019. That involves: 
step1: get 1-to-1 ortholog clusters using 6 fish genomes. Identity default is 25%, default coverage is 50%.
step2: get coding sequences for these clusters
step3: align codons using PRANK
step4: convert alignment to PAML format
step5: run YN00 to get dN, DS
step6: Parse YN00 output to get the data
The output of YN00 anaysis is Guidance_YN00_Out.txt. The final output Ka_Ks_values_Nfur-to-otherfish.txt is generated using parse_values.pl and Guidance_YN00_Out.txt as input.
--Input files: dN-dS_evolutionary_rate/Guidance_YN00_Out.txt
--Output files: dN-dS_evolutionary_rate/Ka_Ks_values_Nfur-to-otherfish.txt

29. get_other_ids.pl
--Description: perl script from Param Singh to get the protein id (XP) and gene id (XM). The final output file Ka_Ks_values_Nfur-to-otherfish_with-XM-XP-Ids.txt is used in Correlomics.py
--Input files: dN-dS_evolutionary_rate/nfurzeri_LongestProteins_ncbi.txt and dN-dS_evolutionary_rate/Ka_Ks_values_Nfur-to-otherfish.txt
