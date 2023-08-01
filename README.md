# ferret_ionocyte
Code repository for interspecies comparison between mouse, ferret, and human
This code repository contains R code for the analysis of 4180 pulmonary rare cells (Ionocyte, neuroendocrine, and tuft) from mouse, human, and ferret. The easiest way to run it is to clone the repo (download ZIP button at top) and open the file 'interspecies_integration_v2.Rmd' in RStudio. Prior to running the code, you should obtain the rare cell feature table for mouse, human, and ferret listed below either through GEO or Human Lung Cell Atlas. Each code chunk can then be run by clicking the run button (top right corner of each chunk).

 ## Human data
 https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293 <br />
 download the full data and then run rare_cells.py to obtain the feature table

 ## Ferret data
 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233654 <br />
 GEO Accession: GSE233654

 ## Mouse data
 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103354 <br />
 GEO Accession: GSE103354, download the RDS file from the supplementary file
