# Read Me
[![DOI](https://zenodo.org/badge/323755412.svg)](https://zenodo.org/badge/latestdoi/323755412)
DOI: 10.5281/zenodo.4536522

This depository contains the data and code for the "Developmental Circadian Disruption Alters Placental Signaling in Mice" paper. 

Important notes:
- Raw sequencing data (fastQ files) and gene count data are available in GEO, accession number GSE169266
- Sample # is labelled with a "10.." prefix for the ID number. For example, sample #4 = ID number 1004.
- Sex of all placental samples was measured by qPCR of fetal tail sample; however, one sample in the RNAseq analysis was found to have mislabelled sex. Sample #9 ("1009") was initially labelled as "Male" but changed to "Female". Likewise, the sex ratio for the dam was changed from 0.6 to 0.7 to reflect the change in sex. 
- Sample #12 ("1012") was found to be an outlier driving all the DESeq2 results and was not included in the RNAseq analysis.

# Data and code

STAR genome alignment and gene count code is available in the "STAR_code.txt" file. 

Figure 1: weight, embryo, and placental outcomes data 
- The "Dam_weights.xlsx" file and "Data_embryo_placenta.xlsx" file contains the data analyzed and graphed in Figure 1. 

Figure 2: placental RNAseq data
- The gene count data is available in the "placenta_genecount.txt" file and corresponding sample informations are available in the "Pheno_data_rnaseq.csv" and "DESeq_mouse.RData" files. The analysis code and code used to produce the graphs is available in the "Github_analysiscode.R" file. The differential expression results for the effect of sex are located in the "Sex_Female_vs_Male.csv" file and the differential expression results for the effect of light treatment are located in the "Treatment_CD_vs_CL.csv" file. 

Figure 3: placental immunofluorescence data
- The raw and averaged immunofluorescence data used in the analysis and graphed in Figure 3 is available in the "IHC_placenta_data.xlsx" file.
