# AERIS_MAP_BEAT_analysis

This folder contains perl and R scripts used for major data analyses in the manuscript "Inflammatory endotype associated airway microbiome in COPD clinical stability and exacerbations – a multi-cohort longitudinal analysis" by Wang et al. The description and usage of the scripts are as below:

1_combat_analysis.R

This script performs Combat analysis for microbiome data, using the OTU table (i.e. L6.txt) as input.

2_comm_type_analysis.R

This script performs community type analysis using two different methods, 1) hclust based on JSD and clValid calculation of Silhouette and 2) PAM clustering based on JSD and calculation of CH index and Silhouette, using the batch-effect adjusted OTU table (i.e. L6_batchrm.txt) as input.

3_check_contam.pl

This script performs taxonomic check for potential contamination based on list from Salter et al. BMC Biology paper, using the OTU table (i.e. L6.txt) as input. The script calls upon a text file 'negative_salter_paper.txt' which should be placed in the same folder with the script.

3_negative_salter_paper.txt

This text file contains list of potential contamination taxa from Salter et al. BMC Biology paper.

4_Haemophilus_rescale.pl

This script rescales the microbiome relative abundance data by downscaling Haemophilus to its average abundance, using the OTU table (i.e. L6.txt) as input.

5_changepoint_calc.pl

This script calculates change-point for microbiome taxa for longitudinal microbiome datasets using PELT algorithm, using the OTU table (i.e. L6.txt) and metadata text file as inputs. The script utilizes an embedded R script (5_PELT.R) which should be placed in the same folder with this script.

5_PELT.R

The R script embedded in 5_changepoint_calc.pl.

6_calculate_OR.pl

This script calculates odds-ratio for the association between microbiome change-points and clinical events of interest (here being the endotypic switches), using the output file from 5_changepoint_calc.pl. The script utilizes an embedded R script (6_calcor.R) which should be placed in the same folder with this script.

6_calcor.R

The R script embedded in 6_calculate_OR.pl.

7_cross_cov_calc.pl

This script calculates cross-covariance between microbiome taxa and clinical measures of interest (here being eosinophilic/neutrophilic counts) for longitudinal samples of each patient, using the OTU table (i.e. L6.txt) and metadata text file as inputs. The script utilizes an embedded R script (7_ccf.R) which should be placed in the same folder with this script.

7_ccf.R

The R script embedded in 7_cross_cov_calc.pl.

