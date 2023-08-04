# KidneyGPS
RShiny web application for navigating genome-wide association study (GWAS) results and their functional and regulatory annotation to guide gene prioritization.

### Description
The KidneyGPS is available at https://kidneygps.ur.de/gps/. It summarizes information on association and annotation of variants and genes mapping to 424 genetic loci identified for glomerular filtration rate based on a GWAS meta-analysis of UK Biobank and CKDGen consortium data (n=1,201,909) [Stanzick et al. Nat. Commun. 2021]. Variants and genes were anotated using various data sources that are described in detail on the website. 

### Dependencies
KidneyGPS depends on the following other R packages:  data.table, shiny, shinyBS, DT, shinyjs, shinyWidgets, htmltools

### Release History
* Version 2.3 (current version)
  * Integration of eGFR association statistics separated by diabetes status for all credible set variants.
* Version 2.2
  * Re-organisation of the GPS tab with extendend filter options.
* Version 2.1
  * separated drug information by indication for kidney diseases according to the ICD-11 codes.
* Version 2.0
  * reanalysis of the 424 loci regarding stability of independent signals. 594 stable independent signals were identified using stepwise conditioning with GCTA. Credible sets were re-calculated for these signals resulting in 35,885 credible set variants.
* Version 1.3.1
  * additional summary numbers in "GPS tab" below GPS table; in minor fixes 
* Version 1.3 
  * Integrated drugability data from Therapeutic Target Database
* Version 1.2
  * additional evidence for human phenotypes, GWAS evidence from Gorski et al. and Winkler et al.
* Version 1.1
  * updated eQTL data
* Verion 1.0 
  * based on results published in Stanzick et al. Nat. Commun. 2021

### Citation
If you are using the app, please cite [Stanzick et al. Nat Commun 2021](https://pubmed.ncbi.nlm.nih.gov/34272381/).   

### Help
Please contact kira.stanzick@ukr.de if you need assistance. 
