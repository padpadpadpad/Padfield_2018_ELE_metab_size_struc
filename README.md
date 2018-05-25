# Analyses of metabolism, temperature and size distributions as published in:

_Padfield et al. (2018) Linking phytoplankton community metabolism to the individual size distribution. Ecology Letters_ [doi: 10.1111/ele.13082](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.13082)

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.1252775.svg)](https://zenodo.org/record/1252775#.WwfLBi_MxBw)

### Outline

This repository contains the final datasets, analyses and figures of the above-mentioned paper. It can recreate the all of the analyses and figures in the main text (Figure 1 to Figure 4).

In addition, the repository contains all of the raw fastq files from the sequencing used in the paper (see below for a description of the pipeline used in order to recreate the phyloseq ready file that is present in `processed_data/`).

Finally, the repository contains R scripts and `.stan` files that allow for Bayesian fitting of the Metabolic Scaling Theory framework presented in the paper.  

### Feedback

- Please report any problems or bugs in the code in the [Issues](https://github.com/padpadpadpad/Iceland_stream_ELE_analyses/issues) tab of the GitHub repository. Alternatively, please email _d.padfield@exeter.ac.uk_.

### Licensing

This code is licensed under GPL-3.

### Running the scripts and analyses

- The project can be `cloned` or for those not familiar with GitHub, a zip file of this project can be downloaded using the "Clone or download" button at the top right of this page.
- Open the R project file in the downloaded folder. [R projects](https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects) automatically assigns the root directory to the directory in which the project resides. Consequently all of the analyses should be runnable without altering paths. These are very easy to open using RStudio. All of the scripts for the analyses can be found in `scripts/`.
- The script `install_packages.R` will install all the packages necessary to run the analyses, split up by each subsequent script
- `adonis_sizestruc.R` contains the analysis of the community composition and plots the size distribution of the communities.
- `MST_mle_model.R` runs the model fitting of the MST framework to the metabolism, size and temperature data of the communities using maximum likelihood.
- `mle_MST_analysis.R` contains the code to analyse the relationsip between temperature-corrected rate and size-corrected biomass for the best fitted model from `MST_mle_model.R`. It also contains the analysis of the raw metabolism values across treatments.
- `MST_stan_model_GP.R` and `MST_stan_model_R.R` contain the R code to run Bayesian model fitting of the MST framework using Stan. The stan models used can be found in `stan_models/`
- `stan_MST_analysis.R` gives the code to analyse the relationsip between temperature-corrected rate and size-corrected biomass for the best fitted model from `MST_stan_model_GP.R` and `MST_stan_model_R.R`.
- All of the data needed to run the analyses are stored in `processed_data/` or `raw_data/`.
- All figures produced by the analysis are saved in `plots/` and are labelled as they are in the main text.

__All analyses are done in R version 3.4.4, on macOS High Sierra 10.13.4. I am unsure whether some older version of R will support all of the packages and whether the analyses will run exactly the same.__

### Accessing the raw sequencing files and running the bioinformatics pipeline

- All of the raw `fastq.gz` files are stored in `sequencing_data/raw_files`. 
- The reference database used for the assignment of taxonomy is present in `sequencing_data/ref_fasta`.
- Sample meta data is stored in `sequencing_data/sequence_file_metadata.csv`.
- To get from raw sequencing data to an object ready for analysis, I used the full-stack [dada2/phyloseq workflow](https://f1000research.com/articles/5-1492/v2) in R. Specifically, I used the `raw_read_processing.R` script within a self-developed pipeline available on [GitHub].(https://github.com/padpadpadpad/AB_dada2_pipeline_R). Reads were trimmed between 25bp and 250bp because of low quality scores outside of these boundaries. With this information and this repository, the processing of the raw files should be possible if desired.
- After the pipeline was ran, we removed any amplicon sequence variant (ASV) that had not been assigned to the Phylum level.
- Please let me know if you have any queries about the sequencing pipeline and I will do my best to help.

### Some information on the Stan models

- Although the Stan models and Bayesian approach does not appear in the paper, I feel like they are a welcome addition to the GitHub reporitory of the manuscript. I am relatively new to Bayesian statistics and writing the models used here was challenging. Consequently, I cannot guarantee that I have adhered to best practice or that the code and models are completely correct. However, they do reassuringly give idential mean estimates to the maximum likelihood approach. And any deviations from best practice or errors were not through a lack of trying. It is this part of the analysis that I would very interested in pursuing with other datasets. I think there is a huge scope to extend the models so please get in touch if you have any queries, ideas or datasets you'd like to try the models on. If I have time I might try write a little more about this approach in a blog post in the near future.
- The models in `stan_models/` fit the equivalent models and give equivalent estimates to the maximum likelihood approach. However, in `stan`, the `generated quantities{}` block allows for the calculation of things from the model fitting process. Using this block, it is possible to calculate mass-corrected biomass and temperature-corrected rate and individual photosynthetic rate from within the model. This gives the benefit of estimating the uncertainty around these predictions, that can be used in new predictions of mass-corrected biomass and temperature-corrected rate.
- I found the model comparisons somewhat problematic. I followed the advice of Andrew Gelman as best as I could, using LOO (Leave One Out) validation and WAIC scores to try and do [model evaluation](http://www.stat.columbia.edu/~gelman/research/unpublished/loo_stan.pdf), but it felt like I was trying to combine frequentist and Bayesian approaches, clinging onto the idea of starting with a complex model and ending with the most parsimonious one.
