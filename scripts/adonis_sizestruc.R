#--------------------------------------------------------------------------#
# analysis of community composition and visualisation of size structure ####
#--------------------------------------------------------------------------#

# This script analyses the community composition of the samples and visualises the size structure of each community
# Recreates Figure 1 of Padfield et al. (2018) Ecology Letters

# clear workspace
mise::mise(vars = TRUE, pkgs = TRUE, console = TRUE, figs = TRUE)

# load packages ####
library(phyloseq)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)

# set.seed
set.seed(42)

# load data ####
ps <- readRDS('processed_data/phyloseq_ready.rds')

# check sample_sums
sample_sums(ps)

# show available ranks in the dataset
rank_names(ps)

# process read data ####
# filter samples under 1000 reads
ps_sub <- subset_samples(ps, sample_sums(ps) > 1000)

# rarefy to an even depth
ps_rare <- rarefy_even_depth(ps_sub)

# filter for just Eukaryotes & Cyanobacteria
ps_aut <- subset_taxa(ps_rare, Kingdom == 'Eukaryota' | Phylum == 'Cyanobacteria')

# check number of samples
nsamples(ps_aut)

# Do PCoA on bray-curtis and plot ####

# bray-curtis ordination
ord1 <- ordinate(ps_aut, method = 'MDS', 'bray')
evals <- ord1$values$Eigenvalues

# plot ordination
plot_ord <- plot_ordination(ps_aut, ord1, color="ancestral", shape = 'treatment') +
  geom_point(size = 3) +
  scale_color_manual('Ancestral', values = c('black', 'red'), labels = c('Cold pond', 'Warm pond')) +
  scale_fill_manual('Ancestral', values = c('black', 'red'), labels = c('Cold pond', 'Warm pond')) +
  scale_shape_discrete('Treatment', labels = c('Ambient Incubator', 'Warm Incubator')) +
  stat_ellipse(aes(col = ancestral, group = ancestral), geom = 'path', type = "t") +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  ylab('PCoA2 [12.7%]') +
  xlab('PCoA1 [57.8%]') +
  theme(legend.position = 'none')

# Run PERMANOVA ####

# calculate distance matrix using bray-curtis dissimilarity
dist_mat <- phyloseq::distance(ps_aut, method = 'bray')

# make a data frame of the sample data
d_samp <- data.frame(sample_data(ps_aut))

# run an Adonis test
mod1 <- vegan::adonis(dist_mat ~ treatment*ancestral, permutations = 9999, data = d_samp)
mod2 <- vegan::adonis(dist_mat ~ treatment + ancestral, permutations = 9999, data = d_samp)

# plot size distribution ####

# read in count data
cd_counts <- readRDS('raw_data/counts.rds') %>%
  separate(., id, c('pond', 'ancestral', 'treatment'), sep = '_', remove = TRUE) %>%
  mutate(., temp = ifelse(treatment == 'A', 16, 20)) %>%
  unite(., id, c(pond, ancestral, treatment, temp), sep = '_', remove = FALSE)

# plot size distribution
plot_size_dist <- ggplot(cd_counts) +
  geom_line(aes(log10(carbon), group = id, col = ancestral), alpha = 0.1, stat = 'density') +
  geom_line(aes(log10(carbon), col = ancestral), stat = 'density') +
  xlab(expression(log[10]~Organism~mass~(µg~C))) +
  ylab(expression(Density)) +
  theme_bw(base_size = 10, base_family = 'Helvetica') +
  theme(strip.text = element_text(hjust = 0, size = 14),
        strip.background = element_blank()) +
  scale_color_manual(values = c('black', 'red'))

# make Figure 1 ####
figure_1 <- gridExtra::grid.arrange(plot_ord + ggtitle('(a) Bray-Curtis dissimilarity'),
                                    plot_size_dist + ggtitle('(b) size distribution') + theme(legend.position = 'none'),
                                    ncol = 2)
ggsave('plots/Figure_1.pdf', figure_1, width = 10, height = 5)
