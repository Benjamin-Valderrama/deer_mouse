#' TO-DO LIST:
#' 
#' 1) Literature that can be used in the discussion:
#' 
#' 1.1) Kynurenine pathway metabolism
#' https://moodle2.units.it/pluginfile.php/246819/mod_resource/content/1/Materiale%20di%20studio%203.pdf


#Housekeeping
options(stringsAsFactors = F)

#Load libraries
library(tidyverse)
library(Tjazi)
library(RColorBrewer)
filling_colors <- brewer.pal(8, "Paired")[c(1:2,7:8)]

#Define working directory
setwd(paste0(getwd(), "/deer_mouse"))

set.seed(12345)
#load data
counts  = read.delim("data/genus_table_from_dada2.csv", sep = ",", row.names = 1)
metadata = read.delim("data/metadata.csv", sep = ",")

#Remove plasma for now as it contains missing values. 
counts = counts[,metadata$Sample_ID]


#' TO DO LIST:
#'
#' 0. Stacked barplot with the composition of the microbiome [DONE]
#' 
#' 1. ALPHA DIV [DONE]
#'     shannon [DONE]
#'     simpson [DONE]
#'     chao1 [DONE]
#'     Test to see diferences [DONE]
#'  
#' 2. BETA DIV  [DONE]
#'     Permanova [DONE]
#'     Plot [DONE]
#'     
#' 3. DIF ABUND [DONE]
#'     GENUS fw_glm(F = ~ Treatment * Nestbuilding + plate) 2[DONE]
#'     GBM [DONE]
#'     GMM [DONE]
#'     
#' 4. Histograms of DEA for each condition by its own (i.e., Treatment and Nestbuilding) on : [DONE]
#'     GENUS [DONE]
#'     GBM [DONE]
#'     GMM [DONE]
#'     
#' 5. Are the strains with Diff abundance in the 1st study also present in the second study? [DONE]
#' 
#' 6. Figures for the paper:
#' 
#' 1) alpha-div, beta-div and stacked barplot [DONE]
#' 2) GBM and GMM -> Need to incorporate the effect size on the right side of the GBM plot
#' 3) Plasma vs Tryptophan



# Presence of 1st study strains -------------------------------------------

#' In the first study, Prevotella & Anaeroplasma were found to be increased in Control, while
#' Desulfovarmiculus, Peptococcus & Holdemanella were found to be increased in OC-like animals
#' 
#' From them, only Anaeroplasma, Peptococcus and Desulfovarmiculus- & Prevotella-related groups were found to be present.
#' However, differences were not statistically significant and by plotting is difficult to see tendencies. 
#' Only the Desulfovibrionaceae (the Desulfovarmiculus-related group) shows a more clear tendency, however, this is inverse
#' to what was reported in the first study.

counts %>% 
  rownames_to_column(var = "Strain") %>% 
  mutate(Strain = word(Strain, 5, 6, sep = "_"),
         Strain = str_replace(Strain, "_", " ")) %>% 
  pivot_longer(cols = starts_with("X"),
               names_to = "Sample_ID",
               values_to = "Counts") %>% 
  left_join(., metadata, by = "Sample_ID") %>% 
  
  # Look at the strains that showed differential abundances in the previous study
  filter(!is.na(Nestbuilding) & str_detect(Strain, "Prevotella|Anaeroplasma|Desulfo|Aestuari|Peptococcus|Holdemanella") & Counts > 0) %>% 
  mutate(Nestbuilding = factor(Nestbuilding, levels = c("NNB", "LNB"))) %>% 
  ggplot(aes(x = Nestbuilding, y = Counts, fill = Nestbuilding)) +
  geom_boxplot() + 
  geom_point() + 
  facet_wrap(.~Strain, scales = "free") + 
  scale_fill_brewer(type = "qual", palette = "Paired", name = "Groups") + 
  theme_minimal() + 
  labs(x = "", y = "Counts") +
  theme(strip.text = element_text(size = 16, face = "bold"),
        
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14),
        
        legend.key.size = unit(0.8, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.6, 0.15))





# Stacked barplot ---------------------------------------------------------
library(RColorBrewer)
library(ggnewscale)

# Transform counts into relative abundances
rel_counts <- apply(counts, 2, function(x){x/sum(x)})

#' Determine the most abundant genera. 
#' 
#' This will be useful if I want to work with the most abundant genera
#' and collapse the other genera into "Others".
#' 
#' FOR THE MOMENT, THIS IS NOT USED IN THE ANALYSIS.
most_abundant_genera <- rel_counts %>% 
        t() %>% 
        as_tibble() %>% 
        pivot_longer(cols = everything(),
                     names_to = "taxa",
                     values_to = "abundance") %>% 
        mutate(taxa = word(taxa, start = 5, end = -1, sep = "_")) %>% 
        group_by(taxa) %>% 
        summarise(mean_abundance = mean(abundance)) %>% 
        filter(!str_detect(taxa, "Unknown")) %>% 
        slice_max(order_by = mean_abundance, n = 14) %>% 
        pull(taxa)

# Colors for the Phylums that appeared in the analysis
my_cols_phylum <- c(brewer.pal(7, "Dark2"),
                    brewer.pal(7, "Set1")[c(2,7,5)],
                    "#7f7f7f")

# Colors to identify the treatments the samples underwent 
my_cols_treatment <- c(brewer.pal(8, "Paired")[c(7,1)])

#' Names of the samples in each group: either escitalopram or water.
#' This will be useful to reorder the samples in the plot with the bars
samples_nnb <- metadata %>% 
        filter(Nestbuilding == "NNB") %>% 
        pull(Sample_ID)

samples_lnb <- metadata %>% 
        filter(Nestbuilding == "LNB") %>% 
        pull(Sample_ID)


#' Stacked barplot with phylum
#' 2: Phylum
#' 3: Class
#' 4: Order
#' 5: Family
#' 6: Genus

# This is to add the column with the Samples_ID to the dataframe
Sample_ID <- rownames(t(rel_counts))


#' Reorder the phylums for the final stacked barplot
phylums <- rel_counts %>%
        t() %>% 
        as_tibble() %>% 
        # Here I added the names of the samples as column
        mutate(Sample_ID = Sample_ID) %>% 
        pivot_longer(cols = !matches("Sample_ID"),
                     names_to = "taxa",
                     values_to = "abundance") %>% 
        
        # Change start to change the taxonomy level used in the plot
        mutate(Phylum = word(taxa, start = 2, end = -1, sep = "_"),
               Phylum = str_replace_all(Phylum, "_.*", "")) %>% 
        pull(Phylum) %>% 
        unique()

phylums <- c("Firmicutes", phylums[phylums != "Firmicutes"])


cols <- c(brewer.pal(12, "Set3")[c(2:4)],
          brewer.pal(12, "Set3")[c(7,8)],
          brewer.pal(8, "Pastel2")[7],
          brewer.pal(12, "Set3")[c(1, 10:12)],
          "#b3b3b3")

#cols[phylums == "Unknown"] = "grey50"


plot_composition <- rel_counts %>%
        t() %>% 
        as_tibble() %>% 
        # Here I added the names of the samples as column
        mutate(Sample_ID = Sample_ID) %>% 
        pivot_longer(cols = !matches("Sample_ID"),
                     names_to = "taxa",
                     values_to = "abundance") %>% 
        
        # Change start to change the taxonomy level used in the plot
        mutate(Phylum = word(taxa, start = 2, end = -1, sep = "_"),
               Phylum = str_replace_all(Phylum, "_.*", "")) %>% 
        inner_join(x = ., y = metadata, by = "Sample_ID") %>% 
        mutate(Sample_ID = factor(x = Sample_ID, 
                                  levels = c(samples_nnb,
                                             samples_lnb)),
               Phylum = factor(x = Phylum,
                               levels = phylums),
               Treatment = factor(x = Treatment, 
                                  levels = c("Water", "Escitalopram"))) %>% 
        
        ggplot(aes(x = Sample_ID, y = abundance, fill = Phylum, color = Phylum)) +
        
        geom_col() +
        
        facet_grid(.~Treatment, scales = "free") +
        
        scale_color_manual(values = cols,
                           breaks = sort(phylums)) +
        scale_fill_manual(values = cols,
                          breaks = sort(phylums)) +
        
        scale_y_continuous(labels = paste0(c(0, 25, 50, 75, 100), "%"), 
                           expand = c(0,0)) +
        scale_x_discrete(expand = c(0,0)) +
        
        new_scale_fill() +
        new_scale_color() +
        
        geom_tile(inherit.aes = FALSE,
                  aes(x = Sample_ID, y = -0.03, 
                      fill = Nestbuilding, color = Nestbuilding), 
                  height = 0.05, show.legend = FALSE) +
        
        scale_color_manual(values = my_cols_treatment) +
        scale_fill_manual(values = my_cols_treatment) +
        
        geom_vline(aes(xintercept = 11.5), linewidth = 1) +
  
        labs(x = "",
             y = "Relative abundance") +
        
        ggtitle("C") +
  
        theme_bw() + 
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.y = element_blank(),
             
              strip.text.x = element_text(size = 16, face = "bold"),
              # strip.background = element_rect(fill = "grey30",
              #                                 color = "grey30"),
              
              legend.title = element_blank(),
              legend.text = element_text(size = 16),
              legend.position = "bottom",
              
              axis.ticks = element_blank(),
              
              axis.text.x = element_blank(),
              axis.text.y = element_text(size = 12),
              axis.title.y = element_text(size = 16),
              
              plot.title = element_text(size = 40, face = "bold", hjust = -0.1))
        


# Alpha diversity ---------------------------------------------------------
library(ggh4x)

alpha_diversity <- get_asymptotic_alpha(species = counts, verbose = FALSE)

alpha_diversity$Nestbuilding <- metadata$Nestbuilding
alpha_diversity$Treatment <- metadata$Treatment

plot_alpha_div <- alpha_diversity %>% 
  pivot_longer(cols = c("Chao1", "Simpson Index", "Shannon Entropy"), 
               values_to = "Values", 
               names_to = "Index") %>% 
  
  # Add a new column with the full treatment name with an enter in between
  mutate(groups_axis = paste0(Nestbuilding,"\n",Treatment),
         groups_axis = factor(groups_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         groups_legend = paste0(Nestbuilding,"\n", Treatment),
         groups_legend = factor(groups_legend, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram"))) %>% 
  
  # Plotting
  ggplot(aes(x = groups_axis,
             y = Values,
             fill = groups_legend)) +
  
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3) +
  
  facet_wrap(~ Index, scales = "free") +
  # ggh4x::facetted_pos_scales(y = list(Index == "Chao1" ~ scale_y_continuous(breaks = c(45,50,55,60,65)),
  #                                     Index == "Shannon Entropy" ~ scale_y_continuous(breaks = seq(1.8, 2.5, length.out = 5)),
  #                                     Index == "Simpson Index" ~ scale_y_continuous(breaks = seq(0.60, 0.80, by = 0.05)))) +
  
  scale_fill_manual(values = filling_colors) +
  
  ylab("") + 
  xlab("") + 
  labs(fill = "Treatment") + 
  
  ggtitle("A") +
    
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.7, units = "in"),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(size = 40, face = "bold", hjust = -0.05))

#ggsave(filename = "results/alpha_div.jpg", plot = last_plot(), width = 20, height = 12, units = "in")


alpha_diversity$plate <- metadata$plate


# Spoiler : Nothing is significant
res_chao1 <-  summary(lm(Chao1 ~ Treatment * Nestbuilding + plate, alpha_diversity))
write(x = "Chao1", file = "results/Alpha_div_tests.txt")
capture.output(knitr::kable(res_chao1$coefficients, digits = 3), file = "results/Alpha_div_tests.txt", append = T)


res_simpson <-  summary(lm(`Simpson Index` ~ Treatment * Nestbuilding + plate, alpha_diversity))
write(x = "\n\nSimpson Index", file = "results/Alpha_div_tests.txt", append = T)
capture.output(knitr::kable(res_simpson$coefficients, digits = 3), file = "results/Alpha_div_tests.txt", append = T)


res_shannon <-  summary(lm(`Shannon Entropy` ~ Treatment * Nestbuilding + plate, alpha_diversity))
write(x = "\n\nShannon Entropy", file = "results/Alpha_div_tests.txt", append = T)
capture.output(knitr::kable(res_shannon$coefficients, digits = 3), file = "results/Alpha_div_tests.txt", append = T)




# Beta diversity ----------------------------------------------------------

genus.exp <- clr_c(counts)

genus.pca = prcomp(t(genus.exp))

pc1 <- round(genus.pca$sdev[1]^2/sum(genus.pca$sdev^2),4) *100
pc2 <- round(genus.pca$sdev[2]^2/sum(genus.pca$sdev^2),4) *100
pc3 <- round(genus.pca$sdev[3]^2/sum(genus.pca$sdev^2),4) *100
pc4 <- round(genus.pca$sdev[4]^2/sum(genus.pca$sdev^2),4) *100


#make data frame for plotting
pca  = data.frame(PC1 = genus.pca$x[,1], 
                  PC2 = genus.pca$x[,2], 
                  PC3 = genus.pca$x[,3], 
                  PC4 = genus.pca$x[,4])

#Include metadata information
pca$ID         = metadata$Larissa_ID
pca$Legend     = factor(paste(metadata$Nestbuilding, metadata$Treatment), levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram") )
pca$Phenotype  = metadata$Nestbuilding
pca$Treatment  = metadata$Treatment  
pca$plate      = metadata$plate

#Plot with ggplot2 
plot_beta_div <- pca %>% 
  mutate(Legend = paste0(Phenotype, "\n", Treatment),
         Legend = factor(x = Legend, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram"))) %>% 
  ggplot(aes(x = PC1, y = PC2, fill = Legend, colour = Legend))+
  stat_ellipse(show.legend = FALSE) + 
  geom_point(shape = 21, size = 6, colour = "black") + 
  xlab(paste("PC1: ", pc1,  "%", sep="")) + 
  ylab(paste("PC2: ", pc2,  "%", sep="")) +
  theme_bw() + 
  ggtitle("B") +
  scale_color_manual(values = filling_colors) +
  scale_fill_manual(values = filling_colors) +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 16),
        
        legend.key.size = unit(0.6, units = "in"),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.position = "bottom",
        
        plot.title = element_text(size = 40, face = "bold", hjust = -0.05))

#ggsave(filename = "results/beta_div.jpg", plot = last_plot(), width = 12, height = 8, units = "in")


# Calculate distance matrix
dist.genus = dist(t(genus.exp), method = "euclidean")

# Perform a PERMANOVA on the euclidean distance matrix to assess overall group differences
vegan::adonis2(dist.genus ~ Treatment * Nestbuilding + plate, data = metadata, method = "euclidean", permutations = 10000)

capture.output(x = vegan::adonis2(dist.genus ~ Treatment * Nestbuilding + plate, data = metadata, method = "euclidean", permutations = 10000),
               file = "results/permanova.txt")





# Figure 1 : Alpha & Beta div + Composition -------------------------------

library(patchwork)

figure1 <- ((plot_alpha_div / plot_beta_div + plot_layout(guides = "keep")) | plot_composition) + 
  plot_layout(guides = "auto",
              widths = unit(c(8,8), c("in", "in")))


ggsave("results/figure1.tiff", plot = figure1, device = "tiff", width = 20, height = 12, dpi = 300, units = "in")




# GBM ---------------------------------------------------------------------

set.seed(19970201)

GBMs <- read.csv(file = "data/GBMs_larissa.csv", row.names = 1, header = T)
GBMs <- GBMs[,c(2:ncol(GBMs))]
GBMs <- column_to_rownames(GBMs, var = "Description")

GMMs <- read.csv(file = "data/GMMs_larissa.csv", row.names = 1, header = T)
GMMs <- GMMs[,c(2:ncol(GMMs))]
GMMs <- column_to_rownames(GMMs, var = "Description")


prep_clr <- function(modules_df){
  
  # Keep only numbers and make sure they are numbers
  modules_counts <- apply(modules_df[3:ncol(modules_df)], c(1,2), function(x) as.numeric(as.character(x)))
  
  # Remove modules with 0 in more than 90% of samples
  n_zeroes_modules <- rowSums(modules_counts == 0)
  modules_clean <- modules_counts[n_zeroes_modules <= round(ncol(modules_counts) * 0.90), ]
  
  # Perform clr 
  clr <- clr_c(modules_clean)
  
}

# Run clr on the dataframes
clr_GBMs <- prep_clr(GBMs)
clr_GMMs <- prep_clr(GMMs)

# GBM & GMM 55 cols
# metadata$Sample_ID has 43

# Keep only the columns that have observations in the metadata
clr_GBMs <- clr_GBMs[, which(colnames(clr_GBMs) %in% metadata$Sample_ID)]
clr_GMMs <- clr_GMMs[, which(colnames(clr_GMMs) %in% metadata$Sample_ID)]


# Remove from the metadata those rows that were deleted because of the 90% threshold filter

# Find which rows are those
waldo::compare(sort(colnames(clr_GBMs)), sort(metadata$Sample_ID))

# Remove them
metadata <- metadata %>% 
  filter(!str_detect(Sample_ID, "X107|X108"))


# Reorder the columns to match the metadata
clr_GBMs <- clr_GBMs[, metadata$Sample_ID]
clr_GMMs <- clr_GMMs[, metadata$Sample_ID]



# GBM glms
GBMs.glm <- fw_glm(x = clr_GBMs,
                  f = ~ Treatment * Nestbuilding + plate,
                  metadata = metadata,
                  adjust.method = "BH")

hist(GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|)`, xlim = c(0, 1), breaks = 20)
hist(GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)




# GMM glms
GMMs.glm <- fw_glm(x = clr_GMMs,
                   f = ~ Treatment * Nestbuilding + plate,
                   metadata = metadata,
                   adjust.method = "BH")

hist(GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|)`, xlim = c(0, 1), breaks = 20)
hist(GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)


# Plot the GBM that show a group effect at q < 0.2 (For interaction between treatment and nestbuilding)
GBM_BH_interaction <- clr_GBMs[GBMs.glm[GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

interesting_gbms <- "Acetate|Butyrate|Propionate|KADC|Nitric|Tryptophan|Quinolinic"


GBM_BH_interaction %>%
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>%
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  #filter(str_detect(name, interesting_gbms)) %>% 
  mutate(name = str_replace(name, ".*ales_", ""),
         Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = paste(Nestbuilding, Group, sep = "\n"),
         Group_legend = paste(Nestbuilding, Group),
         Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  #geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  
  scale_fill_manual(values = filling_colors) +
  scale_color_manual(values = filling_colors) +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups", fill = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.5, units = "in"),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 18, face = "bold"),
        # legend.direction = "horizontal",
        legend.position = "bottom") 

# Saving the plot
#ggsave(filename = "results/GBM_interaction_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")




# Plot the GBM that show a group effect at q < 0.2 (For Nestbuilding, which are the same conditions in the study 1)
GBM_BH_nestbuilding <- clr_GBMs[GBMs.glm[GBMs.glm$`NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]


GBM_BH_nestbuilding %>%
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>% 
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  filter(str_detect(name, interesting_gbms)) %>%
  mutate(name = str_replace(name, ".*ales_", ""),
         Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = factor(Nestbuilding, levels = c("NNB", "LNB"))) %>% 
  #Group_axis = paste(Nestbuilding, Group, sep = "\n"),
  #Group_legend = paste(Nestbuilding, Group),
  #Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
  #Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  #geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 2) +
  
  scale_fill_manual(values = filling_colors[c(2,4)]) +
  scale_color_manual(values = filling_colors[c(2,4)]) +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups", fill = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.8, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        #legend.position = c(0.8, 0.15)
        #legend.position = "none"
        )


#ggsave(filename = "results/GBM_nestbuilding_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")




# Plot the GBM that show a group effect at q < 0.2 (For Nestbuilding, which are the same conditions in the study 1)
GBM_BH_treatment <- clr_GBMs[GBMs.glm[GBMs.glm$`TreatmentWater Pr(>|t|).BH` < 0.2,"feature"],]

GBM_BH_treatment %>%
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>% 
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  filter(str_detect(name, interesting_gbms)) %>%
  mutate(name = str_replace(name, ".*ales_", ""),
         Group_axis = factor(Group, levels = c("Water", "Escitalopram"))) %>% 
  #Group_axis = paste(Nestbuilding, Group, sep = "\n"),
  #Group_legend = paste(Nestbuilding, Group),
  #Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
  #Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  #geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  
  scale_fill_manual(values = c("grey60", "grey40")) +
  scale_color_manual(values = c("grey60", "grey40")) +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups", fill = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))


#ggsave(filename = "results/GBM_escitalopram_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")



# Figure 2 : Gut Brain Modules --------------------------------------------

# Get all the GBM that showed significant results in the interaction
gbm_interaction_data <- GBM_BH_interaction %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Interaction")

# Get all the GBM that showed significant results due to the treatment
gbm_treatment_data <- GBM_BH_treatment %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Escitalopram\nTreatment")

# Get all the GBM that showed significant results due to the basal behavior phenotype
gbm_nestbuilding_data <- GBM_BH_nestbuilding %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Nestbuilding\nBehavior")



gbm_data <- rbind(gbm_interaction_data, gbm_nestbuilding_data, gbm_treatment_data)


gbm_values <- gbm_data %>% 
  filter(str_detect(name, interesting_gbms)) %>% 
  inner_join(., metadata, by = "Sample_ID") %>% 
  mutate(group = paste0(Nestbuilding,"\n",Treatment)) %>% 
  group_by(name, group) %>% 
  mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name) %>% 
  mutate(mean_value_z = scale(mean_value, center = TRUE, scale = TRUE),
         name_fct = str_replace(name, " \\(", "\n\\("),
         name_fct = factor(name_fct, levels = c(
                                                 # Permeability
                                                 "Nitric oxide degradation II\n(NO reductase)",
                                                 "Nitric oxide synthesis II\n(nitrite reductase)",
                                                 "Nitric oxide synthesis I\n(NO synthase)",
                                                 # SCFA
                                                 "Isovaleric acid synthesis II\n(KADC pathway)",
                                                 "Propionate synthesis II",
                                                 "Butyrate synthesis I",
                                                 "Butyrate synthesis II",
                                                 "Acetate synthesis I",
                                                 # Tryptophan metabolism
                                                 "Quinolinic acid degradation",
                                                 "Tryptophan degradation")),
         group = factor(group, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")))
# # Calculate the mean of the score for each GBM within each condition
# gbm_values <- gbm_data %>% 
#   filter(str_detect(name, interesting_gbms)) %>% 
#   mutate(groups = paste0(Group,"\n", Nestbuilding),
#          groups = factor(groups, levels = c("Water\nNNB","Water\nLNB", "Escitalopram\nNNB", "Escitalopram\nLNB")))
#   group_by(name, groups) %>% 
#   summarise(mean_value = mean(value, na.rm = T), .groups = "drop")


# Determine if the interaction or each factor by its own drive significant differences (answer as logical value)
gbm_significance_caterogical <-   gbm_data %>% 
  inner_join(., metadata, by = "Sample_ID") %>% 
  mutate(groups = paste0(Nestbuilding,"\n", Treatment),
         value = TRUE,
         question = factor(x = question, levels = c("Interaction", "Nestbuilding\nBehavior", "Escitalopram\nTreatment")),
         groups = factor(groups, levels = c("NNB\nWater", "LNB\nWater", "NNB\nEscitalopram", "LNB\nEscitalopram"))) %>% 
  select(name, groups, question, value) %>% 
  unique()



significant_gbm_of_interest <- gbm_values %>% pull(name) %>% unique() %>% as.character()

gbm_significance <- GBMs.glm %>% 
  as_tibble() %>% 
  select(feature,
         `TreatmentWater:NestbuildingNNB Estimate`,
         `TreatmentWater Estimate`,
         `NestbuildingNNB Estimate`) %>% 
  filter(feature %in% significant_gbm_of_interest) %>% 
  pivot_longer(cols = !matches("feature"),
               names_to = "question",
               values_to = "estimate") %>% 
  mutate(question = case_when(str_detect(question, ".Water:Nestbuilding.") ~ "Interaction",
                              str_detect(question, "NestbuildingNNB.") ~ "Nestbuilding\nBehavior",
                              TRUE ~ "Escitalopram\nTreatment"),
         question = as.factor(question)) %>% 
  rename("name" = feature) %>% 
  right_join(gbm_significance_caterogical, ., by = c("name", "question")) %>% 
  mutate(value = replace_na(value, FALSE),
         label = ifelse(value, "🞸", ""),
         name_fct = str_replace(name, " \\(", "\n\\("),
         name_fct = factor(name_fct, levels = c(
                                                 # Permeability
                                                 "Nitric oxide degradation II\n(NO reductase)",
                                                 "Nitric oxide synthesis II\n(nitrite reductase)",
                                                 "Nitric oxide synthesis I\n(NO synthase)",
                                                 # SCFA
                                                 "Isovaleric acid synthesis II\n(KADC pathway)",
                                                 "Propionate synthesis II",
                                                 "Butyrate synthesis I",
                                                 "Butyrate synthesis II",
                                                 "Acetate synthesis I",
                                                 # Tryptophan metabolism
                                                 "Quinolinic acid degradation",
                                                 "Tryptophan degradation")),
         group = factor(groups, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")))
  







# Tile plot with the values for each GBM within each condition  
plot_GBM_panelA_fig2 <- (gbm_values %>% 
  ggplot(aes(x = group, y = name_fct, fill = mean_value_z)) +
  geom_tile(color = "black") + 
  scale_fill_gradientn(colors = c("#FDFEFE", "#45B39D", "#0E6655"),
                       breaks = seq(from = -2, to = 2, by = 1),
                       limits = c(-2,2)) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  
  labs(fill = "Z-score") +
  
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.ticks = element_blank(),
        
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.6, "in"),
        legend.key.width = unit(0.3, "in"),
        legend.title = element_text(size = 20),
        
        panel.grid = element_blank(),
        
        plot.title = element_text(size = 40, face = "bold", hjust = -0.2, vjust = -0.15)) +
  
  
  
#' Tile plot with a logical answer to the question: Are the differences seen for each GBM
#' drive by an interaction, or either of the factors included in the analysis.
gbm_significance %>% 
  ggplot(aes(x = question, y = name_fct, fill = estimate)) +
  geom_tile(color = "black") +
  geom_text(aes(label = label), size = 8) +
  
  scale_fill_gradientn(colors = c("#2980B9",  "#A9CCE3", "#FDFEFE", "#F9E79F", "#F1C40F"),
                       labels = seq(-1.5, 1.5, 0.75),
                       breaks = seq(-1.5, 1.5, 0.75),
                       limits = c(-1.8,1.8)) +
  scale_x_discrete(expand = c(0,0), position = "top") +
  scale_y_discrete(expand = c(0,0)) +
  
  labs(fill = "Estimate") +
  
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 16),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        
        panel.grid = element_blank(),
        
        legend.text = element_text(size = 16),
        legend.key.height = unit(0.6, "in"),
        legend.key.width = unit(0.3, "in"),
        legend.title = element_text(size = 20)) +
  
  plot_layout(guides = "collect",
              widths = c(3.5,1.8)))


ggsave(filename = "results/figure2.tiff", device = "tiff", plot = plot_GBM_panelA_fig2,
       width = 17, height = 10, units = "in")



# Figure 2 : GMM ----------------------------------------------------------


GMM_BH_interaction <- clr_GMMs[GMMs.glm[GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]
GMM_BH_nestbuilding <- clr_GMMs[GMMs.glm[GMMs.glm$`NestbuildingNNB Pr(>|t|).BH` < 0.2, "feature"],]
GMM_BH_treatment <- clr_GMMs[GMMs.glm[GMMs.glm$`TreatmentWater Pr(>|t|)` < 0.2, "feature"],]


# Get all the GBM that showed significant results in the interaction
gmm_interaction_data <- GMM_BH_interaction %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Interaction")

# Get all the GBM that showed significant results due to the treatment
gmm_treatment_data <- GMM_BH_treatment %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Escitalopram\nTreatment")

# Get all the GBM that showed significant results due to the basal behavior phenotype
gmm_nestbuilding_data <- GMM_BH_nestbuilding %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% 
  pivot_longer(!matches("Sample_ID")) %>% 
  mutate(question = "Nestbuilding\nBehavior")



gmm_data <- rbind(gmm_interaction_data, gmm_nestbuilding_data, gmm_treatment_data)


#gmm_values <- 
unique(gmm_data$name)
  filter(str_detect(name, interesting_gbms)) %>% 
  inner_join(., metadata, by = "Sample_ID") %>% 
  mutate(group = paste0(Nestbuilding,"\n",Treatment)) %>% 
  group_by(name, group) %>% 
  mutate(mean_value = mean(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(name) %>% 
  mutate(mean_value_z = scale(mean_value, center = TRUE, scale = TRUE),
         name_fct = str_replace(name, " \\(", "\n\\("),
         name_fct = factor(name_fct, levels = c(
           # Permeability
           "Nitric oxide degradation II\n(NO reductase)",
           "Nitric oxide synthesis II\n(nitrite reductase)",
           "Nitric oxide synthesis I\n(NO synthase)",
           # SCFA
           "Isovaleric acid synthesis II\n(KADC pathway)",
           "Propionate synthesis II",
           "Butyrate synthesis I",
           "Butyrate synthesis II",
           "Acetate synthesis I",
           # Tryptophan metabolism
           "Quinolinic acid degradation",
           "Tryptophan degradation")))
# # Calculate the mean of the score for each GBM within each condition
# gbm_values <- gbm_data %>% 
#   filter(str_detect(name, interesting_gbms)) %>% 
#   mutate(groups = paste0(Group,"\n", Nestbuilding),
#          groups = factor(groups, levels = c("Water\nNNB","Water\nLNB", "Escitalopram\nNNB", "Escitalopram\nLNB")))
#   group_by(name, groups) %>% 
#   summarise(mean_value = mean(value, na.rm = T), .groups = "drop")


# Determine if the interaction or each factor by its own drive significant differences (answer as logical value)
gbm_significance_caterogical <-   gbm_data %>% 
  inner_join(., metadata, by = "Sample_ID") %>% 
  mutate(groups = paste0(Nestbuilding,"\n", Treatment),
         value = TRUE,
         question = factor(x = question, levels = c("Interaction", "Nestbuilding\nBehavior", "Escitalopram\nTreatment")),
         groups = factor(groups, levels = c("NNB\nWater", "LNB\nWater", "NNB\nEscitalopram", "LNB\nEscitalopram"))) %>% 
  select(name, groups, question, value) %>% 
  unique()



significant_gbm_of_interest <- gbm_values %>% pull(name) %>% unique() %>% as.character()

gbm_significance <- GBMs.glm %>% 
  as_tibble() %>% 
  select(feature,
         `TreatmentWater:NestbuildingNNB Estimate`,
         `TreatmentWater Estimate`,
         `NestbuildingNNB Estimate`) %>% 
  filter(feature %in% significant_gbm_of_interest) %>% 
  pivot_longer(cols = !matches("feature"),
               names_to = "question",
               values_to = "estimate") %>% 
  mutate(question = case_when(str_detect(question, ".Water:Nestbuilding.") ~ "Interaction",
                              str_detect(question, "NestbuildingNNB.") ~ "Nestbuilding\nBehavior",
                              TRUE ~ "Escitalopram\nTreatment"),
         question = as.factor(question)) %>% 
  rename("name" = feature) %>% 
  right_join(gbm_significance_caterogical, ., by = c("name", "question")) %>% 
  mutate(value = replace_na(value, FALSE),
         label = ifelse(value, "🞸", ""),
         name_fct = str_replace(name, " \\(", "\n\\("),
         name_fct = factor(name_fct, levels = c(
           # Permeability
           "Nitric oxide degradation II\n(NO reductase)",
           "Nitric oxide synthesis II\n(nitrite reductase)",
           "Nitric oxide synthesis I\n(NO synthase)",
           # SCFA
           "Isovaleric acid synthesis II\n(KADC pathway)",
           "Propionate synthesis II",
           "Butyrate synthesis I",
           "Butyrate synthesis II",
           "Acetate synthesis I",
           # Tryptophan metabolism
           "Quinolinic acid degradation",
           "Tryptophan degradation")))








# Tile plot with the values for each GBM within each condition  
plot_GBM_panelA_fig2 <- (gbm_values %>% 
                           ggplot(aes(x = group, y = name_fct, fill = mean_value_z)) +
                           geom_tile(color = "black") + 
                           scale_fill_gradientn(colors = c("#FDFEFE", "#16A085", "#0E6655"),
                                                breaks = seq(from = -2, to = 2, by = 1),
                                                limits = c(-2,2)) +
                           scale_x_discrete(expand = c(0,0), position = "top") +
                           scale_y_discrete(expand = c(0,0)) +
                           
                           labs(fill = "") +
                           
                           #ggtitle("A") +
                           
                           theme_bw() +
                           theme(axis.title = element_blank(),
                                 axis.text.x = element_text(size = 16),
                                 axis.text.y = element_text(size = 16),
                                 axis.ticks = element_blank(),
                                 
                                 legend.text = element_text(size = 16),
                                 legend.key.height = unit(0.6, "in"),
                                 legend.key.width = unit(0.3, "in"),
                                 
                                 panel.grid = element_blank(),
                                 
                                 plot.title = element_text(size = 40, face = "bold", hjust = -0.2, vjust = -0.15)) +
                           
                           
                           
                           #' Tile plot with a logical answer to the question: Are the differences seen for each GBM
                           #' drive by an interaction, or either of the factors included in the analysis.
                           gbm_significance %>% 
                           ggplot(aes(x = question, y = name_fct, fill = estimate)) +
                           geom_tile(color = "black") +
                           geom_text(aes(label = label), size = 8) +
                           
                           scale_fill_gradientn(colors = c("#2980B9",  "#A9CCE3", "#FDFEFE", "#F9E79F", "#F1C40F"),
                                                labels = seq(-1.5, 1.5, 0.75),
                                                breaks = seq(-1.5, 1.5, 0.75),
                                                limits = c(-1.8,1.8)) +
                           scale_x_discrete(expand = c(0,0), position = "top") +
                           scale_y_discrete(expand = c(0,0)) +
                           
                           theme_bw() +
                           theme(axis.text.y = element_blank(),
                                 axis.text.x = element_text(size = 16),
                                 axis.ticks = element_blank(),
                                 axis.title = element_blank(),
                                 
                                 panel.grid = element_blank(),
                                 
                                 legend.text = element_text(size = 16),
                                 legend.key.height = unit(0.6, "in"),
                                 legend.key.width = unit(0.3, "in")) +
                           
                           plot_layout(guides = "collect",
                                       widths = c(3.5,1.8)))








# GMM ---------------------------------------------------------------------


# Plot the GMM that show a group effect at q < 0.2
GMM_BH_interaction <- clr_GMMs[GMMs.glm[GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

GMM_BH_interaction %>%
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>%
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  mutate(name = str_replace(name, ".*ales_", ""),
         Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = paste(Nestbuilding, Group, sep = "\n"),
         Group_legend = paste(Nestbuilding, Group),
         Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 5) +
  
  scale_fill_manual(values = filling_colors) +
  scale_color_manual(values = filling_colors) +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups",
       fill = "Groups") +
  
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) 

# Saving the results
#ggsave(filename = "results/GMM_q02.jpg", plot = last_plot(), width = 12, height = 8, units = "in")



GMM_BH_nestbuilding <- clr_GMMs[GMMs.glm[GMMs.glm$`NestbuildingNNB Pr(>|t|).BH` < 0.2, "feature"],]
                                        
GMM_BH_nestbuilding %>% 
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>% 
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  mutate(name = str_replace(name, ".*ales_", ""),
         Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = factor(Nestbuilding, levels = c("NNB", "LNB"))) %>% 
  #Group_axis = paste(Nestbuilding, Group, sep = "\n"),
  #Group_legend = paste(Nestbuilding, Group),
  #Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
  #Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  #geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  
  scale_fill_manual(values = filling_colors[c(2,4)]) +
  scale_color_manual(values = filling_colors[c(2,4)]) +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups", fill = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.8, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.8, 0.15))




# DEA genus ---------------------------------------------------------------

# Reading the metadata
metadata = read.delim("data/metadata.csv", sep = ",")


# Run the linear models
genus.glm = fw_glm(x = genus.exp,
                   f = ~ Treatment * Nestbuilding + plate,
                   metadata = metadata,
                   adjust.method = "BH")


# Look at the distribution of p-values
hist(genus.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|)`, xlim = c(0, 1), breaks = 20)
hist(genus.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)


# Plot the strains that show a group effect at q < 0.4, because no lower q showed strains
genBH <- genus.exp[genus.glm[genus.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]
#genBH <- genus.exp[genus.glm[genus.glm$`NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

genBH %>%
  t() %>%
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>% 
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  mutate(name = str_replace(name, ".*ales_", ""),
         Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = paste(Nestbuilding, Group, sep = "\n"),
         Group_legend = paste(Nestbuilding, Group),
         Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  
  scale_fill_brewer(type = "qual", palette = "Paired") +
  scale_color_brewer(type = "qual", palette = "Paired") +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups",
       fill = "Groups") +
  
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.6, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.65, 0.2)) 

#ggsave(filename = "results/strains_q04.jpg", plot = last_plot(), width = 12, height = 8, units = "in")


# Histograms --------------------------------------------------------------

# # GBM : Interaction
# png(filename = "results/histograms/GBM_interaction.png",width = 8, height = 6, units = "in", res = 300)
# hist(GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # GBM : Treatment
# png(filename = "results/histograms/GBM_treatment.png", width = 8, height = 6, units = "in", res = 300)
# hist(GBMs.glm$`TreatmentWater Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # GBM : Nestbuilding
# png(filename = "results/histograms/GBM_nestbuild.png", width = 8, height = 6, units = "in", res = 300)
# hist(GBMs.glm$`NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# 
# 
# 
# # GMM : Interaction
# png(filename = "results/histograms/GMM_interaction.png",width = 8, height = 6, units = "in", res = 300)
# hist(GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # GMM : Treatment
# png(filename = "results/histograms/GMM_treatment.png", width = 8, height = 6, units = "in", res = 300)
# hist(GMMs.glm$`TreatmentWater Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # GMM : Nestbuilding
# png(filename = "results/histograms/GMM_nestbuild.png", width = 8, height = 6, units = "in", res = 300)
# hist(GMMs.glm$`NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# 
# 
# 
# # Genus : Interaction
# png(filename = "results/histograms/GENUS_interaction.png",width = 8, height = 6, units = "in", res = 300)
# hist(genus.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # Genus : Treatment
# png(filename = "results/histograms/GENUS_treatment.png", width = 8, height = 6, units = "in", res = 300)
# hist(genus.glm$`TreatmentWater Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()
# # Genus : Nestbuilding
# png(filename = "results/histograms/GENUS_nestbuild.png", width = 8, height = 6, units = "in", res = 300)
# hist(genus.glm$`NestbuildingNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# dev.off()




# Diving into metabolites -------------------------------------------------


# Importing the metabolomic dataset

library(readxl)


# Include all the data but the columns with ratios between bodyparts or molecules
metabolomics <- read_xlsx("data/metabolomic.xlsx", skip = 1, range = "B2:BT50") %>%
  mutate(Phenotype = ifelse(Phenotype == 0, yes = "NNB", no = "LNB"),
         Treatment = ifelse(Treatment == 0, yes = "Water", no = "Escitalopram"),
         Condition = paste(Phenotype, Treatment)) %>% 
  select(where(~ !all(is.na(.))))


# Formatting the dataset to be handled by fw_glm function
metabol.exp <- metabolomics %>% 
  select(!matches("[(]|Condition|Phenotype|Treatment")) %>% 
  mutate(`Animal No` = paste0("X", `Animal No`)) %>% 
  column_to_rownames("Animal No") %>% 
  t()

# metabol.exp[is.na(metabol.exp)] <- 0
# 
# metabol.exp <- metabol.exp %>% clr_c()

  
# Formatting the metadata to be handled by fw_glm function
# metabolomic.xslx has the metadata as well
metadata.metabol <- metabolomics %>% 
  rename(Sample_ID = `Animal No`) %>% 
  mutate(Sample_ID = paste0("X", Sample_ID)) %>% 
  select(matches("Sample_ID|Condition|Phenotype|Treatment"))


# Run lineal models
metabolomics.glm = fw_glm(x = metabol.exp,
                          f = ~ Phenotype * Treatment,
                          metadata = metadata.metabol,
                          adjust.method = "BH")



# Look at the distribution of p-values
hist(metabolomics.glm$`PhenotypeNNB:TreatmentWater Pr(>|t|)`, xlim = c(0, 1), breaks = 20)

hist(metabolomics.glm$`PhenotypeNNB:TreatmentWater Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
hist(metabolomics.glm$`PhenotypeNNB Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
hist(metabolomics.glm$`TreatmentWater Pr(>|t|).BH`, xlim = c(0, 1), breaks = 20)
# After BH correction only Treatment threw significant results


# Plot the strains that show a group effect at q < 0.4, because no lower q showed strains
metabolBH<- metabol.exp[metabolomics.glm[metabolomics.glm$`PhenotypeNNB:TreatmentWater Pr(>|t|).BH` < 0.2,"feature"],]

metabolBH <- metabol.exp[metabolomics.glm[metabolomics.glm$`TreatmentWater Pr(>|t|).BH` < 0.2,"feature"],]

metabolBH %>% 
  t() %>%
  as_tibble() %>%
  add_column(Group = metadata.metabol$Treatment,
             Nestbuilding = metadata.metabol$Phenotype) %>% 
  pivot_longer(!c("Group", "Nestbuilding")) %>% 
  mutate(Group = factor(Group, levels = c("Water", "Escitalopram")),
         Group_axis = paste(Nestbuilding, Group, sep = "\n"),
         Group_legend = paste(Nestbuilding, Group),
         Group_axis = factor(Group_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         Group_legend = factor(Group_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
  ggplot(aes(x = Group_axis,
             y = value,
             color = Group_axis,
             fill = Group_axis)) +
  
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  
  scale_fill_brewer(type = "qual", palette = "Paired") +
  scale_color_brewer(type = "qual", palette = "Paired") +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups",
       fill = "Groups") +
  
  theme_bw() + 
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.83, 0.08))

#ggsave(filename = "results/metabolites_treatment.jpg", plot = last_plot(), width = 18, height = 12, units = "in")




# Correlations between GBM and Metabolites --------------------------------

library(anansi)

GBM_BH <- clr_GBMs[GBMs.glm[GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

#' The names in GBM_BH are in a format different from the one used in metabol.exp. We need to reformat them:
#' 
#' We first get the list of names that was used by Larissa to name the samples
larissa_names <- metadata[metadata$Sample_ID %in% colnames(GBM_BH), "Larissa_ID"]
# We use that list of names to rename the columns in the GBM_BH dataset
colnames(GBM_BH) <- larissa_names

# Change the "X" in the anmes in metabol.exp for a "L"
new_metabol_names <- str_replace(colnames(metabol.exp), "X", "L")
colnames(metabol.exp) <- new_metabol_names

# From metabol.exp we only keep the columns shared between that dataframe and GMB_BH
new_metabol.exp <- metabol.exp[,colnames(metabol.exp) %in% colnames(GBM_BH)]


# Ordering the names and trasposing the dataframes
tableY <- t(new_metabol.exp[,sort(colnames(new_metabol.exp))])
tableX <- t(GBM_BH[,sort(colnames(GBM_BH))])

# Make a dictionary 
dictionary <- set_names(rep(list(rownames(GBM_BH)), length(rownames(metabol.exp))), 
                        str_replace(rownames(metabol.exp), "/ ", "/"))

# Make a dictionary with interactions of interest

metabolites_of_interest <- list("Plasma 5HIAA/5HT", "COLON 5HIAA", "PFC QUIN", "Plasma 5HT")
gbm_of_interest <- list("Tryptophan degradation", 
                      "Quinolinic acid degradation", 
                      c("Acetate synthesis I","Nitric oxide synthesis II (nitrite reductase)", "p-Cresol degradation"),
                      c("Acetate synthesis I", "Glutamate synthesis I", "Isovaleric acid synthesis II (KADC pathway)"))


dictionary_of_interest <- set_names(gbm_of_interest, 
                                    metabolites_of_interest)

# Generating the web
web <- weaveWebFromTables(tableY, cbind(tableX, name = 1), dictionary)
web <- weaveWebFromTables(tableY, tableX, dictionary_of_interest)

anansi_metadata <- metadata[metadata$Larissa_ID %in% rownames(tableX),] %>% 
  mutate(Legend = paste(Nestbuilding, Treatment))


anansi_out <- anansi(web    = web,
                     method = "pearson",
                     groups = anansi_metadata$Legend,
                     adjust.method = "BH",
                     verbose = T)




# Plotting the Rho coefficients for each combination of Metabolites and GBMs
anansiLong <- spinToLong(anansi_output = anansi_out, translate = F)

anansiLong <- anansiLong[anansiLong$model_full_q.values < 0.1,]



# Figure 3 ----------------------------------------------------------------

plot_anansi_interaction <- anansiLong %>% 
  filter(feature_Y == "Plasma 5HT") %>% 
  mutate(type = str_replace(type, " ", "\n"),
         type = factor(x = type, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram", "All")),
         feature_X = str_replace(feature_X, pattern = " \\(", replacement = "\n\\(")) %>% 
  filter(model_disjointed_p.values < 0.05) %>% 
  ggplot(aes(x      = r.values, 
           y      = feature_X, 
           fill   = type, 
           #alpha  = model_disjointed_p.values < 0.05
           )) + 
  
  #Make a vertical dashed red line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  
  #Points show  raw correlation coefficients
  geom_point(shape = 21, size = 5) +
  
  #fix the scales, labels, theme and other layout
  scale_y_discrete(position = "right", limits = rev) +
  scale_x_continuous(limits = c(-1, 1), expand = c(0.01, 0.01)) +
  
  #facet per compound
  ggforce::facet_col(~feature_Y, scales = "free_y", space = "free") +

  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
  
  scale_fill_manual(values = c("NNB\nWater" = filling_colors[1], 
                               "NNB\nEscitalopram"  = filling_colors[2],
                               "LNB\nWater"  = filling_colors[3], 
                               "LNB\nEscitalopram" = filling_colors[4],
                               "All"        = "gray"))+
  ylab("") + 
  xlab("Pearson's rho") +
  
  #ggtitle(label = "B") +
  
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(x = 0.8, units = "in"),
        
        axis.ticks = element_blank(),
        axis.text = element_text(size = 16),
        
        strip.text = element_text(size = 14, face = "bold"),
        
        plot.title = element_text(size = 40, hjust = -0.1))






plot_anansi_nestbuilding <- anansiLong %>% 
  filter(feature_Y != "Plasma 5HT") %>% 
  mutate(type = str_replace(type, " ", "\n"),
         type = factor(x = type, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram", "All")),
         feature_X = str_replace(feature_X, pattern = " \\(", replacement = "\n\\(")) %>% 
  filter(model_disjointed_p.values < 0.05) %>% 
  ggplot(aes(x      = r.values, 
             y      = feature_X, 
             fill   = type, 
             #alpha  = model_disjointed_p.values < 0.05
  )) + 
  
  #Make a vertical dashed red line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  
  #Points show  raw correlation coefficients
  geom_point(shape = 21, size = 5) +
  
  #fix the scales, labels, theme and other layout
  scale_y_discrete(position = "right", limits = rev) +
  scale_x_continuous(limits = c(-1, 1), expand = c(0.01, 0.01)) +
  
  #facet per compound
  ggforce::facet_col(~feature_Y, scales = "free_y", space = "free") +
  
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
  
  scale_fill_manual(values = c("NNB\nWater" = filling_colors[1], 
                               "NNB\nEscitalopram"  = filling_colors[2],
                               "LNB\nWater"  = filling_colors[3], 
                               "LNB\nEscitalopram" = filling_colors[4],
                               "All"        = "gray"))+
  ylab("") + 
  xlab("Pearson's rho") +
  
  #ggtitle(label = "A") +
  
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(x = 0.8, units = "in"),
        
        axis.ticks = element_blank(),
        axis.text = element_text(size = 16),
        
        strip.text = element_text(size = 14, face = "bold"),
        
        plot.title = element_text(size = 40, hjust = -0.1))


#ggsave(filename = "results/all_correlations_rho-values.jpg", height = 12, width = 6, units = "in")

library(patchwork)

figure3<- (plot_anansi_nestbuilding | plot_anansi_interaction) + 
  plot_layout(guides = "collect", 
              widths = unit(c(3,3), c("in", "in")),
              heights = unit(c(5,5), c("in", "in"))) + 
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 40, face = "bold"),
        legend.position = "bottom")

ggsave(plot = figure3, filename = "results/figure3.tiff", device = "tiff", width = 16, height = 8, units = "in")
#



# Plot the actual values to see the trends in the interactions

# This table is the same as tableX but it only inclues pathways related to tryptophan metabolism
tryp_pathways <- c("Quinolinic acid degradation", "Tryptophan degradation")

tableX2 <- tableX[,colnames(tableX) %in% tryp_pathways]
dictionary2 <- set_names(rep(list(tryp_pathways), length(rownames(metabol.exp))), 
                        str_replace(rownames(metabol.exp), "/ ", "/"))


web2 <- weaveWebFromTables(tableY, tableX2, dictionary2)

anansi_metadata <- metadata[metadata$Larissa_ID %in% rownames(tableX),] %>% 
  mutate(Legend = paste(Nestbuilding, Treatment))


anansi_out2 <- anansi(web    = web2,
                     method = "pearson",
                     groups = anansi_metadata$Legend,
                     adjust.method = "BH",
                     verbose = T)

anansi_out2 %>% spinToLong() %>%
  filter(model_disjointed_p.values < 0.05) %>% 
  ggplot(aes(x      = r.values, 
             y      = feature_X, 
             fill   = type, 
             #alpha  = model_disjointed_p.values < 0.05
             )) + 
  
  #Make a vertical dashed red line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", colour = "red")+
  
  #Points show  raw correlation coefficients
  geom_point(shape = 21, size = 3) + 
  
  #facet per compound
  ggforce::facet_col(~feature_Y, space = "free", scales = "free_y") + 
  
  #fix the scales, labels, theme and other layout
  scale_y_discrete(limits = rev, position = "right") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 1/3)) +
  scale_fill_manual(values = c("NNB Water" = "#A6CEE3", 
                               "NNB Escitalopram"  = "#1F78B4", 
                               "LNB Water"  = "#B2DF8A", 
                               "LNB Escitalopram" = "#33A02C",
                               "All"        = "gray"))+
  theme_bw() + 
  ylab("") + 
  xlab("Pearson's rho")


#ggsave("results/tryp_correlations_rho-values.jpg", plot = last_plot(), width = 8, height = 4, units = "in")



#Let's look at all canonical interactions that also have a sufficiently well-fitting model:
outPlots = spinToPlots(anansi_out2,
                       target = anansi_out2@input@web@dictionary &
                         anansi_out2@output@model_results$modelfit@q.values < 0.1 & 
                         anansi_out2@output@model_results$disjointed@p.values < 0.05) 

#load ggplot2 and patchwork for plotting
library(patchwork)

plotted = lapply(outPlots, FUN = function(p){
  
  #Main ggplot call
  ggplot(data = p$data, aes(x = X, y = Y, fill = groups, 
                            #color = groups
                            )) +
    
    #Establish geoms:
    geom_point(shape = 21) +
    geom_smooth(method = "lm") +
    theme_bw() +
    
    #Improve annotation:
    scale_fill_manual(values = c("NNB Water" = "#A6CEE3", 
                                 "NNB Escitalopram"  = "#1F78B4", 
                                 "LNB Water"  = "#B2DF8A", 
                                 "LNB Escitalopram" = "#33A02C",
                                 "All"        = "gray")) +
    
    # scale_color_manual(values = c("NNB Water" = "#A6CEE3", 
    #                              "NNB Escitalopram"  = "#1F78B4", 
    #                              "LNB Water"  = "#B2DF8A", 
    #                              "LNB Escitalopram" = "#33A02C",
    #                              "All"        = "gray")) +
    
    ylab(p$name[1]) +
    xlab(p$name[2]) +
    ggtitle(paste(p$name[1], "vs microbial", p$name[2]))
  
})

# Use patchwork to put all the plots together
wrap_plots(plotted) + plot_layout(guides = 'collect')
#ggsave("results/tryp_pathways_raw_correlations.jpg", plot = last_plot(), width = 10, height = 8, units = "in")


# Call patchwork arrange the plots for serotonin precursor
wrap_plots(plotted[[4]]) + plot_layout(guides = 'collect')
ggsave("results/5HIAA-5THvstryptophan.jpg", plot = last_plot(), width = 6, height = 4, units = "in")


