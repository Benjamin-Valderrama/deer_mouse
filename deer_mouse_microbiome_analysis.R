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

#Define working directory
setwd(paste0(getwd(), "/southafrica_study/"))

set.seed(12345)
#load data
counts  = read.delim("data/genus_table_from_dada2.csv", sep = ",", row.names = 1)
metadata = read.delim("data/metadata.csv", sep = ",")

#Remove plasma for now as it contains missing values. 
counts = counts[,metadata$Sample_ID]


#' TO DO LIST:
#'
#' 0. Stacked barplot with the composition of the microbiome
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
#'     GENUS fw_glm(F = ~ Treatment * Nestbuilding + plate) [DONE]
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
#' 1) alpha-div, beta-div and stacked barplot
#' 2) GBM and GMM
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




# Alpha diversity ---------------------------------------------------------
library(ggh4x)

alpha_diversity <- get_asymptotic_alpha(species = counts, verbose = FALSE)

alpha_diversity$Nestbuilding <- metadata$Nestbuilding
alpha_diversity$Treatment <- metadata$Treatment

alpha_diversity %>% 
  pivot_longer(cols = c("Chao1", "Simpson Index", "Shannon Entropy"), 
               values_to = "Values", 
               names_to = "Index") %>% 
  
  # Add a new column with the full treatment name with an enter in between
  mutate(groups_axis = paste0(Nestbuilding,"\n",Treatment),
         groups_axis = factor(groups_axis, levels = c("NNB\nWater", "NNB\nEscitalopram", "LNB\nWater", "LNB\nEscitalopram")),
         groups_legend = paste(Nestbuilding, Treatment),
         groups_legend = factor(groups_legend, levels = c("NNB Water", "NNB Escitalopram", "LNB Water", "LNB Escitalopram"))) %>% 
  
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
  
  scale_fill_brewer(type = "qual", palette = "Paired") +
  
  ylab("") + 
  xlab("") + 
  labs(fill = "Treatment") + 
  
  theme_bw() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

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
ggplot(data = pca, aes(x = PC1, y = PC2, fill = Legend, colour = Legend))+
  stat_ellipse(show.legend = FALSE) + 
  geom_point(shape = 21, size = 5, colour = "black") + 
  xlab(paste("PC1: ", pc1,  "%", sep="")) + 
  ylab(paste("PC2: ", pc2,  "%", sep="")) +
  theme_bw() + 
  ggtitle("Principal Component Ananlysis of Microbiome Aitchison Distance") +
  scale_color_brewer(type = "qual", palette = "Paired") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  coord_equal()

#ggsave(filename = "results/beta_div.jpg", plot = last_plot(), width = 12, height = 8, units = "in")


# Calculate distance matrix
dist.genus = dist(t(genus.exp), method = "euclidean")

# Perform a PERMANOVA on the euclidean distance matrix to assess overall group differences
vegan::adonis2(dist.genus ~ Treatment * Nestbuilding + plate, data = metadata, method = "euclidean", permutations = 10000)

capture.output(x = vegan::adonis2(dist.genus ~ Treatment * Nestbuilding + plate, data = metadata, method = "euclidean", permutations = 10000),
               file = "results/permanova.txt")





# Functional modules ------------------------------------------------------

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
GBM_BH <- clr_GBMs[GBMs.glm[GBMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

GBM_BH %>%
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
  
  #geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 1/2, coef = 100) +
  geom_point(shape = 21, size = 3, color = "black") +
  
  facet_wrap(~ name, scales = "free_y", ncol = 3) +
  
  scale_fill_brewer(type = "qual", palette = "Paired", name = "Groups") +
  scale_color_brewer(type = "qual", palette = "Paired", name = "Groups") +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12, face="bold"),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.7, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.65, 0.04)) 

# Saving the plot
#ggsave(filename = "results/GBM_interaction_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")




# Plot the GBM that show a group effect at q < 0.2 (For Nestbuilding, which are the same conditions in the study 1)
GBM_BH <- clr_GBMs[GBMs.glm[GBMs.glm$`NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

GBM_BH %>%
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
  
  scale_fill_brewer(type = "qual", palette = "Paired", name = "Groups") +
  scale_color_brewer(type = "qual", palette = "Paired", name = "Groups") +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.8, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"),
        legend.direction = "horizontal",
        legend.position = c(0.8, 0.15))


#ggsave(filename = "results/GBM_nestbuilding_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")




# Plot the GBM that show a group effect at q < 0.2 (For Nestbuilding, which are the same conditions in the study 1)
GBM_BH <- clr_GBMs[GBMs.glm[GBMs.glm$`TreatmentWater Pr(>|t|).BH` < 0.2,"feature"],]

GBM_BH %>%
  t() %>% 
  as.data.frame() %>%
  add_column(Group = metadata$Treatment,
             Nestbuilding = metadata$Nestbuilding, 
             ID = metadata$plate) %>% 
  pivot_longer(!c("Group", "Nestbuilding", "ID")) %>% 
  #filter(Nestbuilding == "LNB") %>% 
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
  
  scale_fill_brewer(type = "qual", palette = "Paired", name = "Groups") +
  scale_color_brewer(type = "qual", palette = "Paired", name = "Groups") +
  
  ylab("") + xlab("") + 
  
  labs(color = "Groups") +
  
  theme_bw() + 
  
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14),
        
        strip.text = element_text(size = 14, face = "bold"), 
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold"))


#ggsave(filename = "results/GBM_escitalopram_q02.jpg", plot = last_plot(), width = 20, height = 12, units = "in")





# Plot the GMM that show a group effect at q < 0.2
GMM_BH <- clr_GMMs[GMMs.glm[GMMs.glm$`TreatmentWater:NestbuildingNNB Pr(>|t|).BH` < 0.2,"feature"],]

GMM_BH %>%
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
        
        legend.key.size = unit(0.4, units = "in"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14, face = "bold")) 

# Saving the results
#ggsave(filename = "results/GMM_q02.jpg", plot = last_plot(), width = 12, height = 8, units = "in")



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

# Generating the web
web <- weaveWebFromTables(tableY, tableX, dictionary)


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

anansiLong$

anansiLong %>% 
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


ggsave(filename = "results/all_correlations_rho-values.jpg", height = 12, width = 6, units = "in")







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


