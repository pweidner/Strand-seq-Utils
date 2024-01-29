#Plotting script for mosaicatcher-pipeline v1.8.4 results
library(tidyverse)
library(ggpubr)
library("RColorBrewer")
library(ggdendro)
library("cowplot")
library("scales")
library(ggforce)
library(ggbeeswarm)
library(cluster)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Set working directory to data/$SAMPLE/
setwd("~/Downloads/")

# Set theme for plots
theme_set(theme(legend.position = "none",
                title = element_text(size = 15, hjust = 0.5), 
                axis.text.x = element_text(size = 15),
                axis.text.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 15, hjust = 0.5),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                panel.border = element_rect(colour = "black", fill = NA, size = 1)))

# SV colors and chroms!
sv_list <-
  c(
    "none",
    "del_h1",
    "del_h2",
    "del_hom",
    "dup_h1",
    "dup_h2",
    "dup_hom",
    "inv_h1",
    "inv_h2",
    "inv_hom",
    "idup_h1",
    "idup_h2",
    "complex"
  )

# SV type colors
colors <-
  structure(
    c(
      "grey",
      "#77AADD",
      "#4477AA",
      "#114477",
      "#CC99BB",
      "#AA4488",
      "#771155",
      "#DDDD77",
      "#AAAA44",
      "#777711",
      "#DDAA77",
      "#AA7744",
      "#774411"
    ),
    names = sv_list
  )

chroms <- paste0("chr", c(1:22,"X","Y"))

################## GRAB RUNTIME ENVIRONMENTAL PARAMETERS ###########################

run <- read_tsv("run_summary.txt")
ash_t <- run[133,1] %>% str_remove("ashleys_threshold: ") %>% as.double()
ash_t <- 0.5
##################ASHLEYS PROBABILITIES / PREDICTIONS ##############################

# Load the labels dataframe!
selection <- read_tsv("labels (3).tsv")

# Probability histogram
p1 <- ggplot(selection) +
  geom_vline(data=selection, aes(xintercept=mean(probability*100)), linetype="dashed", size = 1, color = "gray") +
  geom_vline(aes(xintercept=ash_t*100, color="brown1"), linetype="dashed", size = 1) +
  geom_bar(aes(x=probability*100), fill = "azure2", stat = "bin", position = "identity", bins = 48) +
  stat_bin(aes(x=probability*100, label=..count..), bins = 48, geom = "text", vjust = 1.5) +
  labs(x="Probability", y="Cells") +
  scale_y_continuous(expand = c(0, 0, 0, .5))

# Prediction barplot!
p2 <- ggplot(selection, aes(x=as.factor(prediction))) +
  geom_bar(fill = c("brown1", "darkolivegreen2")) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(mult = c(0,.2))) +
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "black") +
  labs(x="Prediction", y="")

# Arrange the plots!
ggarrange(p1, p2, widths = c(1.75, 0.25), heights = 3)

# Save the pots!
ggsave(filename = "P1869_ASHLEYS_predictions.pdf", width = 10, height = 4)

# Clean up after yourself!
#remove(p1,p2,selection)

##########################################################################

####################### MERGED STATS #####################################

# Load merged stats summary df
stats <- read_tsv("stats/stats-merged.tsv")

# Add more verbose labels for plotting
stats$callset <- as.factor(c("lenient", "stringent"))

# Wrangle the data frame
stats$unique_calls_merged <- as.double(stats$unique_calls_merged)

# Replace all N/A with 0
stats[is.na(stats)] <- 0

# Plot general stats using wrapper!
stats %>%
  select(callset, "Cell Count" = cell_count, "Segments" = segments, "Total calls" = total_calls, "Unique calls" = unique_calls, "Mean SV/cell/mb" = avg_sv_load_per_cell_mb, "Mean SV/cell/mb complex" = avg_sv_load_per_cell_complex_mb, "Mean SCE/cell" = avg_sce_per_cell, "Total SCE" = total_sce) %>% 
  pivot_longer(!callset, names_to = "stats", values_to = "value") %>% 
  ggplot(aes(as.factor(callset), value)) +
  geom_col(width = 0.75, fill = "#92C5DE") +
  facet_wrap(~stats, scales = "free_y", nrow = 2) +
  labs(title = "General stats", x="", y="") +
  theme(axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        plot.title = element_text(size = 15, hjust = 0.5),
        strip.text = element_text(size=15))

# Save general stats plot!
ggsave(filename = "stats_general.pdf", width = 20, height = 13)

# Plot 
stats %>%
  select(callset, "Cell Count" = cell_count, "Segments" = segments, "Total calls" = total_calls, "Unique calls" = unique_calls, "Mean SV/cell/mb" = avg_sv_load_per_cell_mb, "Mean SV/cell/mb complex" = avg_sv_load_per_cell_complex_mb, "Mean SCE/cell" = avg_sce_per_cell, "Total SCE" = total_sce) %>% 
  pivot_longer(!callset, names_to = "stats", values_to = "value") %>%
  ggplot(aes(as.factor(stats))) +
  geom_bar(width = 1) +
  coord_polar(theta = "y", start=0) +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="class", 
       x=NULL, 
       y=NULL, 
       title="Pie Chart of class", 
       caption="Source: mpg")

 
stats %>%
  select(callset, calls_af0to10, calls_af0to10_complex, calls_af10to80, calls_af10to80_complex, calls_af80to100, calls_af80to100_complex) %>% 
  pivot_longer(!callset, names_to = "stats", values_to = "value") %>% 
  ggplot(aes(as.factor(callset), value)) +
  geom_col(width = 0.75, fill = "#92C5DE") +
  theme(legend.position = "none",
        title = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  facet_wrap(~stats,nrow = 1) +
  labs(title = "Calls by AF", x="", y="") 

ggsave(filename = "stats_alls_by_AF.pdf", width = 20, height = 8)        
remove(stats)

#################### MOSAICLASSIFIER ##########################

mosai <- read_tsv("HD2_50kb_v210_stringent_filterTRUE.tsv")

# Replot SV mutational landscape using ggplot to make it more accessable/tweakable!
mosai = mosai %>% 
  mutate(seg = paste(chrom,start,end,sep="_"), cell = gsub("P3070_i", "", cell))

# create matrix of samples x segments
mat = matrix(ncol = length(unique(mosai$seg)),nrow = length(unique(mosai$cell)))
colnames(mat) = sort(unique(mosai$seg))
rownames(mat) = unique(mosai$cell)
# fill matrix with sv call
for(i in 1:ncol(mat)){
  myseg = colnames(mat)[i]
  for(j in 1:nrow(mat)){
    mycell = rownames(mat)[j]
    mycall = NA
    mycall = mosai$sv_call_name[mosai$cell==mycell & mosai$seg == myseg]
    mat[j,i] = ifelse(is_empty(mycall),"NA",mycall)
  }
}

# Calculate dissimilarity matrix using Gower of the cluster package!
gower.dist = mat %>% 
  as.data.frame() %>% 
  mutate_all(as.factor) %>% 
  daisy(metric = c("gower"))

# Cluster using categorical clustering!
divisive.clust <- diana(as.matrix(gower.dist), 
                        diss = T,
                        keep.diss = TRUE,
                        stop.at.k = 5)

# Set cell order for plotting!
cellorder = divisive.clust$order.lab
plot(divisive.clust, main = "Divisive", xax.pretty = TRUE)

#Test medicc output
test <- read_tsv("~/Downloads/P1530_test2_medicc_input_pairwise_distances.tsv")
model <- hclust(dist(test), "ave")
model$labels <- test$sample_id
dhc <- as.dendrogram(model)

data <- dendro_data(dhc, type = "rectangle")
(tree <- ggplot(segment(data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  coord_flip() + 
  scale_y_reverse(expand = c(0, 0.2)) +
  theme_no_axes() +
    theme(panel.border = element_blank()))
p <- plot(hc, hang = -1)

(den <- ggdendrogram(hc, rotate = FALSE, size = 2, label = T) +
    coord_flip() +
    scale_y_reverse())

# Plot SV clustered heatmap
ggplot(mosai %>% mutate(cell = gsub("P1530_i","", cell), chrom = gsub("chr","", chrom)), aes(as.factor(end), factor(cell, levels = cellorder), fill = sv_call_name)) +
  geom_tile(width = 1) +
  facet_row(~as.double(sort(chrom)), scales = "free_x", space = "free",) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = colors) +
  labs(title = "Somatic mutational landscape", x="Genomic Position", y="Cell (n=28)") +
  theme_cowplot() +
  theme_pubclean(base_size = 25) +
  theme(axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        strip.text.x = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
  
n_distinct(mosai$cell)
# Plot clustered heatmap of sv calls!
mos <-  ggplot(mosai, aes(as.factor(end), factor(cell, levels = cellorder), fill = sv_call_name)) +
  geom_tile(width = 1) +
  facet_row(~as.double(sort(chrom)), scales = "free_x", space = "free",) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = colors) +
  labs(title = "Somatic mutational landscape", x="Genomic Position", y="Cell (n=45)") +
  theme_cowplot() +
  theme(axis.ticks.y = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggarrange(tree,mos, widths = c(1,10), heights = c(1,2))

extra <- ggplot(mosai, aes(x=mosai$start,y = mosai$af)) +
  geom_point(alpha = 0) +
  geom_density_2d() +
  labs(title = "",x="",y="") +
  theme()
ggExtra::ggMarginal(extra, type = "densigram")
######################## SV PER CELL BOXPLOT #############################
df <- mosai %>%
  mutate(cell, sv_call_name = gsub("_h1|_h2","_het", sv_call_name)) %>% 
  group_by(cell) %>% 
  count(sv_call_name)
  
ggplot(df, aes(as.factor(sv_call_name), n, fill = sv_call_name)) +
  geom_boxplot(width = 0.66) +
  scale_fill_manual(values = colors) +
  theme(legend.position = "none",
        title = element_text(size = 15), 
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  labs(title = "Average SV Calls per cell", x="", y="")
#ggsave("Avg_SV_per_cell.pdf", height = 4, width =5)

################################ SAMPLE INFO #################################

# Fetch info table 
info <- read.table("OP15.info", header = T)
info2 <- read.table("AGLCD.info", header = T)
# Plot stuff!
p1 <- info %>%
  select(cell, good) %>% 
  pivot_longer(!cell, names_to = "stats", values_to = "value") %>% 
  ggplot(aes(stats,value/1000)) +
  ylim(200, 1200) +
  labs(y="Reads per cell [k]", title = "OP Biomek (P3040)") +
  geom_boxplot(aes(fill = stats)) +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(alpha = 0.6, size = 3)+
  theme(axis.title.y = element_text(angle=90))

# Plot stuff!
p2 <- info2 %>%
  select(cell, good) %>% 
  pivot_longer(!cell, names_to = "stats", values_to = "value") %>% 
  ggplot(aes(stats,value/1000)) +
  ylim(200, 1200) +
  labs(y="Reads per cell [k]", title = "Conventional") +
  geom_boxplot(aes(fill = stats)) +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(alpha = 0.6, size = 3) +
  theme(axis.title.y = element_text(angle = 90))

ggarrange(p1, p2) + theme_pubclean()
################################### PLOIDY ###################################

df <- read_table(file = "Downloads/ploidy_detailled.txt")
colnames(df) <- c("chrom", "start", "end", "logLH-ploidy-1", "logLH-ploidy-2", "3", "4", "5", "6", "ploidy_estimate")
df$chrom <- as.factor(df$chrom)

ggplot(df, aes(start, ploidy_estimate)) +
  geom_col(aes(fill = ifelse(ploidy_estimate <= 2, 'red', "black"))) +
  labs(title="",y="") +
  scale_x_continuous() +
  facet_wrap(~factor(chrom, levels = chroms), ncol = 1) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        title = element_text(size = 15), 
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank())

################################ SAMPLE INFO #################################

# Fetch info table 
info <- read.table("~/Downloads/OP_robot.info", header = T)

# Plot stuff!
info %>%
  select(cell, good) %>% 
  pivot_longer(!cell, names_to = "stats", values_to = "value") %>% 
  ggplot(aes(stats,value/1E6)) +
  geom_boxplot(aes(fill = stats)) +
  scale_fill_brewer(palette = "Set2") +
  geom_jitter(alpha = 0.5)


############
df <- read_table("/Users/pweidne/Documents/Data/Sequencing/AGLCD/MiXCR/alignment.txt")

ggplot(df, aes(fill=condition, y=factor(cell), x=value)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_manual(values = c("gray", "sandybrown")) +
  labs(y="Cell", x = " Successfully aligned reads") +
  theme_pubclean(base_size = 40) +
  theme(axis.text.y = element_text(size = 15))
