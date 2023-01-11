
rm(list = ls())

# https://biologydirect.biomedcentral.com/articles/10.1186/1745-6150-9-3 says 87% coincidence and ~0.3 with proteome

# Load libraries ----------------------------------------------------------

library(readr)
library(dplyr)
library(eulerr)
library(readxl)
library(ggpubr)

# Import ------------------------------------------------------------------

# Data from whole_prot_analysis.R
dist.11 <- read_delim("dist11.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
dist.39 <- read_delim("dist39.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

prots.c <- read_delim("whole_PLT_prot_hsa.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
prots.cm <- read_delim("whole_PLT_prot_mmus.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)


# RNA Weyrich paper -------------------------------------------------------

# https://doi.org/10.1182/blood-2011-03-339705

# https://ashpublications.org/blood/article/121/26/5257/31555/Response-platelet-transcriptome-and-proteome

# Files from the Supplementary

rna2.h <- read_excel("zhe101-sup-tables6.xls") %>% # Human
  filter(RPKM > .3)

rna2.m <- read_excel("zhe101-sup-tables8.xls") %>% # Mouse
  filter(RPKM > .3)


table(rna2.h$Name %in% prots.c$SYMBOL)
table(rna2.m$Name %in% prots.cm$SYMBOL)

# Venn diagrams // Fig. 3A

colors <- c(RColorBrewer::brewer.pal(3, "Pastel1"))

# Overlap between proteins from the human PLT core proteome and human PLT transcripts
rna.prot.venn.human <- list("Human platelet RNA" = rna2.h$Name, "Human platelet protein" = unique(prots.c$SYMBOL))

set.seed(3213)
plot(euler(rna.prot.venn.human),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


# Overlap between proteins from the mouse PLT core proteome and mouse PLT transcripts
rna.prot.venn.mouse <- list("Mouse platelet RNA" = rna2.m$Name, "Mouse platelet protein" = unique(prots.cm$SYMBOL))

set.seed(3213)
plot(euler(rna.prot.venn.mouse),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


# Correlation plots // Fig. 3B

# Human
table(dist.11$`Gene names` %in% rna2.h$Name)

to.plot2 <- dist.11 %>% 
  filter(`Gene names` %in% rna2.h$Name) %>% # Find common genes between proteomics and transcriptomics
  select(`Gene names`,
         median_LFQ) %>% 
  rename(Name = `Gene names`) %>% 
  left_join(rna2.h, by = "Name") %>% # Get all relative abundances
  select(Name,
         median_LFQ,
         RPKM) %>% 
  mutate(log2.RPKM = log2(RPKM)) # Transform to log2 for plotting

# Plot correlation plot with R and p-value
ggscatter(to.plot2, x = "median_LFQ", y = "log2.RPKM",
          add = "reg.line",  # Add regression line
          add.params = list(color = "darkolivegreen", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) + 
  stat_cor(method = "pearson", label.x = 23, label.y = 14) +
  ylab("log2-RPKM (RNA)") +
  xlab("Median LFQ (Protein)") +
  ggtitle("Human")


# Mouse
to.plot3 <- dist.39 %>% 
  filter(SYMBOL %in% rna2.m$Name) %>% # Find common genes between proteomics and transcriptomics
  select(SYMBOL,
         median_LFQ) %>% 
  rename(Name = SYMBOL) %>% 
  left_join(rna2.m, by = "Name") %>% # Get all relative abundances
  select(Name,
         median_LFQ,
         RPKM) %>% 
  mutate(log2.RPKM = log2(RPKM)) # Transform to log2 for plotting

# Plot correlation plot with R and p-value
ggscatter(to.plot3, x = "median_LFQ", y = "log2.RPKM",
          add = "reg.line",  # Add regression line
          add.params = list(color = "darkolivegreen", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 9, label.y = 14) +
  ylab("log2-RPKM (RNA)") +
  xlab("Median LFQ (Protein)") +
  ggtitle("Mouse")
