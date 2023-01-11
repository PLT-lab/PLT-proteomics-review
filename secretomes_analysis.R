
# Load libraries ----------------------------------------------------------

library(conflicted)
library(readr)
library(readxl)
library(stringr)
library(dplyr)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# Load data ---------------------------------------------------------------

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

all.dfs <- read_excel_allsheets("datasets_secretoma.xlsx", tibble = TRUE)

# Prep all datasets

# s2 ----------------------------------------------------------------------

library(ensembldb)
library(EnsDb.Hsapiens.v86)
edb <- EnsDb.Hsapiens.v86

# Filter those with confident = yes

s2 <- all.dfs$s2 %>% 
  filter(`Confidently identified?             (confident when a total of 2 or more peptides were detected)` == "Yes")


df.s2 <- AnnotationDbi::select(edb, keys = unique(s2$Accession), columns = c("UNIPROTID", "SYMBOL", "ENTREZID", "GENEID"), keytype = "PROTEINID") %>% 
  distinct(UNIPROTID, .keep_all = TRUE) %>% 
  rename(UNIPROT = UNIPROTID)


# s3 ----------------------------------------------------------------------

s3 <- all.dfs$s3 %>% 
  select(`Uniprot Accession`,
         Gene) %>% 
  rename(UNIPROT = `Uniprot Accession`,
         SYMBOL = Gene)

df.s3 <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(all.dfs$s3$`Uniprot Accession`), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(s3, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# s4 ----------------------------------------------------------------------

# Mouse

df.s4.m <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(all.dfs$s4$`T: Majority protein IDs`), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# s5 ----------------------------------------------------------------------

s5 <- all.dfs$s5 %>% 
  select(`Gene Names`,
         `Protein IDs`) %>% 
  rename(UNIPROT = `Protein IDs`,
         SYMBOL = `Gene Names`) %>% 
  mutate(UNIPROT = sub("-.*", "", UNIPROT),
         SYMBOL = na_if(SYMBOL, "N/A"))

df.s5 <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(s5$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(s5, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# s6 ----------------------------------------------------------------------

any(is.na(all.dfs$s6$`Protein ID`))
any(is.na(all.dfs$s6$`Gene name`))

s6 <- all.dfs$s6 %>% 
  select(`Protein ID`,
         `Gene name`,
         starts_with("Control_")) %>% 
  na_if(., "NaN") %>% 
  type_convert() %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 10) %>% 
  rename(UNIPROT = `Protein ID`,
         SYMBOL = `Gene name`) %>% 
  mutate(UNIPROT = sub("-.*", "", UNIPROT)) %>% 
  select(UNIPROT,
         SYMBOL)


df.s6 <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(s6$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(s6, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# Whole proteomes ---------------------------------------------------------

# Load reference proteomes - mouse (17127) and human (20385 proteins)
fasta_headers_human <- sub(".*\\|(.*)\\|.*", "\\1", 
                           readLines(paste0(path,"fasta_human_proteome_20221124.fasta"))[grepl("^>", 
                                                                                         readLines(paste0(path,"fasta_human_proteome_20221124.fasta")))])
symbol_human <- sub(".*GN=(.*) PE.*", "\\1", 
                    readLines(paste0(path,"fasta_human_proteome_20221124.fasta"))[grepl("^>", 
                                                                                  readLines(paste0(path,"fasta_human_proteome_20221124.fasta")))])

fasta_headers_mouse <- sub(".*\\|(.*)\\|.*", "\\1", 
                           readLines(paste0(path,"fasta_mouse_proteome_20221124.fasta"))[grepl("^>", 
                                                                                         readLines(paste0(path,"fasta_mouse_proteome_20221124.fasta")))])


symbol_mouse <- sub(".*GN=(.*) PE.*", "\\1", 
                    readLines(paste0(path,"fasta_mouse_proteome_20221124.fasta"))[grepl("^>", 
                                                                                  readLines(paste0(path,"fasta_mouse_proteome_20221124.fasta")))])


# Pass all dfs to list ----------------------------------------------------

library(ggplot2)
library(gprofiler2)

# Human

# Create list with all (human) dfs
l.df.all <- lapply(ls(pattern = "df.s[0-9]$"), function(x) get(x))

result.all <- lapply(l.df.all, "[", , "UNIPROT")

# All proteins detected in all datasets
inter.names.all <- data.frame(nam = Reduce(c, result.all))

# Count occurrence of proteins
res.all <- inter.names.all %>% 
  na.omit() %>% 
  group_by(nam) %>% 
  count(nam) %>%
  filter(nam %in% fasta_headers_human) %>% # Only UNIPROT reviewed proteins
  arrange(desc(n))

# Barplot Supp. Fig. 1B
res.all %>% 
  mutate(color = ifelse(n < 2, "out", "in")) %>%
  ggplot(aes(as.factor(n))) +
  geom_bar(aes(fill = color), color = "black") +
  labs(x = "Number of datasets", y = "Count", title = "") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -.5, colour = "black", size = 3) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_y_continuous(limits = c(0,390), expand = c(0, 0))

# Filter those proteins that appear twice or more

prots.v.all <- res.all %>% 
  filter(n >= 2) # 562



# Annotate
prots.all <- AnnotationDbi::select(org.Hs.eg.db, keys = prots.v.all$nam, columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# Rescue proteins
prot.left.all <- gconvert(query = prots.all$UNIPROT[which(is.na(prots.all$ENSEMBL))], organism = "hsapiens", 
                      target="ENSG", mthreshold = Inf, filter_na = TRUE) %>% 
  rename(ENSEMBL = target,
         UNIPROT = input,
         SYMBOL = name) %>% 
  select(ENSEMBL,
         UNIPROT,
         SYMBOL)

prots.c.all <- prots.all %>% 
  full_join(prot.left.all, by = "UNIPROT") %>% 
  mutate(ENSEMBL = coalesce(!!! select(., matches("ENSEMBL"))),
         SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y"))


# Venn diagrams with whole platelet proteome ------------------------------

# Import PLT human core proteome
library(readr)
whole_PLT_prot_hsa <- read_delim("whole_PLT_prot_hsa.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)


# Overlap between PLT human core proteome and the human secretome
library(eulerr)

dfl.all <- list("Human platelet proteome" = whole_PLT_prot_hsa$UNIPROT, "Platelet secretome" = prots.c.all$UNIPROT)

colors <- c(RColorBrewer::brewer.pal(3, "Pastel1"))

# Fig. 4A
set.seed(3213)
plot(euler(dfl.all),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


# Check overlap with all detected in PLTs, not only those selected

# all <- read_delim("all_prot_detected_whole.txt", 
#                                              delim = "\t", escape_double = FALSE, 
#                                              trim_ws = TRUE)
# 
# table(prots.c.all$UNIPROT %in% all$nam) # Only 2 are not found in any PLT whole proteome dataset
# 
# prots.c.all$SYMBOL[-which(prots.c.all$UNIPROT %in% all$nam)] # IGHV2-70D and CYP7A1
# 
# venn.sec.all <- list("Human platelet proteome" = all$nam, "Platelet secretome" = prots.c.all$UNIPROT)
# 
# set.seed(3213)
# plot(euler(venn.sec.all),
#      quantities = TRUE,
#      fills = list(fill = colors),
#      labels = NULL,
#      legend = list(cex = 1))



# Orthologs ---------------------------------------------------------------

# Human to mouse
ort.human.mouse.all <- gorth(query = prots.all$UNIPROT, source_organism = "hsapiens", 
                         target_organism = "mmusculus", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE)

# Fig. 4B
set.seed(3213)
plot(euler(list("Platelet secretome human" = unique(prots.c.all$SYMBOL), 
                "Platelet secretome mouse.ort" = unique(toupper(ort.human.mouse.all$ortholog_name)))),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))



# Distribution and enrichment - Human -------------------------------------

# Human s5

dist.s5 <- all.dfs$s5 %>% 
  select(`Gene Names`,
         `Protein IDs`,
         Intensity) %>% 
  rename(UNIPROT = `Protein IDs`,
         SYMBOL = `Gene Names`) %>% 
  mutate(UNIPROT = sub("-.*", "", UNIPROT),
         SYMBOL = na_if(SYMBOL, "N/A")) %>% 
  filter(UNIPROT %in% fasta_headers_human) %>% # Filter unreviewed proteins
  distinct(UNIPROT, .keep_all = TRUE) %>% 
  mutate(log2Int = log2(Intensity),
         rank = dense_rank(dplyr::desc(log2Int)))


# Sub-setting based on quantiles

q1.names.5 <- dist.s5 %>% 
  filter(log2Int < summary(log2Int)[2]) %>% 
  pull(UNIPROT)

q13.names.5 <- dist.s5 %>% 
  filter(between(log2Int, summary(log2Int)[2], summary(log2Int)[5])) %>% 
  pull(UNIPROT)

q3.names.5 <- dist.s5 %>% 
  filter(log2Int > summary(log2Int)[5]) %>%
  pull(UNIPROT)


# Enrichent 

# Functions
gos <- function(de, org = "Human"){
  
  # This function does the enrichment for the three GO terms - outputs a list
  
  if (org == "Human"){
    
    org <- org.Hs.eg.db
    
  } else {
    
    org <- org.Mm.eg.db
    
  }
  
  go.cc <- enrichGO(gene = de,
                    OrgDb         = org,
                    keyType       = 'UNIPROT',
                    ont           = "CC",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  
  
  go.bp <- enrichGO(gene = de,
                    OrgDb         = org,
                    keyType       = 'UNIPROT',
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  
  go.mf <- enrichGO(gene = de,
                    OrgDb         = org,
                    keyType       = 'UNIPROT',
                    ont           = "MF",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
  
  go.list <- list(go.cc = go.cc,
                  go.bp = go.bp,
                  go.mf = go.mf)
  
  return(go.list)
  
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# Fig. 5 (Human)

# First quantile
q1.go.5 <- gos(q1.names.5)
q1.go.5.bp <- simplify(q1.go.5$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q1.go.5.bp@result$Description <- firstup(q1.go.5.bp@result$Description)
dotplot(q1.go.5.bp)

# Inter-quantile
q13.go.5 <- gos(q13.names.5)
q13.go.5.bp <- simplify(q13.go.5$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q13.go.5.bp@result$Description <- firstup(q13.go.5.bp@result$Description)
dotplot(q13.go.5.bp)

# Third quantile
q3.go.5 <- gos(q3.names.5)
q3.go.5.bp <- simplify(q3.go.5$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q3.go.5.bp@result$Description <- firstup(q3.go.5.bp@result$Description)
dotplot(q3.go.5.bp)


# Selected proteins to show on the distribution: 

# Filter terms related to PLTs, in each quantile
sel <- q3.go.5.bp@result %>% 
  dplyr::slice(1:10) %>% 
  filter(grepl("Platelet|Coagulation|hemostasis|blood|coagulation", Description)) %>% 
  pull(geneID)

sel2 <- q13.go.5.bp@result %>% 
  dplyr::slice(1:10) %>% 
  filter(grepl("Platelet|Complement", Description)) %>% 
  pull(geneID)

sel3 <- q1.go.5.bp@result %>% 
  dplyr::slice(1:10) %>% 
  filter(grepl("Complement", Description)) %>% 
  pull(geneID)

sels <- c(unlist(strsplit(sel, "/")), unlist(strsplit(sel2, "/")), unlist(strsplit(sel3, "/")))

sel.f <- data.frame(UNIPROT = sels,
                  SYMBOL = dist.s5$SYMBOL[match(sels, dist.s5$UNIPROT)]) %>% 
  distinct(UNIPROT, .keep_all = TRUE)

# From those, the representative proteins are:
sel.prots <- c("VWF", "EGF", "VEGFC", "PDGFA", "SELP", "PF4", "F5", "FERMT3", "TGFB1", "SERPINA4", "PECAM1", "FGB")


# Vertical lines marking the change in quantile

v1 <- dist.s5 %>% 
  filter(UNIPROT %in% q1.names.5) %>% 
  filter(rank == min(rank)) %>% 
  pull(rank)

v2 <- dist.s5 %>% 
  filter(UNIPROT %in% q13.names.5) %>% 
  filter(rank == min(rank)) %>% 
  pull(rank)

vlines_5 <- c(v1, v2)

# Plotting ranked distribution (relative quantification) - Fig. 5 (Human)
dist.s5 %>% 
  ggplot(aes(x = rank, y = log2Int)) +
  geom_vline(xintercept = vlines_5, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.s5$rank <= 5) | (dist.s5$rank >= max(dist.s5$rank) - 5) | (dist.s5$SYMBOL %in% sel.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-MS/MS count") +
  ggrepel::geom_label_repel(data = dist.s5 %>% filter((rank <= 5) | (rank >= max(rank) - 4) | (SYMBOL %in% sel.prots)),
                            aes(label = SYMBOL), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1, label="1st quantile\n511", y = 27), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2, label="3rd quantile\n171", y = 27), colour="darkgrey", angle = 90, size = 3)




# MOUSE -------------------------------------------------------------------

# Rescue proteins
prot.left.m <- gconvert(query = df.s4.m$UNIPROT[which(is.na(df.s4.m$ENSEMBL))], organism = "mmusculus", 
                          target="ENSG", mthreshold = Inf, filter_na = TRUE) %>% 
  rename(ENSEMBL = target,
         UNIPROT = input,
         SYMBOL = name) %>% 
  select(ENSEMBL,
         UNIPROT,
         SYMBOL)

prots.c.m <- df.s4.m %>% 
  full_join(prot.left.m, by = "UNIPROT") %>% 
  mutate(ENSEMBL = coalesce(!!! select(., matches("ENSEMBL"))),
         SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y"))

# Import mouse PLT proteome
whole_PLT_prot_mmus <- read_delim("whole_PLT_prot_mmus.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)

# Overlap between secretome and proteome
dfl.m <- list("Mouse platelet proteome" = whole_PLT_prot_mmus$UNIPROT, "Mouse platelet secretome (Thrombin)" = prots.c.m$UNIPROT)

# Fig. 4A
set.seed(3213)
plot(euler(dfl.m),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))

# Orthologs
ort.mouse.human <- gorth(query = prots.c.m$UNIPROT, source_organism = "mmusculus", 
                             target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE)

# Fig. 4B
set.seed(3213)
plot(euler(list("Platelet secretome mouse" = unique(toupper(prots.c.m$SYMBOL)), 
                "Platelet secretome human.ort" = unique(ort.mouse.human$ortholog_name))),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


# Distribution

# Mouse s4 dataset
# Clean, filter, calculate median protein expression and rank proteins
dist.s4 <- all.dfs$s4 %>% 
  select(`T: Majority protein IDs`,
         `N: MS/MS count`) %>%
  rename(UNIPROT = `T: Majority protein IDs`) %>% 
  left_join(df.s4.m, by = "UNIPROT") %>% 
  distinct(SYMBOL, .keep_all = TRUE) %>% 
  mutate(UNIPROT = sub("-.*", "", UNIPROT),
         log2MS = log2(`N: MS/MS count`),
         rank = dense_rank(dplyr::desc(log2MS))) %>% 
  dplyr::arrange(desc(log2MS)) %>% 
  mutate(rank2 = c(1:350))


# Sub-setting based on quantiles

q1.names.4 <- dist.s4 %>% 
  mutate(UNIPROT = gsub("-.*", "", UNIPROT)) %>%
  filter(log2MS < summary(log2MS)[2]) %>% 
  pull(UNIPROT)

q13.names.4 <- dist.s4 %>% 
  mutate(UNIPROT = gsub("-.*", "", UNIPROT)) %>%
  filter(between(log2MS, summary(log2MS)[2], summary(log2MS)[5])) %>% 
  pull(UNIPROT)

q3.names.4 <- dist.s4 %>% 
  mutate(UNIPROT = gsub("-.*", "", UNIPROT)) %>%
  filter(log2MS > summary(log2MS)[5]) %>%
  pull(UNIPROT)

# Enrichment

# Fig. 5 (Mouse)

# First quantile
q1.go.4 <- gos(q1.names.4, org = "Mouse")
q1.go.4.bp <- simplify(q1.go.4$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q1.go.4.bp@result$Description <- firstup(q1.go.4.bp@result$Description)
dotplot(q1.go.4.bp)

# Inter-quantile
q13.go.4 <- gos(q13.names.4, org = "Mouse")
q13.go.4.bp <- simplify(q13.go.4$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q13.go.4.bp@result$Description <- firstup(q13.go.4.bp@result$Description)
dotplot(q13.go.4.bp)

# Third quantile
q3.go.4 <- gos(q3.names.4, org = "Mouse")
q3.go.4.bp <- simplify(q3.go.4$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q3.go.4.bp@result$Description <- firstup(q3.go.4.bp@result$Description)
dotplot(q3.go.4.bp)


# Selected proteins

selm <- q3.go.4.bp@result %>% 
  dplyr::slice(1:10) %>% 
  filter(grepl("Platelet|Coagulation|Hemostasis|blood|coagulation", Description)) %>% 
  pull(geneID)

selm2 <- q13.go.4.bp@result %>% 
  dplyr::slice(1:10) %>% 
  filter(grepl("Platelet|Complement|Coagulation|Hemostasis|blood|coagulation", Description)) %>% 
  pull(geneID)

selms <- c(unlist(strsplit(selm, "/")), unlist(strsplit(selm2, "/")))

selm.f <- data.frame(UNIPROT = selms,
                    SYMBOL = dist.s4$SYMBOL[match(selms, dist.s4$UNIPROT)]) %>% 
  distinct(UNIPROT, .keep_all = TRUE)

selm.prots <- c("Fgb", "Fga", "Vwf", "F5", "Pf4", "Fermt3", "Selp", "Pdgfb", "Syk", "Gp6", "Serpina10")

# Vertical lines marking the change in quantile

v1m <- dist.s4 %>% 
  filter(UNIPROT %in% q1.names.4) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2m <- dist.s4 %>% 
  filter(UNIPROT %in% q13.names.4) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_4 <- c(v1m, v2m)

# Plotting ranked distribution (relative quantification) - Fig. 5 (Mouse)
dist.s4 %>% 
  ggplot(aes(x = rank2, y = log2MS)) +
  geom_vline(xintercept = vlines_4, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.s4$rank2 <= 5) | (dist.s4$rank2 >= max(dist.s4$rank2) - 5) | (dist.s4$SYMBOL %in% selm.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-MS/MS count") +
  ggrepel::geom_label_repel(data = dist.s4 %>% filter((rank2 <= 5) | (rank2 >= max(rank2) - 4) | (SYMBOL %in% selm.prots)),
                            aes(label = SYMBOL), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1m, label="1st quantile\n264", y = 2), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2m, label="3rd quantile\n89", y = 2), colour="darkgrey", angle = 90, size = 3)





# Output tables -----------------------------------------------------------

# output_tables <- list(all_human_prots = res.all,
#                       selected_human_prots.all = prots.c.all,
#                       selected_mouse_prots = prots.c.m)

# library(openxlsx)
# write.xlsx(output_tables, file = "PLT_secretome_output_tables.xlsx")
