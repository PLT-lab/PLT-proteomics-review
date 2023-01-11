
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

all.dfs <- read_excel_allsheets("all_datasets.xlsx", tibble = TRUE)

# Prep all datasets

# x1 ----------------------------------------------------------------------

all.dfs$x1 <- all.dfs$x1 %>% 
  select(Description,
         Gene_name) %>%  # 251 obs
  mutate(Gene_name = sub('\\s*;.*','', Gene_name))

any(is.na(all.dfs$x1$Description))
any(is.na(all.dfs$x1$Gene_name))

all.dfs$x1$Description[which(is.na(all.dfs$x1$Gene_name))] # Igs

eg.x1.df <- AnnotationDbi::select(org.Hs.eg.db, keys = all.dfs$x1$Gene_name, columns = c("UNIPROT", "ENSEMBL", "ENTREZID"), keytype = "SYMBOL") %>% 
  distinct(SYMBOL, .keep_all = TRUE)


# x2 ----------------------------------------------------------------------

any(is.na(all.dfs$x2$Accession))

# Remove prots not detected
id_rmv_x2 <- all.dfs$x2 %>% 
  select(Accession,
         starts_with("Abundances")) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(Count_NA > 40) %>%
  pull(Accession)

all.dfs$x2 <- all.dfs$x2 %>% 
  filter(!Accession %in% id_rmv_x2) %>% 
  select(Accession,
         `Ensembl Gene ID`,
         `Entrez Gene ID`,
         `Gene Symbol`)

eg.x2.df <- all.dfs$x2 %>% 
  rename(UNIPROT = Accession,
         ENSEMBL = `Ensembl Gene ID`,
         ENTREZID = `Entrez Gene ID`,
         SYMBOL = `Gene Symbol`) %>% 
  mutate(ENSEMBL = sub('\\s*;.*','', ENSEMBL),
         ENSEMBL = sub('\\s*\\..*','', ENSEMBL)) %>% 
  as.data.frame()

# x4 ----------------------------------------------------------------------

eg.x4 <- all.dfs$x4 %>% 
  mutate(Locus = gsub("nxp:NX_", "", Locus),
         Locus = sub("-.*", "", Locus)) %>% 
  filter(!grepl('contaminant', Locus)) %>% 
  pull(Locus)


eg.x4.df <- AnnotationDbi::select(org.Hs.eg.db, keys = eg.x4, columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x6 ----------------------------------------------------------------------

eg.x6 <- all.dfs$x6 %>% 
  select(`Accession...4`) %>% 
  dplyr::rename(Accession = `Accession...4`) %>% 
  mutate(Accession = word(Accession, -2, sep="\\|"))

eg.x6.df <- AnnotationDbi::select(org.Hs.eg.db, keys = eg.x6$Accession, columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE) 


# x7 ----------------------------------------------------------------------

all.dfs$x7 <- all.dfs$x7 %>% 
  filter(!`Potential contaminant` %in% "+",
         !`Only identified by site` %in% "+",
         !Reverse %in% "+") %>% 
  mutate(`Gene names` = sub('\\s*;.*','', `Gene names`),
         `Gene names` = sub("-.*", "", `Gene names`),
         `Protein IDs` = sub('\\s*;.*','', `Protein IDs`)) %>% 
  select(`Protein IDs`,
         `Gene names`,
         ends_with("Non-pregnant)")) # Took only the non-pregnant

x7.names <- all.dfs$x7 %>% 
  select_if(is.character) %>% 
  rename(UNIPROT = `Protein IDs`,
         SYMBOL = `Gene names`)

any(is.na(all.dfs$x7$`Protein IDs`))
any(is.na(all.dfs$x7$`Gene names`))

eg.x7 <- all.dfs$x7 %>% 
  select(-`Gene names`) %>% 
  na_if(0) %>%
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 3) %>% # Filtered out if they were all NAs
  mutate(`Protein IDs` = sub("-.*", "", `Protein IDs`)) %>% 
  pull(`Protein IDs`)

eg.x7.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x7), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(x7.names, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x8 ----------------------------------------------------------------------

eg.x8 <- all.dfs$x8 %>% 
  select(Protein) %>% 
  mutate(Protein = word(Protein, -2, sep="\\|"))

any(is.na(eg.x8$Protein))

eg.x8.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x8$Protein), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x9 ----------------------------------------------------------------------

x9.names <- all.dfs$x9 %>% 
  select(`Gene Name`,
         `Protein ID`) %>% 
  rename(UNIPROT = `Protein ID`,
         SYMBOL = `Gene Name`)

eg.x9.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(all.dfs$x9$`Protein ID`), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(x9.names, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x10 ---------------------------------------------------------------------

eg.x10.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(all.dfs$x10$`Protein accession`), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x11 ---------------------------------------------------------------------

all.dfs$x11 <- all.dfs$x11 %>% 
  select(`Protein IDs`,
         `Gene names`,
         starts_with("LFQ intensity C"),
         `Only identified by site`,
         Reverse,
         `Potential contaminant`)

eg.x11 <- all.dfs$x11 %>% 
  filter(!`Potential contaminant` %in% "+",
         !`Only identified by site` %in% "+",
         !Reverse %in% "+") %>% 
  select(-`Only identified by site`,
         -`Potential contaminant`,
         - Reverse) %>% 
  mutate(`Protein IDs` = sub('\\s*;.*','', `Protein IDs`),
         `Gene names` = sub('\\s*;.*','', `Gene names`)) %>% 
  na_if(0)

any(is.na(eg.x11$`Protein IDs`))
any(is.na(eg.x11$`Gene names`))

eg.x11.v <- eg.x11 %>% 
  select(- `Gene names`) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 6) %>% 
  left_join(eg.x11[,1:2], by = "Protein IDs") %>% 
  select(`Protein IDs`,
         `Gene names`) %>% 
  rename(UNIPROT = `Protein IDs`,
         SYMBOL = `Gene names`)


eg.x11.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x11.v$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x11.v, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)
  

# x12 ---------------------------------------------------------------------

any(is.na(all.dfs$x12$Accession))

eg.x12 <- all.dfs$x12 %>% 
  select(Accession,
         starts_with("C")) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 4) %>% 
  left_join(all.dfs$x12[,1:2], by = "Accession") %>% 
  rename(UNIPROT = Accession,
         SYMBOL = Symbol) %>% 
  select(UNIPROT,
         SYMBOL)

eg.x12.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x12$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x12, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x15 ---------------------------------------------------------------------

any(is.na(all.dfs$x15$Gene.names))

eg.x15 <- all.dfs$x15 %>% 
  filter(Quantified == "Y") %>% # 2871
  select(Gene.names,
         starts_with("CONTROL")) %>% 
  type_convert() %>% 
  mutate(Gene.names = sub('\\s*;.*','', Gene.names)) %>%
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 24) %>%
  pull(Gene.names)

eg.x15.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x15), columns = c("UNIPROT", "ENSEMBL", "ENTREZID"), keytype = "SYMBOL") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x17 ---------------------------------------------------------------------

eg.x17 <- all.dfs$x17 %>% 
  select(`Accession (UNIPROT)`,
         `Gene name`) %>% 
  mutate(`Accession (UNIPROT)` = word(`Accession (UNIPROT)`, -2, sep="\\|")) %>% 
  rename(UNIPROT = `Accession (UNIPROT)`,
         SYMBOL = `Gene name`)

any(is.na(eg.x17$UNIPROT))
any(is.na(eg.x17$SYMBOL))

eg.x17.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x17$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x17, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)

# x21 ---------------------------------------------------------------------

eg.x21 <- all.dfs$x21 %>% 
  filter(!grepl("REVERSED", Name)) %>% 
  select(`Accession #`,
         `Protein Name`) %>% 
  rename(UNIPROT = `Accession #`,
         SYMBOL = `Protein Name`)

eg.x21.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x21$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x21, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x24 ---------------------------------------------------------------------

# Mouse!

eg.x24 <- all.dfs$x24 %>% 
  mutate(SYMBOL = str_match(Reference, "GN=\\s*(.*?)\\s* PE=")[,2]) %>% 
  select(`Uniprot ID`,
         SYMBOL) %>% 
  rename(UNIPROT = `Uniprot ID`)

eg.mx24.df <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(eg.x24$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x24, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x25 ---------------------------------------------------------------------

# There seems to be no NAs in the intensity matrix, so no filtering for the controls is necessary

eg.x25 <- all.dfs$x25 %>% 
  select(`Gene Symbol`,
         `Majority protein IDs`) %>% 
  mutate(`Majority protein IDs` = sub('\\s*;.*','', `Majority protein IDs`),
         `Majority protein IDs` = sub("-.*", "", `Majority protein IDs`)) %>% 
  rename(UNIPROT = `Majority protein IDs`,
         SYMBOL = `Gene Symbol`)


eg.x25.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x25$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x25, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x26 ---------------------------------------------------------------------

# Don't really know which one are the controls, but it seems as if there were some kind of normalization in place
# so maybe the filtering based on NAs is enough

any(is.na(all.dfs$x26$Accession))

eg.x26 <- all.dfs$x26 %>% 
  select(Accession,
         Description,
         ends_with("Count")) %>% 
  mutate(Accession = sub("-.*", "", Accession),
         SYMBOL = str_match(Description, "GN=\\s*(.*?)\\s* PE=")[,2])

eg.x26.v <- eg.x26 %>% 
  select(- SYMBOL,
         - Description) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 8) 

eg.x26.v <- eg.x26.v %>% 
  select(Accession) %>% 
  left_join(eg.x26[,c(1,13)], by = "Accession") %>% 
  rename(UNIPROT = Accession) %>% 
  select(UNIPROT,
         SYMBOL)


eg.x26.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x26.v$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x26.v, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x28 ---------------------------------------------------------------------

# Doesn't seem to be NAs in the intensity matrix

eg.x28 <- all.dfs$x28 %>% 
  select(`Protein IDs`,
         `Gene names`) %>% 
  mutate(`Protein IDs` = sub('\\s*;.*','', `Protein IDs`),
         `Gene names` = sub('\\s*;.*','', `Gene names`),
         `Protein IDs` = sub("-.*", "", `Protein IDs`),
         `Gene names` = sub("-.*", "", `Gene names`)) %>% 
  rename(UNIPROT = `Protein IDs`,
         SYMBOL = `Gene names`)


eg.x28.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x28$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x28, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x29 ---------------------------------------------------------------------

eg.x29 <- all.dfs$x29 %>% 
  select(Entry,
         Description) %>% 
  mutate(Entry = gsub("nxp:NX_", "", Entry),
         Entry = sub("-.*", "", Entry),
         Description = gsub("[[:punct:]]", "", Description),
         SYMBOL = str_match(Description, "Gname\\s*(.*?)\\s* NcbiTaxId9606")[,2]) %>% 
  select(SYMBOL, 
         Entry) %>% 
  rename(UNIPROT = Entry)


eg.x29.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x29$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x29, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x30 ---------------------------------------------------------------------

eg.x30 <- all.dfs$x30 %>% 
  select(UniProtID,
         EntrezID,
         Symbol) %>% 
  rename(UNIPROT = UniProtID,
         SYMBOL = Symbol,
         ENTREZID = EntrezID) %>% 
  mutate(ENTREZID = as.character(ENTREZID))

eg.x30.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x30$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x30, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL"))),
         ENTREZID = coalesce(!!! select(., matches("ENTREZID")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x31 ---------------------------------------------------------------------

eg.x31 <- all.dfs$x31 %>% 
  select(`Uniprot Entry`,
         Genename) %>% 
  rename(UNIPROT = `Uniprot Entry`,
         SYMBOL = Genename)

eg.x31.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x31$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x31, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x33 ---------------------------------------------------------------------

# 114 is the control, here they calculated the ratios

eg.x33 <- all.dfs$x33 %>% 
  select(`Accession #`,
         Name.x) %>% 
  mutate(`Accession #` = word(`Accession #`, -2, sep="\\|"),
         SYMBOL = str_match(Name.x, "GN=\\s*(.*?)\\s* PE=")[,2]) %>% 
  rename(UNIPROT = `Accession #`) %>% 
  filter(!grepl("REVERSED", Name.x)) %>% 
  select(-Name.x)


eg.x33.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x33$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x33, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x35 ---------------------------------------------------------------------

# Mouse // interested in the ultra-purified fraction

eg.x35 <- all.dfs$x35 %>% 
  select(`Accession ID (UniProt)`,
         `Gene names`,
         `LFQ intensity ultra purified`) %>% 
  mutate(`Accession ID (UniProt)` = sub("-.*", "", `Accession ID (UniProt)`),
         `Gene names` = sub("-.*", "", `Gene names`)) %>% 
  rename(UNIPROT = `Accession ID (UniProt)`,
         SYMBOL = `Gene names`) %>% 
  mutate(SYMBOL = sub('\\s*;.*','', SYMBOL)) %>% 
  select_if(is.character)

eg.mx35.df <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(eg.x35$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x35, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x36 ---------------------------------------------------------------------

eg.x36 <- all.dfs$x36 %>% 
  select(`Protein accession`,
         Description) %>% 
  mutate(Description = word(Description, -2, sep="_")) %>% 
  rename(UNIPROT = `Protein accession`,
         SYMBOL = Description)

eg.x36.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x36$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x36, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x37 ---------------------------------------------------------------------

eg.x37 <- all.dfs$x37 %>% 
  select(`Accession No. (SwissProt)`) %>% 
  mutate(`Accession No. (SwissProt)` = sub("(^[^_]+)_.*", "\\1", `Accession No. (SwissProt)`))

eg.x37.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x37$`Accession No. (SwissProt)`), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x39 ---------------------------------------------------------------------

# Mouse

any(is.na(all.dfs$x39$gene_names))

eg.x39 <- all.dfs$x39 %>% 
  select(protein_i_ds,
         starts_with("1_"),
         starts_with("2_"),
         starts_with("3_")) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 5) %>% 
  left_join(all.dfs$x39[,1:2], by = "protein_i_ds") %>% 
  rename(UNIPROT = protein_i_ds,
         SYMBOL = gene_names) %>% 
  mutate(SYMBOL = sub('\\s*;.*','', SYMBOL)) %>% 
  select(UNIPROT,
         SYMBOL)


eg.mx39.df <- AnnotationDbi::select(org.Mm.eg.db, keys = unique(eg.x39$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x39, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x40 ---------------------------------------------------------------------

eg.x40 <- all.dfs$x40 %>% 
  select(Accession,
         Gene) %>% 
  rename(UNIPROT = Accession,
         SYMBOL= Gene)

eg.x40.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x40$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x40, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x41 ---------------------------------------------------------------------

any(is.na(all.dfs$x41$Accession))

eg.x41 <- all.dfs$x41 %>% 
  select(Accession,
         starts_with("Velos"),
         Description) %>% 
  na_if(0) %>% 
  select(!ends_with(" 1"),
         !ends_with(" 5"),
         !ends_with(" 9")) %>% # These are the patients
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 17) %>% 
  mutate(SYMBOL = gsub(".*\\((.*)\\).*", "\\1", Description),
         SYMBOL = gsub("_HUMAN", "", SYMBOL)) %>% 
  select(Accession,
         SYMBOL) %>% 
  rename(UNIPROT = Accession)


eg.x41.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x41$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x41, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x42 ---------------------------------------------------------------------

eg.x42 <- all.dfs$x42 %>% 
  select(Accession,
         `Gene name`,
         starts_with("APB")) %>% 
  mutate(Accession = sub("-.*", "", Accession)) %>% 
  na_if("NaN") %>% 
  type_convert() %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 5) %>% 
  select(Accession,
         `Gene name`) %>% 
  rename(UNIPROT = Accession,
         SYMBOL = `Gene name`)

eg.x42.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x42$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  full_join(eg.x42, by = "UNIPROT") %>% 
  mutate(SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# x43 ---------------------------------------------------------------------

eg.x43 <- read_excel("annotation_x43.xlsx") %>% 
  mutate(`Uniprot ID` = gsub(" ", ",", `Uniprot ID`),
         `Uniprot ID` = gsub("^(.*?),.*", "\\1", `Uniprot ID`),
         `Uniprot ID` = case_when(`Uniprot ID` == "No" ~ NA_character_,
                                  TRUE ~ `Uniprot ID`),
         `Uniprot ID` = gsub("_HUMAN", "", `Uniprot ID`)) %>% 
  select(-`INPUT ID`) %>% 
  rename(UNIPROT = `Uniprot ID`,
         SYMBOL = Symbol,
         ENSEMBL = `Ensembl ID`)

eg.x43.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x43$ENSEMBL), columns = c("UNIPROT", "SYMBOL", "ENTREZID"), keytype = "ENSEMBL") %>% 
  full_join(eg.x43, by = "ENSEMBL") %>% 
  mutate(UNIPROT = coalesce(!!! select(., matches("UNIPROT"))),
         SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y")) %>% 
  distinct(ENSEMBL, .keep_all = TRUE)


# x44 ---------------------------------------------------------------------

eg.x44 <- all.dfs$x44 %>% 
  select(`Uniprot Accession`,
         starts_with("Normalized")) %>% # It's not clear which one are the controls
  mutate(`Uniprot Accession` = gsub("UniProtKB:", "", `Uniprot Accession`)) %>% 
  rename(UNIPROT = `Uniprot Accession`) %>% 
  filter(!grepl("CON_", UNIPROT)) %>%  # Remove contaminants
  na_if(., 0) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 6) # To be on the safe side, filter out only those that are all NAs
  
eg.x44.df <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(eg.x44$UNIPROT), columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)



# Whole proteomes ---------------------------------------------------------

# Load reference proteomes - mouse (17127) and human (20385 proteins)

fasta_headers_human <- sub(".*\\|(.*)\\|.*", "\\1", 
                           readLines("fastas/fasta_human_proteome_20221124.fasta")[grepl("^>", 
                                                                                         readLines("fastas/fasta_human_proteome_20221124.fasta"))])
symbol_human <- sub(".*GN=(.*) PE.*", "\\1", 
                    readLines("fastas/fasta_human_proteome_20221124.fasta")[grepl("^>", 
                                                                                  readLines("fastas/fasta_human_proteome_20221124.fasta"))])

fasta_headers_mouse <- sub(".*\\|(.*)\\|.*", "\\1", 
                           readLines("fastas/fasta_mouse_proteome_20221124.fasta")[grepl("^>", 
                                                                                         readLines("fastas/fasta_mouse_proteome_20221124.fasta"))])


symbol_mouse <- sub(".*GN=(.*) PE.*", "\\1", 
                    readLines("fastas/fasta_mouse_proteome_20221124.fasta")[grepl("^>", 
                                                                                  readLines("fastas/fasta_mouse_proteome_20221124.fasta"))])

# Pass all dfs to list ----------------------------------------------------

library(ggplot2)
library(gprofiler2)

# Human

# Create list with all (human) dfs
l.df <- lapply(ls(pattern = "eg.x[0-9]+.df"), function(x) get(x))

result <- lapply(l.df, "[", , "UNIPROT")

# All proteins detected in all datasets
inter.names <- data.frame(nam = Reduce(c, result)) %>% 
  na.omit()

# Count occurrence of proteins
res <- inter.names %>% 
  group_by(nam) %>% 
  count(nam) %>% 
  filter(nam %in% fasta_headers_human) %>% # Only UNIPROT reviewed proteins
  arrange(desc(n))

# write.table(res, file = "all_prot_detected_whole.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Barplot Supp. Fig. 1A
res %>% 
  mutate(color = ifelse(n < 12, "out", "in")) %>% 
  ggplot(aes(as.factor(n))) +
  geom_bar(aes(fill = color), color = "black") +
  ggbreak::scale_y_break(c(1300, 3500)) +
  ggbreak::scale_y_break(c(570, 1100), expand = c(0, 0)) +
  labs(x = "Number of datasets", y = "Count", title = "") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -.5, colour = "black", size = 3) +
  ylim(0, 3550) +
  scale_fill_brewer(palette = "Pastel1")

# Filter those proteins that appear 12 times or more
prots.v <- res %>% 
  filter(n >= 12) # 2021

# Annotate
prots <- AnnotationDbi::select(org.Hs.eg.db, keys = prots.v$nam, columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


# Rescue proteins
prot.left <- gconvert(query = prots$UNIPROT[which(is.na(prots$ENSEMBL))], organism = "hsapiens", 
         target="ENSG", mthreshold = Inf, filter_na = TRUE) %>% 
  rename(ENSEMBL = target,
         UNIPROT = input,
         SYMBOL = name) %>% 
  select(ENSEMBL,
         UNIPROT,
         SYMBOL) 

prots.c <- prots %>% 
  full_join(prot.left, by = "UNIPROT") %>% 
  mutate(ENSEMBL = coalesce(!!! select(., matches("ENSEMBL"))),
         SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y"))

# Export
# write.table(prots.c, file = "whole_PLT_prot_hsa.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  

# Get mouse orthologs
ort.human.mouse <- gorth(query = prots.c$UNIPROT, source_organism = "hsapiens", 
      target_organism = "mmusculus", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE)




# Mouse

# Create list with all (mouse) dfs
l.dfm <- lapply(ls(pattern = "eg.mx[0-9]+.df"), function(x) get(x))

result.m <- lapply(l.dfm, "[", , "UNIPROT")

# All proteins detected in all datasets
inter.names.m <- data.frame(nam = Reduce(c, result.m))

# Count occurrence of proteins
resm <- inter.names.m %>% 
  na.omit() %>% 
  group_by(nam) %>% 
  count(nam) %>% 
  filter(nam %in% fasta_headers_mouse) %>% # Only UNIPROT reviewed proteins
  arrange(desc(n))

# Barplot Supp. Fig. 1B
resm %>% 
  mutate(color = ifelse(n < 2, "out", "in")) %>% 
  ggplot(aes(as.factor(n))) +
  geom_bar(aes(fill = color), color = "black") +
  labs(x = "Number of datasets", y = "Count", title = "") +
  theme_bw() +
  theme(legend.position = "none") +
  geom_text(aes(label = ..count..), stat = "count", vjust = -.5, colour = "black", size = 3) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_y_continuous(limits = c(0,2800), expand = c(0, 0))

# Filter those proteins that appear twice or more
prots.mv <- resm %>% 
  filter(n >= 2) # 1665 - 59 // 1606

# Annotate
prots.m <- AnnotationDbi::select(org.Mm.eg.db, keys = prots.mv$nam, columns = c("SYMBOL", "ENSEMBL", "ENTREZID"), keytype = "UNIPROT") %>% 
  distinct(UNIPROT, .keep_all = TRUE)


prot.left.m <- gconvert(query = prots.m$UNIPROT[which(is.na(prots.m$ENSEMBL))], organism = "mmusculus", 
                      target = "ENSG", mthreshold = Inf, filter_na = TRUE) %>% 
  rename(ENSEMBL = target,
         UNIPROT = input,
         SYMBOL = name) %>% 
  select(ENSEMBL,
         UNIPROT,
         SYMBOL)

prots.cm <- prots.m %>% 
  full_join(prot.left.m, by = "UNIPROT") %>% 
  mutate(ENSEMBL = coalesce(!!! select(., matches("ENSEMBL"))),
         SYMBOL = coalesce(!!! select(., matches("SYMBOL")))) %>% 
  select(-ends_with(".x"),
         -ends_with(".y"))

# Export
# write.table(prots.cm, file = "whole_PLT_prot_mmus.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Get human orthologs
ort.mouse.human <- gorth(query = prots.cm$UNIPROT, source_organism = "mmusculus",
                         target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE) # 1493



# Venn diagram ------------------------------------------------------------

# Venn diagrams showing overlap between human PLT core proteins and their mouse orthologs

library(eulerr)

ven.prots <- list(Hsa = unique(prots.c$SYMBOL),
                  Mmus.ort = unique(toupper(ort.human.mouse$ortholog_name)))

ven.prots.2 <- list(Hsa.ort = unique(ort.mouse.human$ortholog_name),
                    Mmus = unique(toupper(prots.cm$SYMBOL)))

colors <- c(RColorBrewer::brewer.pal(4, "Set3"))


# Fig. 1B
set.seed(3213)
plot(euler(ven.prots),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


set.seed(3213)
plot(euler(ven.prots.2),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))



# Venn whole proteome -----------------------------------------------------

# Fig. 1A
# Overlap human reference proteome and human PLT core proteome

dfl <- list("Human platelet proteome" = prots.c$UNIPROT, "Whole human proteome" = fasta_headers_human)

colors <- c(RColorBrewer::brewer.pal(3, "Pastel1"))

set.seed(3213)
plot(euler(dfl),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))


# Overlap mouse reference proteome and mouse PLT proteome

dfl.m <- list("Mouse platelet proteome" = prots.cm$UNIPROT, "Whole mouse proteome" = fasta_headers_mouse)

set.seed(3213)
plot(euler(dfl.m),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))



# Whole proteome orthologs

# Fig 1B (upper)
# Human

ort.human.mouse.whole <- gorth(query = fasta_headers_human, source_organism = "hsapiens", 
                               target_organism = "mmusculus", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE) # Get orthologs

dfl.ort.hm <- list("Whole human proteome" = unique(symbol_human), 
                   "Whole mouse proteome, ort" = toupper(unique(ort.human.mouse.whole$ortholog_name)))

set.seed(3213)
plot(euler(dfl.ort.hm),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))

# Mouse

ort.mouse.human.whole <- gorth(query = fasta_headers_mouse, source_organism = "mmusculus", 
                               target_organism = "hsapiens", mthreshold = Inf, filter_na = TRUE) %>% 
  distinct(input, .keep_all = TRUE) # Get orthologs

dfl.ort.mh <- list("Whole mouse proteome" = unique(toupper(symbol_mouse)), 
                   "Whole human proteome, ort" = toupper(unique(ort.mouse.human.whole$ortholog_name)))


set.seed(3213)
plot(euler(dfl.ort.mh),
     quantities = TRUE,
     fills = list(fill = colors),
     labels = NULL,
     legend = list(cex = 1))



# Expression distribution and subsetting - human --------------------------

# Representative PLT proteins
plt.prots <- c("CLEC1B", "SYK", "FERMT3", "ITGB3", "GP6", "PRKCA", "SELP", "ITGB2",
               "PRKG1", "FYN", "MYLK", "AKT2", "AKT1")

# x42 dataset

# Clean, filter, calculate median protein expression and rank proteins
dist.42 <- all.dfs$x42 %>% 
  select(Accession,
         `Gene name`,
         starts_with("APB")) %>% 
  mutate(`Gene name` = gsub("\\..*", "", `Gene name`),
         Accession = gsub("-.*", "", Accession)) %>% 
  distinct(Accession, .keep_all = TRUE) %>% 
  na_if("NaN") %>% 
  type_convert() %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 5) %>% 
  rowwise() %>%
  mutate(median_APB = median(c_across(matches('^APB')), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(rank = dense_rank(dplyr::desc(median_APB))) %>% 
  dplyr::arrange(desc(median_APB)) %>% 
  mutate(rank2 = c(1:1978))


# Sub-setting based on quantiles

q1.names.42 <- dist.42 %>%  # First quantile
  mutate(Accession = gsub("-.*", "", Accession)) %>%
  filter(median_APB < summary(median_APB)[2]) %>% 
  pull(Accession)

q13.names.42 <- dist.42 %>% # Inter-quantile
  mutate(Accession = gsub("-.*", "", Accession)) %>%
  filter(between(median_APB, summary(median_APB)[2], summary(median_APB)[5])) %>% 
  pull(Accession)

q3.names.42 <- dist.42 %>% # Third quantile
  mutate(Accession = gsub("-.*", "", Accession)) %>%
  filter(median_APB > summary(median_APB)[5]) %>%
  pull(Accession)

# Vertical lines marking the change in quantile

v1.42 <- dist.42 %>% # Between the first and inter-quantile
  filter(Accession %in% q1.names.42) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2.42 <- dist.42 %>% # Between inter- and third quantile
  filter(Accession %in% q13.names.42) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_42 <- c(v1.42, v2.42)

# Plotting ranked distribution (relative quantification) - Supp. Fig. 2A

dist.42 %>% 
  ggplot(aes(x = rank, y = median_APB)) +
  geom_vline(xintercept = vlines_42, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.42$rank <= 5) | (dist.42$rank >= max(dist.42$rank) - 5) | (dist.42$`Gene name` %in% plt.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.42 %>% filter((rank <= 5) | (rank >= max(rank) - 4) | (`Gene name` %in% plt.prots)),
                            aes(label = `Gene name`), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1.42, label="1st quantile\n1484", y = 23), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2.42, label="3rd quantile\n496", y = 23), colour="darkgrey", angle = 90, size = 3)




# x11 dataset
# Clean, filter, calculate median protein expression and rank proteins
dist.11 <- eg.x11 %>% 
  select(- `Gene names`) %>% 
  mutate_if(., is.numeric, log2) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA >= 6) %>% 
  left_join(eg.x11[,1:2], by = "Protein IDs") %>% 
  filter(`Protein IDs` %in% fasta_headers_human) %>% # Filter unreviewed proteins
  rowwise() %>%
  mutate(median_LFQ = median(c_across(matches('^LFQ')), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(rank = dense_rank(dplyr::desc(median_LFQ))) %>% 
  dplyr::arrange(desc(median_LFQ)) %>% 
  mutate(rank2 = c(1:1766))

# write.table(dist.11, file = "dist11.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


# Sub-setting based on quantiles

q1.names.11 <- dist.11 %>% # First quantile
  filter(median_LFQ < summary(median_LFQ)[2]) %>% 
  pull(`Protein IDs`)

q13.names.11 <- dist.11 %>% # Inter-quantile
  filter(between(median_LFQ, summary(median_LFQ)[2], summary(median_LFQ)[5])) %>% 
  pull(`Protein IDs`)

q3.names.11 <- dist.11 %>% # Third quantile
  filter(median_LFQ > summary(median_LFQ)[5]) %>%
  pull(`Protein IDs`)

# Vertical lines marking the change in quantile

v1 <- dist.11 %>% 
  filter(`Protein IDs` %in% q1.names.11) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2 <- dist.11 %>% 
  filter(`Protein IDs` %in% q13.names.11) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_11 <- c(v1, v2)


# Plotting ranked distribution (relative quantification) - Fig. 2 (Human)

dist.11 %>% 
  ggplot(aes(x = rank2, y = median_LFQ)) +
  geom_vline(xintercept = vlines_11, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.11$rank2 <= 5) | (dist.11$rank2 >= max(dist.11$rank2) - 5) | (dist.11$`Gene names` %in% plt.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.11 %>% filter((rank2 <= 5) | (rank2 >= max(rank2) - 4) | (`Gene names` %in% plt.prots)),
                            aes(label = `Gene names`), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1, label="1st quantile\n1325", y = 23), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2, label="3rd quantile\n442", y = 23), colour="darkgrey", angle = 90, size = 3)






# x25
# Clean, filter, calculate median protein expression and rank proteins
dist.25 <- all.dfs$x25 %>%
  mutate(`Majority protein IDs` = sub('\\s*;.*','', `Majority protein IDs`),
         `Majority protein IDs` = sub("-.*", "", `Majority protein IDs`)) %>% 
  select(1, 3, 4:6) %>% 
  distinct(`Majority protein IDs`, .keep_all = TRUE) %>% 
  rowwise() %>%
  mutate(median_LFQ = median(c_across(matches('^C')), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(rank = dense_rank(dplyr::desc(median_LFQ))) %>% 
  dplyr::arrange(desc(median_LFQ)) %>% 
  mutate(rank2 = c(1:2504))

# Sub-setting based on quantiles

q1.names.25 <- dist.25 %>% # First quantile
  filter(median_LFQ < summary(median_LFQ)[2]) %>% 
  pull(`Majority protein IDs`)

q13.names.25 <- dist.25 %>% # Inter-quantile
  filter(between(median_LFQ, summary(median_LFQ)[2], summary(median_LFQ)[5])) %>% 
  pull(`Majority protein IDs`)

q3.names.25 <- dist.25 %>% # Third quantile
  filter(median_LFQ > summary(median_LFQ)[5]) %>%
  pull(`Majority protein IDs`)


# Vertical lines marking the change in quantile

v1.25 <- dist.25 %>% 
  filter(`Majority protein IDs` %in% q1.names.25) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2.25 <- dist.25 %>% 
  filter(`Majority protein IDs` %in% q13.names.25) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_25 <- c(v1.25, v2.25)

# Plotting ranked distribution (relative quantification) - Supp. Fig. 2A

dist.25 %>% 
  ggplot(aes(x = rank, y = median_LFQ)) +
  geom_vline(xintercept = vlines_25, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.25$rank <= 5) | (dist.25$rank >= max(dist.25$rank) - 5) | (dist.25$`Gene Symbol` %in% plt.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.25 %>% filter((rank <= 5) | (rank >= max(rank) - 4) | (`Gene Symbol` %in% plt.prots)),
                            aes(label = `Gene Symbol`), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1.25, label="1st quantile\n1879", y = 23), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2.25, label="3rd quantile\n627", y = 23), colour="darkgrey", angle = 90, size = 3)



## x28
# Clean, filter, calculate median protein expression and rank proteins
dist.28 <- all.dfs$x28 %>% 
  select(`Gene names`,
         `Protein IDs`,
         starts_with("LFQ day1-")) %>% 
  mutate(`Gene names` = sub('\\s*;.*','', `Gene names`),
         `Protein IDs` = sub('\\s*;.*','', `Protein IDs`)) %>% 
  rowwise() %>%
  mutate(median_LFQ = median(c_across(matches('^LFQ')), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(rank = dense_rank(dplyr::desc(median_LFQ))) %>% 
  dplyr::arrange(desc(median_LFQ)) %>% 
  mutate(rank2 = c(1:2501))

# Sub-setting based on quantiles

q1.names.28 <- dist.28 %>% # First quantile
  filter(median_LFQ < summary(median_LFQ)[2]) %>% 
  pull(`Protein IDs`)

q13.names.28 <- dist.28 %>% # Inter-quantile
  filter(between(median_LFQ, summary(median_LFQ)[2], summary(median_LFQ)[5])) %>% 
  pull(`Protein IDs`)

q3.names.28 <- dist.28 %>% # Third quantile
  filter(median_LFQ > summary(median_LFQ)[5]) %>%
  pull(`Protein IDs`)


# Vertical lines marking the change in quantile

v1.28 <- dist.28 %>% 
  filter(`Protein IDs` %in% q1.names.28) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2.28 <- dist.28 %>% 
  filter(`Protein IDs` %in% q13.names.28) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_28 <- c(v1.28, v2.28)

# Plotting ranked distribution (relative quantification) - Supp. Fig. 2A

dist.28 %>% 
  ggplot(aes(x = rank, y = median_LFQ)) +
  geom_vline(xintercept = vlines_28, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.28$rank <= 5) | (dist.28$rank >= max(dist.28$rank) - 5) | (dist.28$`Gene names` %in% plt.prots), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.28 %>% filter((rank <= 5) | (rank >= max(rank) - 4) | (`Gene names` %in% plt.prots)),
                            aes(label = `Gene names`), 
                            max.overlaps = Inf,
                            min.segment.length = 0,
                            seed = 42,
                            nudge_x = .15,
                            box.padding = 0.5,
                            nudge_y = 1,
                            segment.curvature = -0.1,
                            segment.ncp = 3,
                            segment.angle = 20) +
  geom_text(aes(x=v1.28, label="1st quantile\n1877", y = 24), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2.28, label="3rd quantile\n626", y = 24), colour="darkgrey", angle = 90, size = 3)





# Expression distribution and subsetting - mouse --------------------------

# Representative PLT proteins
plt.prots.m <- stringr::str_to_title(plt.prots)


# x39
# Clean, filter, calculate median protein expression and rank proteins
cont <- read_excel("list_contaminants_mouse_Mann.xlsx") %>% # Contaminant list
  select(`Gene names`,
         `Protein IDs`) %>% 
  mutate(`Gene names` = sub('\\s*;.*','', `Gene names`),
         `Protein IDs` = sub('\\s*;.*','', `Protein IDs`))

dist.39 <- all.dfs$x39 %>% 
  select(protein_i_ds,
         starts_with("1_"),
         starts_with("2_"),
         starts_with("3_")) %>% 
  rowwise() %>%
  mutate(Count_NA = sum(is.na(cur_data()))) %>%
  ungroup %>% 
  filter(!Count_NA == 5) %>% 
  left_join(all.dfs$x39[,1:2], by = "protein_i_ds") %>% 
  rename(UNIPROT = protein_i_ds,
         SYMBOL = gene_names) %>% 
  mutate(SYMBOL = sub('\\s*;.*','', SYMBOL)) %>% 
  filter(UNIPROT %in% fasta_headers_mouse) %>% # Filter unreviewed proteins
  rowwise() %>%
  mutate(median_LFQ = median(c_across(matches('dag')), na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(cont = ifelse((SYMBOL %in% cont$`Gene names`) | 
                         (UNIPROT %in% cont$`Protein IDs`) | (SYMBOL %in% c("Hbb-b2", "Hbb-b1", "Hba")), "+", NA)) %>% # Remove contaminants
  filter(!cont %in% "+") %>% 
  mutate(rank = dense_rank(dplyr::desc(median_LFQ))) %>% 
  dplyr::arrange(desc(median_LFQ)) %>% 
  mutate(rank2 = c(1:1714))

# write.table(dist.39, file = "dist39.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Sub-setting based on quantiles

q1.names.39 <- dist.39 %>% 
  filter(median_LFQ < summary(median_LFQ)[2]) %>% 
  pull(UNIPROT)

q13.names.39 <- dist.39 %>% 
  filter(between(median_LFQ, summary(median_LFQ)[2], summary(median_LFQ)[5])) %>% 
  pull(UNIPROT)

q3.names.39 <- dist.39 %>% 
  filter(median_LFQ > summary(median_LFQ)[5]) %>%
  pull(UNIPROT)

# Vertical lines marking the change in quantile

v1.39 <- dist.39 %>% 
  filter(UNIPROT %in% q1.names.39) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2.39 <- dist.39 %>% 
  filter(UNIPROT %in% q13.names.39) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_39 <- c(v1.39, v2.39)


# Plotting ranked distribution (relative quantification) - Fig. 2 (Mouse)

dist.39 %>% 
  ggplot(aes(x = rank2, y = median_LFQ)) +
  geom_vline(xintercept = vlines_39, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.39$rank2 <= 5) | (dist.39$rank2 >= max(dist.39$rank2) - 5) | (dist.39$SYMBOL %in% plt.prots.m), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.39 %>% filter((rank2 <= 5) | (rank2 >= max(rank2) - 4) | (SYMBOL %in% plt.prots.m)),
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
  geom_text(aes(x=v1.39, label="1st quantile\n1286", y = 12), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2.39, label="3rd quantile\n430", y = 12), colour="darkgrey", angle = 90, size = 3)





# x35
# Clean, filter, calculate median protein expression and rank proteins
dist.35 <- all.dfs$x35 %>% 
  select(`Accession ID (UniProt)`,
         `Gene names`,
         `LFQ intensity ultra purified`) %>% 
  mutate(`Accession ID (UniProt)` = sub("-.*", "", `Accession ID (UniProt)`),
         `Gene names` = sub("-.*", "", `Gene names`)) %>% 
  rename(UNIPROT = `Accession ID (UniProt)`,
         SYMBOL = `Gene names`) %>% 
  mutate(SYMBOL = sub('\\s*;.*','', SYMBOL)) %>% 
  distinct(UNIPROT, .keep_all = TRUE) %>% 
  rowwise() %>%
  mutate(median_LFQ = `LFQ intensity ultra purified`) %>% 
  ungroup() %>% 
  mutate(rank = dense_rank(dplyr::desc(median_LFQ))) %>% 
  dplyr::arrange(desc(median_LFQ)) %>% 
  mutate(rank2 = c(1:4345))

# Sub-setting based on quantiles

q1.names.35 <- dist.35 %>% 
  filter(median_LFQ < summary(median_LFQ)[2]) %>% 
  pull(UNIPROT)

q13.names.35 <- dist.35 %>% 
  filter(between(median_LFQ, summary(median_LFQ)[2], summary(median_LFQ)[5])) %>% 
  pull(UNIPROT)

q3.names.35 <- dist.35 %>% 
  filter(median_LFQ > summary(median_LFQ)[5]) %>%
  pull(UNIPROT)

# Vertical lines marking the change in quantile

v1.35 <- dist.35 %>% 
  filter(UNIPROT %in% q1.names.35) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

v2.35 <- dist.35 %>% 
  filter(UNIPROT %in% q13.names.35) %>% 
  filter(rank2 == min(rank2)) %>% 
  pull(rank2)

vlines_35 <- c(v1.35, v2.35)

# Plotting ranked distribution (relative quantification) - Supp. Fig. 2B

dist.35 %>% 
  ggplot(aes(x = rank, y = median_LFQ)) +
  geom_vline(xintercept = vlines_35, alpha = .2, linetype = 2) +
  geom_point(color = ifelse((dist.35$rank <= 5) | (dist.35$rank >= max(dist.35$rank) - 5) | (dist.35$SYMBOL %in% plt.prots.m), "red", "black")) +
  theme_classic() +
  labs(x = "Rank", y = "Log2-LFQ intensities") +
  ggrepel::geom_label_repel(data = dist.35 %>% filter((rank <= 5) | (rank >= max(rank) - 4) | (SYMBOL %in% plt.prots.m)),
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
  geom_text(aes(x=v1.35, label="1st quantile\n3260", y = 33), colour="darkgrey", angle = 90, size = 3) +
  geom_text(aes(x=v2.35, label="3rd quantile\n1087", y = 33), colour="darkgrey", angle = 90, size = 3)



# GO enrichment - HUMAN ---------------------------------------------------

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
  
} # Enrichment GO terms

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
} # Transform to uppercase the first letter of description

# x11 dataset

# Fig. 2 (Human)

# First quantile
q1.go.11 <- gos(q1.names.11)
q1.go.11.bp <- simplify(q1.go.11$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q1.go.11.bp@result$Description <- firstup(q1.go.11.bp@result$Description)
dotplot(q1.go.11.bp)

# Inter-quantile
q13.go.11 <- gos(q13.names.11)
q13.go.11.bp <- simplify(q13.go.11$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q13.go.11.bp@result$Description <- firstup(q13.go.11.bp@result$Description)
dotplot(q13.go.11.bp)

# Third quantile
q3.go.11 <- gos(q3.names.11)
q3.go.11.bp <- simplify(q3.go.11$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q3.go.11.bp@result$Description <- firstup(q3.go.11.bp@result$Description)
dotplot(q3.go.11.bp)



# GO enrichment - MOUSE ---------------------------------------------------

# x39 dataset

# Fig. 2 (Mouse)

# First quantile
q1.go.39 <- gos(q1.names.39, org = "Mouse")
q1.go.39.bp <- simplify(q1.go.39$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q1.go.39.bp@result$Description <- firstup(q1.go.39.bp@result$Description)
dotplot(q1.go.39.bp)

# Inter-quantile
q13.go.39 <- gos(q13.names.39, org = "Mouse")
q13.go.39.bp <- simplify(q13.go.39$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q13.go.39.bp@result$Description <- firstup(q13.go.39.bp@result$Description)
dotplot(q13.go.39.bp)

# Third quantile
q3.go.39 <- gos(q3.names.39, org = "Mouse")
q3.go.39.bp <- simplify(q3.go.39$go.bp, cutoff = 0.7, by = "p.adjust", select_fun = min)

q3.go.39.bp@result$Description <- firstup(q3.go.39.bp@result$Description)
dotplot(q3.go.39.bp)



# Output tables -----------------------------------------------------------

# output_tables <- list(all_human_prots = res,
#                       selected_human_prots = prots.c,
#                       all_mouse_prots = resm,
#                       selected_mouse_prots = prots.cm)

# library(openxlsx)
# write.xlsx(output_tables, file = "PLT_lysate_output_tables.xlsx")
