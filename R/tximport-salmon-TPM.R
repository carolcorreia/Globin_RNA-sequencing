## Globin analysis
## Summarising gene-level TPM estimates after Salmon
## 22/05/17

setwd("/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon")
img_dir <- getwd()

load("tximport-salmon-TPM.RData")

source("https://bioconductor.org/biocLite.R")
biocLite("rtracklayer", type = "source")
biocLite("tximport")
biocLite("GenomicFeatures", type = "source")

install.packages("rjson")

library(plyr)
library(tidyverse)
library(ggrepel)
library(magrittr)
library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(tximport)
library(rjson)
library(reshape2)
library(biomaRt)

# Summary of salmon stats
h_salmon_stats <- read.table("salmon_summary_human.txt",
                             header = TRUE)

h_salmon_stats


p_salmon_stats <- read.table("salmon_summary_pig.txt",
                             header = TRUE)

p_salmon_stats


c_salmon_stats <- read.table("salmon_summary_cattle.txt",
                             header = TRUE)

c_salmon_stats

# Clean names of samples in stats tables
h_salmon_stats$File_name %<>%
    str_replace("_salmon_quant.log", "") %>% 
    str_replace("NGD_(FALSE|TRUE)_", "") %>% 
    str_replace("GD_(FALSE|TRUE)_", "") %>% 
    str_replace("Subj10", "Hsa_P10_D") %>% 
    str_replace("Subj11", "Hsa_P11_D") %>% 
    str_replace("Subj12", "Hsa_P12_D") %>% 
    str_replace("Subj13", "Hsa_S01_U") %>% 
    str_replace("Subj14", "Hsa_S02_U") %>% 
    str_replace("Subj15", "Hsa_S03_U") %>% 
    str_replace("Subj16", "Hsa_S04_U") %>% 
    str_replace("Subj17", "Hsa_S05_U") %>% 
    str_replace("Subj18", "Hsa_S06_U") %>% 
    str_replace("Subj19", "Hsa_P19_U") %>% 
    str_replace("Subj20", "Hsa_P20_U") %>% 
    str_replace("Subj21", "Hsa_P21_U") %>% 
    str_replace("Subj22", "Hsa_P22_U") %>% 
    str_replace("Subj23", "Hsa_P23_U") %>% 
    str_replace("Subj24", "Hsa_P24_U") %>% 
    str_replace("Subj1", "Hsa_S01_D") %>% 
    str_replace("Subj2", "Hsa_S02_D") %>% 
    str_replace("Subj3", "Hsa_S03_D") %>% 
    str_replace("Subj4", "Hsa_S04_D") %>% 
    str_replace("Subj5", "Hsa_S05_D") %>% 
    str_replace("Subj6", "Hsa_S06_D") %>% 
    str_replace("Subj7", "Hsa_P07_D") %>% 
    str_replace("Subj8", "Hsa_P08_D") %>% 
    str_replace("Subj9", "Hsa_P09_D") 

p_salmon_stats$File_name %<>%
    str_replace("_salmon_quant.log", "") %>% 
    str_replace("7197", "Ssc_01") %>% 
    str_replace("7199", "Ssc_02") %>% 
    str_replace("7210", "Ssc_03") %>% 
    str_replace("7312", "Ssc_04") %>% 
    str_replace("7349", "Ssc_05") %>%
    str_replace("7413", "Ssc_06") %>%
    str_replace("7437", "Ssc_07") %>%
    str_replace("7439", "Ssc_08") %>%
    str_replace("7467", "Ssc_09") %>%
    str_replace("7468", "Ssc_10") %>%
    str_replace("7472", "Ssc_11") %>%
    str_replace("7474", "Ssc_12") %>%
    str_replace("(CT|C)", "_U") %>% 
    str_replace("(GD|T)", "_D")


c_salmon_stats$File_name %<>%
    str_replace("_salmon_quant.log", "") %>% 
    str_replace("A", "") %>% 
    str_replace("_W-1_F", "") %>%
    str_replace("6511", "Bta_01_U") %>%
    str_replace("6514", "Bta_02_U") %>%
    str_replace("6520", "Bta_03_U") %>%
    str_replace("6522", "Bta_04_U") %>%
    str_replace("6526", "Bta_05_U") %>%
    str_replace("6635", "Bta_06_U") %>%
    str_replace("6636", "Bta_07_U") %>%
    str_replace("6637", "Bta_08_U") %>%
    str_replace("6644", "Bta_09_U") %>%
    str_replace("6698", "Bta_10_U")

h_salmon_stats$Library_type %<>% 
    str_replace("IU", "Inward Unstranded")
p_salmon_stats$Library_type %<>% 
    str_replace("IU", "Inward Unstranded")
c_salmon_stats$Library_type %<>% 
    str_replace("ISF", "Inward Stranded Forward")

# Export renamed salmon stats
write.table(h_salmon_stats,
            file = "h_salmon_stats.txt",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)

write.table(p_salmon_stats,
            file = "p_salmon_stats.txt",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)

write.table(c_salmon_stats,
            file = "c_salmon_stats.txt",
            sep = "\t",
            col.names = TRUE,
            quote = FALSE)

# Get file names
path_human <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/human_TPM"
files_human <- list.files(path        = path_human,
                          pattern     = "^(GD|NGD)",
                          all.files   = TRUE,
                          full.names  = FALSE,
                          recursive   = FALSE,
                          ignore.case = FALSE)


files_human


path_pig <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/pig_TPM"
files_pig <- list.files(path        = path_pig,
                        pattern     = "^7",
                        all.files   = TRUE,
                        full.names  = FALSE,
                        recursive   = FALSE,
                        ignore.case = FALSE)

files_pig


path_cattle <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/cattle_TPM"
files_cattle <- list.files(path        = path_cattle,
                           pattern     = "^A",
                           all.files   = TRUE,
                           full.names  = FALSE,
                           recursive   = FALSE,
                           ignore.case = FALSE)

files_cattle

# Get full path for each file
human_samples <- (file.path(path_human, files_human))
names(human_samples) <- paste0(files_human)
human_samples

pig_samples <- (file.path(path_pig, files_pig))
names(pig_samples) <- paste0(files_pig)
pig_samples

cattle_samples <- (file.path(path_cattle, files_cattle))
names(cattle_samples) <- paste0(files_cattle)
cattle_samples

# Create annotation DB from GTF files
# ftp://ftp.ensembl.org/pub/release-88/gtf/homo_sapiens/Homo_sapiens.GRCh38.88.gtf.gz

# ftp://ftp.ensembl.org/pub/release-88/gtf/sus_scrofa/Sus_scrofa.Sscrofa10.2.88.gtf.gz

# ftp://ftp.ensembl.org/pub/release-88/gtf/bos_taurus/Bos_taurus.UMD3.1.88.gtf.gz

human_gtf <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/GTF_ensembl/Homo_sapiens.GRCh38.88.gtf"
humanDB <- makeTxDbFromGFF(file = human_gtf,
                           format = "gtf",
                           dataSource = "Ensembl GTF Homo_sapiens.GRCh38.88",
                           organism = "Homo sapiens")

humanDB


pig_gtf <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/GTF_ensembl/Sus_scrofa.Sscrofa10.2.88.gtf"
pigDB <- makeTxDbFromGFF(file = pig_gtf,
                         format = "gtf",
                         dataSource = "Ensembl GTF Sus_scrofa.Sscrofa10.2.88",
                         organism = "Sus scrofa")

pigDB


cattle_gtf <- "/Users/ccorreia/Dropbox/CSF/Animal_Genomics/Globin/salmon/GTF_ensembl/Bos_taurus.UMD3.1.88.gtf"
cattleDB <- makeTxDbFromGFF(file = cattle_gtf,
                            format = "gtf",
                            dataSource = "Ensembl GTF Bos_taurus.UMD3.1.88",
                            organism = "Bos taurus")

cattleDB

# Create dataframe of transcripts and gene IDs (it has to be in this order)
keytypes(humanDB)
H_keys <- keys(humanDB, keytype = "GENEID")
H_tx2gene <- AnnotationDbi::select(humanDB,
                              keys = H_keys,
                              keytype = "GENEID",
                              columns = "TXNAME")

H_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(H_tx2gene)


keytypes(pigDB)
P_keys <- keys(pigDB, keytype = "GENEID")
P_tx2gene <- AnnotationDbi::select(pigDB,
                                   keys = P_keys,
                                   keytype = "GENEID",
                                   columns = "TXNAME")

P_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(P_tx2gene)

keytypes(cattleDB)
C_keys <- keys(cattleDB, keytype = "GENEID")
C_tx2gene <- AnnotationDbi::select(cattleDB,
                                   keys = C_keys,
                                   keytype = "GENEID",
                                   columns = "TXNAME")

C_tx2gene %<>% 
    dplyr::select(TXNAME, GENEID)

head(C_tx2gene)

# Import transcript-level estimates summarized to the gene-level
# (it's the default when txOut = FALSE or not provided)
# (ignoreTxVersion = TRUE otherwise the transcripts names in the GTF file 
# won't match the ones in salmon's quant.sf)
human_txi <- tximport(human_samples,
                      type = "salmon",
                      tx2gene = H_tx2gene,
                      ignoreTxVersion = TRUE)

names(human_txi)
head(human_txi$abundance)


pig_txi <- tximport(pig_samples,
                    type = "salmon",
                    tx2gene = P_tx2gene,
                    ignoreTxVersion = TRUE)

names(pig_txi)
head(pig_txi$abundance)


cattle_txi <- tximport(cattle_samples,
                       type = "salmon",
                       tx2gene = C_tx2gene,
                       ignoreTxVersion = TRUE)

names(cattle_txi)
head(cattle_txi$abundance)


# Get datasets names from ensembl biomart
listMarts(mart = NULL, host="www.ensembl.org", path="/biomart/martservice",
          port = 80)
mart_datasets <- useMart("ENSEMBL_MART_ENSEMBL")
listDatasets(mart_datasets)

# Get gene symbols
human_TPM <- as.data.frame(human_txi$abundance)
mart_human <- useMart("ENSEMBL_MART_ENSEMBL",
                      dataset = "hsapiens_gene_ensembl")
filters_human <- listFilters(mart_human)
filters_human

h_genes <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = rownames(human_TPM),
                 mart = mart_human)

head(h_genes)
#anns2 <- anns[match(rownames(dds), anns[, 1]), ]
#rownames(anns2) <- rownames(dds)

"ENSG00000188536" %in% rownames(human_txi$abundance) #HBA2
"ENSG00000206172" %in% rownames(human_txi$abundance) #HBA1
"ENSG00000244734" %in% rownames(human_txi$abundance) #HBB



pig_TPM <- as.data.frame(pig_txi$abundance)
mart_pig <- useMart("ENSEMBL_MART_ENSEMBL",
                    dataset = "sscrofa_gene_ensembl")
filters_pig <- listFilters(mart_pig)
filters_pig

p_genes <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = rownames(pig_TPM),
                 mart = mart_pig)

head(p_genes) # very few returned gene symbols

"ENSSSCG00000007978" %in% rownames(pig_txi$abundance) #HBA1
"ENSSSCG00000014725" %in% rownames(pig_txi$abundance) #HBB
"HBA" %in% p_genes$hgnc_symbol
"HBA" %in% p_genes$external_gene_name
"HBA1" %in% p_genes$hgnc_symbol
"HBA1" %in% p_genes$external_gene_name
"HBB" %in% p_genes$hgnc_symbol
"HBB" %in% p_genes$external_gene_name



cattle_TPM <- as.data.frame(cattle_txi$abundance)
mart_cattle <- useMart("ENSEMBL_MART_ENSEMBL",
                       dataset = "btaurus_gene_ensembl")
filters_cattle <- listFilters(mart_cattle)
filters_cattle

c_genes <- getBM(attributes = c("ensembl_gene_id",
                                "hgnc_symbol",
                                "external_gene_name"),
                 filters = "ensembl_gene_id",
                 values = rownames(cattle_TPM),
                 mart = mart_cattle)

head(c_genes) 

"ENSBTAG00000026417" %in% rownames(cattle_txi$abundance) #HBA
"ENSBTAG00000026418" %in% rownames(cattle_txi$abundance) #HBA
"ENSBTAG00000038748" %in% rownames(cattle_txi$abundance) #HBB

# Remove non-expressed genes
h_TPM_nozeros <- human_TPM[rowSums(human_TPM) > 0, ]
dim(h_TPM_nozeros)
dim(human_TPM)

p_TPM_nozeros <- pig_TPM[rowSums(pig_TPM) > 0, ]
dim(p_TPM_nozeros)
dim(pig_TPM)

c_TPM_nozeros <- cattle_TPM[rowSums(cattle_TPM) > 0, ]
dim(c_TPM_nozeros)
dim(cattle_TPM)

# Remove lowly expressed genes (< 1 TPM)
h_TPM_filt <- h_TPM_nozeros[rowSums(h_TPM_nozeros >= 1) >= 12, ]
dim(h_TPM_filt)
dim(h_TPM_nozeros)

p_TPM_filt <- p_TPM_nozeros[rowSums(p_TPM_nozeros >= 1) >= 12, ]
dim(p_TPM_filt)
dim(p_TPM_nozeros)

c_TPM_filt <- c_TPM_nozeros[rowSums(c_TPM_nozeros >= 1) >= 10, ]
dim(c_TPM_filt)
dim(c_TPM_nozeros)


# Tidy nozeros TPM dfs
h_TPM_nozeros %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(h_TPM_nozeros)

p_TPM_nozeros %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(p_TPM_nozeros)

c_TPM_nozeros %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(c_TPM_nozeros)

# Tidy filt TPM dfs
h_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(h_TPM_filt)

p_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(p_TPM_filt)

c_TPM_filt %<>% 
    rownames_to_column(var = "Ensembl_gene_ID") %>% 
    melt(value.name = "TPM") %>% 
    dplyr::rename(sample = variable)

head(c_TPM_filt)

# Add plotting labels to nozeros df
h_TPM_nozeros$labels <- h_TPM_nozeros$sample

h_TPM_nozeros$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("NGD_(FALSE|TRUE)_", "") %>% 
    str_replace("GD_(FALSE|TRUE)_", "") %>% 
    str_replace("Subj10", "Hsa_P10_D") %>% 
    str_replace("Subj11", "Hsa_P11_D") %>% 
    str_replace("Subj12", "Hsa_P12_D") %>% 
    str_replace("Subj13", "Hsa_S01_U") %>% 
    str_replace("Subj14", "Hsa_S02_U") %>% 
    str_replace("Subj15", "Hsa_S03_U") %>% 
    str_replace("Subj16", "Hsa_S04_U") %>% 
    str_replace("Subj17", "Hsa_S05_U") %>% 
    str_replace("Subj18", "Hsa_S06_U") %>% 
    str_replace("Subj19", "Hsa_P19_U") %>% 
    str_replace("Subj20", "Hsa_P20_U") %>% 
    str_replace("Subj21", "Hsa_P21_U") %>% 
    str_replace("Subj22", "Hsa_P22_U") %>% 
    str_replace("Subj23", "Hsa_P23_U") %>% 
    str_replace("Subj24", "Hsa_P24_U") %>% 
    str_replace("Subj1", "Hsa_S01_D") %>% 
    str_replace("Subj2", "Hsa_S02_D") %>% 
    str_replace("Subj3", "Hsa_S03_D") %>% 
    str_replace("Subj4", "Hsa_S04_D") %>% 
    str_replace("Subj5", "Hsa_S05_D") %>% 
    str_replace("Subj6", "Hsa_S06_D") %>% 
    str_replace("Subj7", "Hsa_P07_D") %>% 
    str_replace("Subj8", "Hsa_P08_D") %>% 
    str_replace("Subj9", "Hsa_P09_D") 

unique(h_TPM_nozeros$labels)


p_TPM_nozeros$labels <- p_TPM_nozeros$sample

p_TPM_nozeros$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("A", "") %>% 
    str_replace("7197", "Ssc_01") %>% 
    str_replace("7199", "Ssc_02") %>% 
    str_replace("7210", "Ssc_03") %>% 
    str_replace("7312", "Ssc_04") %>% 
    str_replace("7349", "Ssc_05") %>%
    str_replace("7413", "Ssc_06") %>%
    str_replace("7437", "Ssc_07") %>%
    str_replace("7439", "Ssc_08") %>%
    str_replace("7467", "Ssc_09") %>%
    str_replace("7468", "Ssc_10") %>%
    str_replace("7472", "Ssc_11") %>%
    str_replace("7474", "Ssc_12") %>%
    str_replace("(CT|C)", "_U") %>% 
    str_replace("(GD|T)", "_D")

unique(p_TPM_nozeros$labels)


c_TPM_nozeros$labels <- c_TPM_nozeros$sample

c_TPM_nozeros$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("A", "") %>% 
    str_replace("_W-1_F", "") %>%
    str_replace("6511", "Bta_01_U") %>%
    str_replace("6514", "Bta_02_U") %>%
    str_replace("6520", "Bta_03_U") %>%
    str_replace("6522", "Bta_04_U") %>%
    str_replace("6526", "Bta_05_U") %>%
    str_replace("6635", "Bta_06_U") %>%
    str_replace("6636", "Bta_07_U") %>%
    str_replace("6637", "Bta_08_U") %>%
    str_replace("6644", "Bta_09_U") %>%
    str_replace("6698", "Bta_10_U")

unique(c_TPM_nozeros$labels)


# Add plotting labels to filt df
h_TPM_filt$labels <- h_TPM_filt$sample

h_TPM_filt$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("NGD_(FALSE|TRUE)_", "") %>% 
    str_replace("GD_(FALSE|TRUE)_", "") %>% 
    str_replace("Subj10", "Hsa_P10_D") %>% 
    str_replace("Subj11", "Hsa_P11_D") %>% 
    str_replace("Subj12", "Hsa_P12_D") %>% 
    str_replace("Subj13", "Hsa_S01_U") %>% 
    str_replace("Subj14", "Hsa_S02_U") %>% 
    str_replace("Subj15", "Hsa_S03_U") %>% 
    str_replace("Subj16", "Hsa_S04_U") %>% 
    str_replace("Subj17", "Hsa_S05_U") %>% 
    str_replace("Subj18", "Hsa_S06_U") %>% 
    str_replace("Subj19", "Hsa_P19_U") %>% 
    str_replace("Subj20", "Hsa_P20_U") %>% 
    str_replace("Subj21", "Hsa_P21_U") %>% 
    str_replace("Subj22", "Hsa_P22_U") %>% 
    str_replace("Subj23", "Hsa_P23_U") %>% 
    str_replace("Subj24", "Hsa_P24_U") %>% 
    str_replace("Subj1", "Hsa_S01_D") %>% 
    str_replace("Subj2", "Hsa_S02_D") %>% 
    str_replace("Subj3", "Hsa_S03_D") %>% 
    str_replace("Subj4", "Hsa_S04_D") %>% 
    str_replace("Subj5", "Hsa_S05_D") %>% 
    str_replace("Subj6", "Hsa_S06_D") %>% 
    str_replace("Subj7", "Hsa_P07_D") %>% 
    str_replace("Subj8", "Hsa_P08_D") %>% 
    str_replace("Subj9", "Hsa_P09_D") 

unique(h_TPM_filt$labels)


p_TPM_filt$labels <- p_TPM_filt$sample

p_TPM_filt$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("A", "") %>% 
    str_replace("7197", "Ssc_01") %>% 
    str_replace("7199", "Ssc_02") %>% 
    str_replace("7210", "Ssc_03") %>% 
    str_replace("7312", "Ssc_04") %>% 
    str_replace("7349", "Ssc_05") %>%
    str_replace("7413", "Ssc_06") %>%
    str_replace("7437", "Ssc_07") %>%
    str_replace("7439", "Ssc_08") %>%
    str_replace("7467", "Ssc_09") %>%
    str_replace("7468", "Ssc_10") %>%
    str_replace("7472", "Ssc_11") %>%
    str_replace("7474", "Ssc_12") %>%
    str_replace("(CT|C)", "_U") %>% 
    str_replace("(GD|T)", "_D")

unique(p_TPM_filt$labels)


c_TPM_filt$labels <- c_TPM_filt$sample

c_TPM_filt$labels %<>%
    str_replace("_quant.sf", "") %>% 
    str_replace("A", "") %>% 
    str_replace("_W-1_F", "") %>%
    str_replace("6511", "Bta_01_U") %>%
    str_replace("6514", "Bta_02_U") %>%
    str_replace("6520", "Bta_03_U") %>%
    str_replace("6522", "Bta_04_U") %>%
    str_replace("6526", "Bta_05_U") %>%
    str_replace("6635", "Bta_06_U") %>%
    str_replace("6636", "Bta_07_U") %>%
    str_replace("6637", "Bta_08_U") %>%
    str_replace("6644", "Bta_09_U") %>%
    str_replace("6698", "Bta_10_U")

unique(c_TPM_filt$labels)

# Convert labels into factors
h_TPM_nozeros$labels %<>% 
    factor(levels = c("Hsa_S01_U", "Hsa_S02_U", "Hsa_S03_U", "Hsa_S04_U",
                      "Hsa_S05_U", "Hsa_S06_U", "Hsa_P19_U", "Hsa_P20_U",
                      "Hsa_P21_U", "Hsa_P22_U", "Hsa_P23_U", "Hsa_P24_U",
                      "Hsa_S01_D", "Hsa_S02_D", "Hsa_S03_D", "Hsa_S04_D",
                      "Hsa_S05_D", "Hsa_S06_D", "Hsa_P07_D", "Hsa_P08_D",
                      "Hsa_P09_D", "Hsa_P10_D", "Hsa_P11_D", "Hsa_P12_D"))
levels(h_TPM_nozeros$labels)

h_TPM_filt$labels %<>% 
    factor(levels = c("Hsa_S01_U", "Hsa_S02_U", "Hsa_S03_U", "Hsa_S04_U",
                      "Hsa_S05_U", "Hsa_S06_U", "Hsa_P19_U", "Hsa_P20_U",
                      "Hsa_P21_U", "Hsa_P22_U", "Hsa_P23_U", "Hsa_P24_U",
                      "Hsa_S01_D", "Hsa_S02_D", "Hsa_S03_D", "Hsa_S04_D",
                      "Hsa_S05_D", "Hsa_S06_D", "Hsa_P07_D", "Hsa_P08_D",
                      "Hsa_P09_D", "Hsa_P10_D", "Hsa_P11_D", "Hsa_P12_D"))
levels(h_TPM_filt$labels)




p_TPM_nozeros$labels %<>% 
    factor(levels = c("Ssc_01_U", "Ssc_02_U", "Ssc_03_U", "Ssc_04_U",
                      "Ssc_05_U", "Ssc_06_U", "Ssc_07_U", "Ssc_08_U",
                      "Ssc_09_U", "Ssc_10_U", "Ssc_11_U", "Ssc_12_U",
                      "Ssc_01_D", "Ssc_02_D", "Ssc_03_D", "Ssc_04_D",
                      "Ssc_05_D", "Ssc_06_D", "Ssc_07_D", "Ssc_08_D",
                      "Ssc_09_D", "Ssc_10_D", "Ssc_11_D", "Ssc_12_D"))
levels(p_TPM_nozeros$labels)

p_TPM_filt$labels %<>% 
    factor(levels = c("Ssc_01_U", "Ssc_02_U", "Ssc_03_U", "Ssc_04_U",
                      "Ssc_05_U", "Ssc_06_U", "Ssc_07_U", "Ssc_08_U",
                      "Ssc_09_U", "Ssc_10_U", "Ssc_11_U", "Ssc_12_U",
                      "Ssc_01_D", "Ssc_02_D", "Ssc_03_D", "Ssc_04_D",
                      "Ssc_05_D", "Ssc_06_D", "Ssc_07_D", "Ssc_08_D",
                      "Ssc_09_D", "Ssc_10_D", "Ssc_11_D", "Ssc_12_D"))
levels(p_TPM_filt$labels)



c_TPM_nozeros$labels %<>% 
    factor(levels = c("Bta_01_U", "Bta_02_U", "Bta_03_U", "Bta_04_U", "Bta_05_U",
                      "Bta_06_U", "Bta_07_U", "Bta_08_U", "Bta_09_U", "Bta_10_U"))
levels(c_TPM_nozeros$labels)

c_TPM_filt$labels %<>% 
    factor(levels = c("Bta_01_U", "Bta_02_U", "Bta_03_U", "Bta_04_U", "Bta_05_U",
                      "Bta_06_U", "Bta_07_U", "Bta_08_U", "Bta_09_U", "Bta_10_U"))
levels(c_TPM_filt$labels)


# Add treatment column
h_TPM_nozeros$treatment <- h_TPM_nozeros$labels
h_TPM_nozeros$treatment %<>% 
    str_replace("Hsa_(S|P)\\d\\d_", "")
unique(h_TPM_nozeros$treatment)

h_TPM_filt$treatment <- h_TPM_filt$labels
h_TPM_filt$treatment %<>% 
    str_replace("Hsa_(S|P)\\d\\d_", "")
unique(h_TPM_filt$treatment)


p_TPM_nozeros$treatment <- p_TPM_nozeros$labels
p_TPM_nozeros$treatment %<>% 
    str_replace("Ssc_\\d\\d_", "")
unique(p_TPM_nozeros$treatment)

p_TPM_filt$treatment <- p_TPM_filt$labels
p_TPM_filt$treatment %<>% 
    str_replace("Ssc_\\d\\d_", "")
unique(p_TPM_filt$treatment)


c_TPM_nozeros$treatment <- c_TPM_nozeros$labels
c_TPM_nozeros$treatment %<>% 
    str_replace("Bta_\\d\\d_", "")
unique(c_TPM_nozeros$treatment)

c_TPM_filt$treatment <- c_TPM_filt$labels
c_TPM_filt$treatment %<>% 
    str_replace("Bta_\\d\\d_", "")
unique(c_TPM_filt$treatment)

# Total number of genes by treatment after tximport summarisation and zeros removed
h_TPM_nozeros %>%
    dplyr::filter(treatment == "U") %>% 
    plyr::count("Ensembl_gene_ID") %>% 
    dim()

h_TPM_nozeros %>%
    dplyr::filter(treatment == "D") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

unique(h_TPM_nozeros$Ensembl_gene_ID) %>% 
    length()

p_TPM_nozeros %>%
    dplyr::filter(treatment == "U") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

p_TPM_nozeros %>%
    dplyr::filter(treatment == "D") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

unique(p_TPM_nozeros$Ensembl_gene_ID) %>% 
    length()

unique(c_TPM_nozeros$Ensembl_gene_ID) %>% 
    length()

# Total number of genes by treatment after tximport summarisation and TPM >= 1

h_TPM_filt %>%
    dplyr::filter(treatment == "U") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

h_TPM_filt %>%
    dplyr::filter(treatment == "D") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

unique(h_TPM_filt$Ensembl_gene_ID) %>% 
    length()

p_TPM_filt %>%
    dplyr::filter(treatment == "U") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

p_TPM_filt %>%
    dplyr::filter(treatment == "D") %>% 
    dplyr::count(Ensembl_gene_ID) %>% 
    dim()

unique(p_TPM_filt$Ensembl_gene_ID) %>% 
    length()

unique(c_TPM_filt$Ensembl_gene_ID) %>% 
    length()


# Bind gene-level TPM filt dfs
A_TPM_filt <- dplyr::bind_rows(h_TPM_filt, p_TPM_filt)
A_TPM_filt <- dplyr::bind_rows(A_TPM_filt, c_TPM_filt)
dim(A_TPM_filt)
dim(h_TPM_filt)
dim(p_TPM_filt)
dim(c_TPM_filt)

# Add species column A_TPM_filt
A_TPM_filt$species <- A_TPM_filt$Ensembl_gene_ID
A_TPM_filt$species %<>%
    stringr::str_replace("ENSG0.+", "Human") %>% 
    stringr::str_replace("ENSSSCG0.+", "Porcine") %>% 
    stringr::str_replace("ENSBTAG0.+", "Bovine") 

# Treatment factors allG_TPM_filt
A_TPM_filt$treatment %<>% factor(levels = c("U", "D"))

A_TPM_filt$treatment %<>%
    stringr::str_replace("U", "Undepleted") %>% 
    stringr::str_replace("D", "Globin depleted") %>%
    factor(levels = c("Undepleted", "Globin depleted"))

# Density plot of gene-level TPM after filtering (FIGURE 1)
A_TPM_filt_plot <- ggplot(A_TPM_filt) +
                        geom_density(aes(log(TPM + 1),
                                         group = treatment,
                                         fill = treatment),
                                     alpha = 0.5) +
                        scale_fill_manual("Treatment",
                                          values = c("#af8dc3", "#7fbf7b")) +
                        facet_grid(. ~ species) +
                        theme_bw() +
                        ylab("Density of gene-level TPM estimates") +
                        xlab(expression(paste(log[2], "(TPM + 1)")))

A_TPM_filt_plot

ggsave("A_TPM_filt_plot.png",
       plot      = A_TPM_filt_plot,
       limitsize = FALSE,
       dpi       = 300,
       path      = img_dir)

# Proportion of globin within all genes after filtering
h_TPM_filt %>% # Sum total TPM per gene, within each treat
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    dplyr::summarise(sum_TPM = sum(TPM)) -> summary_h_TPM_filt

summary_h_TPM_filt %>% # Sum all TPMs by treat
    dplyr::filter(treatment == "U") -> U_summary_h_TPM_filt
colSums(as.data.frame(U_summary_h_TPM_filt$sum_TPM))

U_summary_h_TPM_filt %>% # Sum globin TPMs by treat
    dplyr::filter(Ensembl_gene_ID %in% c("ENSG00000206172",
                                         "ENSG00000188536",
                                         "ENSG00000244734"))

summary_h_TPM_filt %>%
    dplyr::filter(treatment == "D") -> D_summary_h_TPM_filt
colSums(as.data.frame(D_summary_h_TPM_filt$sum_TPM))

D_summary_h_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSG00000206172",
                                         "ENSG00000188536",
                                         "ENSG00000244734"))


p_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    dplyr::summarise(sum_TPM = sum(TPM)) -> summary_p_TPM_filt

summary_p_TPM_filt %>%
    dplyr::filter(treatment == "U") -> U_summary_p_TPM_filt
colSums(as.data.frame(U_summary_p_TPM_filt$sum_TPM))

U_summary_p_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSSSCG00000007978",
                                         "ENSSSCG00000014725"))

summary_p_TPM_filt %>%
    dplyr::filter(treatment == "D") -> D_summary_p_TPM_filt
colSums(as.data.frame(D_summary_p_TPM_filt$sum_TPM))

D_summary_p_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSSSCG00000007978",
                                         "ENSSSCG00000014725"))



c_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    summarise(sum_TPM = sum(TPM)) -> summary_c_TPM_filt

summary_c_TPM_filt %>%
    dplyr::filter(treatment == "U") -> U_summary_c_TPM_filt
colSums(as.data.frame(U_summary_c_TPM_filt$sum_TPM))

U_summary_c_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSBTAG00000026417",
                                         "ENSBTAG00000026418",
                                         "ENSBTAG00000038748"))

# Dot plot of globin gene-level filtered TPM (FIGURE 2)

U_summary_c_TPM_filt
ggplot(U_summary_c_TPM_filt, aes(x = size, fill = type)) +
    geom_dotplot(method="histodot", stackgroups = TRUE)





# Subset globin genes by Ensembl ID
h_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSG00000206172",
                                         "ENSG00000188536",
                                         "ENSG00000244734")) -> Hg_TPM_filt
unique(Hg_TPM_filt$Ensembl_gene_ID)
head(Hg_TPM_filt)


p_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSSSCG00000007978",
                                         "ENSSSCG00000014725")) -> Pg_TPM_filt
unique(Pg_TPM_filt$Ensembl_gene_ID)
head(Pg_TPM_filt)


c_TPM_filt %>% 
    dplyr::filter(Ensembl_gene_ID %in% c("ENSBTAG00000026417",
                                         "ENSBTAG00000026418",
                                         "ENSBTAG00000038748")) -> Cg_TPM_filt
unique(Cg_TPM_filt$Ensembl_gene_ID)
head(Cg_TPM_filt)

# Bind globin TPM filt dfs
allG_TPM_filt <- dplyr::bind_rows(Hg_TPM_filt, Pg_TPM_filt)
allG_TPM_filt <- dplyr::bind_rows(allG_TPM_filt, Cg_TPM_filt)
dim(allG_TPM_filt)
dim(Hg_TPM_filt)
dim(Pg_TPM_filt)
dim(Cg_TPM_filt)

# Add species column allG_TPM_filt
allG_TPM_filt$species <- allG_TPM_filt$Ensembl_gene_ID
allG_TPM_filt$species %<>%
    str_replace("ENSG0.+", "Human") %>% 
    str_replace("ENSSSCG0.+", "Porcine") %>% 
    str_replace("ENSBTAG0.+", "Bovine") 

# Add gene symbol allG_TPM_filt
allG_TPM_filt$gene_symbol <- allG_TPM_filt$Ensembl_gene_ID
allG_TPM_filt$gene_symbol %<>%
    str_replace("ENSG00000206172", "HBA1") %>% 
    str_replace("ENSG00000188536", "HBA2") %>% 
    str_replace("ENSG00000244734", "HBB") %>% 
    str_replace("ENSSSCG00000007978", "HBA1") %>% 
    str_replace("ENSSSCG00000014725", "HBB") %>%
    str_replace("ENSBTAG00000026417", "HBA1") %>% 
    str_replace("ENSBTAG00000026418", "HBA2") %>% 
    str_replace("ENSBTAG00000038748", "HBB")

head(allG_TPM_filt)

# Treatment factors allG_TPM_filt
allG_TPM_filt$treatment %<>% factor(levels = c("U", "D"))

allG_TPM_filt$treatment %<>%
    str_replace("U", "Undepleted") %>% 
    str_replace("D", "Globin depleted") %>%
    factor(levels = c("Undepleted", "Globin depleted"))

# Get means for each globin from filt dfs
h_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    summarise(mean(TPM)) %>%
    dplyr::filter(Ensembl_gene_ID %in% c("ENSG00000206172",
                                         "ENSG00000188536",
                                         "ENSG00000244734")) %>% 
    dplyr::rename(meanTPM = `mean(TPM)`) -> Hg_mean_filt

p_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    summarise(mean(TPM)) %>%
    dplyr::filter(Ensembl_gene_ID %in% c("ENSSSCG00000007978",
                                         "ENSSSCG00000014725")) %>% 
    dplyr::rename(meanTPM = `mean(TPM)`) -> Pg_mean_filt

c_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment) %>% 
    summarise(mean(TPM)) %>%
    dplyr::filter(Ensembl_gene_ID %in% c("ENSBTAG00000026417",
                                         "ENSBTAG00000026418",
                                         "ENSBTAG00000038748")) %>% 
    dplyr::rename(meanTPM = `mean(TPM)`) -> Cg_mean_filt

# Add mean TPM value for each globin within treatment
head(allG_TPM_filt)

allG_TPM_filt$meanTPM <- c(0)

head(allG_TPM_filt)


allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000188536" &
          allG_TPM_filt$treatment == "Globin depleted"), "meanTPM"] <- 76187.7564
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000188536" &
          allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 372700.7252
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000206172" &
          allG_TPM_filt$treatment == "Globin depleted"), "meanTPM"] <- 108598.0196
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000206172" &
          allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 345426.1538
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000244734" &
          allG_TPM_filt$treatment == "Globin depleted"), "meanTPM"] <- 21547.4487
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSG00000244734" &
              allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 107306.2230


allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSSSCG00000007978" &
          allG_TPM_filt$treatment == "Globin depleted"), "meanTPM"] <- 97326.7448
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSSSCG00000007978" &
              allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 337118.3320
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSSSCG00000014725" &
          allG_TPM_filt$treatment == "Globin depleted"), "meanTPM"] <- 134952.6925
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSSSCG00000014725" &
              allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 465716.4975


allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSBTAG00000026417" &
          allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 244.7277
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSBTAG00000026418" &
              allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 244.7441
allG_TPM_filt[(allG_TPM_filt$Ensembl_gene_ID == "ENSBTAG00000038748" &
              allG_TPM_filt$treatment == "Undepleted"), "meanTPM"] <- 394.2828


unique(allG_TPM_filt$meanTPM)
allG_TPM_filt %>% 
    dplyr::group_by(Ensembl_gene_ID, treatment, meanTPM) %>% 
    dplyr::summarise()



# Plot allG_TPM_filt (FIGURE 3)
allG_TPM_filt_plot <- ggplot(allG_TPM_filt) +
                        geom_jitter(aes(gene_symbol,
                                        log2(TPM + 1),
                                        colour = treatment),
                                    size = 3,
                                    alpha = 0.7) +
                        scale_colour_manual("Treatment",
                                            values = c("#af8dc3", "#7fbf7b")) +
                        scale_y_continuous(limits = c(0, 20)) +
                        geom_point(aes(gene_symbol,
                                       log2(meanTPM + 1),
                                       shape = treatment),
                                   colour = "black",
                                   size = 3) +
                        scale_shape_manual(expression(paste("Mean ", log[2], "(TPM + 1)")),
                        values = c(17, 15)) +
                        #labs(shape = expression(paste("Mean ", log[2], "(TPM + 1)"))) +
                        facet_grid(. ~ species) +
                        theme_bw() +
                        theme(axis.text.x = element_text(face = "italic")) +
                        xlab(NULL) +
                        ylab(expression(paste(log[2], "(TPM + 1)")))

allG_TPM_filt_plot

ggsave("allG_TPM_filt_plot.png",
       plot      = allG_TPM_filt_plot,
       limitsize = FALSE,
       dpi       = 300,
       path      = img_dir)
#
#guides(shape = guide_legend("Treatment"),
#       colour = guide_legend("Treatment")) +


# SD
sd(Table$Size)

save.image(file = "tximport-salmon-TPM.RData")


devtools::session_info()









