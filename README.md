# Gene-Set-Enrichment-Analysis-clusterProfiler-
Gene set enrichment analysis using clusterProfiler
setwd("/Users/ogunbawoadebisi/Downloads/meta_analyse")

#To install clusterProfiler package
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("clusterProfiler")
#To reproduce examples in this document, you need to install several extra packages:
install.packages(c("forcats", "ggplot2", "ggnewscale", "ggupset"))
BiocManager::install(c("org.Hs.eg.db", "enrichplot",
                       "ChIPseeker", "TxDb.Hsapiens.UCSC.hg19.knownGene"))
#Bioinformatics tools that depends on clusterProfiler for different topics,

db <- utils::available.packages(repo=BiocManager::repositories())
pkgs <- tools::package_dependencies('clusterProfiler', db=db,
                                    sort(pkgs)
  
library(AnnotationDbi)
library(clusterProfiler)
library(AnnotationHub)
library(readxl)
library("org.Hs.eg.db")
library(enrichplot)
## geneList for GSEA examples
data(geneList, package="DOSE")

geneList <-read_xlsx("Prioritised gene.xlsx")
ego <- enrichGO(geneList$GeneId, OrgDb = "org.Hs.eg.db", ont="BP", readable=TRUE)
ego
enrichbp <- ego@result
library(openxlsx)
write.xlsx(enrichbp, "enrichbp.xlsx")
# Extract enriched GO terms
enriched_term_bp <- as.data.frame(enrichbp)[, c("ID", "Description")]
# Print enriched GO terms with descriptions
print(enriched_term_bp)
# Retrieve phenotype information for enriched GO terms
phenotype_bp <- bitr(enriched_term_bp$ID, fromType = "GO", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")


#for MF 
ego2 <- enrichGO(geneList$GeneId, OrgDb = "org.Hs.eg.db", ont="MF", readable=TRUE)
ego2
enrichmf <- ego2@result
library(openxlsx)
write.xlsx(enrichmf, "enrichmf.xlsx")
# Extract enriched GO terms
enriched_term_mf <- as.data.frame(enrichmf)[, c("ID", "Description")]
# Print enriched GO terms with descriptions
print(enriched_term_mf)

#for CC
ego3 <- enrichGO(geneList$GeneId, OrgDb = "org.Hs.eg.db", ont="CC", readable=TRUE)
ego3
enrichcc <- ego3@result
library(openxlsx)
write.xlsx(enrichcc, "enrichcc.xlsx")
# Extract enriched GO terms
enriched_term_cc <- as.data.frame(enrichcc)[, c("ID", "Description")]
# Print enriched GO terms with descriptions
print(enriched_term_cc)

