################################################################################
##  This code is inteded to analyze differential expression from RNAseq mouse
## tissue data. Data was mapped to mm10 with star, then sorted with samtools
## sort -n, inexed with samtools index.  All files were analyzed with
## picard-tools CollectRnaseqMetrics.jar.  Gene counts were performed with
## HTSeq using -m insertion-strict -s 'reverse'
##
##  Each library is a pool of RNA from 4 mice.  Embryo and Placenta data was
## sequenced across two lanes and data was combined for this analysis within
## the DESeq2 pipeline.  FASTQ and BAM submission provide the complete data
## as a single file.
##
##  MGI symbols were pulled from Ensembl using biomaRt to the latest annotation.
## these may be different than those published, albeit the ensembl_gene_id
## should remain constant.
##
## Rodrigo GULARTE MERIDA
## Unit of Animal Genomics
## GIGA -- Research 
## Université de Liège
################################################################################

## ----setup--------------------------------------------------------------------
source("vignettes/main.R")

tmpDir <- "tmp/02.rnaseq_deg"
tabDir <- "tables/02.rnaseq_deg"

sapply(c(tmpDir, tabDir), usethis::use_directory)

cols <- brewer.pal(5, "Set1")[c(4:5,3,1,2)]
names(cols) <- strains

## ----libraries----------------------------------------------------------------
## library(GenomicFeatures)
library(DESeq2)
library(biomaRt)

## ----read_data_and_design-----------------------------------------------------
htcount <- read.csv(gzfile("data-raw/RNAseq_RawCounts.csv.gz"),
                        row.names = 1L)

cssDesign <- read.csv("data-raw/cssDesign.csv",
                      row.names = 1L)

cssDesign$strain <- factor(cssDesign$strain,
                           levels = c("B6.A-15BC1",
                                      "B6.A-17BC1", "B6.A-19BC1",
                                      "B6.A-XF1", "B6.CTRL"),
                           labels = strains,
                           ordered = FALSE)

## ----DESeq2_Pipeline----------------------------------------------------------
## Creating the class: DESeqDataSet objects for the complete raw count data
## from HTSeq
## ddsFCT : ddsFullCountTable : complete raw count data from HTSeq,
##  with metadata, and experimental design.  Formula/Design is `~ strain`, as
## it will be separated by tissue at a later step
## dds : collapsed data set.  This data adds both technical replicates
##  from embryo and placenta thus from 70 to 50 libraries as desired
ddsFCT <- DESeqDataSetFromMatrix(countData = htcount,
                                 colData = cssDesign,
                                 design = ~ strain)
ddsCR <- collapseReplicates(ddsFCT, groupby = ddsFCT$libraryName)

## RE-orders the levels of the strain factos, such that B6.CTRL is the first
## one. recomended in the beginner's guide to DESeq2
ddsCR$strain <- relevel(ddsCR$strain, "B6.C")
as.data.frame(colData(ddsCR))

## ----tissue_specific_objects--------------------------------------------------
## Separating the collapsed data set into individual tissues,
## Removing non-used levels from the tissue column
## Confirming that B6.CTRL is the first level as the first level is used
##  for comparison 
dds <- lapply(tissues, function(TSS) ddsCR[ , ddsCR$tissue == TSS ])
sapply(1:5, function(i) dds[[i]]$tissue <<- factor(dds[[i]]$tissue, tissues))
sapply(1:5, function(i) dds[[i]]$tissue <<- droplevels(dds[[i]]$tissue))
sapply(1:5, function(i) dds[[i]]$strain <<- relevel(dds[[i]]$strain, "B6.C"))
names(dds) <- tissues

## ----read_countnor_malization-------------------------------------------------
## Running the DESeq pipeline in one step : this function :
##  1. estimates size factors w/ `estimateSizeFactors`
##  2. estimates dispersions w/ `estimateDispertions`
##  3. runs GLM w/ `nbinomWaldTest`
invisible(sapply(tissues, function(i)
    dds[[i]] <<- DESeq(dds[[i]], fitType = "parametric")))
htcount.norm <- lapply(dds, counts, normalized = TRUE)

## ----tissue_specific_differential_gene_expression_analysis--------------------
## Getting results by tissue in a per strain vs CTRL basis
pituitary.de <- lapply(as.character(unique(strains))[1:4], function(STR)
    results(dds[["pituitary"]], contrast = c("strain", STR, "B6.C")))
liver.de <- lapply(as.character(unique(strains))[1:4], function(STR)
    results(dds[["liver"]], contrast = c("strain", STR, "B6.C")))
heart.de <- lapply(as.character(unique(strains))[1:4], function(STR)
    results(dds[["heart"]], contrast = c("strain", STR, "B6.C")))
embryo.de <- lapply(as.character(unique(strains))[1:4], function(STR)
    results(dds[["embryo"]], contrast = c("strain", STR, "B6.C")))
placenta.de <- lapply(as.character(unique(strains))[1:4], function(STR)
    results(dds[["placenta"]], contrast = c("strain", STR, "B6.C")))

## ----mgi_symbol_annotation----------------------------------------------------
## adding MGI gene symbols to each results objects
##  uses biomaRt package
## Note: thesaurus will change based to most recent annotation
gtf.genes <- sapply(strsplit(rownames(pituitary.de[[1]]),
                             split = "\\+"), "[", 1)
mmu <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
thesaurus <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
              filters = "ensembl_gene_id",
              values = gtf.genes,
              mart = mmu)
idx <- match(gtf.genes, thesaurus$ensembl_gene_id)


invisible(
sapply(1:4, function(s) {
  pituitary.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  liver.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  heart.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  embryo.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
  placenta.de[[s]]$mgi.symbol <<- thesaurus$mgi_symbol[idx]
})
)

## ----adding_comparison_tissue_and_id_metadata_to_each_object------------------
## merging all *.de to data.frame
for(i in 1:4) {
  pituitary.de[[i]]$comparison <-
      gsub(".*strain", "",
           attr(pituitary.de[[i]], "elementMetadata")[2,2])
  pituitary.de[[i]]$tissue <- "pituitary"
  pituitary.de[[i]]$id <- rownames(pituitary.de[[i]])
}
pituitary.df <- do.call(rbind, pituitary.de)
  
for(i in 1:4) {
  liver.de[[i]]$comparison <- gsub(".*strain", "",
                              attr(liver.de[[i]], "elementMetadata")[2,2])
  liver.de[[i]]$tissue <- "liver"
  liver.de[[i]]$id <- rownames(liver.de[[i]])
}
liver.df <- do.call(rbind, liver.de)


for(i in 1:4) {
  heart.de[[i]]$comparison <- gsub(".*strain", "",
                               attr(heart.de[[i]], "elementMetadata")[2,2])
  heart.de[[i]]$tissue <- "heart"
  heart.de[[i]]$id <- rownames(heart.de[[i]])
}
heart.df <- do.call(rbind, heart.de)

for(i in 1:4) {
  embryo.de[[i]]$comparison <- gsub(".*strain", "",
                              attr(embryo.de[[i]], "elementMetadata")[2,2])
  embryo.de[[i]]$tissue <- "embryo"
  embryo.de[[i]]$id <- rownames(embryo.de[[i]])
}
embryo.df <- do.call(rbind, embryo.de)


for(i in 1:4) {
  placenta.de[[i]]$comparison <- gsub(".*strain", "",
                           attr(placenta.de[[i]], "elementMetadata")[2,2])
  placenta.de[[i]]$tissue <- "placenta"
  placenta.de[[i]]$id <- rownames(placenta.de[[i]])
}
placenta.df <- do.call(rbind, placenta.de)

## ----combining_all_tissue_specific_analysis_to_a_single_df--------------------
## merge all
cssDE <- rbind(pituitary.df, liver.df, heart.df, embryo.df, placenta.df)
#% write.csv(cssDE, file = file.path(tabDir, "cssDE.csv"),
#%           quote = FALSE, row.names = FALSE)

## ----subset_of_significant_deg________________________________________________
## selected significant DEG manually 
#% sigDE <- read.csv("Significant_cssDE.csv")
sigDE <- cssDE[!is.na(cssDE$padj) & cssDE$padj <= 0.05, ] 

## save to data
usethis::use_data(cssDE, sigDE, overwrite = TRUE)

#% ## 
#% all.genes <- list(pituitary.df, liver.df, heart.df, embryo.df, placenta.df)
#% names(all.genes) <- tissues
#% 
#% ##
#% sig.genes <- lapply(tissues, function(i)
#%     as.character(sigDE$id[sigDE$tissue == i]))
#% names(sig.genes) <- tissues
#% 
#% ## 
#% sig.genes <- lapply(tissues, function(SG)
#%     all.genes[[SG]][as.character(all.genes[[SG]]$id) %in% sig.genes[[SG]]$id ])
#% sig.genes <- do.call(rbind, data.frame(sig.genes))
#% 
#% ## write.csv(sig.genes, "ForExp_Matrix.csv")

