#' =============================================================================
#'  The code in this file prepare the complete experimental backcross,
#' RNA-seq, and qPCR validation
#' 
#' Pre- and post- processed cross objects 
#' - css: original backcross of adult 60 day males + body composition
#' - em: embryo backcross, 13dpc embryo weight + placenta
#' - rbc: replicate backcross of adult 60 day males + body composition
#' - N3: persistence cross of adult 60d males + body composition
#'
#' RNA-seq
#' - htcount: raw counts RNA-seq
#' - cssDesign: experimental design for DESEq2
#'
#' q RT-PCR
#' expTR: qPCR Expression of Technical Replicates with BioMarkHD
#' expSS: qPCR expression of biological littermates
#' expVD: qPCR expression of biological littermates Validation 1
#' expVD2: qPCR expression of biological littermates Validation 2 (Mid1 B6.A-19)
#'
#' 
#' Rodrigo GULARTE MERIDA
#' Unit of Animal Genomics
#' GIGA -- Research
#' Université de Liège
#' Belgique
#' =============================================================================


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
## library(tgv_2023)
library(qtl)
library(tidyverse)

set.seed(2011)

## tmp files
tmpDir <- "tmp/00.dataset/"
tabDir <- "tables/00.dataset/"

sapply(list(tmpDir, tabDir), usethis::use_directory)

## ----error handling-----------------------------------------------------------
options(error = dump.frames("css", to.file = TRUE), width = 160,
        max.print = 1000)

## ----data---------------------------------------------------------------------
## ---- phenotype abreviations
abv <- read.delim("data-raw/abv.dat", row.names = 1)
usethis::use_data(abv, overwrite = TRUE)

## ---- full pedigree
#' tgv.pedigree <- read.pedigree(TRUE, path = "data-raw/pedigree.csv")
#' tgv.pedigree <- tgv.pedigree %>% get.wm.age %>% estimate.growth.rates %>% get.parental.age
tgv.pedigree <- read.pedigree(FALSE, path = "data-raw/pedigree.csv")  %>%
    get.wm.age %>%
    estimate.growth.rates %>%
    get.parental.age

## ----original_cross-----------------------------------------------------------
## check 8028 and 8029 real genotype, make a commit when it is corrected
pre.css <- read.cross("csvs", genfile = "data-raw/css_geno.csv",
                  phefile = "data-raw/css_pheno.csv",
                  genotypes = c("BB", "AB"), alleles = c("B", "A"),
                  na.strings = c("=", "-"))
## recode experimental generations
pre.css$pheno$generation = gsub("BC1", "N2", 
                                gsub("BC2", "N3", pre.css$pheno$generation))
pre.css$pheno$population = "Discovery Backcross (N2)"
pre.css$pheno$population.year = 2012
pre.css$pheno$cohort.name = "Cohort I (N2)"                       
    
css <- jittermap(pre.css)
pmap1 <- pull.map(css)


## ----genetic_map_cleanup------------------------------------------------------
## Removing animals excluded from the analysis
css.excl <- read.table("data-raw/excluded.animals.txt", header = TRUE, sep = "\t")
css.excl <- rbind(css.excl, read.table("data-raw/rbc.excluded.animals.txt",
                               header = TRUE, sep = "\t"))
css.excl$id <- paste("-", css.excl$id, sep = "")
css <- css[,css.excl$id]

## Remove markers with errors and markers outside the congenic region
dropped <- read.table("data-raw/bad.markers.txt",
                      sep = "\t", header = TRUE,
                      colClasses = rep("character", 2))
                      
dropped <- rbind(dropped, read.table("data-raw/em.bad.markers.txt",
                                     sep = "\t",  header = TRUE,
                                     colClasses = rep("character", 2)))

dropped <- unique(dropped)
css <- drop.markers(css, markers = dropped$marker)

## css Y markers
css <- drop.markers(css, markers = c("rs48098632", "rs46643293", "rs49888853"))

pmap2 <- pull.map(css)


## Plotting physical maps
pdf(file.path(tmpDir, "CSS_RBC_PhysMap.pdf"),
    width = 11.7, height = 8.3, paper = "a4r",
    pointsize = 20)
par(mar = c(4, 4, 1.5, 1)+.1)
plot.map(pmap1, pmap2, horizontal = TRUE,
  xlab = "Location (Mb)",
  main = "",  #"B6.A-(N) Markers Genotyped",
  cex.axis = 3, cex.lab = 4)
dev.off()

## ----CSS_subsets--------------------------------------------------------------
## Separating cross into three crosses, one for each chromosome.
( c15 <- css[, css$pheno$strain == "B6.A-15 BC1"] )
( c17 <- css[, css$pheno$strain == "B6.A-17 BC1"] )
( c19 <- css[, css$pheno$strain == "B6.A-19 BC1"] )
( cX <-  css[, css$pheno$strain == "B6.A-X F1"] )

## reading in list of sorted animals based on recombination end point
#sorted.mouse.ids <- readLines("data-raw/id.sorted.txt")[-1]
sorted.mouse.ids <- c(readLines("data-raw/id.sorted.txt")[-1],
                      readLines("data-raw/rbc.sorted.txt")[-1])

## ----genetic_map_processing---------------------------------------------------
## estimate recombination frequency  
cat("Estimating Recombination Frequency\n")
css <- est.rf(css)
c15 <- est.rf(c15)
c17 <- est.rf(c17)
c19 <- est.rf(c19)

# estimate linkage map
cat("Estimating Genetic Map\n")
kmap <- est.map(css, error.prob = 0.0001, map.function = "kosambi")
kmap15 <- est.map(c15, error.prob = 0.0001, map.function = "kosambi")
kmap17 <- est.map(c17, error.prob = 0.0001, map.function = "kosambi")
kmap19 <- est.map(c19, error.prob = 0.0001, map.function = "kosambi")
kmapX <- est.map(cX, error.prob = 0.0001, map.function = "kosambi")

pdf(file.path(tmpDir, "CSS_RBC_maps_phy+gen.pdf"))
plot.map(css, kmap)
plot.map(c15, kmap15, chr = 15)
plot.map(c17, kmap17, chr = 17)
plot.map(c19, kmap19, chr = 19)
plot.map(cX, kmapX, chr = "X")
dev.off()

# replace physical map with a genetic map
cat("replacing map\n")
css <- replace.map(css, kmap)
c15 <- replace.map(c15, kmap15)
c17 <- replace.map(c17, kmap17)
c19 <- replace.map(c19, kmap19)
cX <- replace.map(cX, kmapX)


## ----genotype_processing------------------------------------------------------
# calculate genotype probabilities to detect genotyping errors
cat("Calculating Error LOD\n")
css <- calc.errorlod(css, error.prob=0.001)
c15 <- calc.errorlod(c15, error.prob=0.001)
c17 <- calc.errorlod(c17, error.prob=0.001)
c19 <- calc.errorlod(c19, error.prob=0.001)
cX <- calc.errorlod(cX,  error.prob=0.001)

# calculate genotype probabilities, simulate and impute
n.draws = 500
cat("Calculating Genotype Probabilities\n")
css <- calc.genoprob(css, step = 0.01, error.prob=0.001, map.function = "kosambi")
css <- sim.geno(css, step=0.1, n.draws=n.draws, error.prob=0.001)
css <- fill.geno(css, method = "no_dbl_XO", map.function = "kosambi")
##
c15 <- calc.genoprob(c15, step = 0.01, error.prob=0.001, map.function = "kosambi")
c15 <- sim.geno(c15, step=0.1, n.draws=n.draws, error.prob=0.001)
c15 <- fill.geno(c15, method = "no_dbl_XO", map.function = "kosambi")
##
c17 <- calc.genoprob(c17, step = 0.01, error.prob=0.001, map.function = "kosambi")
c17 <- sim.geno(c17, step=0.1, n.draws=n.draws, error.prob=0.001)
c17 <- fill.geno(c17, method = "no_dbl_XO", map.function = "kosambi")
##
c19 <- calc.genoprob(c19, step = 0.01, error.prob=0.001, map.function = "kosambi")
c19 <- sim.geno(c19, step=0.1, n.draws=n.draws, error.prob=0.001)
c19 <- fill.geno(c19, method = "no_dbl_XO", map.function = "kosambi")
##
cX <- calc.genoprob(cX, step = 0.01, error.prob=0.001, map.function = "kosambi")
cX <- sim.geno(cX, step=0.1, n.draws=n.draws, error.prob=0.001)
cX <- fill.geno(cX, method = "no_dbl_XO", map.function = "kosambi")


## ----subset_non_recombinant_and_recombinant_mice------------------------------
## Recombinant class files are:
## B6/B6 isogenic: "B6_NR.txt" "emB6_NR.txt" "rbc_B6_NR.txt"
## B6/A non recombinants: "Ht_NR.txt" "emHt_NR.txt" "rbc_Ht_NR.txt" 
## B6/_ recombinants: "recombinants.txt"  "em_recombinants.txt"
##                    "rbc_recombinants.txt"
css.b6.nr <- read.table(
    "data-raw/B6_NR.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)
css.ht.nr <- read.table(
    "data-raw/Ht_NR.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)
css.recombinants <- read.table(
    "data-raw/recombinants.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)

## use %in% to keep the css order
## replace list with cross
css.b6.nr <- css[, as.character(css.b6.nr$id)]
css.ht.nr <- css[, as.character(css.ht.nr$id)]
css.recombinants <- css[, as.character(css.recombinants$id)]

png(file.path(tmpDir, "CSS_Genotype_Images_%03d.png"),
    width = 11.7, height = 8.3,
    units = "in", pointsize = 20, res = 500)
geno.image2(css[, sorted.mouse.ids])
geno.image2(css.b6.nr[, sorted.mouse.ids])
geno.image2(css.ht.nr[, sorted.mouse.ids])
geno.image2(css.recombinants[, sorted.mouse.ids])
dev.off()

write.cross(css.b6.nr, "csvs", filestem = file.path(tabDir, "CSS_B6_NR"))
write.cross(css.ht.nr, "csvs", filestem = file.path(tabDir, "CSS_Ht_NR"))
write.cross(css.recombinants, "csvs", filestem = file.path(tabDir, "CSS_REC"))

usethis::use_data(tgv.pedigree,
                  sorted.mouse.ids, 
                  pre.css, css, c15, c17, c19, cX,
                  overwrite = TRUE)
usethis::use_data(css.b6.nr, css.ht.nr, css.recombinants, overwrite = TRUE)

## ----Discovery_and_RBC_Cross_Summaries----------------------------------------
## Plot Cross Summary
pdf(file.path(tmpDir, "CSS_RBC_Cross_Summary.pdf"),
    width = 20, height = 20)
plot(css)
plot(c15, pheno.col = names(c15$pheno)[-grep("dam4.*", names(c15$pheno))])
plot(c17, pheno.col = names(c17$pheno)[-grep("dam4.*", names(c17$pheno))])
plot(c19, pheno.col = names(c19$pheno)[-grep("dam4.*", names(c19$pheno))])
plot(cX, pheno.col = names(cX$pheno)[-grep("dam4.*", names(cX$pheno))])
dev.off()


geno.maps <- read.delim("data-raw/genotype.panels.map.txt")
mm <- reshape2::melt(geno.maps, value.name = "true.false",
                     variable.name = "panel",
                     id.vars = c("id", "chr", "pos"))
mm$true.false[mm$true.false == FALSE] <- NA

panel.colors <- c("TRUE" = "#C42021", "FALSE" = "gray85")

pdf(file.path(tmpDir, "GenotypeMap_ggplot.pdf"),
    width = 85/25.4, height = 210/25.4)
ggplot(mm[mm$chr %in% c(15, 17, 19, "X", "Y"),],
       aes(panel, pos, fill = true.false)) +
    geom_tile(height = 0.6, width = 0.85) +
    scale_fill_manual(values = panel.colors) +
    facet_wrap(~ chr, ncol = 5, scales = "free_y") +
    xlab("Genotyping Panel") + ylab("Position (Mb)") +
    coord_cartesian(expand = FALSE) +
    theme_bw() +
    theme(legend.position = "bottom")
dev.off()

pdf(file.path(tmpDir, "Genotype_Map_all_Chromosomes_ggplot.pdf"),
    width = 85/25.4, height = 210/25.4)
mm %>%
    mutate(panel = factor(panel, c("P2", "P1"))) %>%
    ggplot(aes(pos, panel, fill = true.false)) +
    geom_tile(height = 0.4, width = 0.85) +
    scale_fill_manual(values = panel.colors) +
    facet_wrap(chr ~ ., ncol = 1) +
    xlab("Position (Mb)") +
    ylab("Genotyping Panel") +
    theme_bw() +
    theme(legend.position = "bottom") 
dev.off()

usethis::use_data(geno.maps, overwrite = TRUE)

## ----Embryo_Cross_Data--------------------------------------------------------
##  This code is ment to select the non recombinant mice and perform an
## analysis of variance of the B6 non recombinant mice from the different
## strains and the controls

## read cross data
## embryos 1..99 were test embryos, only use >= 100
pre.em <- read.cross("csvs", genfile = "data-raw/em_geno.csv",
                 phefile = "data-raw/em_pheno.csv",
                 genotypes = c("BB", "AB"),
                 alleles = c("B", "A"),
                 na.strings = c("=", "-", NA))
## recode experimental generations
pre.em$pheno$generation = "N2"
pre.em$pheno$population = "Embryo Backcross (N2)"
pre.em$pheno$population.year = 2012
pre.em$pheno$cohort.name = "Cohort II (N2)"
em <- jittermap(pre.em)

## exclude animals such as founders, F1 breeders, controls and animals with dubious genotyeps
em.excl <- read.table("data-raw/excluded.embryos.txt", header = TRUE, sep = "\t")
em.excl$id <- paste("-", em.excl$id, sep = "")
em <- em[, as.character(em.excl$id)]

## dropping markers with errors
## dropped <- read.table("data-raw/em.bad.markers.txt", sep = "\t",  header = TRUE, colClasses = rep("character", 2))
em <- drop.markers(em, markers = dropped$marker)
em <- drop.nullmarkers(em)

## selecting non-recombinant animals and generating a cross object from them
em.b6.nr <- read.table("data-raw/emB6_NR.txt", header = TRUE, sep = "\t", colClasses = c("numeric", "factor", "integer", "numeric"), fill = TRUE)
em.b6.nr <- em[, as.character(em.b6.nr$id)]

## selecting non-recombinants Heterozygous and generating a cross object from them
em.ht.nr <- read.table("data-raw/emHt_NR.txt", header = TRUE, sep = "\t", colClasses = c("numeric", "factor", "integer", "numeric"), fill = TRUE)
em.ht.nr <- em[,  as.character(em.ht.nr$id)]

## selecting recombinants and all other non-excluded, and generating a cross object from them
em.recombinants <- read.table("data-raw/em_recombinants.txt", header = TRUE, sep = "\t", colClasses = c("numeric", "factor", "integer", "numeric"), fill = TRUE)

if(! any(em.recombinants$id %in% em.b6.nr$pheno$id) &
   ! any(em.recombinants$id %in% em.ht.nr$pheno$id)) {
    em.recombinants <- em[, (em$pheno$id %in% em.recombinants$id)]
    message("ALL OK")
} else {
    warning("There are B6/B6 Non-Recombinants or B6/A individuals in the list of recombinant mice")
}


## reading in list of sorted animals
em.sorted.ids <- read.delim("data-raw/em.sorted.txt", colClasses = "character")
## em.sht <- readLines("data-raw/em.sorted.het.txt")

em.mf.sorted.ids <- em$pheno[, c("id", "sex")]
em.mf.sorted.ids$id <- factor(em.mf.sorted.ids$id, em.sorted.ids$id)
em.mf.sorted.ids <- em.mf.sorted.ids %>% arrange(sex, id) %>%
    mutate(id = as.character(id))

## calculating recombination frequencies
em <- est.rf(em)
em <- calc.errorlod(em)
em <- fill.geno(em, method = "no_dbl_XO", map.function = "kosambi")
em <- calc.genoprob(em, step = 0.05, map.function = "kosambi")

##
png(file.path(tmpDir, "EM_Genotype_image%03d.png"),
    width = 11.7, height = 8.3, units = "in",
    pointsize = 20, res = 300)
## main = "EM Genotypes for All Individuals")
geno.image2(em[, em.mf.sorted.ids$id], main = "")
## main = "EM Genotypes for b6/b6 Non-Recombinant Individuals")
geno.image2(em.b6.nr, main = "")
## main = "EM Genotypes for a/b6 Non-Recombinant Individuals")
geno.image2(em.ht.nr, main = "")
dev.off()

pdf(file.path(tmpDir, "EM_Genotype_image_pub.pdf"),
    width = 11.7, height = 8.3,
    paper = "a4r", pointsize = 20)
## main = "EM Genotypes for All Individuals")
geno.image2(em[, em.sorted.ids$id], main = "")
dev.off()

write.cross(em.b6.nr, "csvs", filestem = file.path(tabDir, "EM_B6_NR"))
write.cross(em.ht.nr, "csvs", filestem = file.path(tabDir, "EM_Ht_NR"))
write.cross(em.recombinants, "csvs", filestem = file.path(tabDir, "EM_REC"))

usethis::use_data(pre.em, em, em.b6.nr, em.ht.nr, em.recombinants,
                  em.mf.sorted.ids, overwrite = TRUE)

## ----Replicate_Back_Cross-----------------------------------------------------
pre.rbc <- read.cross("csvs", genfile = "data-raw/css_rbc_geno.csv",
                  phefile = "data-raw/css_pheno.csv", 
                  genotypes = c("BB", "AB"),
                  alleles = c("B", "A"),
                  na.strings = c("=", "-", NA))
## recode experimental generations
pre.rbc$pheno$generation = gsub("BC1", "N2",
                                gsub("BC2", "N3",pre.rbc$pheno$generation))
pre.rbc$pheno$population = "Replicate Backcross (N2)"
pre.rbc$pheno$population.year = 2013
pre.rbc$pheno$cohort.name = "Cohort III (N2)"
( rbc <- jittermap(pre.rbc) )

##
rbc.excl <- read.delim("data-raw/excluded.animals.txt")
rbc.excl$id <- paste("-", rbc.excl$id, sep = "")
( rbc <- rbc[, rbc.excl$id] )
##
rbc <- drop.markers(rbc, markers = dropped$marker)
( rbc <- drop.markers(rbc, markers = c("rs48098632", "rs46643293", "rs49888853")) )
##
( rbc <- drop.nullmarkers(rbc) )

rbc.b6.nr <- read.table(
    "data-raw/rbc_B6_NR.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)

rbc.ht.nr <- read.table(
    "data-raw/rbc_Ht_NR.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)

rbc.recombinants <- read.table(
    "data-raw/rbc_recombinants.txt",
    header = TRUE, sep = "\t",
    colClasses = c("numeric", "factor", "integer", "numeric"),
    fill = TRUE)


if(length(intersect(rbc.b6.nr$id, rbc.ht.nr)) == 0 &
   length(intersect(rbc.b6.nr$id, rbc.recombinants)) == 0 &
   length(intersect(rbc.ht.nr$id, rbc.recombinants)) == 0 ) {
    ## use %in% to keep the rbc order
    rbc.b6.nr <- rbc[, as.character(rbc.b6.nr$id)]
    rbc.ht.nr <- rbc[, as.character(rbc.ht.nr$id)]
    rbc.recombinants <- rbc[, as.character(rbc.recombinants$id)]
} else {
    warning("repeated ids exist in genotype class lists")
}

png(file.path(tmpDir, "RBC_Genotype_Images_%03d.png"),
    width = 11.7, height = 8.3,
    units = "in", pointsize = 20, res = 500)
geno.image2(rbc[, sorted.mouse.ids])
geno.image2(rbc.b6.nr[, sorted.mouse.ids])
geno.image2(rbc.ht.nr[, sorted.mouse.ids])
geno.image2(rbc.recombinants[, sorted.mouse.ids])
dev.off()

write.cross(rbc.b6.nr, "csvs", filestem = file.path(tabDir, "RBC_B6_NR"))
write.cross(rbc.ht.nr, "csvs", filestem = file.path(tabDir, "RBC_Ht_NR"))
write.cross(rbc.recombinants, "csvs", filestem = file.path(tabDir, "RBC_REC"))

usethis::use_data(pre.rbc, rbc, rbc.b6.nr, rbc.ht.nr, rbc.recombinants,
                  overwrite = TRUE)

## ----Persistance_Back_Cross(N3)-----------------------------------------------
pre.css.n3 <- read.cross("csvs", genfile = "../data/css_geno.csv",
                     phefile = "../data/css_N3_pheno.csv",
                     genotypes = c("BB", "AB"),
                     alleles = c("B", "A"), na.strings = c("=", "-", NA))
## recode experimental generations
pre.css.n3$pheno$generation <- gsub("BC1", "N2",
                                    gsub("BC2", "N3", pre.css.n3$pheno$generation))
pre.css.n3$pheno$population = "N3"
pre.css.n3$pheno$population.year = 2013
pre.css.n3$pheno$cohort.name = "Cohort IV (N3)"
css.n3 <- jittermap(pre.css.n3)

css.n3 <- drop.markers(css.n3, markers = dropped$marker)
css.n3 <- drop.nullmarkers(css.n3)

png(file.path(tmpDir, "CSS_N3_Genotype_Images_%03d.png"),
    width = 11.7, height = 8.3,
    units = "in", pointsize = 20, res = 500)
geno.image2(css.n3)
dev.off()

write.cross(css.n3, "csvs", filestem = file.path(tabDir, "CSS_N3"))

usethis::use_data(pre.css.n3, css.n3, overwrite = TRUE)


## ----RNAseq_cssDE-------------------------------------------------------------
htcount <- read.csv(gzfile("data-raw/RNAseq_RawCounts.csv.gz"), row.names = 1L)
cssDesign <- read.csv("data-raw/cssDesign.csv", row.names = 1L)

cssDesign$strain <- factor(cssDesign$strain,
                           levels = c("B6.A-15BC1",
                                      "B6.A-17BC1", "B6.A-19BC1",
                                      "B6.A-XF1", "B6.CTRL"),
                           labels = strains,
                           ordered = FALSE)

usethis::use_data(htcount, cssDesign, overwrite = TRUE)


## ----qRT_PCR_Validation-------------------------------------------------------
expTR <- read.csv(gzfile("data-raw/expression_technical_replicate.csv.gz"))
expSS <- read.csv(gzfile("data-raw/expression_strain_specific.csv.gz"))
expVD <- read.csv(gzfile("data-raw/expression_validation.csv.gz"))
expVD2 <- read.csv(gzfile("data-raw/expression_validation_2.csv.gz"))

usethis::use_data(expTR, expSS, expVD, expVD2, overwrite = TRUE)


## __EOF__
