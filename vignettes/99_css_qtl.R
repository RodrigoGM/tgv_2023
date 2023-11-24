## =============================================================================
##   The code contained in this file is to perform a QTL analysis of Body Weight,
## body composition
##
## Rodrigo GULARTE MERIDA
## GIGA - Research
## Universite de Liege
## Aveneu de l'Hopital 1 (B34)
## 4031 Sart-Tilman, Liège
## Belgique
## =============================================================================

## ----error_handling-----------------------------------------------------------
options(error = dump.frames("css", to.file = TRUE), width = 160,
        max.print = 1000)

## ----libraries----------------------------------------------------------------
##library(tgv_2023)
devtools::load_all()

## ----parameters---------------------------------------------------------------
force.preproc <- FALSE
## force.analysis <- FALSE
run.permutations <- TRUE

## Setting up phenotypes
phenos <- c("bw30", "bw40", "bw50", "bw60", "na", "nt", "femur", "ecw", "liver",
            "spleen", "kidney", "heart", "testis", "gfp", "rfp", "mfp", "ffp",
            "brain", "muscle", "tail", "tf", "g1", "g2", "g3", "pwg")

bw.phenos <- c("bw40", "bw50", "bw60", "g1", "g2", "g3", "pwg")
mm.phenos <- c("na", "nt", "femur", "tail")
bc.phenos <- c("ecw", "liver", "spleen", "kidney", "heart", "testis", "brain",
               "muscle", "gfp", "rfp", "mfp", "ffp", "tf")

## ----setup_covariates---------------------------------------------------------
BW30.15 <- c15$pheno$bw30
BW30.17 <- c17$pheno$bw30
BW30.19 <- c19$pheno$bw30

SAC.15 <- c15$pheno$bw60
SAC.17 <- c17$pheno$bw60
SAC.19 <- c19$pheno$bw60

## ----perform_linkage_analyses-------------------------------------------------
b30.15 <- scanone(c15, pheno.col = bw.phenos, method = "em")
b30.17 <- scanone(c17, pheno.col = bw.phenos, method = "em")
b30.19 <- scanone(c19, pheno.col = bw.phenos, method = "em")
##
bw.15 <-  scanone(c15, pheno.col = bw.phenos, method = "em", addcovar = BW30.15)
bw.17 <-  scanone(c17, pheno.col = bw.phenos, method = "em", addcovar = BW30.17)
bw.19 <-  scanone(c19, pheno.col = bw.phenos, method = "em", addcovar = BW30.19)
##
mm.15 <- scanone(c15, pheno.col = mm.phenos, method = "em", addcovar = SAC.15)
mm.17 <- scanone(c17, pheno.col = mm.phenos, method = "em", addcovar = SAC.17)
mm.19 <- scanone(c19, pheno.col = mm.phenos, method = "em", addcovar = SAC.19)
##
bc.15 <- scanone(c15, pheno.col = bc.phenos, method = "em", addcovar = SAC.15)
bc.17 <- scanone(c17, pheno.col = bc.phenos, method = "em", addcovar = SAC.17)
bc.19 <- scanone(c19, pheno.col = bc.phenos, method = "em", addcovar = SAC.19)

## ----permutations-------------------------------------------------------------
if(run.permutations) {
    nc <- 8
    bw.15.p <-  scanone(c15, pheno.col = bw.phenos, method = "em",
                        addcovar = BW30.15, n.perm = 1000, n.cluster = nc)
    bw.17.p <-  scanone(c17, pheno.col = bw.phenos, method = "em",
                        addcovar = BW30.17, n.perm = 1000, n.cluster = nc)
    bw.19.p <-  scanone(c19, pheno.col = bw.phenos, method = "em",
                        addcovar = BW30.19, n.perm = 1000, n.cluster = nc)
    ##
    mm.15.p <- scanone(c15, pheno.col = mm.phenos, method = "em",
                       addcovar = SAC.15, n.perm = 1000, n.cluster = nc)
    mm.17.p <- scanone(c17, pheno.col = mm.phenos, method = "em",
                       addcovar = SAC.17, n.perm = 1000, n.cluster = nc)
    mm.19.p <- scanone(c19, pheno.col = mm.phenos, method = "em",
                       addcovar = SAC.19, n.perm = 1000, n.cluster = nc)
    ##
    bc.15.p <- scanone(c15, pheno.col = bc.phenos, method = "em",
                       addcovar = SAC.15, n.perm = 1000, n.cluster = nc)
    bc.17.p <- scanone(c17, pheno.col = bc.phenos, method = "em",
                       addcovar = SAC.17, n.perm = 1000, n.cluster = nc)
    bc.19.p <- scanone(c19, pheno.col = bc.phenos, method = "em",
                       addcovar = SAC.19, n.perm = 1000, n.cluster = nc)
    ##
    usethis::use_data(bw.15.p, bw.17.p, bw.19.p,
                      mm.15.p, mm.17.p, mm.19.p,
                      bc.15.p, bc.17.p, bc.19.p,
                      overwrite = TRUE)
} else {
    print("permutations not ran, loading last savedrun")
    data(bw.15.p, bw.17.p, bw.19.p,
         mm.15.p, mm.17.p, mm.19.p,
         bc.15.p, bc.17.p, bc.19.p)
}


## ----preliminary_plots--------------------------------------------------------
abv <- read.delim("data-raw/abv.dat", row.names = 1)
cols = brewer.pal(9, "Set1")

## ----permutation_thresholds---------------------------------------------------
qtl.plot.min3 <- function(scanone, lodcolumn, ...) {
    ylim <- c(0, 
              ifelse(max(lodcolumn) < 4, 4, max(lodcolumn)))
    plot(scanone, lodcolumn = lodcolumn,
         ylim = ylim,
         xaxs = "i", yaxs = "i", ...)
}


## ----body weight phenotypes
layout(mat = matrix(c(1:7,0), nrow = 2, byrow = TRUE))
## chr 15
sapply(bw.phenos, function(PHEN) {
    qtl.plot.min3(bw.15,
                  lodcolumn = grep(paste0("^", PHEN, "$"), bw.phenos),
                  chr = 15, xlab = "Chr 15 : Map Position (cM)")
})
## chr17
sapply(bw.phenos, function(PHEN)
    qtl.plot.min3(bw.17,
                  lodcolumn = grep(paste0("^", PHEN, "$"), bw.phenos),
                  chr = 17, xlab = "Chr 17 : Map Position (cM)"))
## chr19
sapply(bw.phenos, function(PHEN)
    qtl.plot.min3(bw.19,
                  lodcolumn = grep(paste0("^", PHEN, "$"), bw.phenos),
                  chr = 19, xlab = "Chr 19 : Map Position (cM)"))

## ----length (mm) phenotypes
layout(mat = matrix(c(1:4, rep(0, 4)), nrow = 2, byrow = TRUE))
## chr 15
lapply(mm.phenos, function(PHEN)
    qtl.plot.min3(mm.15,
                  lodcolumn = grep(paste0("^", PHEN, "$"), mm.phenos),
                  chr = 15, xlab = "Chr 15 : Map Position (cM)"))
## chr 17
lapply(mm.phenos, function(PHEN)
    qtl.plot.min3(mm.17,
                  lodcolumn = grep(paste0("^", PHEN, "$"), mm.phenos),
                  chr = 17, xlab = "Chr 17 : Map Position (cM)"))
## chr 19
lapply(mm.phenos, function(PHEN)
    qtl.plot.min3(mm.19,
                  lodcolumn = grep(paste0("^", PHEN, "$"), mm.phenos),
                  chr = 19, xlab = "Chr 19 : Map Position (cM)"))

## ----body composition phenotypes
layout(mat = matrix(c(1:13, 0, 0), nrow = 3, byrow = TRUE))
## chr 15
lapply(bc.phenos, function(PHEN)
    qtl.plot.min3(bc.15,
                  lodcolumn = grep(paste0("^", PHEN, "$"), bc.phenos),
                  chr = 15, xlab = "Chr 15 : Map Position (cM)"))
## chr 17
lapply(bc.phenos, function(PHEN)
    qtl.plot.min3(bc.17, lodcolumn = grep(paste0("^", PHEN, "$"), bc.phenos),
         chr = 17, xlab = "Chr 17 : Map Position (cM)"))
## chr 19
lapply(bc.phenos, function(PHEN)
    qtl.plot.min3(bc.19, lodcolumn = grep(paste0("^", PHEN, "$"), bc.phenos),
         chr = 19, xlab = "Chr 19 : Map Position (cM)"))


## epistatic QTL models
n.draws = 500
cat("Calculating Genotype Probabilities\n")
c15 <- sim.geno(c15, step=0.1, n.draws=n.draws, error.prob=0.001)
c17 <- sim.geno(c17, step=0.1, n.draws=n.draws, error.prob=0.001)
c19 <- sim.geno(c19, step=0.1, n.draws=n.draws, error.prob=0.001)
cX <- sim.geno(cX, step=0.1, n.draws=n.draws, error.prob=0.001)


## multiple QTL analysis
bw30.stq <- lapply(list(c15, c17, c19), stepwiseqtl, pheno.col = bw.phenos, max.qtl = 5, method = "imp", keeplodprofile = TRUE, keeptrace = TRUE)
##
mm.stq <- lapply(list(c15, c17, c19), stepwiseqtl, pheno.col = mm.phenos, max.qtl = 5, method = "imp", keeplodprofile = TRUE, keeptrace = TRUE)
##
bc.stq <- lapply(list(c15, c17, c19), stepwiseqtl, pheno.col = bc.phenos, max.qtl = 5, method = "imp", keeplodprofile = TRUE, keeptrace = TRUE)

usethis::use_data(c15, c17, c19, overwrite = TRUE)
usethis::use_data(bw30.stq, mm.stq, bc.stq, overwrite = TRUE)

## __EOF__


