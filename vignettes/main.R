## main

## libraries
library(qtl)
library(RColorBrewer)
library(tidyverse)
library(gplots)
library(lme4)
library(multcomp)


devtools::load_all()

## ----parameters---------------------------------------------------------------

## Setting up phenotypes
phenos <- c("bw30", "bw40", "bw50", "bw60", "na", "nt", "femur", "ecw", "liver",
            "spleen", "kidney", "heart", "testis", "gfp", "rfp", "mfp", "ffp",
            "brain", "muscle", "tail", "tf", "g1", "g2", "g3", "pwg")

bw.phenos <- c("bw40", "bw50", "bw60", "g1", "g2", "g3", "pwg")
mm.phenos <- c("na", "nt", "femur", "tail")
bc.phenos <- c("ecw", "liver", "spleen", "kidney", "heart", "testis", "brain",
               "muscle", "gfp", "rfp", "mfp", "ffp", "tf")

## strain orders
strains <-  c("B6.A-15 BC1", "B6.A-17 BC1", "B6.A-19 BC1", "B6.A-X F1", "B6.C")

## population orders
popuation.names <- c("Discovery Population (N2)", "Replicate Backcross (N2)",
                     "N3")
population.years <- c("2012", "2013")

## tissues
tissues <- c("pituitary", "liver", "heart", "embryo", "placenta")
Tissues <- c("Pituitary", "Liver", "Heart", "Embryo", "Placenta")

## ---- color maps
