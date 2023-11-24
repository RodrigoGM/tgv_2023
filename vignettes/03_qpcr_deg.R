### Summary Plots For Vd experession se

## libraries
library(ggplot2)
library(RColorBrewer)
library(lattice)
#library(GGally)

## ----data---------------------------------------------------------------------
data(expTR, expSS, expVD, expVD2)

##
expVD2$tissue = "Pituitary"

expTR$ppln2 <- expTR$ppln
expVD$ppln2 <- paste("VD", expVD$ppln, sep = ":")
expVD2$ppln2 <- paste("VD2", expVD2$ppln, sep = ":")
expVD2$strain <- gsub("BC1", "N2", expVD2$strain)

## adult tissues:
expPT <- rbind(expTR[expTR$tissue == "Pituitary", ],
               expVD[expVD$tissue == "Pituitary", ],
               expVD2[expVD2$tissue == "Pituitary", ])
expHT <- rbind(expTR[expTR$tissue == "Heart", ],
               expVD[expVD$tissue == "Heart", ],
               expVD2[expVD2$tissue == "Heart", ])
expLV <- rbind(expTR[expTR$tissue == "Liver", ],
               expVD[expVD$tissue == "Liver", ],
               expVD2[expVD2$tissue == "Liver", ])
## embryo tissues:
expEM <- rbind(expTR[expTR$tissue == "Embryo", ],
               expVD[expVD$tissue == "Embryo", ],
               expVD2[expVD2$tissue == "Embryo", ])
expPL <- rbind(expTR[expTR$tissue == "Placenta", ],
               expVD[expVD$tissue == "Placenta", ],
               expVD2[expVD2$tissue == "Placenta", ])

## extra data variables
tissues <- c("Heart", "Liver", "Pituitary", "Embryo", "Placenta")
hkg <- c("ENSMUSG00000001525", "ENSMUSG00000014767",
         "ENSMUSG00000025534", "ENSMUSG00000025630")

targets <- lapply(tissues, function(x) system(paste("cat ", paste(x, "targets", sep = "."), "| sed '/^#/d'"), intern = TRUE ))
names(targets) <- tissues

## Setting up factors to organize the plot
## expPT$tissue <- factor(expPT$tissue, levels = tissues)

expPT$cohort <- gsub(" .*", "", expPT$strain)
expPT$strainExt <- paste(expPT$strain, " (", expPT$ppln, ")", sep = "")
expPT$strainExt <- factor(expPT$strainExt,
                          levels = c("B6.A-15 N2 (TechRep)","B6.A-15 N2 (Bio)","B6.A-15 N3 (N3)",
                              "B6.A-17 N2 (TechRep)","B6.A-17 N2 (Bio)","B6.A-17 N2 (RBC)","B6.A-17 N3 (N3)",
                              "B6.A-19 N2 (TechRep)","B6.A-19 N2 (Bio)","B6.A-19 N2 (RBC)","B6.A-19 N3 (N3)",
                              "B6.A-X F1 (TechRep)","B6.A-X F1 (Bio)",
                              "B6.C (TechRep)","B6.C (Bio)","B6.C (RBC)","B6.C (N3)"))

expPT$ppln <- factor(expPT$ppln, levels = c("TechRep", "Bio", "RBC", "N3"))
expPT$ppln2 <- factor(expPT$ppln2, levels = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"))

expPT$strain <- factor(expPT$strain, levels = c("B6.A-15 N2", "B6.A-15 N3", "B6.A-17 N2", "B6.A-17 N3", "B6.A-19 N2", "B6.A-19 N3", "B6.A-X F1", "B6.C"))

# expPT$Detector <- factor(expPT$Detector, levels =  c("ENSMUSG00000035299-1", "ENSMUSG00000035299-2", "ENSMUSG00000035299-3", "ENSMUSG00000028957-1", "ENSMUSG00000028957-2", "ENSMUSG00000022389", "ENSMUSG00000034591", "ENSMUSG00000015357", "ENSMUSG00000020538", "ENSMUSG00000035686", "ENSMUSG00000038550", "ENSMUSG00000091971", "ENSMUSG00000090877", "ENSMUSG00000063889-1", "ENSMUSG00000063889-2", "ENSMUSG00000064380", "ENSMUSG00000001506", "ENSMUSG00000028763"))

expPT <- expPT[order(expPT$tissue, expPT$ppln, expPT$strain, expPT$Detector), ]

## Heart Data
expHT$cohort <- gsub(" .*", "", expHT$strain)
expHT$strainExt <- paste(expHT$strain, " (", expHT$ppln, ")", sep = "")
expHT$strainExt <- factor(expHT$strainExt,
                          levels = c("B6.A-15 N2 (TechRep)","B6.A-15 N2 (Bio)","B6.A-15 N3 (N3)",
                              "B6.A-17 N2 (TechRep)","B6.A-17 N2 (Bio)","B6.A-17 N2 (RBC)","B6.A-17 N3 (N3)",
                              "B6.A-19 N2 (TechRep)","B6.A-19 N2 (Bio)","B6.A-19 N2 (RBC)","B6.A-19 N3 (N3)",
                              "B6.A-X F1 (TechRep)","B6.A-X F1 (Bio)",
                              "B6.C (TechRep)","B6.C (Bio)","B6.C (RBC)","B6.C (N3)"))

expHT$ppln <- factor(expHT$ppln, levels = c("TechRep", "Bio", "RBC", "N3"))
expHT$ppln2 <- factor(expHT$ppln2, levels = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"))

expHT$strain <- factor(expHT$strain, levels = c("B6.A-15 N2", "B6.A-15 N3", "B6.A-17 N2", "B6.A-17 N3", "B6.A-19 N2", "B6.A-19 N3", "B6.A-X F1", "B6.C"))

##expHT$Detector <- factor(expHT$Detector, levels =  c("ENSMUSG00000035299-1", "ENSMUSG00000035299-2", "ENSMUSG00000035299-3", "ENSMUSG00000028957-1", "ENSMUSG00000028957-2", "ENSMUSG00000022389", "ENSMUSG00000034591", "ENSMUSG00000015357", "ENSMUSG00000020538", "ENSMUSG00000035686", "ENSMUSG00000038550", "ENSMUSG00000091971", "ENSMUSG00000090877", "ENSMUSG00000063889-1", "ENSMUSG00000063889-2", "ENSMUSG00000064380", "ENSMUSG00000001506", "ENSMUSG00000028763"))

expHT <- expHT[order(expHT$tissue, expHT$ppln, expHT$strain, expHT$Detector), ]

## Liver Data
expLV$cohort <- gsub(" .*", "", expLV$strain)
expLV$strainExt <- paste(expLV$strain, " (", expLV$ppln, ")", sep = "")
expLV$strainExt <- factor(expLV$strainExt,
                          levels = c("B6.A-15 N2 (TechRep)","B6.A-15 N2 (Bio)","B6.A-15 N3 (N3)",
                              "B6.A-17 N2 (TechRep)","B6.A-17 N2 (Bio)","B6.A-17 N2 (RBC)","B6.A-17 N3 (N3)",
                              "B6.A-19 N2 (TechRep)","B6.A-19 N2 (Bio)","B6.A-19 N2 (RBC)","B6.A-19 N3 (N3)",
                              "B6.A-X F1 (TechRep)","B6.A-X F1 (Bio)",
                              "B6.C (TechRep)","B6.C (Bio)","B6.C (RBC)","B6.C (N3)"))

expLV$ppln <- factor(expLV$ppln, levels = c("TechRep", "Bio", "RBC", "N3"))
expLV$ppln2 <- factor(expLV$ppln2, levels = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"))

expLV$strain <- factor(expLV$strain, levels = c("B6.A-15 N2", "B6.A-15 N3", "B6.A-17 N2", "B6.A-17 N3", "B6.A-19 N2", "B6.A-19 N3", "B6.A-X F1", "B6.C"))

## expLV$Detector <- factor(expLV$Detector, levels =  c("ENSMUSG00000035299-1", "ENSMUSG00000035299-2", "ENSMUSG00000035299-3", "ENSMUSG00000028957-1", "ENSMUSG00000028957-2", "ENSMUSG00000022389", "ENSMUSG00000034591", "ENSMUSG00000015357", "ENSMUSG00000020538", "ENSMUSG00000035686", "ENSMUSG00000038550", "ENSMUSG00000091971", "ENSMUSG00000090877", "ENSMUSG00000063889-1", "ENSMUSG00000063889-2", "ENSMUSG00000064380", "ENSMUSG00000001506", "ENSMUSG00000028763"))

expLV <- expLV[order(expLV$tissue, expLV$ppln, expLV$strain, expLV$Detector), ]

#save(list = c(ls(pattern = "exp.*")), file = "TR+VD+VD2.ddCt.data.rda")


## color selection
## load("TR+VD+VD2.ddCt.data.rda")
strainColors <- brewer.pal(5, "Set1")[c(4, 4, 5, 5, 3, 3, 1, 2)]
names(strainColors) <- c("B6.A-15 N2", "B6.A-15 N3", "B6.A-17 N2", "B6.A-17 N3", "B6.A-19 N2", "B6.A-19 N3", "B6.A-X F1", "B6.C")
### strainColors <- c("B6.A-15 N2" = strainColors[1], "B6.A-15 N3" = strainColors[1], "B6.A-17 N2" = strainColors[3], "B6.A-17 N3" = strainColors[3], "B6.A-19 N2" = strainColors[5], "B6.A-19 N3" = strainColors[5], "B6.A-X F1" = strainColors[7], "B6.C" = strainColors[8])
##strainColScale <- scale_color_manual(values = strainColors)


give.n <- function(x){
    return(c(y = min(x)-0.5, label = length(x)))
}

ggExprs <- ggplot(expSS, aes(strain, exprs))
ggPT <- ggplot(expPT, aes(strain, exprs))

## pdf("TR+VD+VD2_BoxPlot_Exprs_A0.pdf", 1189/25.4, 841/25.4)
ggExprs + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + strainColScale +
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
## dev.off()


## pdf("TR+VD+VD2_BoxPlot_PT_Exprs_A0.pdf", 1189/25.4, 841/25.4)
ggPT + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + strainColScale +
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
## dev.off()


## Removing population outliers
## id     reason
## 5869   [RNA] 6.46 ng/uL
## Top outliers 
## 1252 452 5617 5646 5966 6028 6080 6107 6127 6159 6166 6178 6280 6311 6313 6350 6355 6363 6410 
## 5638 5638 6379 1292 5182 6302 5187 5182 5201 5273 5638 5184 6379 5203 5638 5187 5187 6379 5187 5273

expPT$exprs[expPT$id == 5869] <- NA
expPT$exprs[expPT$id == 7172] <- NA
expPT$exprs[expPT$id == 7171] <- NA
expPT$exprs[expPT$id == 7274] <- NA
expPT$exprs[expPT$id == 7194] <- NA
expPT$exprs[expPT$id == 7195] <- NA

expPT$exprs[expPT$id == 5273 & expPT$Detector == "ENSMUSG00000028763"] <- NA
expPT$exprs[expPT$id %in% c(493, 5160, 5929) & expPT$Detector %in% c("ENSMUSG00000090877", "ENSMUSG00000091971")] <- NA
expPT$exprs[expPT$id %in% c(6236, 7584, 7857, 5930, 6102, 7856) & expPT$Detector %in% c("ENSMUSG00000063889-2")] <- NA
expPT$exprs[expPT$id %in% c(5664, 6402) & expPT$Detector %in% c("ENSMUSG00000063889-2")] <- NA

ggPT <- ggplot(expPT, aes(strain, exprs, color = strain, fill = strain))

ggPT + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + strainColScale +
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")


## Mid1
mid1PT <- expPT[ expPT$Detector %in% c("ENSMUSG00000035299-1", "ENSMUSG00000035299-2", "ENSMUSG00000035299-3") & expPT$strain %in% c("B6.A-19 N2", "B6.A-19 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 12, 49))
colnames(nullMatrix) <- names(expPT)
nullMatrix$strain <- rep(c("B6.A-19 N2", "B6.C"), each = 6)
nullMatrix$Detector <- rep(c("ENSMUSG00000035299-1", "ENSMUSG00000035299-2", "ENSMUSG00000035299-3"), 4)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 6)
nullMatrix$exprs <- 0

mid1PT <- rbind(mid1PT, nullMatrix)

ggMid1 <- ggplot(mid1PT, aes(strain, exprs, fill = strain))
Mid1.exprs.pvals <- data.frame(x = c(1.5, 1.5), y = c(5, 5), strain = c(NA, NA), ppln2 = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"), pval = c("p < 0.0001", "p < 0.0001", "p < 0.0001", "p < 0.0001", "NA", "NA"))

pdf("Mid1_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggMid1 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) ## + strainColScale +
        xlab("Isogenic Provenance (Cohort)") + ylab(label = expression(Expression~~2^(-~Delta~Delta~Ct))) +
        stat_summary(fun.data = give.n, geom = "text", col = "black") +
        geom_text(aes(x, y, label = pval), dat = Mid1.exprs.pvals, col = "black") +
        geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge", alpha = 2/10, dotsize = .2)
dev.off()
##        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  

## Crem
cremPT <- expPT[expPT$Detector %in% c("ENSMUSG00000063889-1", "ENSMUSG00000063889-2") & expPT$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]
## cremPT$exprs[cremPT$exprs > 20] <- NA  ## outlier
nullMatrix <- data.frame(matrix(NA, 6, 49))
colnames(nullMatrix) <- names(expPT)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 3)
nullMatrix$Detector <- rep("ENSMUSG00000063889-0", 6)
nullMatrix$ppln2 <- rep(c("VD:RBC", "VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

cremPT <- rbind(cremPT, nullMatrix)

ggCrem <- ggplot(cremPT, aes(strain, exprs, color = strain, fill = strain))
Crem.exprs.pvals <- data.frame(x = c(1.5, 1.5), y = c(3, 3), strain = c(NA, NA), ppln2 = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"), pval = c("p = 0.0283", "p < 0.0001", "NA", "p > 0.1", "p = 0.0627", "NA"))

pdf("Crem_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggCrem + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) +
        stat_summary(fun.data = give.n, geom = "text", col = "black") +
        geom_dotplot(colour = "black", binaxis = "y", stackdir = "center", position = "dodge", alpha = 4/10, dotsize = .4) +
        geom_text(aes(x, y, label = pval), dat = Crem.exprs.pvals, col = "black")
dev.off()



##        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  

## SnoRNA Gm26449 - ENSMUSG00000064380

gmPT <- expPT[expPT$Detector == "ENSMUSG00000064380" & expPT$strain %in% c("B6.A-17 N2", "B6.A-17 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 6, 49))
colnames(nullMatrix) <- names(expPT)
nullMatrix$strain <- rep(c("B6.A-17 N2", "B6.C"), each = 3)
nullMatrix$Detector <- rep(c("ENSMUSG00000064380-0", "ENSMUSG00000064380-9"), 3)
nullMatrix$ppln2 <- rep(c("VD:RBC", "VD2:RBC"), 3)
nullMatrix$exprs <- 0

gmPT <- rbind(gmPT, nullMatrix)



ggGm <- ggplot(gmPT, aes(strain, exprs, fill = strain))
Gm.exprs.pvals <- data.frame(x = c(1.5, 1.5), y = c(5, 5), strain = c(NA, NA), ppln2 = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"), pval = c("p = 0.0232", "p = 0.0046", "p > 0.1 ", "p > 0.1", "p > 0.1", "p > 0.1"))

pdf("Gm26448_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggGm + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
        scale_fill_manual(values = strainColors) + 
        xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) +
        stat_summary(fun.data = give.n, geom = "text", col = "black") + 
        geom_text(aes(x, y, label = pval), dat = Gm.exprs.pvals, col = "black")  +
        geom_dotplot(colour = "black", binaxis = "y", stackdir = "center", position = "dodge", alpha = 4/10, dotsize = .8)
dev.off()
##        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  


## Per3 -- ENSMUSG00000028957-1 ENSMUSG00000028957-2
perHT <- expHT[expHT$Detector %in% c("ENSMUSG00000028957-1", "ENSMUSG00000028957-2") & expHT$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 6, 49))
colnames(nullMatrix) <- names(expHT)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 3)
nullMatrix$Detector <- rep(c("ENSMUSG00000028957-0"), 6)
nullMatrix$ppln2 <- rep(c("VD:RCB", "VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

perHT <- rbind(perHT, nullMatrix)

ggPer3 <- ggplot(perHT, aes(strain, exprs, color = strain, fill = strain))
Per3.exprs.pvals <- data.frame(x = c(1.5, 1.5), y = c(2.2, 2.2), strain = c(NA, NA), ppln2 = c("TechRep", "VD:Bio", "VD:RBC", "VD:N3", "VD2:Bio", "VD2:RBC"), pval = c("p = 0.0232", "p = 0.0046", "p > 0.1 ", "p > 0.1", "p > 0.1", "p > 0.1"))

pdf("Per3_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggPer3 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) +
            stat_summary(fun.data = give.n, geom = "text", col = "black") +
                geom_text(aes(x, y, label = pval), dat = Per3.exprs.pvals, col = "black")
dev.off()
#        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  


## Tef -- ENSMUSG00000022389
tefHT <- expHT[expHT$Detector %in% c("ENSMUSG00000022389") & expHT$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expHT)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

tefHT <- rbind(tefHT, nullMatrix)

ggTef3 <- ggplot(tefHT, aes(strain, exprs, color = strain, fill = strain))

pdf("Tef3_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggTef3 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Slc41a2 -- ENSMUSG00000034591
slcLV <- expLV[expLV$Detector %in% c("ENSMUSG00000034591") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

slcLV <- rbind(slcLV, nullMatrix)

ggSlc41 <- ggplot(slcLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Slc41a2_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggSlc41 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Clpx -- ENSMUSG00000015357
clpxLV <- expLV[expLV$Detector %in% c("ENSMUSG00000015357") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

clpxLV <- rbind(clpxLV, nullMatrix)

ggClpx <- ggplot(clpxLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Clpx_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggClpx + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()

## Srebf1 -- ENSMUSG00000020538
srebLV <- expLV[expLV$Detector %in% c("ENSMUSG00000020538") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

srebLV <- rbind(srebLV, nullMatrix)

ggSreb <- ggplot(srebLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Srebf1_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggSreb + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Thrsp -- ENSMUSG00000035686
thrLV <- expLV[expLV$Detector %in% c("ENSMUSG00000035686") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

thrLV <- rbind(thrLV, nullMatrix)

ggThr <- ggplot(thrLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Thrsp_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggThr + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Gm129 -- ENSMUSG00000038550
gm129LV <- expLV[expLV$Detector %in% c("ENSMUSG00000038550") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

gm129LV <- rbind(gm129LV, nullMatrix)

ggGm129 <- ggplot(gm129LV, aes(strain, exprs, color = strain, fill = strain))

pdf("Gm129_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggGm129 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Hspa1a -- ENSMUSG00000091971
hspaLV <- expLV[expLV$Detector %in% c("ENSMUSG00000091971") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

hspaLV <- rbind(hspaLV, nullMatrix)

ggHspa <- ggplot(hspaLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Hspa1a_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggHspa + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Hspa1b -- ENSMUSG00000090877
hspbLV <- expLV[expLV$Detector %in% c("ENSMUSG00000090877") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

hspbLV <- rbind(hspbLV, nullMatrix)

ggHspb <- ggplot(hspbLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Hspa1b_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggHspb + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Col1a1 -- ENSMUSG00000001506
col1LV <- expLV[expLV$Detector %in% c("ENSMUSG00000001506") & expLV$strain %in% c("B6.A-19 N2", "B6.A-19 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-19 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

col1LV <- rbind(col1LV, nullMatrix)

ggCol1 <- ggplot(col1LV, aes(strain, exprs, color = strain, fill = strain))

pdf("Col1a1_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggCol1 + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()


## Hspg2 -- ENSMUSG00000028763
hspgLV <- expLV[expLV$Detector %in% c("ENSMUSG00000028763") & expLV$strain %in% c("B6.A-15 N2", "B6.A-15 N3", "B6.C"), ]

nullMatrix <- data.frame(matrix(NA, 4, 49))
colnames(nullMatrix) <- names(expLV)
nullMatrix$strain <- rep(c("B6.A-15 N2", "B6.C"), each = 2)
nullMatrix$Detector <- rep(c("BLANK1", "BLANK2"), 2)
nullMatrix$ppln2 <- rep(c("VD2:Bio", "VD2:RBC"), 2)
nullMatrix$exprs <- 0

hspgLV <- rbind(hspgLV, nullMatrix)

ggHspg <- ggplot(hspgLV, aes(strain, exprs, color = strain, fill = strain))

pdf("Hspg2_TR_VD_VD2.pdf", width = 315/25.4, height = 210/25.4)
ggHspg + geom_boxplot(color = "black") + facet_grid(Detector ~ ppln2, scales = "free") +
    scale_fill_manual(values = strainColors) + strainColScale +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  
            xlab("Isogenic Provenance (Cohort)") + ylab(expression(Expression~~2^(-~Delta~Delta~Ct))) + stat_summary(fun.data = give.n, geom = "text", col = "black")
dev.off()

