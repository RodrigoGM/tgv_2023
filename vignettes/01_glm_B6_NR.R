## =============================================================================
##  This code is to compare body weight and body composition phenotypes
## from non-recombinant C57BL/6J mice derived from the different
## backcross populations.
##
## Rodrigo GULARTE MERIDA
## GIGA - Research
## Universite de Liege
## Aveneu de l'Hopital 1 (B34)
## 4031 Sart-Tilman, Liège
## Belgique
## =============================================================================


## ----libraries----------------------------------------------------------------
source("vignettes/main.R")

## ----parameters---------------------------------------------------------------
tmpDir <- "tmp/01.glm_B6_NR"
tabDir <- "tables/01.glm_B6_NR"
sapply(c(tmpDir, tabDir), usethis::use_directory)

## ----pull_phenotype_data_from_B6.B6_non_recombinants--------------------------
b6.nr <- list(
    "css" = css.b6.nr$pheno,
    "rbc" = rbc.b6.nr$pheno,
    "n3" = css.n3$pheno
)

strain.relevel <- c("B6.C", "B6.A-15 BC1", "B6.A-17 BC1", "B6.A-19 BC1", "B6.A-X F1")
b6.nr$css$strain <- factor(b6.nr$css$strain, strain.relevel)
b6.nr$rbc$strain <- factor(b6.nr$rbc$strain, strain.relevel)
b6.nr$n3$strain <- factor(b6.nr$n3$strain, strain.relevel)

## -----------------------------------------------------------------------------

## mean and sd
sapply(strains, function(STR)
    apply(b6.nr$css[b6.nr$css$strain == STR, phenos], 2, mean, na.rm = TRUE))
sapply(strains, function(STR)
    apply(b6.nr$css[b6.nr$css$strain == STR, phenos], 2, sd, na.rm = TRUE))
sapply(strains, function(STR)
    apply(b6.nr$rbc[b6.nr$rbc$strain == STR, phenos], 2, mean, na.rm = TRUE))
sapply(strains, function(STR)
    apply(b6.nr$rbc[b6.nr$rbc$strain == STR, phenos], 2, sd, na.rm = TRUE))

                                        # mean and se
options(digits = 2)
mse <- lapply(strains, function(STR) {
    do.call(rbind, (
        apply(b6.nr$css[b6.nr$css$strain == STR, phenos], 2, mean_se, 1.96)
    ))
})
names(mse) <- strains

sapply(strains, function(x) {
    mse[[x]]$se <<- mse[[x]]$y - mse[[x]]$ymin
    paste(format(mse[[x]]$y, digits = 2),
          format(mse[[x]]$se, digits = 2), sep = "+-")
})

rmse <- lapply(strains, function(STR) {
    do.call(rbind, (
        apply(b6.nr$rbc[b6.nr$rbc$strain == STR, phenos], 2, mean_se, 1.96)
    ))
})

names(rmse) <- strains
sapply(strains, function(x) {
    rmse[[x]]$se <<- rmse[[x]]$y - rmse[[x]]$ymin
    paste(format(rmse[[x]]$y, digits = 2),
          format(rmse[[x]]$se, digits = 2), sep = "+-")
})


## ----glm_isogenic_vs_b6_ctrl_no_covariates_model------------------------------
## discovery backcross
options(digits = 9)
j30.out <- glm(bw30 ~ strain, data = b6.nr$css)
j30.comp <- glht(j30.out, linfct = mcp(strain = "Tukey"),
                 alternative = "two.sided")

bw.out <- apply(b6.nr$css[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$css)
})
names(bw.out) <- bw.phenos
bw.comp <- lapply(bw.out, glht, linfct = mcp(strain = "Tukey"),
                  alternative = "two.sided")

mm.out <- apply(b6.nr$css[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$css)}
    )
names(mm.out) <- mm.phenos
mm.comp<- lapply(mm.out, glht, linfct = mcp(strain = "Tukey"),
                 alternative = "two.sided")

bc.out <- apply(b6.nr$css[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$css)
})
names(bc.out) <- bc.phenos
bc.comp<- lapply(bc.out, glht, linfct = mcp(strain = "Tukey"),
                 alternative = "two.sided")

bb.smry <- lapply(c(list("bw30" = j30.comp), bw.comp, mm.comp, bc.comp), summary)

bb.smry.df <- lapply(phenos, function(i) {
    xx <- do.call(cbind, bb.smry[[i]]$test)[, -c(1:2,7)]
    xx <- as.data.frame(xx) %>% rownames_to_column(var = "comparison")
    xx$tissue <- i
    return(xx)
}) %>% do.call(rbind, .)


## replicate backcross
rj30.out <- glm(bw30 ~ strain, data = b6.nr$rbc)
rj30.comp <- glht(rj30.out, linfct = mcp(strain = "Tukey"),
                  alternative = "two.sided")

rbw.out <- apply(b6.nr$rbc[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$rbc)
})
names(rbw.out) <- bw.phenos
rbw.comp <- lapply(rbw.out, glht, linfct = mcp(strain = "Tukey"),
                   alternative = "two.sided")


rmm.out <- apply(b6.nr$rbc[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$rbc)
})
names(rmm.out) <- mm.phenos
rmm.comp<- lapply(rmm.out, glht, linfct = mcp(strain = "Tukey"),
                  alternative = "two.sided")

rbc.out <- apply(b6.nr$rbc[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain, data = b6.nr$rbc)}
    )
names(rbc.out) <- bc.phenos
rbc.comp<- lapply(rbc.out, glht, linfct = mcp(strain = "Tukey"),
                  alternative = "two.sided")


## 
rbb.smry <- lapply(c(list("bw30" = rj30.comp), rbw.comp, rmm.comp, rbc.comp), summary)
rbb.smry.df <- lapply(phenos, function(i) {
    xx <- do.call(cbind, rbb.smry[[i]]$test)[, -c(1:2,7)]
    xx <- as.data.frame(xx) %>% rownames_to_column(var = "comparison")
    xx$tissue <- i
    return(xx)
}) %>% do.call(rbind, .)

bb.smry.df$cohort <- "Discovery Backcross (N2)"
rbb.smry.df$cohort <- "Replicate Backcross (N2)"

rbind(bb.smry.df, b6.nr$rbc.smry.df) %>%
    xlsx::write.xlsx(file = file.path(tabDir, "b6_comparisons_effects_sd_pval_2023.xlsx"),
                     sheetName = "BodyComp SMRY", row.names = FALSE)


## ----testing_other_models-----------------------------------------------------
## having Litter size as a covariate
j30.ls.out <- glm(bw30 ~ strain + ls, data = b6.nr$css)
j30.ls.comp <- glht(j30.out, linfct = mcp(strain = "Tukey"),
                    alternative = "two.sided")

bw.ls.out <- apply(b6.nr$css[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$css)
})
names(bw.ls.out) <- bw.phenos
bw.ls.comp <- lapply(bw.ls.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

mm.ls.out <- apply(b6.nr$rbc[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$rbc)
})
names(mm.ls.out) <- mm.phenos
mm.ls.comp<- lapply(mm.ls.out, glht, linfct = mcp(strain = "Tukey"),
                    alternative = "two.sided")

bc.ls.out <- apply(b6.nr$css[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$css)
})
names(bc.ls.out) <- bc.phenos
bc.ls.comp <- lapply(bc.ls.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")


## having multiple dams (md) as a covariate
j30.md.out <- glm(bw30 ~ strain + md, data = b6.nr$css)
j30.md.comp <- glht(j30.out, linfct = mcp(strain = "Tukey"),
                    alternative = "two.sided")

bw.md.out <- apply(b6.nr$css[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + md, data = b6.nr$css)
})
names(bw.md.out) <- bw.phenos
bw.md.comp <- lapply(bw.md.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

mm.md.out <- apply(b6.nr$css[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + md, data = b6.nr$css)
})
names(mm.md.out) <- mm.phenos
mm.md.comp<- lapply(mm.md.out, glht, linfct = mcp(strain = "Tukey"),
                    alternative = "two.sided")

bc.md.out <- apply(b6.nr$css[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + md, data = b6.nr$css)
})
names(bc.md.out) <- bc.phenos
bc.md.comp <- lapply(bc.md.out, glht, linfct = mcp(strain = "Tukey"), 
                     alternative = "two.sided")


## having dam.age as a covariate
## estimate as median age of all dams in cage @ litter DOB
## check in old laptop or `litters` to see if where I have it
j30.da.out <- glm(bw30 ~ strain + dam.age, data = b6.nr$css)

j30.da.comp <- glht(j30.out, linfct = mcp(strain = "Tukey"),,
                    alternative
                    alternative = "two.sided")

bw.da.out <- apply(b6.nr$css[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + dam.age, data = b6.nr$css)
})
names(bw.da.out) <- bw.phenos
bw.da.comp <- lapply(bw.da.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

mm.da.out <- apply(b6.nr$css[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + dam.age, data = b6.nr$css)
})
names(mm.da.out) <- mm.phenos
mm.da.comp <- lapply(mm.da.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

bc.da.out <- apply(b6.nr$css[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + dam.age, data = b6.nr$css)
})
names(bc.da.out) <- bc.phenos
bc.da.comp <- lapply(bc.da.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")


## All covariates// Full model
j30.f.out <- glm(bw30 ~ strain + ls + sp + md, data = b6.nr$css)
j30.f.comp <- glht(j30.out, linfct = mcp(strain = "Tukey"),
                   alternative = "two.sided")

bw.f.out <- apply(b6.nr$css[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls + sp + md, data = b6.nr$css)
})
names(bw.f.out) <- bw.phenos
bw.f.comp <- lapply(bw.f.out, glht, linfct = mcp$1)
})(strain = "Tukey"),
alternative = "two.sided")

mm.f.out <- apply(b6.nr$css[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls + sp + md, data = b6.nr$css)
})
names(mm.f.out) <- mm.phenos
mm.f.comp<- lapply(mm.f.out, glht, linfct = mcp(strain = "Tukey"),
                   alternative = "two.sided")

bc.f.out <- apply(b6.nr$css[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls + sp + md, data = b6.nr$css)
})
names(bc.f.out) <- bc.phenos
bc.f.comp <- lapply(bc.f.out, glht, linfct = mcp(strain = "Tukey"),
                    alternative = "two.sided")




## ----Replicate_Backcross_Covariate_Testing-----------------------------------
## only on litter size
rj30.ls.out <- glm(bw30 ~ strain, data = b6.nr$rbc)
rj30.ls.comp <- glht(rj30.ls.out, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

rbw.ls.out <- apply(b6.nr$rbc[, bw.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$rbc)
})
names(rbw.ls.out) <- bw.phenos
rbw.ls.comp <- lapply(rbw.ls.out, glht, linfct = mcp(strain = "Tukey"),
                      alternative = "two.sided")

rmm.ls.out <- apply(b6.nr$rbc[, mm.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$rbc)
})
names(rmm.ls.out) <- mm.phenos
rmm.ls.comp<- lapply(rmm.ls.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")

rbc.ls.out <- apply(b6.nr$rbc[, bc.phenos], 2, function(PHEN) {
    glm(PHEN ~ strain + ls, data = b6.nr$rbc)
})
names(rbc.ls.out) <- bc.phenos
rbc.ls.comp<- lapply(rbc.ls.out, glht, linfct = mcp(strain = "Tukey"),
                     alternative = "two.sided")


## ----Summary_Plots:_MEANS----------------------------------------------------
colp <- brewer.pal(10, "Paired")[c(9,10,7,8,3,4,5,6,1,2)]
new.axis.labels <- c("B6.A-15 N2 B/B", "B6.A-17 N2 B/B", 
                     "B6.A-19 N2 B/B", "B6.A-X N2 B/B","B6 CTRL")

pdf(file.path(tmpDir, "CSS_B6NR_means.pdf"), width = 10, height = 6.5)
par(mar = c(8.3, 4, 2, 2)+.1, cex.lab = 1.4, cex.axis = 1.4, pch = 15)
##
layout(mat = matrix(c(1:8), byrow = TRUE, ncol = 4))
lapply(c("bw30", bw.phenos), function(PHEN) {
    plotmeans(b6.nr$css[, PHEN] ~ b6.nr$css$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(2, 10, by= 2)])
    axis(1, at = 1:5, labels = new.axis.labels, las = 2)
})
##
layout(mat = matrix(c(1:4, rep(0, 4)), byrow = TRUE, ncol = 4))
lapply(mm.phenos, function(PHEN) {
    plotmeans(b6.nr$css[, PHEN] ~ b6.nr$css$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(2, 10, by= 2)])
    axis(1, at = 1:5, labels = new.axis.labels, las = 2)
})
##
layout(mat = matrix(c(1:8), byrow = TRUE, ncol = 4))
lapply(bc.phenos[1:8], function(PHEN) {
    plotmeans(b6.nr$css[, PHEN] ~ b6.nr$css$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(2, 10, by= 2)])
    axis(1, at = 1:5, labels = new.axis.labels, las = 2)
})
##
layout(mat = matrix(c(1:5, rep(0, 3)), byrow = TRUE, ncol = 4))
lapply(bc.phenos[9:13], function(PHEN) {
    plotmeans(b6.nr$css[, PHEN] ~ b6.nr$css$strain,
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(2, 10, by= 2)])
    axis(1, at = 1:5, labels = new.axis.labels, las = 2)
})
dev.off()


pdf(file.path(tmpDir, "RBC_B6_NR_means.pdf"), width = 10, height = 6.5)
par(mar = c(8.3, 4, 2, 2)+.1, cex.lab = 1.4, cex.axis = 1.4, pch = 15)
##
layout(mat = matrix(c(1:8), byrow = TRUE, ncol = 4))
lapply(c("bw30", bw.phenos), function(PHEN) {
    plotmeans(b6.nr$rbc[, PHEN] ~ b6.nr$rbc$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(4, 10, by= 2)])
    axis(1, at = 1:4, labels = new.axis.labels[2:5], las = 2)
})
##
layout(mat = matrix(c(1:4, rep(0, 4)), byrow = TRUE, ncol = 4))
lapply(mm.phenos, function(PHEN) {
    plotmeans(b6.nr$rbc[, PHEN] ~ b6.nr$rbc$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(4, 10, by= 2)])
    axis(1, at = 1:4, labels = new.axis.labels[2:5], las = 2)
})
##
layout(mat = matrix(c(1:8), byrow = TRUE, ncol = 4))
lapply(bc.phenos[1:8], function(PHEN) {
    plotmeans(b6.nr$rbc[, PHEN] ~ b6.nr$rbc$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(4, 10, by= 2)])
    axis(1, at = 1:4, labels = new.axis.labels[2:5], las = 2)
})
##
layout(mat = matrix(c(1:5, rep(0, 3)), byrow = TRUE, ncol = 4))
lapply(bc.phenos[9:13], function(PHEN) {
    plotmeans(b6.nr$rbc[, PHEN] ~ b6.nr$rbc$strain, 
              p = 0.95, xlab = "", ylab = abv[PHEN, "full.name"],
              connect = FALSE, cex = 3.2, xaxt = "n",
              col = colp[seq(4, 10, by= 2)])
    axis(1, at = 1:4, labels = new.axis.labels[2:5], las = 2)
})
dev.off()




####################################################################
##  Body weights over time
####################################################################

## data
library(lattice)
wt <- read.csv("weights.csv", na.strings = c("", "-", NA, "."))
wt <- wt[!is.na(wt$weight),]


pdf("../figures/CSS_F0_F1_BC1_vs_time.pdf", width = 10, height = 4)
layout(mat = matrix(1:3, nrow = 1))
par(mar = c(5.5, 5, 1.5, 1)+.1, cex.axis = 1.4, cex.lab = 1.5)
with(wt[as.character(wt$strain) == "B6.A-15" & wt$sex == "M", ], 
     plotmeans(weight ~ age, ylim = c(10, 30), col = colp[2], lwd = 4, pch = 15, n.label = FALSE,
               ylab = "Weight (g)", xlab = "Age (d)", main = expression(F[0])))
with(wt[as.character(wt$strain) == "B6.A-17" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[4], 
               lwd = 4, pch = 16, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-19" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[6], 
               lwd = 4, pch = 17, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-X" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[8], 
               lwd = 4, pch = 18, n.label = FALSE, xaxt = "n"))
with(wt[wt$id <= 6415 & as.character(wt$strain) == "B6" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[10], 
               lwd = 4, pch = 19, n.label = FALSE, xaxt = "n"))
abline(h = seq(10, 30, 5), col = "gray50")
legend("topleft", 
       legend = c("B6.A - 15", "B6.A - 17", "B6.A - 19", "B6.A - X", "B6"), 
       col = colp[c(2, 4, 6, 8, 10)], lwd = 2, pch = 15:19, 
       bty = "n", bg = "white")

par(mar = c(5.5, 5, 1.5, 1)+.1, cex.axis = 1.4, cex.lab = 1.5)
with(wt[wt$id <= 6415 & 
	as.character(wt$strain) == "B6.A-15 F1" & wt$sex == "M", ],
     plotmeans(weight ~ age, ylim = c(10, 30), col = colp[2],
               lwd = 4, pch = 15, n.label = FALSE,
               ylab = "Weight (g)", xlab = "Age (d)", main = expression(F[1])))
with(wt[as.character(wt$strain) == "B6.A-17 F1" & wt$sex == "M", ],
     plotmeans(weight ~ age, add = TRUE, col = colp[4], 
               lwd = 4, pch = 16, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-19 F1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[6], 
               lwd = 4, pch = 17, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-X F1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[8], 
               lwd = 4, pch = 18, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[10], 
               lwd = 4, pch = 19, n.label = FALSE, xaxt = "n"))
abline(h = seq(10, 30, 5), col = "gray50")

par(mar = c(5.5, 5, 1.5, 1)+.1, cex.axis = 1.4, cex.lab = 1.5)
with(wt[as.character(wt$strain) == "B6.A-15 BC1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, ylim = c(10, 30), col = colp[2], 
               lwd = 4, pch = 15, n.label = FALSE,
               ylab = "Weight (g)", xlab = "Age (d)", 
               main = expression(BC[1])))
with(wt[as.character(wt$strain) == "B6.A-17 BC1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[4], 
               lwd = 4, pch = 16, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-19 BC1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[6], 
               lwd = 4, pch = 17, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.A-X F1" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[8], 
               lwd = 4, pch = 18, n.label = FALSE, xaxt = "n"))
with(wt[as.character(wt$strain) == "B6.C" & wt$sex == "M", ], 
     plotmeans(weight ~ age, add = TRUE, col = colp[10], 
               lwd = 4, pch = 19, n.label = FALSE, xaxt = "n"))
abline(h = seq(10, 30, 5), col = "gray50")
dev.off()

##lit.size <- lm(weight ~ ls*age, data = wt)
##weights.ls <- wt$weights + mean(wt$ls)


pdf(file.path(tmpDir, "CSS_Comparisons.pdf"), width = 10, height = 10)
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.comp$bw40, main = expression(bw40==strain))
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.ls.comp$bw40, main = expression(bw40==strain+litter~size), 
     col = "red")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.sp.comp$bw40, main = expression(bw40==strain+sire~presence), 
     col = "blue")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.md.comp$bw40, main = expression(bw40==strain+multiple~dams), 
     col = "green3")
par(mar = c(6, 12, 4, 1)+.1)
plot.new() ## plot(bw.da.comp$bw40, main = expression(bw40==strain+dam~age), col = "purple2")
##
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.comp$bw50, main = expression(bw50==strain))
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.ls.comp$bw50, main = expression(bw50==strain+ls), col = "red")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.sp.comp$bw50, main = expression(bw50==strain+sire~presence), col = "blue")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.md.comp$bw50, main = expression(bw50==strain+multiple~dams), col = "green3")
par(mar = c(6, 12, 4, 1)+.1)
plot.new() ## plot(bw.da.comp$bw50, main = expression(bw50==strain+dam~age), col = "purple2")
##
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.comp$bw60, main = expression(bw60==strain))
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.ls.comp$bw60, main = expression(bw60==strain+ls), col = "red")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.sp.comp$bw60, main = expression(bw60==strain+sire~presence), col = "blue")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.md.comp$bw60, main = expression(bw60==strain+multiple~dams), col = "green3")
par(mar = c(6, 12, 4, 1)+.1)
plot.new() ## plot(bw.da.comp$bw60, main = expression(bw60==strain+dam~age), col = "purple2")
##
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.comp$pwg, main = expression(pwg==strain))
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.ls.comp$pwg, main = expression(pwg==strain+ls), col = "red")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.sp.comp$pwg, main = expression(pwg==strain+sire~presence), col = "blue")
par(mar = c(6, 12, 4, 1)+.1)
plot(bw.md.comp$pwg, main = expression(pwg==strain+multiple~dams), col = "green3")
par(mar = c(6, 12, 4, 1)+.1)
plot.new() # plot(bw.da.comp$pwg, main = expression(pwg==strain+dam~age), col = "purple2")
dev.off()

#################################################################################
##  scatter plot of t-value
#################################################################################
j30.at <- merge(summary(j30.comp)$test$tstat, 
                summary(rj30.comp)$test$tstat, by = 0, all.x = TRUE)

bw.t <- sapply(lapply(bw.comp, summary), function(x) x$test$tstat)
bw.rt <- sapply(lapply(rbw.comp, summary), function(x) x$test$tstat)
bw.at <- merge(bw.t, bw.rt, by = 0, all.x = TRUE)

mm.t <- sapply(lapply(mm.comp, summary), function(x) x$test$tstat)
mm.rt <- sapply(lapply(rmm.comp, summary), function(x) x$test$tstat)
mm.at <- merge(mm.t, mm.rt, by = 0, all.x = TRUE)

bc.t <- sapply(lapply(bc.comp, summary), function(x) x$test$tstat)
bc.rt <- sapply(lapply(rbc.comp, summary), function(x) x$test$tstat)
bc.at <- merge(bc.t, bc.rt, by = 0, all.x = TRUE)

gl <- c(-4, -2, 0, 2, 4)
lw <- c(1, 1, 1.3, 1, 1)
cll <- "gray"
lt <- c(2, 2, 1, 2, 2)
xylim <- c(-5, 5)

pdf(file.path(tmpDir, "CSS_vs_RBC_z-value_scatter.pdf"),
    width = 297/25.4, height = 2*(297/4)/25.4, paper = "a4r",
    useDingbats = FALSE)
##
layout(mat = matrix(c(1:8), nrow=2, byrow=TRUE))
par(mar = c(5, 5, 1, 1)+.1)
plot(j30.at$x, j30.at$y,
     xlab = paste(abv["bw30", "full.name"], "\nOriginal Cross"),
     ylab = paste(abv["bw30", "full.name"], "\nReplicate Cross"),
     xlim = xylim, ylim = xylim,
     pch = 19)
abline(v = gl, h = gl, col = cll, lwd = lw, lty = lt)
text(j30.at$x, j30.at$y, 
     labels = j30.at$Row.names, pos = 2, cex = 0.6)
sapply(bw.phenos, function(i) {
    plot(x = bw.at[, paste(i, ".x", sep = "")], 
         y = bw.at[,paste(i, ".y", sep = "")],
         xlab = paste(abv[i, "full.name"], "\nOriginal Cross"),
         ylab = paste(abv[i, "full.name"], "\nReplicate Cross"),
         xlim = xylim, ylim = xylim,
         pch = 19)
    abline(v = gl, h = gl, col = cll, lwd = lw, lty = lt)
    text(x = bw.at[, paste(i, ".x", sep = "")], 
         y = bw.at[,paste(i, ".y", sep = "")],
         labels = bw.at$Row.names, pos = 2, cex = 0.6)
})
##
layout(mat = matrix(c(1:4, rep(0, 4)), nrow = 2, byrow = TRUE))
sapply(mm.phenos, function(i) {
    plot(x = mm.at[, paste(i, ".x", sep = "")], 
         y = mm.at[,paste(i, ".y", sep = "")],
         xlab = paste(abv[i, "full.name"], "\nOriginal Cross"),
         ylab = paste(abv[i, "full.name"], "\nReplicate Cross"),
         xlim = xylim, ylim = xylim,
         pch = 19)
    abline(v = gl, h = gl, col = cll, lwd = lw, lty = lt)
    text(x = mm.at[, paste(i, ".x", sep = "")], 
         y = mm.at[,paste(i, ".y", sep = "")],
         labels = mm.at$Row.names, pos = 2, cex = 0.6)
})
##
layout(mat = matrix(c(1:8), nrow=2, byrow=TRUE))
sapply(bc.phenos, function(i) {
    plot(x = bc.at[, paste(i, ".x", sep = "")], 
         y = bc.at[,paste(i, ".y", sep = "")],
         xlab = paste(abv[i, "full.name"], "\nOriginal Cross"),
         ylab = paste(abv[i, "full.name"], "\nReplicate Cross"),
         xlim = xylim, ylim = xylim,
         pch = 19)
    abline(v = gl, h = gl, col = cll, lwd = lw, lty = lt)
    text(x = bc.at[, paste(i, ".x", sep = "")], 
         y = bc.at[,paste(i, ".y", sep = "")],
         labels = bc.at$Row.names, pos = 2, cex = 0.6)
})
dev.off()


pdf(file.path(tmpDir, "CSS_N2_BB_Strain_GLM_plots.pdf"), 
    width = 210/25.4, height = 210/25.4)
layout(mat = matrix(1:4, nrow = 2, byrow= TRUE))
sapply(bw.phenos, function(PHE) {
    plot(bw.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"), 
          side = 3, line = 3, at = 0, font = 2)
})
sapply(mm.phenos, function(PHE) {
    plot(mm.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"),
          side = 3, line = 3, at = 0, font = 2)
})
sapply(bc.phenos, function(PHE) {
    plot(bc.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"), 
          side = 3, line = 3, at = 0, font = 2)
})
dev.off()

pdf(file.path(tmpDir, "RBC_N2_BB_Strain_GLM_plots.pdf"), 
    width = 210/25.4, height = 210/25.4)
layout(mat = matrix(1:4, nrow = 2, byrow= TRUE))
sapply(bw.phenos, function(PHE) {
    plot(rbw.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"), 
          side = 3, line = 3, at = 0, font = 2)
})
sapply(mm.phenos, function(PHE) {
    plot(rmm.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"), 
          side = 3, line = 3, at = 0, font = 2)
})
sapply(bc.phenos, function(PHE) {
    plot(rbc.out[[PHE]], cex = 0.8)
    mtext(text = paste(abv[PHE, "full.name"], "~ strain"), 
          side = 3, line = 3, at = 0, font = 2)
})
dev.off()


## ----scatter_plot_of_t.value--------------------------------------------------
## merge Discovery and Replicate backcross t.stat
#% j30.at <- merge(summary(j30.comp)$test$tstat,
#%                 summary(rj30.comp)$test$tstat, by = 0, all.x = TRUE)
#% bw.t <- sapply(lapply(bw.comp, summary), function(x) x$test$tstat)
#% bw.rt <- sapply(lapply(rbw.comp, summary), function(x) x$test$tstat)
#% bw.at <- merge(bw.t, bw.rt, by = 0, all.x = TRUE)
#% mm.t <- sapply(lapply(mm.comp, summary), function(x) x$test$tstat)
#% mm.rt <- sapply(lapply(rmm.comp, summary), function(x) x$test$tstat)
#% mm.at <- merge(mm.t, mm.rt, by = 0, all.x = TRUE)
#% bc.t <- sapply(lapply(bc.comp, summary), function(x) x$test$tstat)
#% bc.rt <- sapply(lapply(rbc.comp, summary), function(x) x$test$tstat)
#% bc.at <- merge(bc.t, bc.rt, by = 0, all.x = TRUE)
#% all.t <- merge(j30.at, bw.at, by = "Row.names")
#% all.t <- merge(all.t, mm.at, by = "Row.names")
#% all.t <- merge(all.t, bc.at, by = "Row.names")
#% 
#% j30.t.ht <- merge(summary(j30.comp)$test$tstat, 
#%                   summary(rj30.comp)$test$tstat, by = 0, all.x = TRUE)

## write to file, then annotate in excel to match e.g.
## note: N2 = Discovery Backcross (N2), N2.2 = Replicate Backcross (N2)
#>  head(tval) %>% knitr::kable(.)
#> |comparison                       |trait                   |geno           |    N2| N2.2|
#> |:--------------------------------|:-----------------------|:--------------|-----:|----:|
#> |B6.A-17 N2 B/B vs B6.A-15 N2 B/B |Body Weight at 30 d (g) |B6/B6 vs B6/B6 | -1.38|   NA|
#> |B6.A-19 N2 B/B vs B6.A-15 N2 B/B |Body Weight at 30 d (g) |B6/B6 vs B6/B6 | -0.43|   NA|
#> |B6.A-19 N2 B/B vs B6.A-17 N2 B/B |Body Weight at 30 d (g) |B6/B6 vs B6/B6 |  0.91| 1.88|
#> |B6.A-X F1 B/B vs B6.A-15 N2 B/B  |Body Weight at 30 d (g) |B6/B6 vs B6/B6 | -1.82|   NA|
#> |B6.A-X F1 B/B vs B6.A-17 N2 B/B  |Body Weight at 30 d (g) |B6/B6 vs B6/B6 | -0.48| 3.09|
#> |B6.A-X F1 B/B vs B6.A-19 N2 B/B  |Body Weight at 30 d (g) |B6/B6 vs B6/B6 | -1.35| 0.26|
#% write.table(melt(all.t), file = "B6NR_t-values.txt",
#%             quote = FALSE, row.names = FALSE, sep = "\t")
#% 
#% tval <- read.delim("tables/01.glm_B6_NR/B6NR_t-values.txt")
#% tval <- read.delim("tables/01.glm_B6_NR/B6NR_t-values.txt")
data(tval)

comp.bb <- unique(grep("B6.C", tval$comparison, value = TRUE))
comp.ab <- c("B6.A-15 N2 A/B vs B6.A 15 N2 B/B", 
             "B6.A-17 N2 A/B vs B6.A 17 N2 B/B", 
             "B6.A-19 N2 A/B vs B6.A 19 N2 B/B")
trait.bb <- c("Body Weight at 30 d (g)", 
              "Body Weight at 40 d (g)", "Body Weight at 50 d (g)", 
              "Body Weight at 60 d (g)", "Carcass Weight (g)", 
              "Kidney (g)", "Heart (g)", "Brain Weight (g)")

pdf("Selected_traits_t-value_scatter.pdf", width = 210/25.4, height = 120/25.4)
ggplot(tval[tval$comparison %in% c(comp.bb, comp.ab) & 
            tval$trait %in% trait.bb ,], 
       aes(N2, N2.2, fill = comparison, colour = comparison, shape = geno)) +
    geom_point() +
    facet_wrap( ~ trait, ncol = 4) +
    xlim(-5, 5) + ylim(-5, 5) +
    xlab("Discovery Backcross (N2)") + 
    ylab("Replicate Backcross (N2)") + 
    theme_bw() +
    theme(legend.position="bottom", text = element_text(size = 9))
dev.off()


pdf("All_traits_t-value_scatter.pdf", width = 210/25.4, height = 210/25.4)
ggplot(tval[tval$comparison %in% c(comp.bb, comp.ab), ], 
       aes(N2, N2.2, fill = comparison, colour = comparison, shape = geno)) + 
    geom_point() +
    facet_wrap( ~ trait, ncol = 5) + 
    xlim(-5, 5) + ylim(-5, 5) +
    xlab("Discovery Backcross (N2)") + 
    ylab("Replicate Backcross (N2)") + 
    theme_bw() +
    theme(legend.position = "bottom", text = element_text(size = 9))
dev.off()

##ggplot(tval[tval$comparison %in% comp.bb,], 
##      aes(N2, N2.2, fill = trait, colour = trait, shape = comparison)) + 
##  geom_point() + xlim(-5, 5) + ylim(-5, 5)


## __EOF__
