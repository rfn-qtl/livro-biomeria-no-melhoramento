##################### Chapter 29 - Script examples of using R #######################

# In this chapter we will no follow necessarily the same order or structure of threvious ones. On the other hand, we will provide a sequence of analysis usually applied to make inferences about a breeding dataset. 


###########################################################################################
# Example 1 - Mixed models equations, experimental designs, spatial analysis, components of variance, heritability, and selection index 
###########################################################################################

#grain yield in the well-watered experiment at Murony (MAIZE DATASET)
m13w <- read.csv("mur13w.csv")
head(m13w)
tail(m13w)
summary(m13w)

# Set variables required to model spatial variation
m13w$Col <- as.factor(m13w$Column)
m13w$Row <- as.factor(m13w$Row)
m13w$Rep <- as.factor(m13w$Replicate)
m13w$block <- as.factor(m13w$block)

# Spatial analysis adjustments
m13w$Colnum <- as.numeric(m13w$Column)
m13w$Rownum <- as.numeric(m13w$Row)
nrow <- max(m13w$Rownum)
ncol <- max(m13w$Colnum)
nseg.row <- nrow
nseg.col <- ncol

# plotting maps of agricultural field experiments that are laid out in grids
library(desplot)
desplot(m13w, grain.yield ~ Col * Row,
        out1 = Replicate,
        out2 = block)

# running different experimental designs via mixed model equations
library("SpATS")
# Perform mixed models with Genotypes as fixed and random to estimate BLUEs and BLUPs
# As it is possible see, deffferent models can be accomodated as random or fixed effects
# Even spatial analysis
# Fist, let's run genotypes as fixed
yld_hybrid_fixed <- SpATS(response = "grain.yield", 
                          fixed = ~ Rep, 
                          random = ~ Col + Row, 
                          spatial = ~ PSANOVA(Colnum, Rownum, nseg = c(nseg.col, nseg.row)), 
                          genotype = "Variety_ID", 
                          genotype.as.random = FALSE, 
                          data = m13w)
# obtaining the spatial trends - from raw data to BLUES and BLUPS
plot(yld_hybrid_fixed) 

# Now, run as random
yld_hybrid_random <- SpATS(response = "grain.yield", 
                           fixed = ~ Rep, 
                           random = ~ Col + Row, 
                           spatial = ~ PSANOVA(Colnum, Rownum,nseg = c(nseg.col, nseg.row)),
                           genotype = "Variety_ID", 
                           genotype.as.random = TRUE, 
                           weights = NULL,
                           family = gaussian(),
                           data = m13w)
# spatial plots and trends
plot.SpATS(yld_hybrid_random)

# Deviance and LRT - to test a significance of a factor, first we need to reduce the model (remove the factor), and run again
# First, let's save the deviance of the full model ((i.e., - 2 times the restricted log-likelihood))
#  Deviance via R-package JOPS
library(JOPS)
yld_hybrid_random.full  <- SpATS.nogeno(response = "grain.yield", 
                          fixed = ~ Rep, 
                          random = ~ Col + Row + Variety_ID, 
                          spatial = ~ PSANOVA(Colnum, Rownum, nseg = c(nseg.col, nseg.row)),
                          weights = NULL,
                          family = gaussian(),
                          control = list(maxit = 25, tolerance = 1e-05,
                          monitoring = 0, update.psi = FALSE),
                          data = m13w)

yld_hybrid_random.full$deviance[1]

yld_hybrid_random.red  <- SpATS.nogeno(response = "grain.yield", 
                                        fixed = ~ Rep, 
                                        random = ~ Col + Row, 
                                        spatial = ~ PSANOVA(Colnum, Rownum, nseg = c(nseg.col, nseg.row)),
                                        weights = NULL,
                                        family = gaussian(),
                                        control = list(maxit = 25, tolerance = 1e-05,
                                                       monitoring = 0, update.psi = FALSE),
                                        data = m13w)
yld_hybrid_random.red$deviance[1]

# The LRT test is done by the absolute difference between the models  
yld_hybrid_random.full$deviance[1] - yld_hybrid_random.red$deviance[1]
# Then, see the probability using chi-square distribution with 1 DF
# in this case, the row effect was not significant
pchisq(abs(yld_hybrid_random.full$deviance[1] - yld_hybrid_random.red$deviance[1]), 1, lower.tail = F) 
# see some standard values
pchisq(3.65, 1, lower.tail = F) # significant 5%
pchisq(6.63, 1, lower.tail = F) # significant 1%

# Components of variance can be obtained via
yld_hybrid_random$var.comp

# to estimate the residuals, for this specific model we can use the following equations
sigma.res <- yld_hybrid_random$residuals
sigma.res[is.na(sigma.res)] <- 0 # eliminate NA
sigma.res <- var(sigma.res)

# to obtain the heritability via the package function we can use 
getHeritability(yld_hybrid_random)
# or manually via
yld_hybrid_random$var.comp[1] / sum(yld_hybrid_random$var.comp[c(1,8)], sigma.res)
# Broad-sense heritability based on Cullis method
Vg <- yld_hybrid_random$var.comp["Variety_ID"]
ng <- length(unique(m13w$Variety_ID))
C11_g <- yld_hybrid_random$vcov$C11_inv
dim(C11_g)
trC11_g <-sum(diag(C11_g))
av2 <- 2/ng * (trC11_g - (sum(C11_g) - trC11_g) / ng-1) # mean var of a difference between genotypic BLUPS
(H2.Cullis <- 1 - av2 / (2 * Vg))

# Estimate BLUPs for Grain Yield
blups_m13w <- predict.SpATS(yld_hybrid_random, which = "Variety_ID")
# Reliability
(rel <- mean(1 - blups_m13w$standard.errors^2 / yld_hybrid_random$var.comp[1]))
# weights for ID's - adjust residual for further analysis
vcov.mme <- yld_hybrid_random$vcov$C11_inv
w <- diag(vcov.mme)

# Estimate BLUEs
blues_m13w <- predict.SpATS(yld_hybrid_fixed, which = "Variety_ID")

# Estimate BLUPs for anthesis
ant_hybrid_random <- SpATS(response = "anthesis", 
                           fixed = ~ Rep, 
                           random = ~ Col + Row, 
                           spatial = ~ PSANOVA(Colnum, Rownum,nseg = c(nseg.col, nseg.row)),
                           genotype = "Variety_ID", 
                           genotype.as.random = TRUE, 
                           weights = NULL,
                           family = gaussian(),
                           data = m13w)

# Now, estimate BLUPs for  anthesis
blups_m13w.ant <- predict.SpATS(ant_hybrid_random, which = "Variety_ID")

########### Selection index
# phenotypic correlation between trait
cor(m13w[,c("anthesis", "grain.yield")], use = "pairwise", method = "pearson")
# phenotypic covariance between trait
P <- cov(m13w[,c("anthesis", "grain.yield")], use = "pairwise", method = "pearson")

# genotypic correlation between traits
cor(blups_m13w[,"predicted.values"], blups_m13w.ant["predicted.values"], method = "pearson")
# genotypic covariance between trait
G <- cov(cbind(blups_m13w[,"predicted.values"], blups_m13w.ant["predicted.values"]))

# Smith-Hazel
# define the economic weights per trait
ecoW <- c(1, 2) # in this case two times more for grain yield
# then, the selection weights per trait 
b <- solve(as.matrix(P)) %*% as.matrix(G) %*% as.matrix(ecoW)
# Finally, the vector of SI per genotype
ISsh <- as.matrix(m13w[,c("anthesis", "grain.yield")]) %*% b  

# Pasek-Baker
# define the desired genetic gains per traits in genetic standard deviations
desired <- c(1, 2)
G.scaled <- cov(scale(cbind(blups_m13w[,"predicted.values"], blups_m13w.ant["predicted.values"])))
b <- solve(as.matrix(G.scaled)) %*% as.matrix(desired)
ISpb <- as.matrix(m13w[,c("anthesis", "grain.yield")]) %*% b

# correlation between indices
cor(ISsh, ISpb, use = "pairwise", method = "pearson")


###########################################################################################
# Example 2 - MET, AMMI, GGE-Biplot, adaptability, and stability
###########################################################################################

library(sommer)
data(DT_example)
DT <- DT_example

str(DT_example)
head(DT_example)

# Fitting genotype by environment models
fitMET <- mmer(Yield ~ Env,
               random= ~ Block + Year + Name + Name:Env,
               rcov= ~ units,
               data=DT, 
               verbose = FALSE)
summary(fitMET)$varcomp
blocks <- length(unique(DT_example$Block))
loc <- length(unique(DT_example$Env))

# Broad-sense heritability
vpredict(fitMET, h2 ~ V3 / ( V3 + V4/loc + V5/(loc*blocks) ) ) # trials level
vpredict(fitMET, h2 ~ V3 / ( V3 + V4 + V5 ) ) # plot level

# H2 Cullis
Vg <- fitMET$sigma$Name
ng <- length(unique(DT_example$Name))
C22_g <- fitMET$PevU$Name$Yield
trC22_g <-sum(diag(C22_g))
av2 <- 2/ng * (trC22_g - (sum(C22_g) - trC22_g) / ng-1) # mean var of a difference between genotypic BLUPS
(H2.Cullis <- 1 - av2 / (2 * Vg))

# weights for ID's - adjust residual for further analysis
w <- diag(C22_g)
 
# predicting the BLUP - main effect
BLUPs <- predict.mmer(object = fitMET, classify = "Name")
BLUPs$pvals
# reliability
mean(1 - BLUPs$pvals$standard.error^2 / Vg)

# predicting the BLUP per enviroment
BLUPS.env <- predict.mmer(object = fitMET, classify = c("Name","Env"))
BLUPS.env$pvals


# running another option, considering an diagonal model  (DG)
fitMET.DG <- mmer(Yield ~ Env + Block,
               random= ~ Year + vsr(usr(Env), Name),
               rcov= ~ units,
               data=DT, 
               verbose = FALSE)
summary(fitMET.DG)$varcomp


# running another option, considering an diagonal model  (DG)
fitMET.DG <- mmer(Yield ~ Env + Block,
                  random= ~ Year + vsr(dsr(Env), Name),
                  rcov= ~ units,
                  data=DT, 
                  verbose = FALSE)
summary(fitMET.DG)$varcomp


# MET: unstructured model (US)
fitMET.US <- mmer(Yield ~ Env + Block,
                  random= ~ Year + vsr(usr(Env), Name),
                  rcov= ~ units,
                  data=DT, 
                  verbose = FALSE)
summary(fitMET.US)$varcomp

# Comparing the models - LRT, AIC, and BIC
lrt12 <- anova(fitMET, fitMET.DG)
lrt23 <- anova(fitMET.DG, fitMET.US)


# Stability, daptability, AMMI/GGE – Biplot via statgenGxE package 
# Bart-Jan van Rossum
# remotes::install_github("Biometris/statgenGxE", ref = "develop", dependencies = TRUE)
# https://cran.r-project.org/web/packages/statgenGxE/vignettes/statgenGxE.html
library(statgenGxE)
library(sommer)
data(DT_example)
DT <- DT_example
head(DT)

## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = DT, genotype = "Name", trial = "Env")

## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "Yield")
summary(dropsVarComp)

## Extract variance components.
vc(dropsVarComp)

## Compute heritability.
herit(dropsVarComp)

## Plot the results of the fitted model.
plot(dropsVarComp)

## Predictions of the genotype main effect.
predGeno <- predict(dropsVarComp)
head(predGeno)

## predictions at the level of genotype x scenarioFull.
predGenoTrial <- predict(dropsVarComp, predictLevel = "trial")
head(predGenoTrial)


# Finlay-Wilkinson Analysis
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "Yield")
summary(dropsFW)
## Create scatter plot for Finlay Wilkinson analysis.
## Color genotypes by geneticGroup.
plot(dropsFW, plotType = "scatter", colorGenoBy = "Year")
## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line", colorGenoBy = "Year")


# AMMI Analysis
data(dropsPheno)
## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID", trial = "Experiment")
## Run gxeAmmi for grain.yield.
dropsAm <- gxeAmmi(TD = dropsTD, trait = "grain.yield")
summary(dropsAm)
## Run gxeAmmi. Let algorithm determine number of principal components.
dropsAm2 <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = NULL)
summary(dropsAm2)
## Run gxeAmmi with three principal components.
## Exclude genotypes 11430 and A3.
dropsAm3 <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = 3, 
                    excludeGeno = c("11430", "A3"))
## Run gxeAmmi per year in the data.
dropsAmYear <- gxeAmmi(TD = dropsTD, trait = "grain.yield", byYear = TRUE)
## Create an AMMI1 plot.
plot(dropsAm, plotType = "AMMI1")
## Create an AMMI2 biplot with symmetric scaling.
plot(dropsAm, scale = 0.5, plotType = "AMMI2")
## Create an AMMI2 biplot.
## Color genotypes based on variable geneticGroup. Use custom colors.
## Color environments based on variable scenarioFull.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", 
     colorGenoBy = "geneticGroup", colGeno = c("red", "blue", "green", "yellow"),
     colorEnvBy = "scenarioFull")
## Create an AMMI2 biplot with convex hull around the genotypes.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", plotConvHull = TRUE, colorEnvBy = "scenarioFull")
## Create an AMMI2 biplot.
## Align environment Mur13W with the positive x-axis.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", colorEnvBy = "scenarioFull",
     rotatePC = "Mur13W")


# GGE Analysis
## Run gxeGGE with default settings.
dropsGGE <- gxeGGE(TD = dropsTD, trait = "grain.yield")
summary(dropsGGE) 
## Create a GGE2 biplot.
plot(dropsGGE, plotType = "GGE2")
## Compute mega environments.
dropsMegaEnv <- gxeMegaEnv(TD = dropsTD, trait = "grain.yield")
## Summarize results.
summary(dropsMegaEnv)

# Stability measures
## Compute stability measures for dropsTD.
dropsStab <- gxeStability(TD = dropsTD, trait = "grain.yield")
## In the summary print the top two percent of the genotypes.
summary(dropsStab, pctGeno = 2)
## Create plots for the different stability measures.
## Color genotypes by geneticGroup.
plot(dropsStab, colorGenoBy = "geneticGroup")


###########################################################################################
# Example 3 - Path analysis, PCA, Cluster analysis, FA, 
###########################################################################################

# Path Analysis
library(lavaan)
library(semPlot)
library(OpenMx)
library(GGally)
library(corrplot)

#loading data - EX1
pheno <- readRDS("pheno")
head(pheno)
cor1 <- cor(pheno[, c("MPS", "SRA",  "NAE")], use = "pairwise")
corrplot(cor1, method = 'square')
#Fit Confirmatory Factor Analysis Models
fit1 = cfa("MPS ~ SRA + NAE", data = pheno)
summary(fit1, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
# graph
semPaths(fit1, 'std', layout = 'circle')

# EX2 - two regressions
data2 <- mtcars
head(data2)
model2 = 'mpg ~ hp + gear + cyl + disp + carb + am + wt 
hp ~ cyl + disp + carb'
fit2 <- cfa(model2, data = data2)
summary(fit2, fit.measures = TRUE, standardized = TRUE, rsquare = TRUE)
semPaths(fit2, 'std', 'est', curveAdjacent = TRUE, style = "OpenMx", layout = "tree2")


# PCA - Principal components analysis
# if (!requireNamespace('BiocManager', quietly = TRUE))
#   install.packages('BiocManager')
# 
# BiocManager::install('PCAtools')
library(PCAtools)
data2 <- mtcars
head(data2)

# the analysis, removing 0% of the variables based on variance
p <- pca(t(data2), center = T, removeVar = 0.0)

# identify the best number of PCA to retain
elbow <- findElbowPoint(p$variance)
elbow
# Thus
library(ggplot2)
screeplot(p,
          components = getComponents(p, 1:length(p$loadings)),
          vline = elbow) +
  geom_label(aes(x = elbow + 1, y = 50,
                 label = 'elbow', vjust = -1, size = 8))

which(cumsum(p$variance) > 98)[1]

# outputs
head(p$loadings)
head(p$rotated)
head(p$variance)
head(p$sdev)
# Other plots
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)
pairsplot(p)
# Determine the variables that drive variation among each PC
plotloadings(p, labSize = 2)

# now Remove this % of variables based on low variance. 
p <- pca(t(data2), center = T, rank = elbow, removeVar = round((ncol(mtcars) - elbow) / ncol(mtcars), 1))
# outputs
head(p$loadings)
head(p$rotated)
head(p$variance)
head(p$sdev)

# A scree and other plots
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)


# Clustering
pca <- prcomp(data2, center = F, scale. = F) 
var.pca <- summary(pca)$importance[, 1:11] 
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization

# or we can use this function to determine the optimal number of clusters
set.seed(123)
dim(data2)
fviz_nbclust(scale(data2), kmeans, method = "silhouette")
fviz_nbclust(scale(data2), kmeans, method = "gap_stat")
fviz_nbclust(scale(data2), kmeans, method = "wss")

# K-Means Clustering
k <- kmeans(scale(data2), centers = 2, nstart = 25)
fviz_cluster(k, data = scale(data2))
k$size
sum(k$size)
k$cluster


# Canonical Correlation Analysis
data2 <- mtcars
x.var <- c("mpg", "hp", "gear", "am", "wt", "drat") 
y.var <- c("cyl", "disp", "carb", "qsec", "vs")

library(yacca)
cca.fit <- cca(data2[,x.var], data2[ ,y.var])
cca.regfit <- cca(data2[,x.var], data2[ ,y.var], reg.param=1) # Some minimal regularization
#View the results
cca.fit
summary(cca.fit)
plot(cca.fit)
cca.regfit
# F Test for Canonical Correlations Using Rao’s Approximation
F.test.cca(cca.fit)
helio.plot(cca.fit)


# Factor Analysis
library(psych)
# Kaiser provided the following values for interpreting the results:
# * 0.00 to 0.49 unacceptable
# * 0.50 to 0.59 miserable
# * 0.60 to 0.69 mediocre
# * 0.70 to 0.79 middling
# * 0.80 to 0.89 meritorious
# * 0.90 to 1.00 marvelous
KMO(data2)
# The results (MSA=0.49) are unacceptable. I am going to deal with that by eliminating all low-contribution variables.
data2 <- data2[, KMO(data2)$MSAi>0.50] # Get rid of all variables with MSA < 0.50
# Bartlett’s test compares the correlation matrix to an identity matrix (a matrix filled with zeroes).
cortest.bartlett(data2)
# Determine Number of Factors to Extract
ev <- eigen(cor(data2)) # get eigenvalues
ev$values
scree(data2, pc=FALSE)
fa.parallel(data2, fa="fa")

# Extract (and rotate) factors
Nfacs <- 2  # This is for four factors. You can change this as needed.
fit <- factanal(data2, Nfacs, rotation="promax")
print(fit, digits=2, cutoff=0.3, sort=TRUE)

# plot Factor 1 by Factor 2
load <- fit$loadings[,1:2]
plot(load, type="n") # set up plot
text(load, labels = names(data2), cex=.7)


# You can also visualize the factor model to ease your interpretation. Though, as you can see, this does not work very well when you have a lot of variables.
loads <- fit$loadings
fa.diagram(loads)
dim(fit$loadings)
head(fit$loadings)

# Consistency: Cronbach’s Alpha
alpha(data2, check.keys = TRUE)


###########################################################################################
# Example 4 - Mating designs, Wald test, and comparing models via AIC, BIC, and LRT
###########################################################################################
# Mating designs (Diallel)
library(sommer)

#Loading and adjusting the file for Full diallel (F1s + parents)
pheno <- readRDS("pheno")
pheno$female <- as.character(pheno$female)
pheno$female[is.na(pheno$female)] <- pheno$gid[is.na(pheno$female)]
pheno$female <- as.factor(pheno$female)
pheno$male <- as.character(pheno$male)
pheno$male[is.na(pheno$male)] <- pheno$gid[is.na(pheno$male)]
pheno$male <- as.factor(pheno$male)
head(pheno)

# Running the model
mix.FD <- mmer(MPS ~ rep + N, 
            random = ~female + male + male:female, 
            data = pheno)

# wald test
wald.test(b = mix.FD$Beta$Estimate, Sigma = mix.FD$VarBeta, Terms = 2) # rep
wald.test(b = mix.FD$Beta$Estimate, Sigma = mix.FD$VarBeta, Terms = 3) # N

# components of variance
summary(mix.FD)$varcomp

# Narrow-sense heritability 
vpredict(mix.FD, h2 ~ V1 / ( V1 + V4 ) ) # for females
vpredict(mix.FD, h2 ~ V2 / ( V2 + V4 ) ) # for males
vpredict(mix.FD, h2 ~ V3 / ( V3 + V4 ) ) # for SCA

# obtaining the breeding values (GCA and SCA)
females.FD <- predict.mmer(object = mix.FD, classify = "female")
females.FD$pvals
males.FD <- predict.mmer(object = mix.FD, classify = "male")
males.FD$pvals
SCA.FD <- predict.mmer(object = mix.FD, classify = c("male","female"))
SCA.FD$pvals


#Loading and adjusting the file for Partial diallel (only F1s)
pheno <- readRDS("pheno")
pheno <- pheno[!is.na(pheno$female),]
head(pheno)

# Running the model
mix.PD <- mmer(MPS ~ rep + N, 
               random = ~female + male + male:female, 
               data = pheno)

# wald test
wald.test(b = mix.PD$Beta$Estimate, Sigma = mix.PD$VarBeta, Terms = 2) # rep
wald.test(b = mix.PD$Beta$Estimate, Sigma = mix.PD$VarBeta, Terms = 3) # N

# components of variance
summary(mix.PD)$varcomp

# Narrow-sense heritability 
vpredict(mix.PD, h2 ~ V1 / ( V1 + V4 ) ) # for females
vpredict(mix.PD, h2 ~ V2 / ( V2 + V4 ) ) # for males
vpredict(mix.PD, h2 ~ V3 / ( V3 + V4 ) ) # for SCA

# obtaining the breeding values (GCA and SCA)
females.PD <- predict.mmer(object = mix.PD, classify = "female")
females.PD$pvals
males.PD <- predict.mmer(object = mix.PD, classify = "male")
males.PD$pvals
SCA.PD <- predict.mmer(object = mix.PD, classify = c("male","female"))
SCA.PD$pvals

# Comparing the models - LRT, AIC, and BIC
lrt <- anova(mix.FD, mix.PD)

############################ the end ##########################