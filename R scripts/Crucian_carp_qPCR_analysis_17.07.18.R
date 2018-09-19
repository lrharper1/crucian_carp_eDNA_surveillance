#' ---
#' Title: "Crucian carp (*Carassius carassius*) environmental DNA (eDNA) analysis"
#' Author: "Lynsey Rebecca Harper"
#' Date: "17th July 2018"
#' ---
#' 
#' 
#' A hydrolysis probe qPCR assay was developed for targeted detection of eDNA
#' from the non-native crucian carp, which also receives conservation management 
#' as a Biodiversity Action Plan species in the UK. Water samples from 20 ponds in 
#' North Norfolk, Eastern England, were screened for the species. Fyke net
#' and historical data on presence-absence and density of crucian carp is 
#' available for these ponds. Method sensitivity was compared and eDNA analysis
#' evaluated as a non-invasive alternative to distribution and abundance 
#' estimation for this species.
#' 
#' 
#' ## Prepare working environment
#' 
#' Clear R memory, set working directory and load required packages. Then,
#' load functions for calculating kappa coefficient, plotting model 
#' residuals and testing model fit.
#' 

## Clear memory
rm(list=ls())

## Direct R to folder containing packages on University workstation
#.libPaths(c("C:\\R\\rlib", .libPaths("rlib")))

## set working directory to the location of the script
# install.packages("rstudioapi") # first time for each computer
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## Check working directory
getwd()

## Load required packages
p <- c("ggplot2","munsell","lazyeval","grid","gridExtra","lme4","glmmADMB",
       "coda","MASS","car","scales","AICcmodavg","xtable","gtools",
       "xlsx","reshape2","dplyr","plyr","arm","RVAideMemoire","permute",
       "ResourceSelection","bbmle","RColorBrewer", "MuMIn", "rpart",
       "LMERConvenienceFunctions","ggmap","mapproj","geosphere","jpeg",
       "proto","rjson","RgoogleMaps","maps","labeling","ggsn","png",
       "coin","modeltools","mvtnorm","pROC","car")
new.packages <- p[!(p %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "https://cran.ma.imperial.ac.uk/",
                                          dependencies=TRUE)

lapply(p, require, character.only = TRUE)

## Load custom functions
f <- c("CheckResidsFunction.R", "OverdispersalFunction.R", 
       "CheckConvergenceFunction.R", "HighstatLibV6.R",
       "MyLibrary.R", "glmmML_pres_function.R")

lapply(f, source)


#'
#' To ensure reproducibility, print details about the version of R being
#' used for analysis.
#' 

sessionInfo()


#'
#' ## qPCR dataset
#' 
#' Examine data.
#' 

density_df <- read.csv("../Pond data/Density_analysis.csv", header=TRUE)

summary(density_df)
head(density_df)
names(density_df)
str(density_df)

## Replace NA values with 0
density_df[is.na(density_df)] <- 0



#'---
#'
#' ## 1) Plot species distribution
#' 

## Check coordinate points plot ok with a ggplot
gtest <- ggplot(density_df, aes(x=Longitude,y=Latitude)) + geom_point(aes(size=Copy_Number))
gtest

## Now we want to apply them to a map. First get the map by downloading
## the Google road map for Norfolk. The code below will fetch the
## map from Google, centred on a specific latitude and longitude.
## A road map is just one option - google if you want a different 
## option.

## Find latitudes and longitudes to be centre of distribution map
mean(density_df$Latitude)
mean(density_df$Longitude)

## Install developmental version of 'ggmap' package
install.packages("devtools", dependencies = TRUE)
devtools::install_github("dkahle/ggmap", force=TRUE)
library(ggmap)

## Get Norfolk map
Norfolk <- get_map(location = c(lon = 1.054849, lat = 52.82313), 
                   zoom = 10, 
                   maptype = 'terrain',
                   color = 'bw')
ggmap(Norfolk)

## Now plot ponds on the map.
## Here, there are points at each site, coloured by presence-absence of 
## crucian carp
g1 <- ggmap(Norfolk, extent = "device", crop = T)
g1 <- g1 + geom_point(data = density_df, 
                      aes(x=Longitude, y=Latitude, fill=Stocked), 
                      size=3, pch=22)
g1 <- g1 + scale_fill_manual(name = "Crucian carp",
                             breaks = c("N","Y"),
                             labels = c("Absent","Present"),
                             values = c("grey30","white"))
g1 <- g1 + theme(legend.justification=c(0,0), 
                 legend.position=c(0.75,0.85),
                 legend.text = element_text(size=20, angle=0),
                 legend.title = element_text(size=20, angle=0),
                 legend.background = element_rect(fill="white",
                                                  size=0.5, 
                                                  linetype="solid",
                                                  colour ="black"),
                 legend.key = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA, size=1))
g1 <- g1 + scalebar(x.min = 0.75549, x.max = 1.343829,
                    y.min = 52.67828, y.max = 52.91184,
                    dd2km = TRUE, model = 'WGS84', 
                    location = "topleft", anchor = c(x=0.7,y=53.05), 
                    dist = 10, height = 0.02,
                    st.dist = 0.05, st.bottom = FALSE, st.size = 5)
g1



#' ---
#' 
#' ## 2) Crucian carp presence-absence
#'
#' First, compare agreement between eDNA analysis and fyke netting for 
#' crucian carp detection
#' 

## Calculate proportion of negative and positive detections by each 
## method
net_table <- table(density_df$Stocked)
prop <- prop.table(net_table) 
net_prop <- data.frame(prop)

eDNA_table <- table(density_df$eDNA)
prop <- prop.table(eDNA_table) 
eDNA_prop <- data.frame(prop)

## Create new dataframe of proportions
method_prop <- cbind(net_prop, eDNA_prop[,2])
colnames(method_prop) <- c("Crucian_carp","Netting","eDNA")

## Use melt in reshape package to show positive and negative proportion 
## for each method and enable plotting proportion by method and coloured 
## by carp detection
method_prop <- melt(method_prop, id.vars='Crucian_carp')
colnames(method_prop) <- c("Crucian_carp","Method","Proportion")

## Add frequencies of positive and negative ponds by each method in new 
## column
method_prop$Frequency <- c(10,10,11,9)

## Plot of proportion of ponds with and without crucian carp as 
## determined by netting and eDNA
p1 <- ggplot(method_prop[order(method_prop$Crucian_carp,decreasing=T),],
             aes(x=Method, y=Proportion, fill=Crucian_carp))
p1 <- p1 + geom_bar(stat="identity", colour="black")
p1 <- p1 + geom_text(aes(label=Frequency), size=8, position=position_fill(), vjust=1.5)
p1 <- p1 + scale_y_continuous(labels = percent_format(),
                              expand = c(0, 0), limits = c(0, 1.05))
p1 <- p1 + labs(x="Method", y="Ponds surveyed")
p1 <- p1 + scale_fill_manual(name="Crucian carp",
                             values=c("grey50","deepskyblue2"),
                             breaks=c("N","Y"),
                             labels=c("Absent",
                                      "Present"))
p1 <- p1 + theme(panel.background = element_rect(fill = 'white'),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 plot.title = element_text(face="bold"),
                 text = element_text(size=26))
p1



#' ---
#' 
#' ## 2) Biological replicate variation
#' 
#' Five biological replicates were taken for each pond. Examine the 
#' variation in detection and copy number amongst replicates.
#' 

## Import data
rep_df <- read.csv("../Pond data/Replicates_analysis.csv", header=TRUE)

## Examine data
summary(rep_df)
head(rep_df)
names(rep_df)
str(rep_df)

## Make replicate, or sample, a factor
rep_df$Replicate <- as.factor(rep_df$Replicate)
str(rep_df)

## Subset data for only ponds containing crucian carp
crucian <- subset(rep_df, rep_df$Crucian == "Y")

## Order levels of qPCR by crucian carp positive then crucian carp negative
crucian$qPCR <- relevel(crucian$qPCR, 'Y')

## Plot number of positive replicates for crucian carp across sites
p2a <- ggplot(crucian, aes(x = Site, fill = qPCR)) 
p2a <- p2a + geom_bar(width=0.5) + ggtitle('(a)')
p2a <- p2a + labs(x="", y="Number of biological replicates")
p2a <- p2a + scale_y_continuous(expand = c(0,0), limits = c(0,5), breaks = seq(0, 5, by = 1))
p2a <- p2a + scale_fill_manual(values=c("deepskyblue2","grey50"),
                               name="Crucian carp",
                               breaks=c("Y","N"),
                               labels=c("Positive",
                                        "Negative"))
p2a <- p2a + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20),
                   legend.position="bottom")
p2a

## Plot variation in copy number by site
p2b <- ggplot(crucian, aes(x = Site, y = copy_no))
p2b <- p2b + ggtitle('(b)')
p2b <- p2b + geom_jitter(aes(x=Site, y=copy_no),
                         colour="deepskyblue2", 
                         shape=19, cex=2, width=0.3,
                         show.legend=FALSE)
p2b <- p2b + geom_boxplot(outlier.colour="black", 
                          outlier.shape=19, 
                          outlier.size=2,
                          outlier.alpha=1,
                          show.legend=FALSE,
                          alpha = 0.5)
p2b <- p2b + labs(x="Site", y=expression(paste("DNA copy number (copies/",mu,"l)")))
p2b <- p2b + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20))
p2b

## Produce plots as one figure
library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(p2a), ggplotGrob(p2b), size = "last"))


## What caused variation in copy number of individual water samples?
## Variables include: 
## - volume of water filtered
## - number of filter papers used
## - eDNA concentration after DNA extraction
## - Water conditions
## 

## Explore data and potential relationships with DNA copy number
hist(crucian$copy_no)
plot(crucian$copy_no ~ crucian$Vol_filtered)
plot(crucian$copy_no ~ crucian$No_filters)
plot(crucian$copy_no ~ crucian$Qubit)
plot(crucian$copy_no ~ crucian$Sediment)
plot(crucian$copy_no ~ crucian$Algal)
plot(crucian$copy_no ~ crucian$Duckweed)

## The distribution of the response variable, DNA copy number, is 
## heavily skewed to the left which indicates our data contains a
## substantial number of 0 values. This could create problems for
## statistical modeling. 
## DNA copy number appears to have relationships with qubit 
## concentration and presence of duckweed in water samples.


## First, we need to identify co-varying independent variables
## Examine matrix of Spearman's rank correlations for dataset as no 
## assumptions are made about linearity in a relationship between two 
## variables
cor(crucian[, c(6:8)], method="spearman")

## Pairwise plot but with strength of collinearity indicated text size
## Booth et al. (1994) and Zuur et al. (2009) suggest that correlations
## between pairs of variables with magnitudes +/-0.5 indicate high
## collinearity, although others say +/-0.3
plot(crucian[, c(6:8,11:13)], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## This suggests volume of water filtered and the number of filters used
## are collinear. There is also strong collinearity between algal content of 
## samples with volume filtered and number of filters. Volume filtered 
## and number of filters are minorly collinear with most other variables.

## Fit model with all variables and calculate VIF values to further 
## assess collinearity
glmBioRep <- glm(copy_no ~ Vol_filtered + No_filters + Qubit + Duckweed 
                 + Sediment + Algal,
                 family = poisson,
                 data = crucian)

## Use vif() function in car package to calculate VIFs
vif(glmBioRep)
sqrt(vif(glmBioRep)) > 2

## Error produced about aliased coefficients which suggests there is 
## perfect multicollinearity in dataset

## Remove number of filters and see if that resolves collinearity
glmBioRep <- glm(copy_no ~ Vol_filtered + Qubit + Duckweed 
                 + Sediment + Algal,
                 family = poisson,
                 data = crucian)

vif(glmBioRep)
sqrt(vif(glmBioRep)) > 2

## VIF values less than 3 thus removal of number of filters has resolved any
## collinearity.


## Classification trees allow detailed investigation into the relative
## importance of explanatory variables
f1 <- formula(copy_no ~ Vol_filtered + Qubit + Algal + Sediment + Duckweed)
tree <- rpart(f1, data = crucian, method = "class",
                  minsplit=5, cp = 0.001)

par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(tree,uniform=TRUE, margin=0.1)
text(tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))

## Volume filtered appears to be strongest predictor of DNA copy 
## number, followed by qubit concentration, algal, sediment and 
## duckweed content of water.


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))

## Tree indicates 2 explanatory variables will give optimal model fit.


## Fit nested GLMM to data to examine variation in DNA copy number across
## sample replicates and account for spatial autocorrelation in dataset.
CNmodel <- glmer(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit
                 + Sediment + Duckweed + Algal,
                 family = "poisson",
                 data = crucian)

summary(CNmodel)
anova(CNmodel)
drop1(CNmodel, test = "Chi")   # for binomial/integer y models, 
                               # statistics for significance of each term
display(CNmodel)
se.ranef(CNmodel)              # levels of random factor centred around 0


## Manually test for overdispersion
419.7/42
1-pchisq(419.7, df=42)

## Model appears to be overdispersed.
## This happens when the observed and predicted variances of the 
## dependent variable are the same. Generally, more covariates, 
## zero-inflated distribution, random terms, or alternative link 
## function are needed to counter overdispersion. 

## However, model summary() can give an inflated value for residual 
## deviance. Consequently, the classic overdispersion statistic
## (= residual deviance/residual df) can be unreliable for GLMM.
## Use customised function to check for overdispersion
overdisp_fun(CNmodel)

## Model not overdispersed
## Ratio reported is the same as the overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(crucian$copy_no ~ fitted(CNmodel))

## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(CNmodel)

## 61.81% variance explained by fixed and random effects
## 44.29% variance explained by fixed effects only

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(crucian$copy_no, fitted(CNmodel))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data


## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
chkres(CNmodel)
sresid <- resid(CNmodel, type = "pearson")
shapiro.test(sresid)   # P = 5.437e-08

## Q-Q plot and heteroscedasticity plot
mcp.fnc(CNmodel)

## Substantial deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.
## Discernible pattern can be seen with data points clustered on one side
## Also getting some negative predicted values therefore there appears
## to be heteroscedascity in model. Plot standardised residuals against
## each independent variable to identify source of heterogeneity i.e.
## independent variable that is non-linearly associated with y
plot(sresid ~ crucian$Qubit)     
plot(sresid ~ crucian$Duckweed)  # heteroscedastic
plot(sresid ~ crucian$Sediment)  # heteroscedastic    
plot(sresid ~ crucian$Algal)

## Assumption 3: no collinearity
## All collinearity was removed during preliminary examination of
## data.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(crucian$Vol_filtered)
hi <- (2*2)/50
plot(leverage(crucian$Vol_filtered), type = "h")
abline(0.08, 0, lty = 2)
points(leverage(crucian$Vol_filtered))

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(crucian$Qubit)
hi <- (2*2)/50
plot(leverage(crucian$Qubit), type = "h")
abline(0.08, 0, lty = 2)
points(leverage(crucian$Qubit))


## All model validation checks indicate the model is not a good fit to 
## the data. If assumption of normal residuals is violated, transformation
## of the dependent variable is sometimes recommended. However, this has
## also been advised against as it causes effects to be seen in data that 
## may not really be there. If assumption of homoscedasticity is violated,
## changing the distribution family of the model may improve fitting. In
## the case of count data, quasi-Poisson or negative binomial families
## are recommended. We will try a negative binomial distribution.
CNmodel.nb <- glmer.nb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                       + Sediment + Duckweed + Algal,
                       data = crucian)

summary(CNmodel.nb)
anova(CNmodel.nb)
drop1(CNmodel.nb, test = "Chi")   # for binomial/integer y models, 
                                  # statistics for significance of each term

display(CNmodel.nb)
se.ranef(CNmodel.nb)             # levels of random factor centred around 0


## Test for overdispersion
419.7/41
1-pchisq(419.7, df=41)    # overdispersed...
overdisp_fun(CNmodel.nb)  # not overdispersed

## Plot the fitted data against the observed data
plot(crucian$copy_no ~ fitted(CNmodel.nb))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(crucian$copy_no, fitted(CNmodel.nb))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
chkres(CNmodel.nb)
sresid <- resid(CNmodel.nb, type = "pearson")
shapiro.test(sresid)   # P = 5.467e-08

## Q-Q plot and heteroscedasticity plot
mcp.fnc(CNmodel.nb)

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.
## Discernible pattern can be seen with data points clustered on one side
## Also getting some negative predicted values therefore there appears
## to be heteroscedascity in model. Plot standardised residuals against
## each independent variable to identify source of heterogeneity i.e.
## independent variable that is non-linearly associated with y
plot(sresid ~ crucian$Qubit)   
plot(sresid ~ crucian$Duckweed) # heteroscedastic
plot(sresid ~ crucian$Sediment) # heteroscedastic      
plot(sresid ~ crucian$Algal)

## Assumption 3: no collinearity
## All collinearity removed earlier.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(crucian$Vol_filtered)
hi <- (2*2)/50
plot(leverage(crucian$Vol_filtered), type = "h")
abline(0.08, 0, lty = 2)
points(leverage(crucian$Vol_filtered))

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(crucian$Qubit)
hi <- (2*2)/50
plot(leverage(crucian$Qubit), type = "h")
abline(0.08, 0, lty = 2)
points(leverage(crucian$Qubit))


## All model validation checks indicate that model with negative binomial
## distribution is not a good fit to the data either. There are a lot of 
## zeros in our data which may be causing problems for model fitting. 
## Therefore, we will try zero-inflated Poisson and negative binomial
## distributions to see if they improve model fit.
## We can also model our data using "hurdle models" which model all zeros 
## together rather than treating zeros as true or false zeros. Try both 
## Poisson and negative binomial.
## Use glmmADMB package to fit zero-inflated models. We will also refit
## the Poisson and negative binomial models using glmmadmb for comparison
## with zero-inflated.
CNmodel <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                    + Sediment + Duckweed + Algal,
                    data = crucian,
                    family = "poisson")

CNmodel.nb <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                       + Sediment + Duckweed + Algal,
                       data = crucian,
                       family = "nbinom")

CNmodel.zip <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                        + Sediment + Duckweed + Algal,
                        data = crucian,
                        zeroInflation = TRUE,
                        family = "poisson")

CNmodel.zinb <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                         + Sediment + Duckweed + Algal,
                         data = crucian,
                         zeroInflation = TRUE,
                         family = "nbinom")

CNmodel.hzip <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                         + Sediment + Duckweed + Algal,
                         data = crucian,
                         family = "truncpoiss")
 
CNmodel.hzinb <- glmmadmb(copy_no ~ (1|Site/Replicate) + Vol_filtered + Qubit 
                          + Sediment + Duckweed + Algal,
                          data = crucian,
                          family = "truncnbinom")

## Model summaries:
summary(CNmodel)
summary(CNmodel.nb)
summary(CNmodel.zip)   # did not fit
summary(CNmodel.zinb)  # did not fit
summary(CNmodel.hzip)  # did not fit
summary(CNmodel.hzinb) 
## Hurdle models did not fit to data

## Estimates of significance for each model
drop1(CNmodel, test = "Chi")       # did not estimate properly
drop1(CNmodel.nb, test = "Chi")    # duckweed and algal significant
drop1(CNmodel.hzinb, test = "Chi") # did not estimate properly

## Test models for overdispersion
overdisp_fun(CNmodel)       # overdispersed
overdisp_fun(CNmodel.nb)    # not overdispersed
overdisp_fun(CNmodel.hzinb) # NA values produced

## Plot the fitted data against the observed data
plot(crucian$copy_no ~ fitted(CNmodel))       # poor fit
plot(crucian$copy_no ~ fitted(CNmodel.nb))    # reasonable fit
plot(crucian$copy_no ~ fitted(CNmodel.hzinb)) # poor fit

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(crucian$copy_no, fitted(CNmodel))       # poor fit
hoslem.test(crucian$copy_no, fitted(CNmodel.nb))    # good fit
hoslem.test(crucian$copy_no, fitted(CNmodel.hzinb)) # good fit?

## Perform model validation checks to see if models are good fit to 
## data and making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
sresid.p <- resid(CNmodel)
hist(sresid.p, freq=F)              # non-normal distribution
qqnorm(sresid.p, cex=1.8, pch=20)   
qqline(sresid.p, lty=2, lwd=2)      # non-normal, massive tail
shapiro.test(sresid.p)              # data not normal

sresid.nb <- resid(CNmodel.nb)
hist(sresid.nb, freq=F)              # residuals look better
qqnorm(sresid.nb, cex=1.8, pch=20)   # again, residuals look more normal
qqline(sresid.nb, lty=2, lwd=2)
shapiro.test(sresid.nb)              # still significant but larger p-value

sresid.hzinb <- resid(CNmodel.hzinb)
hist(sresid.hzinb, freq=F)              # could not plot
qqnorm(sresid.hzinb, cex=1.8, pch=20)   # could not plot
qqline(sresid.hzinb, lty=2, lwd=2)
shapiro.test(sresid.hzinb)              # could not calculate

## Compare AIC of models fitted
library(bbmle)
AICtab(CNmodel, CNmodel.nb, CNmodel.hzinb)

## AIC indicates hurdle zero-inflated negative binomial model is best
## fit to data but model validation checks indicate that the negative binomial
## model is the best fit.


## There is an alternative package for fitting zero-inflated models,
## called glmmTMB.
## The ziformula argument, or zi, describes how the probability of an 
## extra zero (i.e. structural zero) will vary with predictors.
## We might assume that zero-inflation will at least vary according to 
## site, replicate, and concentration of extracted DNA.
## The zero-inflation probability is always modeled with logit-link to 
## keep it between 0 and 1.
library(glmmTMB)
altCNmodel <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                      + Qubit + Sediment + Duckweed + Algal,
                      data = crucian,
                      family = poisson)

altCNmodel.nb <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                         + Qubit + Sediment + Duckweed + Algal,
                         data = crucian,
                         family = nbinom2)

altCNmodel.zip <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                          + Qubit + Sediment + Duckweed + Algal,
                          data = crucian,
                          zi = ~ (1|Site/Replicate) + Qubit,
                          family = poisson)

altCNmodel.zinb <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                           + Qubit + Sediment + Duckweed + Algal,
                           data = crucian,
                           zi = ~ (1|Site/Replicate) + Qubit,
                           family = nbinom2)

altCNmodel.hzip <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                           + Qubit + Sediment + Duckweed + Algal,
                           data = crucian,
                           ziformula = ~ (1|Site/Replicate) + Qubit,
                           family = truncated_poisson)

altCNmodel.hzinb <- glmmTMB(copy_no ~ (1|Site/Replicate) + Vol_filtered
                            + Qubit + Sediment + Duckweed + Algal,
                            data = crucian,
                            ziformula = ~ (1|Site/Replicate) + Qubit,
                            family = truncated_nbinom2)

## Model summaries:
summary(altCNmodel)
summary(altCNmodel.nb)
summary(altCNmodel.zip)
summary(altCNmodel.zinb)  
summary(altCNmodel.hzip)
summary(altCNmodel.hzinb)   # did not fit properly

## Estimates of significance for each model
drop1(altCNmodel, test = "Chi")       # qubit, sediment and duckweed significant
drop1(altCNmodel.nb, test = "Chi")    # duckweed and algal significant
drop1(altCNmodel.zip, test = "Chi")   # all significant bar duckweed, convergence problem
drop1(altCNmodel.zinb, test = "Chi")  # all significant bar duckweed, convergence problem
drop1(altCNmodel.hzip, test = "Chi")  # volume of water significant
drop1(altCNmodel.hzinb, test = "Chi") # volume of water significant, convergence problem

## Zero-inflated negative binomial model did not fit correctly - 
## possibly due to overparameterisation, i.e. zero observations may not 
## depend on variables specified but another variable missing from model. 
## Models with convergence problems (non-positive definite Hessian 
## matrices) should be excluded from further consideration in
## general.

## Plot the fitted data against the observed data
plot(crucian$copy_no ~ fitted(altCNmodel))       # good fit
plot(crucian$copy_no ~ fitted(altCNmodel.nb))    # reasonable fit
plot(crucian$copy_no ~ fitted(altCNmodel.hzip))  # good fit

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(crucian$copy_no, fitted(altCNmodel))       # good fit
hoslem.test(crucian$copy_no, fitted(altCNmodel.nb))    # good fit
hoslem.test(crucian$copy_no, fitted(altCNmodel.hzip))  # good fit

## Perform model validation checks to see if models are good fit to 
## data and making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
sresid.altp <- resid(altCNmodel)
hist(sresid.altp, freq=F)              # normal distribution
qqnorm(sresid.altp, cex=1.8, pch=20)   
qqline(sresid.altp, lty=2, lwd=2)      # mostly normal, slight tail
shapiro.test(sresid.altp)              # data normal

sresid.altnb <- resid(altCNmodel.nb)
hist(sresid.altnb, freq=F)             # normal distribution
qqnorm(sresid.altnb, cex=1.8, pch=20)  # not normal, tails at both ends
qqline(sresid.altnb, lty=2, lwd=2)
shapiro.test(sresid.altnb)             # data not normal

sresid.althzip <- resid(altCNmodel.hzip)
hist(sresid.althzip, freq=F)              # semi-normal distribution
qqnorm(sresid.althzip, cex=1.8, pch=20)   # not normal, tails at both ends
qqline(sresid.althzip, lty=2, lwd=2)
shapiro.test(sresid.althzip)              # data not normal

## Compare AIC of all models fitted using glmmTMB
AICtab(altCNmodel, altCNmodel.nb, altCNmodel.hzip)

## AIC indicates that the hurdle zero-inflated Poisson model is best, whereas
## model validation checks suggest Poisson model is best fit to data

## Final model:
summary(altCNmodel)
drop1(altCNmodel, test = "Chi")


## PLOT MODEL FIT
## Obtain predicted values for full data set: crucian
## se.fit = TRUE will obtain the standard error for each of these 
## predictions
fit <- data.frame(predict(altCNmodel, newdata=crucian, se.fit = TRUE))

## Create new data set with fitted values and original data for plots 1 
## and 2
box.dat <- cbind(crucian, fit)

## Plot DNA copy number against Qubit concentration
## Create a range of Qubit concentration values which increase by 2.285
## to predict DNA copy number as concentration increases
range <- seq(from=min(crucian$Qubit), 
             to=max(crucian$Qubit), by=2.285)

## Create new data frame where only concentration changes
## Factors set as one level and covariates set to mean
d1 <- data.frame(Site=rep('SKEY1', length(range)),
                 Replicate=rep('1', length(range)),
                 Duckweed=rep("Y", length(range)),
                 Sediment=rep("Y", length(range)),
                 Algal=rep("Y", length(range)),
                 Vol_filtered=rep(0.925, length(range)),
                 Qubit=range)

## Get predictions for this new dataset
fit <- predict(altCNmodel, newdata=d1, se.fit=TRUE) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## New data frame for plotting
dat.conc <- cbind(d1, fit, ciu, cil)

## PLOT 3A: Qubit concentration and DNA copy number
p3a <- ggplot() + ggtitle('(a)')
p3a <- p3a + geom_jitter(aes(x=crucian$Qubit, y=crucian$copy_no), colour="black", cex=1.5)
p3a <- p3a + geom_line(aes(x=dat.conc$Qubit, y=dat.conc$fit), size = 1)
p3a <- p3a + geom_ribbon(aes(x=dat.conc$Qubit, ymin = dat.conc$cil, ymax = dat.conc$ciu), alpha = 0.2)
p3a <- p3a + scale_colour_gradient(limits=c(0, 500), low="gray30", high="gray80")
p3a <- p3a + coord_cartesian(xlim=c(0,125),ylim=c(0,500))
p3a <- p3a + scale_x_continuous(breaks=seq(0,125,25))
p3a <- p3a + labs(x=expression(paste("Qubit concentration (ng/",mu,"l)")),
                  y=expression(paste("DNA copy number (copies/",mu,"l)")))
p3a <- p3a + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20))
p3a


## PLOT 3B: Duckweed content and DNA copy number
p3b <- ggplot(box.dat) + ggtitle('(b)')
p3b <- p3b + geom_jitter(aes(x=Duckweed, y=copy_no), colour="black", cex=1.5, width=0.2)
p3b <- p3b + geom_boxplot(aes(x=Duckweed, y=fit), outlier.colour = "red", alpha=0.8, show.legend=FALSE)
p3b <- p3b + labs(x="Duckweed present", y="")
p3b <- p3b + scale_x_discrete(labels = c("No", "Yes"))
p3b <- p3b + coord_cartesian(ylim=c(0,500))
p3b <- p3b + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20))
p3b


## PLOT 4C: Sediment and DNA copy number
p3c <- ggplot(box.dat) + ggtitle('(c)')
p3c <- p3c + geom_jitter(aes(x=Sediment, y=copy_no), colour="black", cex=1.5, width=0.2)
p3c <- p3c + geom_boxplot(aes(x=Sediment, y=fit), outlier.colour = "red", alpha=0.8, show.legend=FALSE)
p3c <- p3c + labs(x="Sediment present", y="")
p3c <- p3c + scale_x_discrete(labels = c("No", "Yes"))
p3c <- p3c + coord_cartesian(ylim=c(0,500))
p3c <- p3c + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20))
p3c


## Both plots
grid.newpage()
grid.draw(cbind(ggplotGrob(p3a), ggplotGrob(p3b), ggplotGrob(p3c), size = "last"))




#' ---
#' 
#' 3) Influences of external factors on eDNA detection in ponds
#' 
#' Examine variation in crucian carp copy number attributable to density,
#' physicochemical properties of ponds and surrounding habitat. These
#' variables may increase habitat suitability of ponds for crucian carp
#' and density must be included as well to account for this.
#' 
#' Variables include: 
#' - CPUE
#' - Conductivity
#' - Alkalinity
#' - Temperature
#' - pH
#' - Surface oxygen
#' - submerged macrophyte (%)
#' - shading (%)
#' 

metadata <- read.csv("../Pond data/CrucianCarpDeterminants.csv", header=TRUE)

summary(metadata)
head(metadata)
names(metadata)
str(metadata)

## Exploratory plots
hist(metadata$Copy_Number)
plot(metadata$Copy_Number ~ metadata$CPUE)
plot(metadata$Copy_Number ~ metadata$Conductivity)
plot(metadata$Copy_Number ~ metadata$Alkalinity)
plot(metadata$Copy_Number ~ metadata$pH)
plot(metadata$Copy_Number ~ metadata$Temperature)
plot(metadata$Copy_Number ~ metadata$O2_surface)
plot(metadata$Copy_Number ~ metadata$perc_sub_macrophyte)
plot(metadata$Copy_Number ~ metadata$per_shading)


## The distribution of the response variable, DNA copy number, is 
## heavily skewed to the left which indicates our data contains a
## substantial number of 0 values. This could create problems for
## statistical modeling. 
## DNA copy number may have relationships with CPUE, temperature, 
## surface oxygen, conductivity and submerged macrophytes.

## First, we identify co-varying independent variables
## Examine matrix of Spearman's rank correlations for dataset as no 
## assumptions are made about linearity in a relationship between two 
## variables
cor(metadata[, 3:11], method="spearman")

## Pairwise plot but with strength of collinearity indicated text size
## Booth et al. (1994) and Zuur et al. (2009) suggest that correlations
## between pairs of variables with magnitudes +/-0.5 indicate high
## collinearity, although others say +/-0.3
plot(metadata[, 3:11], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## There is strong collinearity between:
## - no. individuals and CPUE
## - conductivity and alkalinity
## - no. individuals and CPUE with temperature
## - no. individuals, CPUE, pH, temperature, and shading with surface 
##   oxygen
## - pH and shading

## Including only the variables that are expected to be most relevant to
## DNA copy number
plot(metadata[, c(4:5,7:8,10)], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)


## Now check the variance inflation factors (VIFs) of variables to 
## assess the extent of remaining collinearity.

## Fit all explanatory variables to the data, then calculate 
## VIFs for each variable from the resulting model. Variables
## with VIF > 3 are cause for concern.
glmCN <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
             + pH + Temperature + O2_surface 
             + perc_sub_macrophyte + per_shading,
             family = poisson,
             data = metadata)

## Use vif() function in car package to calculate VIFs
vif(glmCN)
sqrt(vif(glmCN)) > 2


## Extreme multicollinearity present in dataset
## Drop surface oxygen as it has highest VIF value

glmCN <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
             + pH + Temperature + perc_sub_macrophyte 
             + per_shading,
             family = poisson,
             data = metadata)
vif(glmCN)
sqrt(vif(glmCN)) > 2


## Drop percentage of shading as it is variable we are least interested
## in with the next highest VIF value

glmCN <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
             + pH + Temperature + perc_sub_macrophyte,
             family = poisson,
             data = metadata)
vif(glmCN)
sqrt(vif(glmCN)) > 2


## All VIFs now below 10, which would be acceptable by some statisticians.
## However, Zuur et al. (2009) recommend removing variables until all
## variables have VIF<3.
## Remove Alkalinity as it is variable we are least interested in with 
## the next highest VIF value.

glmCN <- glm(Copy_Number ~ CPUE + Conductivity + pH + Temperature 
             + perc_sub_macrophyte,
             family = poisson,
             data = metadata)
vif(glmCN)
sqrt(vif(glmCN)) > 2

## All VIFs now below 5, which is more acceptable although CPUE and 
## temperature have VIF>3.
## CPUE is more likely to influence eDNA concentration than temperature
## thus exclude temperature from the model

glmCN <- glm(Copy_Number ~ CPUE + Conductivity + pH 
             + perc_sub_macrophyte,
             family = poisson,
             data = metadata)
vif(glmCN)
sqrt(vif(glmCN)) > 2


## Classification trees allow detailed investigation into the relative
## importance of explanatory variables
f1 <- formula(Copy_Number ~ CPUE + Conductivity 
              + pH + perc_sub_macrophyte)

crucian_tree <- rpart(f1, data = metadata, method = "class",
                      minsplit=5, cp = 0.001)

par(xpd = NA, mar = c(1.5, 1.5, 1.5, 1.5))
plot(crucian_tree,uniform=TRUE, margin=0.1)
text(crucian_tree, use.n = T, cex = 1.0)
par(xpd = F, mar = c(4.5, 4.5, 0.5, 0.5))

## Conductivity appears to be strongest predictor of DNA copy number,
## followed by CPUE.


## Cross-validation (pruning tree)
par(mar = c(4.5, 4.5, 4.5, 0.5))
plotcp(crucian_tree)
par(mar = c(4.5, 4.5, 0.5, 0.5))

## A tree of size of 1 or 3 is optimal i.e. 3 explanatory variables. 
## The best tree will be subjective.


## Apply GLMM to examine variation in DNA copy number accounting for 
## spatial autocorrelation in dataset by modeling site as a random
## factor.

envm <- glmer(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
              + perc_sub_macrophyte,
              family = "poisson",
              data = metadata)
summary(envm)
anova(envm)
drop1(envm, test = "Chi") # for binomial/integer y models, 
                          # statistics for significance of each term
display(envm)
se.ranef(envm) 

## Overdispersion test (chi-square)
85.4/4
1-pchisq(85.4, df=4)

## Model appears to be overdispersed.
## This happens when the observed and predicted variances of the 
## dependent variable are the same. Generally, more covariates, 
## zero-inflated distribution, random terms, or alternative link 
## function are needed to counter overdispersion. 

## However, model summary() can give an inflated value for residual 
## deviance. Consequently, the classic overdispersion statistic
## (= residual deviance/residual df) can be unreliable for GLMM.
## Use customised function to check for overdispersion.
overdisp_fun(envm)   

## Model not overdispersed
## Ratio reported is equivalent to the overdispersion parameter
## 1 = no overdispersion, >2 = excessive overdispersion

## Plot the fitted data against the observed data
plot(metadata$Copy_Number ~ fitted(envm))

## Get R-squared value
## marginal R-squared = proportion of variance in response variable explained
## by fixed variables only
## conditional R-squared = proportion of variance in response explained by
## fixed and random variables
r.squaredGLMM(envm)

## 79.09% variance explained by fixed and random effects
## 79.09% variance explained by fixed effects only

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(metadata$Copy_Number, fitted(envm))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data


## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
chkres(envm)
sresid <- resid(envm, type = "pearson")
shapiro.test(sresid)  # P = 0.0319

## Q-Q plot and heteroscedasticity plot
mcp.fnc(envm)

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.
## Discernible pattern can be seen with data points clustered on one side
## Also getting some negative predicted values therefore there appears
## to be heteroscedascity in model. Plot standardised residuals against
## each independent variable to identify source of heterogeneity i.e.
## independent variable that is non-linearly associated with y.
plot(sresid ~ metadata$CPUE)
plot(sresid ~ metadata$Conductivity)
plot(sresid ~ metadata$pH)
plot(sresid ~ metadata$perc_sub_macrophyte)

## Possibly some heteroscedascity from macrophyte cover and conductivity

## Assumption 3: no collinearity
## All collinearity was removed during preliminary examination of
## data.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$CPUE)
hi <- (2*2)/10
plot(leverage(metadata$CPUE), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$CPUE)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$Conductivity)
hi <- (2*2)/10
plot(leverage(metadata$Conductivity), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$Conductivity)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$pH)
hi <- (2*2)/10
plot(leverage(metadata$pH), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$pH)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$perc_sub_macrophyte)
hi <- (2*2)/10
plot(leverage(metadata$perc_sub_macrophyte), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$perc_sub_macrophyte)


## All model validation checks indicate the model is not a good fit to 
## the data. If assumption of normal residuals is violated, transformation
## of the dependent variable is sometimes recommended. However, this has
## also been advised against as it causes effects to be seen in data that 
## may not really be there. If assumption of homoscedasticity is violated,
## changing the distribution family of the model may improve fitting. In
## the case of count data, quasi-Poisson or negative binomial families
## are recommended. We will try a negative binomial distribution.
envm.nb <- glmer.nb(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
                    + perc_sub_macrophyte,
                    data = metadata)

summary(envm.nb)
anova(envm.nb)
drop1(envm.nb, test = "Chi") # for binomial/integer y models, 
                             # statistics for significance of each term
display(envm.nb)
se.ranef(envm.nb) 

## Test for overdispersion
85.4/3
1-pchisq(85.4, df=3)   # overdispersed...
overdisp_fun(envm.nb)  # not overdispersed...

## Plot the fitted data against the observed data
plot(metadata$Copy_Number ~ fitted(envm.nb))

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(metadata$Copy_Number, fitted(envm.nb))

## Model supposedly fits well as p-value is not significant 
## i.e. no significant difference between the model and the observed data

## Perform model validation checks to ensure model is good fit to data and 
## making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
chkres(envm.nb)
sresid <- resid(envm.nb, type = "pearson")
shapiro.test(sresid) # P = 0.0327, same as Poisson

## Q-Q plot and heteroscedasticity plot
mcp.fnc(envm.nb)

## Some deviation from normality as residuals not normally distributed
## therefore model may not be that reliable.
## Discernible pattern can be seen with data points clustered on one side
## Also getting some negative predicted values therefore there appears
## to be heteroscedascity in model. Plot standardised residuals against
## each independent variable to identify source of heterogeneity i.e.
## independent variable that is non-linearly associated with y
plot(sresid ~ metadata$CPUE)
plot(sresid ~ metadata$Conductivity)
plot(sresid ~ metadata$pH)
plot(sresid ~ metadata$perc_sub_macrophyte)
## Possibly some heteroscedascity from macrophyte cover and conductivity

## Assumption 3: no collinearity
## All collinearity was removed during preliminary examination of
## data.

## Assumption 4: no serial auto-correlation
## can arise if there are unaccounted for relationships in data set or
## if there is confounding effect of time or space
## Use graphical approach with Auto-Correlation Function (ACF)
acf(sresid, main = "Auto-correlation plot")
## Significant autocorrelation at first time lag

## Assumption 5: model not biased by influential observations
## Data point influential if hi<2p/n
leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$CPUE)
hi <- (2*2)/10
plot(leverage(metadata$CPUE), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$CPUE)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$Conductivity)
hi <- (2*2)/10
plot(leverage(metadata$Conductivity), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$Conductivity)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$pH)
hi <- (2*2)/10
plot(leverage(metadata$pH), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$pH)

leverage <- function(x) {1/length(x) + (x-mean(x))^2/sum((x-mean(x))^2)}
leverage(metadata$perc_sub_macrophyte)
hi <- (2*2)/10
plot(leverage(metadata$perc_sub_macrophyte), type = "h")
abline(0.4, 0, lty = 2)
points(metadata$perc_sub_macrophyte)


## All model validation checks indicate that model with negative binomial
## distribution is not a good fit to the data either. There are a lot of 
## zeros in our data which may be causing problems for model fitting. 
## Therefore, we will try zero-inflated Poisson and negative binomial
## distributions to see if they improve model fit.
## We can also model our data using "hurdle models" which model all zeros 
## together rather than treating zeros as true or false zeros. Try both 
## Poisson and negative binomial.
## Use glmmADMB package to fit zero-inflated models. We will also refit
## the Poisson and negative binomial models using glmmadmb for comparison
## with zero-inflated.
envm <- glmmadmb(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
                 + perc_sub_macrophyte,
                 data = metadata,
                 family = "poisson")

envm.nb <- glmmadmb(Copy_Number ~ (1|Site) +  CPUE + pH + Conductivity 
                    + perc_sub_macrophyte,
                    data = metadata,
                    family = "nbinom")

envm.zip <- glmmadmb(Copy_Number ~ (1|Site) +  CPUE + pH + Conductivity 
                     + perc_sub_macrophyte,
                     data = metadata,
                     zeroInflation = TRUE,
                     family = "poisson")

envm.zinb <- glmmadmb(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
                      + perc_sub_macrophyte,
                      data = metadata,
                      zeroInflation = TRUE,
                      family = "nbinom")

envm.hzip <- glmmadmb(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
                      + perc_sub_macrophyte,
                      data = metadata,
                      family = "truncpoiss")

envm.hzinb <- glmmadmb(Copy_Number ~ (1|Site) + CPUE + pH + Conductivity 
                       + perc_sub_macrophyte,
                       data = metadata,
                       family = "truncnbinom")


## Model summaries:
summary(envm)       
summary(envm.nb)
summary(envm.zip)
summary(envm.zinb)
summary(envm.hzip)
summary(envm.hzinb)

## Estimates of significance for each model
drop1(envm, test = "Chi")       # CPUE and pH significant
drop1(envm.nb, test = "Chi")    # CPUE, cond and macrophyte cover significant
drop1(envm.zip, test = "Chi")   # did not fit properly
drop1(envm.zinb, test = "Chi")  # did not fit properly
drop1(envm.hzip, test = "Chi")  # cond significant
drop1(envm.hzinb, test = "Chi") # cond and macrophyte cover significant

## Test models for overdispersion
overdisp_fun(envm)       # overdispersed
overdisp_fun(envm.nb)    # not overdispersed
overdisp_fun(envm.zip)   # overdispersed
overdisp_fun(envm.zinb)  # overdispersed
overdisp_fun(envm.hzip)  # NA values
overdisp_fun(envm.hzinb) # NA values

## Plot the fitted data against the observed data
plot(metadata$Copy_Number ~ fitted(envm))       # poor fit
plot(metadata$Copy_Number ~ fitted(envm.nb))    # good fit
plot(metadata$Copy_Number ~ fitted(envm.zip))   # poor fit
plot(metadata$Copy_Number ~ fitted(envm.zinb))  # poor fit
plot(metadata$Copy_Number ~ fitted(envm.hzip))  # good fit
plot(metadata$Copy_Number ~ fitted(envm.hzinb)) # poor fit

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(metadata$Copy_Number, fitted(envm))       # poor fit
hoslem.test(metadata$Copy_Number, fitted(envm.nb))    # good fit
hoslem.test(metadata$Copy_Number, fitted(envm.zip))   # poor fit
hoslem.test(metadata$Copy_Number, fitted(envm.zinb))  # poor fit
hoslem.test(metadata$Copy_Number, fitted(envm.hzip))  # good fit
hoslem.test(metadata$Copy_Number, fitted(envm.hzinb)) # good fit?

## Perform model validation checks to see if models are good fit to 
## data and making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
sresid.p <- resid(envm)
hist(sresid.p, freq=F)              # non-normal distribution
qqnorm(sresid.p, cex=1.8, pch=20)   
qqline(sresid.p, lty=2, lwd=2)      # non-normal, massive tail
shapiro.test(sresid.p)              # data not normal

sresid.nb <- resid(envm.nb)
hist(sresid.nb, freq=F)              # normal distribution
qqnorm(sresid.nb, cex=1.8, pch=20)   # again, residuals look normal
qqline(sresid.nb, lty=2, lwd=2)
shapiro.test(sresid.nb)              # data normal

sresid.zip <- resid(envm.zip)
hist(sresid.zip, freq=F)              # non-normal distribution
qqnorm(sresid.zip, cex=1.8, pch=20)   # non-normal, massive tail
qqline(sresid.zip, lty=2, lwd=2)
shapiro.test(sresid.zip)              # data not normal

sresid.zinb <- resid(envm.zinb)
hist(sresid.zinb, freq=F)              # non-normal distribution
qqnorm(sresid.zinb, cex=1.8, pch=20)   # non-normal, massive tail
qqline(sresid.zinb, lty=2, lwd=2)
shapiro.test(sresid.zinb)              # data not normal

## PROBLEMS PERFORMING CHECKS FOR HURDLE MODELS

## Compare AIC of models fitted (don't include hurdle models due to 
## problems with fitting and estimation)
AICtab(envm, envm.nb, envm.zip, envm.zinb)

## AIC and model validation checks indicate that negative binomal model 
## is best fit to data.


## There is an alternative package for fitting zero-inflated models,
## called glmmTMB.
## The ziformula argument, or zi, describes how the probability of an 
## extra zero (i.e. structural zero) will vary with predictors.
## We might assume that zero-inflation will at least vary according to 
## site, replicate, and concentration of extracted DNA.
## The zero-inflation probability is always modeled with logit-link to 
## keep it between 0 and 1.
library(glmmTMB)

altenvm <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                   + Conductivity + perc_sub_macrophyte,
                   data = metadata,
                   family = poisson)

altenvm.nb <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                      + Conductivity + perc_sub_macrophyte,
                      data = metadata,
                      family = nbinom2)

altenvm.zip <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                       + Conductivity + perc_sub_macrophyte,
                       data = metadata,
                       zi = ~ (1|Site),
                       family = poisson)

altenvm.zinb <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                        + Conductivity + perc_sub_macrophyte,
                        data = metadata,
                        zi = ~ (1|Site),
                        family = nbinom2)

altenvm.hzip <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                        + Conductivity + perc_sub_macrophyte,
                        data = metadata,
                        ziformula = ~ (1|Site),
                        family = truncated_poisson)

altenvm.hzinb <- glmmTMB(Copy_Number ~ (1|Site) + CPUE + pH
                         + Conductivity + perc_sub_macrophyte,
                         data = metadata,
                         ziformula = ~ (1|Site),
                         family = truncated_nbinom2)

## Model summaries:
summary(altenvm)
summary(altenvm.nb)
summary(altenvm.zip)
summary(altenvm.zinb)
summary(altenvm.hzip)
summary(altenvm.hzinb)

## Estimates of significance for each model
drop1(altenvm, test = "Chi")       # CPUE and cond significant
drop1(altenvm.nb, test = "Chi")    # CPUE, cond, and macrophyte significant
drop1(altenvm.zip, test = "Chi")   # problems estimating
drop1(altenvm.zinb, test = "Chi")  # convergence and estimation problems
drop1(altenvm.hzip, test = "Chi")  # CPUE and macrophyte significant
drop1(altenvm.hzinb, test = "Chi") # CPUE and macrophyte significant

## Plot the fitted data against the observed data
plot(metadata$Copy_Number ~ fitted(altenvm))       # good fit
plot(metadata$Copy_Number ~ fitted(altenvm.nb))    # good fit
plot(metadata$Copy_Number ~ fitted(altenvm.hzip))  # good fit
plot(metadata$Copy_Number ~ fitted(altenvm.hzinb)) # good fit

## Hosmer and Lemeshow Goodness of Fit Test
hoslem.test(metadata$Copy_Number, fitted(altenvm))       # good fit
hoslem.test(metadata$Copy_Number, fitted(altenvm.nb))    # good fit
hoslem.test(metadata$Copy_Number, fitted(altenvm.hzip))  # good fit
hoslem.test(metadata$Copy_Number, fitted(altenvm.hzinb)) # good fit

## Perform model validation checks to see if models are good fit to 
## data and making reliable predictions
## Assumption 1: residuals are normally distributed
## Assumption 2: no heteroscedascity
sresid.p <- resid(altenvm)
hist(sresid.p, freq=F)              # normal distribution
qqnorm(sresid.p, cex=1.8, pch=20)   
qqline(sresid.p, lty=2, lwd=2)      # normal distribution
shapiro.test(sresid.p)              # data normal

sresid.nb <- resid(altenvm.nb)
hist(sresid.nb, freq=F)              # non-normal distribution
qqnorm(sresid.nb, cex=1.8, pch=20)   # not normal, large tail
qqline(sresid.nb, lty=2, lwd=2)
shapiro.test(sresid.nb)              # data not normal

sresid.hzip <- resid(altenvm.hzip)
hist(sresid.hzip, freq=F)              # normal distribution
qqnorm(sresid.hzip, cex=1.8, pch=20)   # normal distribution
qqline(sresid.hzip, lty=2, lwd=2)
shapiro.test(sresid.hzip)              # data normal

sresid.hzinb <- resid(altenvm.hzinb)
hist(sresid.hzinb, freq=F)              # non-normal distribution
qqnorm(sresid.hzinb, cex=1.8, pch=20)   # not normal, large tail
qqline(sresid.hzinb, lty=2, lwd=2)
shapiro.test(sresid.hzinb)              # data not normal

## Compare AIC of models fitted
AICtab(altenvm, altenvm.nb, altenvm.hzip, altenvm.hzinb)

## AIC indicates the zero-inflated Poisson hurdle model is the best fit 
## to the data but model validation checks indicate that the Poisson model
## is best.


## Final model:
summary(altenvm)
drop1(altenvm, test = "Chi")


## PLOT MODEL FIT
## Obtain predicted values for the full data set: crucian
## se.fit = TRUE will obtain the standard error for each of these 
## predictions

## PLOT 4A: CPUE and DNA copy number
## Create a range of CPUE values which increase by 12.6
## to predict copy number as CPUE increases
range <- seq(from=min(metadata$CPUE), 
             to=max(metadata$CPUE), by=12.6)

## Create new data frame where only temperature changes
## Factors are set to first level in dataset
## Continuous variables are set as the mean values in this study
d2 <- data.frame(Site=rep('RAIL', length(range)),
                 CPUE = range,
                 pH = mean(metadata$pH),
                 Conductivity = mean(metadata$Conductivity),
                 perc_sub_macrophyte = mean(metadata$perc_sub_macrophyte))

## Get predictions for this new dataset
fit <- data.frame(predict(altenvm, newdata=d2, se.fit = TRUE,
                          na.action=na.exclude, type="response")) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## Create new data frame
dat.cpue <- cbind(d2, fit, ciu, cil)

## Plot showing relationship between temperature and DNA copy number
## with predicted values from model
p4a <- ggplot() + ggtitle('(a)')
p4a <- p4a + geom_jitter(aes(x=metadata$CPUE, y=metadata$Copy_Number), colour="black", cex=3, height=0.7, width=0.2)
p4a <- p4a + geom_line(aes(x=dat.cpue$CPUE, y=dat.cpue$fit), size=1)
p4a <- p4a + geom_ribbon(aes(x=dat.cpue$CPUE, ymin = dat.cpue$cil, ymax=dat.cpue$ciu), alpha=0.2)
p4a <- p4a + coord_cartesian(xlim=c(0,150), ylim=c(0,250))
p4a <- p4a + labs(x=expression("CPUE estimate"), 
                  y=expression(paste("DNA copy number (copies/",mu,"l)")))
p4a <- p4a + theme(panel.background = element_rect(fill = 'white'),
                 plot.title = element_text(face="bold", hjust = 0),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 legend.title = element_blank(),
                 text = element_text(size=20),
                 legend.position="none")
p4a


## PLOT 4B: Conductivity and DNA copy number
## Create a range of conductivity values which increase by 54
## to predict copy number as conductivity increases
range <- seq(from=min(metadata$Conductivity), 
             to=max(metadata$Conductivity), by=54)

## Create new data frame where only conductivity changes
## Factors are set to first level in dataset
## Continuous variables are set as the mean values in this study
d3 <- data.frame(Site=rep('RAIL', length(range)),
                 CPUE = mean(metadata$CPUE),
                 pH = mean(metadata$pH),
                 Conductivity = range,
                 perc_sub_macrophyte = mean(metadata$perc_sub_macrophyte))

## Get predictions for this new dataset
fit <- data.frame(predict(altenvm, newdata=d3, se.fit = TRUE,
                          na.action=na.exclude, type="response")) 

## Calculate upper 95% CIs from the SEs of the predictions
ciu <- (fit$fit+1.96*fit$se.fit) 

## Calculate lower 95% CIs from the SEs of the predictions
cil <- (fit$fit-1.96*fit$se.fit)

## Create new data frame
dat.cond <- cbind(d3, fit, ciu, cil)

## Plot showing relationship between conductivity and DNA copy number
## with predicted values from model
p4b <- ggplot() + ggtitle('(b)')
p4b <- p4b + geom_jitter(aes(x=metadata$Conductivity, y=metadata$Copy_Number), colour="black", cex=3, height=0.7, width=0.2)
p4b <- p4b + geom_line(aes(x=dat.cond$Conductivity, y=dat.cond$fit), size = 1)
p4b <- p4b + geom_ribbon(aes(x=dat.cond$Conductivity, ymin = dat.cond$cil, ymax = dat.cond$ciu), alpha = 0.2)
p4b <- p4b + coord_cartesian(xlim=c(0,800), ylim=c(0,250))
p4b <- p4b + labs(x=expression(paste("Conductivity (",mu,"s/cm)")),
                  y="")
p4b <- p4b + theme(panel.background = element_rect(fill = 'white'),
                   plot.title = element_text(face="bold", hjust = 0),
                   axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                   axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                   axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                   axis.text.x = element_text(colour="black"),
                   axis.text.y = element_text(colour="black"),
                   legend.title = element_blank(),
                   text = element_text(size=20),
                   legend.position="none")
p4b


## All plots
grid.newpage()
grid.draw(cbind(ggplotGrob(p4a), 
                ggplotGrob(p4b), size = "last"))




#' ---
#' 
#' 4) Occupancy modelling
#' 

## Examine the environmental metadata to be used and check there
## is no collinearity in data set

## View pairwise plots of variables
pairs(metadata[,4:11])

## Below produces a pairwise plot but the strength of collinearity 
## is indicated by the size of the text
## Excessive collinearity generally considered present if r > 0.5.
plot(metadata[,4:11], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## There is strong collinearity between:
## - conductivity and alkalinity
## - no. individuals and CPUE with temperature
## - no. individuals, CPUE, pH, temperature, and shading with surface 
##   oxygen
## - pH and shading

## However, we must account for the influence of crucian carp density 
## (CPUE) on eDNA detection probability. We also need to account for 
## effect of pH on eDNA detection probability at sample level.

## Remove alkalinity, surface oxygen, and percentage of shading
## Examine collinearity
plot(metadata[, c(4:5,7:8,10)], 
     lower.panel = panel.cor, 
     upper.panel = panel.smooth2)

## Removal of these variables would seem to resolve collinearity
## as best we can.


## Now check the variance inflation factors (VIFs) of variables to 
## further assess the extent of collinearity.
## Fit all explanatory variables to the data, then calculate 
## VIFs for each variable from the resulting model. Variables
## with VIF > 3 are cause for concern.
glmEnv <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
              + pH + Temperature + O2_surface
              + perc_sub_macrophyte + per_shading,
              family = poisson,
              data = metadata)

## Use vif() function in car package to calculate VIFs
vif(glmEnv)
sqrt(vif(glmEnv)) > 2

## Collinearity is indeed a problem
## Remove variables in order of highest of VIF value, re-run model and 
## reassess VIF values. This is the strategy recommended by Zuur et al.
## (2009, 2010).
## Remove percentage of surface oxygen first
glmEnv <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
              + pH + Temperature + perc_sub_macrophyte + per_shading,
              family = poisson,
              data = metadata)

vif(glmEnv)
sqrt(vif(glmEnv)) > 2

## VIF values now much lower but collinearity still a problem
## CPUE now has highest VIF value but is most likely to influence
## eDNA detection, thus we will remove percentage of shading which
## has the next highest VIF value.
glmEnv <- glm(Copy_Number ~ CPUE + Conductivity + Alkalinity
              + pH + Temperature + perc_sub_macrophyte,
              family = poisson,
              data = metadata)

vif(glmEnv)
sqrt(vif(glmEnv)) > 2

## Now, remove alkalinity as it is the least informative variable for 
## eDNA detection probability.
glmEnv <- glm(Copy_Number ~ CPUE + Conductivity + pH + Temperature 
              + perc_sub_macrophyte,
              family = poisson,
              data = metadata)

vif(glmEnv)
sqrt(vif(glmEnv)) > 2

## Now only CPUE and temperature are displaying evidence of collinearity
## (VIF>3).

## Remove water temperature as CPUE is more biologically intuitive, i.e.
## more fish will result in higher detection probability. However, it
## is likely that CPUE and temperature have a combined effect on eDNA
## detectability
glmEnv <- glm(Copy_Number ~ CPUE + Conductivity + pH + perc_sub_macrophyte,
              family = poisson,
              data = metadata)

vif(glmEnv)
sqrt(vif(glmEnv)) > 2


## Now, import data in format required by eDNAoccupancy package:
## qPCR results
crucianDetectionData <- read.csv("../Pond data/crucianDetectionData.csv", header=TRUE)

## Data checks
summary(crucianDetectionData)
head(crucianDetectionData)
names(crucianDetectionData)
str(crucianDetectionData)

## Environmental metadata
crucianSurveyData <- read.csv("../Pond data/crucianSurveyData.csv", header=TRUE)

## Data checks
summary(crucianSurveyData)
head(crucianSurveyData)
names(crucianSurveyData)
str(crucianSurveyData)

## Remove environmental variables previously identified as being 
## collinear - alkalinity, percentages of surface oxygen and shading
crucianSurveyData <- crucianSurveyData[,-c(4,6:7,9)]


## Load eDNAoccupancy package
#install.packages("eDNAoccupancy_0.2.0.tar.gz", repos=NULL, type="source")
library(eDNAoccupancy)

## Compute detection matrices y and K are computed from the detection data.
crucianDetections = occData(crucianDetectionData, 
                            siteColName = 'site',
                            sampleColName = 'sample')

## number of detections per sample
head(crucianDetections$y)

## number of qPCR replicates per sample
head(crucianDetections$K)

## Center and scale numeric-valued covariate measurements
crucianSurveyData.sc = scaleData(crucianSurveyData)

## Calculate standard deviation around mean of each variable
sd(crucianSurveyData.sc$CPUE, na.rm=TRUE)
sd(crucianSurveyData.sc$cond, na.rm=TRUE)
sd(crucianSurveyData.sc$pH, na.rm=TRUE)
sd(crucianSurveyData.sc$macrophyte, na.rm=TRUE)


## Fit a multi-scale occupancy model that evaluates eDNA availability
## as a condition of covariates the sample and qPCR replication levels 
## of the sampling hierarchy
## set.seed can be used to ensure reproducibility

## We expect macrophyte cover to influence whether crucian carp are 
## present in ponds, due to prey availability, cover from predation,
## O2 availability.

## At the sample level, we expect fish density to affect whether
## we detect eDNA, where more fish means more chance of detection.
## pH affects eDNA degradation, which is slower at higher pH thus
## we would expect more detection at higher pH.
## Temperature influences eDNA through the organism, causing their
## metabolic activity to increase and resulting in more DNA shed
## so we expect to detect eDNA more often at higher temperature.
## Conductivity gives an indication of dissolved salts and inorganic
## materials. It may tell us about inhibitory compounds so we would
## expect less eDNA detection as conductivity increases.
## Macrophyte cover may reduce UV exposure, and in heavily vegetated 
## ponds close to terrestrialisation create anoxic conditions, which 
## both slow eDNA degradation.

## At the qPCR replicate level, we only expect density and conductivity
## to influence eDNA detection. Fish density is ecpected to increase
## eDNA concentration and therefore chance of eDNA amplification.
## Conductivity may represent inhibitory substances so we would expect
## higher conductivity to result in failed qPCR reactions and no
## detection.

## Create list of variables on which model selection should be
## performed. Variables will be tried in various combinations
## at the sample and replicate level.
varlist <- c("CPUE","pH","cond","macrophyte","CPUE.rep","cond.rep")

## For ease of running a conditional loop, we need to keep
## thing simple by duplicating the CPUE and conductivity 
## columns with new names.
crucianSurveyData.sc$CPUE.rep <- crucianSurveyData.sc$CPUE
crucianSurveyData.sc$cond.rep <- crucianSurveyData.sc$cond

## seed a list of model combinations
combos <- NULL

## Create loop to create all combinations of all lengths
for (n in 0:6){
  test <- combn(varlist, n, simplify = F)
  combos <- append(combos, test)
}

## There are 64 possible model combinations - now write a loop 
## over them. We'll need two short vectors to describe which 
## variables are sample level and which are qPCR replicate level
samp.level <- c("CPUE","pH","cond","macrophyte")
rep.level <- c("CPUE.rep","cond.rep")

## Seed the output data frame (run every time)
model.summary <- data.frame(model.number = numeric(),  
                            model.vars = character(),
                            PPLC = numeric(),
                            WAIC = numeric())

## Now run for loop to perform all model combinations
for (n in 1:length(combos)){
  print(paste0("Starting model number ",n))
  model.number <- n  # for output
  vars <- combos[[n]]    # pick out the vars to be used
  model.vars <- paste(vars, collapse = ', ') # for output
  samp.vars.vec <- vars[which(vars %in% samp.level)]  # pick out which are sample-level...
  rep.vars.vec <- vars[which(vars %in% rep.level)]    # and which are replicate-level
  
  # prepare lists of variables - or replacements if no variables of that type
  if (length(samp.vars.vec) > 0){  # if using any sample variable, write it out...
    samp.vars <- paste("~ ", paste(samp.vars.vec, collapse = ' + '))
  } else {                         # otherwise write "~ 1"
    samp.vars <- "~ 1"
  }
  
  if (length(rep.vars.vec) > 0){   # as above for replicate variables
    rep.vars <- paste("~ ", paste(rep.vars.vec, collapse = ' + '))
  } else {
    rep.vars <- "~ 1"
  }
  
  # construct the model
  testm = occModel(formulaSite = ~ 1,
                   formulaSiteAndSample = samp.vars,
                   formulaReplicate = rep.vars,
                   siteData = crucianSurveyData.sc,
                   detectionMats=crucianDetections,
                   niter=11000,
                   niterInterval=2000,
                   siteColName = 'site')
  
  
  ## Evaluate fit of sample model
  PPLC <- posteriorPredictiveLoss(testm, burnin=1000)[[1]]
  WAIC <- WAIC(testm, burnin=1000)[[1]]
  
  # prepare output    
  out <- cbind(model.number, model.vars, PPLC, WAIC)
  model.summary <- rbind(model.summary,out)
  
}

## View data frame
summary(model.summary) 

## Need to turn the two output columns back into numeric variables
model.summary$PPLC <- as.numeric(as.character(model.summary$PPLC))
model.summary$WAIC <- as.numeric(as.character(model.summary$WAIC))

## Save data frame as a csv file 
write.csv(model.summary, "OccupancyModelTestSummary.csv")


## After inspecting the results, it would seem that the occupancy
## model that best fits the data includes macrophyte cover as a 
## covariate of site occurrence, no covariates of sample occurrence,
## and CPUE and conductivity as covariates of detection probability
## in qPCR replicates.


## Fit best model:
## set.seed can be used to ensure reproducibility
bestm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ 1,
                 formulaReplicate = ~ CPUE + cond,
                 siteData = crucianSurveyData.sc,
                 detectionMats=crucianDetections,
                 niter=11000,
                 niterInterval=2000,
                 siteColName = 'site')

posteriorSummary(bestm, burnin=1000, mcError=TRUE)

## Assess whether the Markov chain used to compute these estimates appears
## to have converged using trace plots of the parameters
plotTrace(bestm, c('beta..Intercept.', 'alpha..Intercept.', 'delta..Intercept.'), 
          burnin=1000)

## Autocorrelation plots of the parameters
plotACF(bestm, c('beta..Intercept.', 'alpha..Intercept.', 'delta..Intercept.'), 
        burnin=1000)

## Estimate posterior summaries of derived parameters. The probability of 
## eDNA occurrence in ponds was assumed to vary with macrophyte cover (psi), 
## and the conditional probability of eDNA detection was assumed to 
## be constant in samples (theta) but vary as a function of CPUE and 
## conductivity in qPCR replicates (p). 
## The posterior medians of these derived parameters are estimated as follows:
psi = posteriorSummaryOfSiteOccupancy(bestm, burnin=1000)
theta = posteriorSummaryOfSampleOccupancy(bestm, burnin=1000)
p = posteriorSummaryOfDetection(bestm, burnin=1000)
cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

## Evaluate fit of best model
posteriorPredictiveLoss(bestm, burnin=1000)
WAIC(bestm, burnin=1000)
posteriorSummaryOfAUC(bestm, burnin=1000)


## Compare the 'best' model to the null model for the entire sampling
## hierarchy:
nullm = occModel(formulaSite = ~ 1,
                 formulaSiteAndSample = ~ 1,
                 formulaReplicate = ~ 1,
                 siteData = crucianSurveyData.sc,
                 detectionMats=crucianDetections,
                 niter=11000,
                 niterInterval=2000,
                 siteColName = 'site')
   
posteriorSummary(nullm, burnin=1000, mcError=TRUE)

plotTrace(nullm, c('beta..Intercept.', 'alpha..Intercept.', 'delta..Intercept.'), 
          burnin=1000)

plotACF(nullm, c('beta..Intercept.', 'alpha..Intercept.', 'delta..Intercept.'), 
        burnin=1000)

## Estimate posterior summaries of derived parameters. The probability of 
## eDNA occurrence in ponds was assumed to be constant (psi), and the 
## conditional probability of eDNA detection was assumed to be constant
## in samples (theta) and qPCR replicates (p). 
## The posterior medians of these derived parameters are estimated as follows:
psi = posteriorSummaryOfSiteOccupancy(nullm, burnin=1000)
theta = posteriorSummaryOfSampleOccupancy(nullm, burnin=1000)
p = posteriorSummaryOfDetection(nullm, burnin=1000)
cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

## Evaluate fit of null model
posteriorPredictiveLoss(nullm, burnin=1000)
WAIC(nullm, burnin=1000)
posteriorSummaryOfAUC(nullm, burnin=1000)


## Set seed to ensure reproducibility of results for final model
set.seed(3726)

## Final model:
finalm = occModel(formulaSite = ~ 1,
                  formulaSiteAndSample = ~ 1,
                  formulaReplicate = ~ CPUE + cond,
                  siteData = crucianSurveyData.sc,
                  detectionMats=crucianDetections,
                  niter=11000,
                  niterInterval=2000,
                  siteColName = 'site')

posteriorSummary(finalm, burnin=1000, mcError=TRUE)

psi = posteriorSummaryOfSiteOccupancy(finalm, burnin=1000)
theta = posteriorSummaryOfSampleOccupancy(finalm, burnin=1000)
p = posteriorSummaryOfDetection(finalm, burnin=1000)
cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1])

## Round results to 2 decimal places
round(cbind(psi=psi$median, theta=theta$median[,1], p=p$median[,1]), 2)

## Obtain 95% CIs for each level of model
CI <- round(cbind(psi$lower, psi$upper, 
                  theta$lower[,1], theta$upper[,1],
                  p$lower[,1], p$upper[,1]), 2)
CI


## Plot effects of covariate on eDNA detection probability based on
## best model.
## Make new data frame for plotting:
plotting <- data.frame(crucianSurveyData[,c(2:3)],
                       p$median[,1], p$lower[,1], p$upper[,1])

## Change column names
names(plotting) <- c("CPUE","cond",
                     "p_median","p_lower","p_upper")

## Plot detection probability in qPCR replicates:
p5a <- ggplot(plotting, aes(x=CPUE, y=p_median)) + ggtitle("(a)")
p5a <- p5a + geom_point(pch=16, cex=2, colour="black")
p5a <- p5a + geom_errorbar(aes(ymin=p_lower,ymax=p_upper), width=3)
p5a <- p5a + coord_cartesian(xlim=c(0,140), ylim=c(0,1))
p5a <- p5a + scale_x_continuous(breaks=seq(0,140,20))
p5a <- p5a + scale_y_continuous(breaks=seq(0,1,0.1))
p5a <- p5a + labs(x="CPUE estimate", 
                  y="Probability of eDNA detection")
p5a <- p5a + theme(panel.background = element_rect(fill = 'white'),
                 plot.title = element_text(face="bold", hjust = 0),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 legend.title = element_blank(),
                 text = element_text(size=20),
                 legend.position="none")
p5a

p5b <- ggplot(plotting, aes(x=cond, y=p_median)) + ggtitle("(b)")
p5b <- p5b + geom_point(pch=16, cex=2, colour="black")
p5b <- p5b + geom_errorbar(aes(ymin=p_lower,ymax=p_upper), width=3)
p5b <- p5b + coord_cartesian(xlim=c(0,800), ylim=c(0,1))
p5b <- p5b + scale_x_continuous(breaks=seq(0,800,100))
p5b <- p5b + scale_y_continuous(breaks=seq(0,1,0.1))
p5b <- p5b + labs(x=expression(paste("Conductivity (",mu,"s/cm)")),
                  y="")
p5b <- p5b + theme(panel.background = element_rect(fill = 'white'),
                 plot.title = element_text(face="bold", hjust = 0),
                 axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                 axis.title.x = element_text(margin = unit(c(8, 0, 0, 0), "mm")),
                 axis.title.y = element_text(margin = unit(c(0, 5, 0, 0), "mm")),
                 axis.text.x = element_text(colour="black"),
                 axis.text.y = element_text(colour="black"),
                 legend.title = element_blank(),
                 text = element_text(size=20),
                 legend.position="none")
p5b

## Arrange all plots
## Load 'ggpubr' package
library(ggpubr)
ggarrange(p5a,p5b, ncol=2, nrow=1) 

