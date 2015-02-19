
library(nlme)                           # lme()
library(qvalue)                         # qvalue()

ctf8Dat <- read.csv(jPaste(whereAmI, "data/ctf8_dataforRick.csv"))
str(ctf8Dat)
## 'data.frame':	6144 obs. of  19 variables:
## $ OrfName: Factor w/ 4791 levels "YAL002W","YAL004W",..: 4791 3578 2859 3851 ..
##  $ Plate1 : num  0.904 1.103 1 1.234 1.149 ...
##  $ Plate2 : num  0.842 0.97 0.941 1.116 1.026 ...
##  $ Plate3 : num  0.967 1.15 1.001 1.261 1.07 ...
##  $ Plate4 : num  1.02 1.02 1.05 1.08 1.05 ...
##  $ Plate5 : num  0.982 0.972 0.993 0.956 0.904 ...
##  $ Plate6 : num  0.867 0.937 0.992 0.898 0.947 ...
##  $ Plate7 : num  0 1.19 1 1.18 1.1 ...
##  $ Plate8 : num  0 0.972 0.861 0.949 0.962 ...
##  $ Plate9 : num  0 1.02 0.813 1 0.896 ...
##  $ Plate10: num  0 0.978 1.006 1.03 1.003 ...
##  $ Plate11: num  0 1.06 1.03 1.1 1.01 ...
##  $ Plate12: num  0 1.005 0.997 1.026 1.034 ...
##  $ Plate13: num  1 1.06 0.865 1.099 0.915 ...
##  $ Plate14: num  0.877 0.922 0.845 1 0.946 ...
##  $ Plate15: num  1.051 1.128 0.906 1.241 0.936 ...
##  $ Plate16: num  0.851 1.016 0.935 0.906 1.001 ...
##  $ Plate17: num  0.831 1.006 0.929 0.91 0.966 ...
##  $ Plate18: num  0.746 0.828 0.817 0.741 0.918 ...

## renaming the ORF factor
names(ctf8Dat)[1] <- "ORF"              # I like this better
orfTab <- with(ctf8Dat, sort(table(ORF)))
(orfFreqTab <- table(orfTab))
##    1    2    3    4    5    6    8   63  245 
## 4353  144   51  198   12   29    2    1    1 

## write some tabulation results to file
foo <- as.data.frame(orfFreqTab)
orfSpotSummary <-
  data.frame(tReps = as.numeric(as.character(foo$orfTab)),
             nOrfs = foo$Freq)
orfSpotSummary$totSpots <- with(orfSpotSummary, tReps * nOrfs)
rownames(orfSpotSummary) <- NULL
orfSpotSummary

jWriteTable(rbind(orfSpotSummary, colSums(orfSpotSummary)),
            jPaste(whereAmI, "results/orfSporSummary-raw.txt"),
            sep = "\t") 

## info on which is plate is what w.r.t. biological replicates and
## treatment
pDat <- read.csv(jPaste(whereAmI, "data/ctf8_plate_info.csv"))
str(pDat)
## 'data.frame':	18 obs. of  3 variables:
##  $ PlateNo             : int  1 2 3 4 5 6 7 8 9 10 ...
##  $ Biological.Replicate: int  1 1 1 2 2 2 3 3 3 4 ...
##  $ Expt.Control        : Factor w/ 2 levels "C","E": 1 1 1 2 2 2 1 1 1 2 ...

## harmonize ctf8Dat and pDat in terms of "plate"
## in anticipation of using merge() below
names(ctf8Dat)
plateLevels <- jPaste("Plate", 1:nrow(pDat))
pDat$plate <- factor(jPaste("Plate", pDat$PlateNo),
                     levels = plateLevels)

names(pDat)
names(pDat) <- c("plateNo", "bRep", "tx", "plate")

## fit linear mixed model to data for each ORF
## this will take a couple minutes, ~4 mins
ctf8Res <-
  data.frame(ORF = levels(ctf8Dat$ORF),
             t(with(ctf8Dat,
                    sapply(seq_len(nlevels(ORF)), function(i) {
                      jOrf <- levels(ORF)[i]
                      jDat <- data.frame(y = unlist(subset(ctf8Dat, subset = ORF == jOrf,
                                           select = -ORF)),
                                         plate = factor(rep(plateLevels, each = orfTab[jOrf]),
                                           levels = plateLevels))
                      jDat <- merge(jDat, pDat, sort = FALSE)
                      if(sum(jDat$y > 0)) {
                        lmeFit <- lme(y ~ tx, jDat, random = ~1 | bRep)
                        lmeTable <- summary(lmeFit)$tTable
                        foo <- c(nSpots = nrow(jDat),
                                 nSpotsGood = sum(jDat$y > 0),
                                 propSpotsGood = mean(jDat$y > 0),
                                 txEffEst = lmeTable["txE", "Value"],
                                 txEffSe = lmeTable["txE", "Std.Error"],
                                 txEffPval = lmeTable["txE", "p-value"],
                                 bRepSd = sqrt(getVarCov(lmeFit)),
                                 sigma = lmeFit$sigma)
                      } else {
                        foo <- c(nSpots = nrow(jDat),
                                 nSpotsGood = sum(jDat$y > 0),
                                 propSpotsGood = mean(jDat$y > 0),
                                 txEffEst = NA, txEffSe = NA,
                                 txEffPval = NA, bRepSd = NA, sigma = NA)
                      }
                      return(foo)
                    }))))
               
str(ctf8Res)
## data.frame':	4791 obs. of  9 variables:
##  $ ORF          : Factor w/ 4791 levels "YAL002W","YAL004W",..: 1 2 3 4 5 6 7 ..
##  $ nSpots       : num  18 18 18 18 18 18 18 18 18 18 ...
##  $ nSpotsGood   : num  18 18 18 18 18 18 18 18 18 18 ...
##  $ propSpotsGood: num  1 1 1 1 1 1 1 1 1 1 ...
##  $ txEffEst     : num  -0.0199 0.0446 -0.1477 -0.0701 0.0458 ...
##  $ txEffSe      : num  0.0724 0.0334 0.0887 0.0338 0.0337 ...
##  $ txEffPval    : num  0.797 0.253 0.171 0.106 0.245 ...
##  $ bRepSd       : num  0.0788 0.0315 0.1052 0.0381 0.0376 ...
##  $ sigma        : num  0.0707 0.0453 0.047 0.0279 0.0296 ...

summary(ctf8Res$txEffPval)
##    Min.  1st Qu.   Median     Mean  3rd Qu.     Max.     NA's 
##  0.0000   0.1629   0.3825   0.4267   0.6753   0.9997 167.0000 

## get q-values
lmeOK <- !is.na(ctf8Res$txEffPval)
table(lmeOK)
## FALSE  TRUE 
##   167  4624 
qValRes <- qvalue(ctf8Res$txEffPval[lmeOK])
qValRes$pi0                             # 0.752723
1 - qValRes$pi0                         # 0.2472769
ctf8Res$txEffQval <- NA
ctf8Res$txEffQval[lmeOK] <- qValRes$qvalues

# write results to a file
jWriteTable(ctf8Res,
            jPaste(whereAmI, "results/ctf8-lmeResults.txt"))
