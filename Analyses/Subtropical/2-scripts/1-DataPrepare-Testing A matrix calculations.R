#############################################
# Title: 1-DataPrepare.R
# Date: 18 March 2025
# Author: Ky L Mathews
# Collaborator: Katie O'Connor
# Description: Preparing data for A matrix analysis
    # - Read in phenotypic data
    # - Obtain vector of genotypes
    # - Extract pedigree from KMD
    # - Generate A matrices for 2:6 generations and ALL
    # - Extract traits of interest, yldpp, avfrtwt, cyldpp and cavgfrtwt
    # - Output data to .RData for analysis


## Setup ----
rm(list = ls())
require(tidyverse)
require(pedicure)
require(asreml)
require(dwreml)
require(AGHmatrix)
source("G:/Delivery/R&DDel/HortForestSc/Horticulture/STRWBERY/BREED/Katie/Strawberry Breeding Program/KMDR/ASBP_KMD_login.R")
source("../../Analyses/00-functions/Pedigree4-Ngenerations.R")

## Import data ----
load("G:/Delivery/R&DDel/HortForestSc/Horticulture/STRWBERY/BREED/BS22000/2023-24/1 Subtropical 2023-24/1 Maroochy/MRF Trials 2024/Analyses/Input/s2024-clean-spatial.RData")
d0 <- s2024_spatial

head(d0)

## Extract genotype names ----
gkeep.name <- levels(d0$GKeep)


## Pedigree from KMDR ----
kmd.login <- ASBP_KMD_login(katmandoo_url = "https://lavml093/strawberry") #mathewsk, see Edge for password
ped0 <- KMDR::export_pedigree_info(url = kmd.login$url,
                                   token = kmd.login$token,
                                   genotypes = gkeep.name,
                                   columns=list("GenotypeName","GenotypeAliasName","FemaleParentName",
                                                "MaleParentName","FamilyGroup"),
                                   purdy_level = 0,
                                   generations = 100, # will give x generations back
                                   include_founder_order=TRUE,
                                   set_blanks_to_NA = TRUE,
                                   check_alias = TRUE, # TRUE means will check the Alias column as well
                                   use_captions = FALSE)

# ped0 <- read.xlsx("G:/Delivery/R&DDel/HortForestSc/Horticulture/STRWBERY/BREED/Pedigree/Full pedigree list.xlsx",
#                   sheet = "Pedigree")
names(ped0)
names(ped0) <- gsub("GenotypeName", "Genotype", names(ped0))

#check that all genotypes are in the pedigree file
unique(gkeep.name %in% ped0$Genotype) # All good - yay!!

# Checking for selfs
ped1 <- unique(ped0[ped0$Genotype!=0, c("Genotype", "FemaleParentName", "MaleParentName")])
nrow(subset(ped1, FemaleParentName==MaleParentName & FemaleParentName !=0))  #0 # no selfed lines


## Obtain pedigree files for different levels
#KLM: Note we could use KMDR::export_pedigree_info and change the number of generations
## generations in KMD function = ngen-1 in ped.Ngen().
## I've used generation  = 100 in KMD function to get the full pedigree, and then subset it here
## to know for sure that all pedigree files are from the same source of information
## Discussed with KO and agreed that extracting once from KMD is the way forward.

ped.ls <- list()
for(i in 1:6) ped.ls[[i]] <-  ped.Ngen(ped.full=ped1, gen2select=gkeep.name, me = "GenotypeName", mum = "FemaleParentName",
                                       dad = "MaleParentName", ngen = i)
ped.ls[[7]] <- ped1
names(ped.ls) <- c(paste0("A", 2:7), "AF")
unlist(lapply(ped.ls, nrow))

ped2 <- checkPed(ped1)
head(ped2)
tail(ped2)
## Generate Ainverse matrices --------------------------------------------------
inb.founder <- 1-0.5^1

x <- ped.ls[[3]]

#See if you can put in Genotypes names instead of row and column

ainv.ls <- lapply(ped.ls, function(x){x$self <- 0
                                      # x$coi <- inb.founder #need coi for new founders
                                      # x$coi[!is.na(ped$self)] <- NA
                                      x$coi <- as.numeric(NA)
                                      x.ped <- checkPed(x, self = "self", coi="coi",f = inb.founder, verbose = TRUE)
                                      a.ainv <- asreml::ainverse(x.ped, fgen = list("self", inb.founder))
                                      p.nrm <- pedicure::nrm(x.ped, self = "self", coi = "coi", f = inb.founder, inverse = TRUE)
                                      agh.2 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 2)))
                                      agh.8 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 8)))
                                      dimnames(a.ainv)[[2]] <- dimnames(p.nrm)[[2]] <- dimnames(agh.2)[[2]] <- dimnames(agh.8)[[2]] <- c("row", "column", "value")
                                      yy <- list("ped"=x.ped, "a.ainv" = a.ainv, 
                                                 "p.nrm" = p.nrm, "agh.2" = agh.2,
                                                 "agh.8" = agh.8)
                                      return(yy)
                                      })
## Tests to send to BC
# a.ainv <- ainv.ls[[4]]$a.ainv
# p.nrm <- ainv.ls[[4]]$p.nrm
# mean(attr(a.ainv, 'inbreeding')[gkeep.name])
# mean(attr(p.nrm, "F")[gkeep.name])
# 
# head(a.ainv)
# head(p.nrm)
# tail(a.ainv)
# tail(p.nrm)

#Now add in ainv for AIsweep
ainv.ls[[8]] <- ainv.ls[[7]]
ainv.ls[[8]]$a.ainv <- AIsweep(ainv.ls[[8]]$a.ainv, keep = gkeep.name)
ainv.ls[[8]]$p.nrm <- AIsweep(ainv.ls[[8]]$p.nrm, keep = gkeep.name)
ainv.ls[[8]]$agh.2 <- AIsweep(ainv.ls[[8]]$agh.2, keep = gkeep.name)
ainv.ls[[8]]$agh.8 <- AIsweep(ainv.ls[[8]]$agh.8, keep = gkeep.name)

# KO getting this error 250319 for ainv.ls[[8]]$p.nrm <- AIsweep(ainv.ls[[8]]$p.nrm, keep = gkeep.name)
# Error in y[(x[, 2] - 1) * nrow + x[, 1]] <- x[, 3] : 
#   only 0's may be mixed with negative subscripts

names(ainv.ls) <- c(paste0("A", 2:7), "AF", "AI")

unlist(lapply(ainv.ls, function(x) mean(attr(x$a.ainv,'inbreeding')[gkeep.name]))) #coefficient of inbreeding not returned for AIsweep
unlist(lapply(ainv.ls, function(x) mean(attr(x$p.nrm,'F')[gkeep.name])))

## Compare the values from the different algorithms
lapply(ainv.ls, function(x) {c(dim(x$a.ainv), dim(x$p.nrm), dim(x$agh.2), dim(x$agh.8))})
# Curiously the ag.2 returns 1 less row for A4 does it for agh.8 too  A5 and above
x <- ped.ls[[4]] #names(ped.ls)
x$self <- 0
ped <- x
save("ped", file = "4bc.RData")
x$coi <- inb.founder
x.ped <- checkPed(x, coi="coi",f = inb.founder, verbose = TRUE)
x.ped$self <- as.numeric(x.ped$self)
a.ainv <- asreml::ainverse(x.ped, fgen = list("self", inb.founder))
dimnames(a.ainv)[[2]] <- c("row", "column", "value")
a.gen <- attr(a.ainv,"rowNames")
length(a.gen)
head(a.ainv)
a.ainv <- as.data.frame(a.ainv)
a.ainv$row <- a.gen[a.ainv$row]
a.ainv$column <- a.gen[a.ainv$column]
dim(a.ainv)

p.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder, inverse = TRUE, zap=1e-10)
dimnames(p.nrm)[[2]] <- c("row", "column", "value")
p.gen <- attr(p.nrm, "rowNames")
p.nrm <- as.data.frame(p.nrm)
p.nrm$row <- p.gen[p.nrm$row]
p.nrm$column <- p.gen[p.nrm$column]
head(p.nrm)
dim(p.nrm)

agh.2 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 2)), tol = 1e-14)  # this makes it to dim(agh.2) = 2682
dimnames(agh.2)[[2]] <- c("row", "column", "value")
a2.gen <- attr(agh.2, "rowNames")
agh.2 <- as.data.frame(agh.2)
agh.2$row <- a2.gen[agh.2$row]
agh.2$column <- a2.gen[agh.2$column]
head(agh.2)
dim(agh.2)

agh.8 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 8)), tol = 1e-10)
dimnames(agh.8)[[2]] <- c("row", "column", "value")
a8.gen <- attr(agh.8, "rowNames")
agh.8 <- as.data.frame(agh.8)
# agh.8$row[a8.gen=="2015-068"]
# agh.8$column[a8.gen=="Red_Rhapsody"]
agh.8$row <- a8.gen[agh.8$row]
agh.8$column <- a8.gen[agh.8$column]
head(agh.8)
dim(agh.8)

x1 <- merge(as.data.frame(a.ainv), as.data.frame(p.nrm), by = c("row", "column"), all.x = TRUE)
x1 <- merge(as.data.frame(x1), as.data.frame(agh.2), by = c("row", "column"), all.x = TRUE)
x1 <- merge(as.data.frame(x1), as.data.frame(agh.8), by = c("row", "column"), all.x = TRUE)
x1 <- x1[order(x1$row, x1$column),]
names(x1) <- c("row", "column", "a.ainv", "p.nrm", "agh.2", "agh.8")

subset(x1, is.na(agh.2)==TRUE)
subset(x1, is.na(agh.8)==TRUE)
head(x1)
subset(x1, row=="2015-068" & column=="Red_Rhapsody")

# row       column a.ainv p.nrm agh.2        agh.8
# # 283 2015-068 Red_Rhapsody      0     0    NA 0.0001900691
# row       column a.ainv p.nrm agh.2        agh.8
# 288 2015-068 Red_Rhapsody      0     0    NA 0.0002005186

#KLM2025-05-12: shows that agh.2 is now 0 not NA.
#      row       column a.ainv p.nrm        agh.2        agh.8
# 2015-068 Red_Rhapsody      0     0 1.191516e-15 0.0001900691


subset(ped2, Genotype %in% c("2015-068", "Red_Rhapsody"))
# Genotype      FemaleParentName   MaleParentName
# Red_Rhapsody         2005-063 Suncoast_Delight
# 2015-068         Red_Rhapsody          Tristan

subset(ped2, FemaleParentName %in% "Red_Rhapsody")
subset(ped2, FemaleParentName %in% "2015-068" |MaleParentName %in% "2015-068")
subset(ped2, Genotype %in% "2015-216")

#HHMMM.. what are A values...

A.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder, inverse = FALSE)
A.nrm[c("2015-068", "Red_Rhapsody", "2015-216"), c("2015-068", "Red_Rhapsody", "2015-216")]
A.nrm[A.nrm==0.5]
dim(A.nrm)
Ai.nrmsp <- mat2sparse(as.matrix(solve(A.nrm)))
Ai.nrmtrip <- mat2triplet(as.matrix(solve(A.nrm)))
dim(Ai.nrmsp)
length(Ai.nrmtrip$i)
length(unique(Ai.nrmtrip$i))
head(cbind(Ai.nrmtrip$i, Ai.nrmtrip$j, Ai.nrmtrip$x))
dim(as.matrix(solve(A.nrm)))
kk <- as.matrix(solve(A.nrm))
dim(kk)
kk <- lower.tri(kk,diag=TRUE)
table(kk)
head(Ai.nrmsp)
kk[kk==TRUE]# <- 1
kk[lower.tri(kk)==FALSE] <- 0
kk[1:5, 1:5]
sum(kk)
head(Ai.nrmsp)

k1 <- 0
for(i in 1:845) k1 <- sum(i, k1)

dim(Ai.nrmsp)
Ai.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder, inverse = TRUE)
dim(Ai.nrm)


A.2 <- Amatrix(x.ped[, 1:3], ploidy = 2)
Ai.2 <- solve(A.2)
Ai.2sp <- mat2sparse(Ai.2)

dim(A.2)
dim(Ai.2)
dim(Ai.2sp)

#c.f A.nrm
A.nrm[c("2015-068", "Red_Rhapsody", "2015-216"), c("2015-068", "Red_Rhapsody", "2015-216")]
A.2[c("2015-068", "Red_Rhapsody", "2015-216"), c("2015-068", "Red_Rhapsody", "2015-216")]

A.nrm[1:10, 1:10]
A.2[1:10, 1:10]

A.nrmsp <- mat2sparse(as.matrix(A.nrm))
A.2sp <- mat2sparse(A.2)

head(A.nrmsp)
head(A.2sp)
dim(A.nrmsp)
dim(A.2sp)

a.cf <- cbind.data.frame(A.nrmsp, A.2sp)
head(a.cf)
names(a.cv) <- c("row.n","col.n", "val.n", "row.2", "col.2", "val.2")
subset(a.cv, row.n!=row.2) #0
subset(a.cv, col.n!=col.2) #0
subset(a.cv, val.n!=val.2) #0, so everything matches before solve...

Ai.nrm <- solve(as.matrix(A.nrm))
Ai.2 <- solve(A.2)
dim(Ai.nrm)
dim(Ai.2)
table(Ai.nrm==Ai.2) #Ai is the same for both

Ai.nrmsp <- mat2sparse(Ai.nrm)
Ai.2sp <- mat2sparse(Ai.2)
dim(Ai.nrmsp)
dim(Ai.2sp)

ai.cf <- cbind.data.frame(Ai.nrmsp, Ai.2sp)
head(ai.cf)
names(ai.cf) <- c("row.n","col.n", "val.n", "row.2", "col.2", "val.2")
subset(ai.cf, row.n!=row.2) #0
subset(ai.cf, col.n!=col.2) #0
subset(ai.cf, val.n!=val.2) #0, so everything matches before solve...
subset(ai.cf, row.n==170)

Aip.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder, inverse = TRUE)
dim(Aip.nrm)

ai.cf <- merge(Ai.nrmsp, Aip.nrm, by = c("row", "col"), all.y = TRUE)
names(ai.cf)
subset(ai.cf, val != value)
subset(ai.cf, is.na(val)==TRUE)
# row col val value
# 203  76  37  NA     0
Ai.nrm[76, 37] #0
Aip.nrm.df <- data.frame(Aip.nrm)
subset(Aip.nrm.df, value==0)



A.8 <- Amatrix(x.ped[, 1:3], ploidy = 8)
A.8[c("2015-068", "Red_Rhapsody", "2015-216"), c("2015-068", "Red_Rhapsody", "2015-216")]



lapply(ainv.ls, function(x) {lapply(x, function(x1) names(x1))})
lapply(ainv.ls, function(x) {lapply(x, function(x1) is.matrix(x1))})

x <- ainv.ls[[1]]
cf0 <- lapply(ainv.ls, function(x) {x1 <- merge(as.data.frame(x$a.ainv), as.data.frame(x$p.nrm), by = c("row", "column"), all.x = TRUE)
                                    x1 <- merge(as.data.frame(x1), as.data.frame(x$agh.2), by = c("row", "column"), all.x = TRUE)
                                    x1 <- merge(as.data.frame(x1), as.data.frame(x$agh.8), by = c("row", "column"), all.x = TRUE)
                                    x1 <- x1[order(x1$row, x1$column),]
                                    names(x1) <- c("row", "column", "a.ainv", "p.nrm", "agh.2", "agh.8")
                                    x1})
cf0[[1]]
str(cf0[[1]])
cf0.cor <- lapply(cf0, function(x) cor(x[,3:6], use = "pairwise"))

cf0.diff <- lapply(cf0, function(x) {x$ai.pn <- x$a.ainv - x$p.nrm
                                     x$ai.a2 <- x$a.ainv - x$agh.2
                                     x$ai.a8 <- x$a.ainv - x$agh.8
                                     x$pn.a2 <- x$p.nrm - x$agh.2
                                     x$pn.a8 <- x$p.nrm - x$agh.8
                                     x$a2.a8 <- x$agh.2 - x$agh.8
                                    return(x)})
lapply(cf0.diff, function(x) summary(x$ai.pn))
lapply(cf0.diff, function(x) summary(x$ai.a2))
lapply(cf0.diff, function(x) summary(x$ai.a8))
lapply(cf0.diff, function(x) summary(x$pn.a2))
lapply(cf0.diff, function(x) summary(x$pn.a8))
lapply(cf0.diff, function(x) summary(x$a2.a8))

#Check where agh.2 has one less record than other calculations
lapply(cf0.diff, function(x) x[is.na(x$agh.2)==TRUE,])
x <- cf0.diff[[3]]
head(x)
subset(x, is.na(pn.a2)==TRUE)


#p.rnm = agh.2 always.
#A.inv is always different
#agh.8 is different from agh.2 for all generation numbers

## Subset phenotypic data for traits of interest -------------------------------
unique(d0$Trait)
yldtraits <- levels(d0$Trait)[grep("yldpp|avfrtwt", levels(d0$Trait))]
trait2keep <- yldtraits[grep("_w|cavfrtwt_endAug", yldtraits)]
d1 <- data.frame(droplevels(subset(d0, Trait %in% trait2keep)))

names(d1)

## Output data for analyses ----------------------------------------------------

save(list = c("d1", "ped.ls","ainv.ls"), file = paste0("3-RData/", Sys.Date(), "-Data4Analysis.RData"))


#####################
### END OF SCRIPT ###
#####################