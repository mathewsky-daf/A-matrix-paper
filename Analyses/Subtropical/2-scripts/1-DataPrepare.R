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
require(dwrPlus)
require(AGHmatrix)
source("G:/Delivery/R&DDel/HortForestSc/Horticulture/STRWBERY/BREED/Katie/Strawberry Breeding Program/KMDR/ASBP_KMD_login.R")
source("../../Analyses/00-functions/Pedigree4-Ngenerations.R")
#source("../../Analyses/00-functions/mat2sparse.R")

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
inb.founder <- 1e-06

#See if you can put in Genotypes names instead of row and column
ainv.ls <- lapply(ped.ls, function(x){x$self <- 0
                                      x$coi <- inb.founder #need coi for new founders
                                      x$coi[!is.na(x$self)] <- NA
                                      x.ped <- checkPed(x, self = "self", coi="coi",f = inb.founder, verbose = TRUE)
                                      a.ainv <- asreml::ainverse(x.ped, fgen = list("self", inb.founder))
                                      p.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder, inverse = TRUE, zap=1e-10)
                                      agh.2 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 2)), tol = 1e-10)
                                      dgh.2 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 2, dominance = TRUE)), tol = 1e-10)
                                      agh.8 <- mat2sparse(solve(Amatrix(x.ped[, 1:3], ploidy = 8)), tol = 1e-10)
                                      dimnames(a.ainv)[[2]] <- dimnames(p.nrm)[[2]] <- dimnames(agh.2)[[2]] <- dimnames(agh.8)[[2]] <- c("row", "column", "value")
                                      yy <- list("ped"=x.ped, "a.ainv" = a.ainv, 
                                                 "p.nrm" = p.nrm, "agh.2" = agh.2,
                                                 "dgh.2" = dgh.2, "agh.8" = agh.8)
                                      return(yy)
                                      })

unlist(lapply(ainv.ls, function(x) mean(attr(x$a.ainv,'inbreeding')[gkeep.name]))) #coefficient of inbreeding not returned for AIsweep
unlist(lapply(ainv.ls, function(x) mean(attr(x$p.nrm,'F')[gkeep.name])))

## Compare the values from the different algorithms
lapply(ainv.ls, function(x) {c(dim(x$a.ainv), dim(x$p.nrm), dim(x$agh.2), dim(x$dgh.2), dim(x$agh.8))})
# Curiously the ag.2 returns 1 less row for A4 above - investigating this with BC and DB
#KLM2025-05-12: we thought it as about the tolerance level for mat2sparse being different for nrm, so we've made them the same in the above
# but the problem remains- although it is now fixed for agh.8, except for A4...

#KLM2025-05-12: my current thought is to leave the calculation as is, we can only manipulate it by changing the tolerance and that might be dataset dependent... 

## Subset phenotypic data for traits of interest -------------------------------
unique(d0$Trait)
yldtraits <- levels(d0$Trait)[grep("yldpp|avfrtwt", levels(d0$Trait))]
trait2keep <- yldtraits[grep("_w|_endAug", yldtraits)]

#Remove columns we don't need
names2rm <- c("SiteName", "TrialTypeName", "TrialNumber", "Genotype.orig", "WorkingPlantdate", "Age", "Measurement_Age")


d1 <- data.frame(droplevels(subset(d0, Trait %in% trait2keep, select = names(d0[!names(d0) %in% names2rm]))))

d1$HarvestWeek <- factor(d1$HarvestWeek, levels = c(paste0("0", 4:9), 10:22, "endAug"))
levels(d1$HarvestWeek)
names(d1) <- gsub("Trait", "TH", names(d1))
d1$Trait <- factor(unlist(lapply(strsplit(as.character(d1$TH), "_"), function(x) x[1])))
d1 <- d1[order(d1$Trait, d1$Year, d1$HarvestWeek, d1$Plot),]
head(d1)

d1$TYH <- paste(d1$Trait, d1$Year, d1$HarvestWeek, sep = "_")
d1$TYH <- factor(d1$TYH, levels = unique(d1$TYH))
unique(d1$Trait)
unique(d1$TYH)
levels(d1$TYH)

d1$GDrop <- as.character(d1$GDrop)
d1$GDrop[d1$Genotype == "Phenomenal"] <- "Phenomenal"
d1$GKeep[d1$Genotype == "Phenomenal"] <- NA
d1$GDrop <- as.factor(d1$GDrop)
d1$GKeep <- factor(d1$GKeep)
d1$GDrop <- factor(d1$GDrop)

levels(d1$GDrop)

## Output data for analyses ----------------------------------------------------

save(list = c("d1", "ped.ls","ainv.ls"), file = paste0("3-RData/", Sys.Date(), "-Data4Analysis.RData"))


#####################
### END OF SCRIPT ###
#####################