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



## Generate Ainverse matrices --------------------------------------------------
inb.founder <- 1-0.5^1
ainv.ls <- lapply(ped.ls, function(x){x$self <- 0
                                      x.ped <- checkPed(x, self = "self", f = inb.founder, verbose = TRUE)
                                      a.ainv <- asreml::ainverse(x.ped, fgen = list("self", inb.founder))
                                      p.nrm <- pedicure::nrm(x.ped, self = "self", coi = NULL, f = inb.founder)
                                      yy <- list("ped"=x.ped, "a.ainv" = a.ainv, "p.nrm" = p.nrm)
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

names(ainv.ls) <- c(paste0("A", 2:7), "AF", "AI")

unlist(lapply(ainv.ls, function(x) mean(attr(x$a.ainv,'inbreeding')[gkeep.name]))) #coefficient of inbreeding not returned for AIsweep
unlist(lapply(ainv.ls, function(x) mean(attr(x$p.nrm,'F')[gkeep.name])))

## Subset phenotypic data for traits of interest -------------------------------
unique(d0$Trait)
yldtraits <- levels(d0$Trait)[grep("yldpp|avfrtwt", levels(d0$Trait))]
trait2keep <- yldtraits[grep("_w|cavfrtwt_endAug|cyldppNA_endSep", yldtraits)]
d1 <- data.frame(droplevels(subset(d0, Trait %in% trait2keep)))

names(d1)

## Output data for analyses ----------------------------------------------------

save(list = c("d1", "ped.ls","ainv.ls"), file = paste0("3-RData/", Sys.Date(), "-Data4Analysis.RData"))


#####################
### END OF SCRIPT ###
#####################