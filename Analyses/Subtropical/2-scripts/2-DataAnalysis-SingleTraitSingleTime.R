#############################################
# Title: 2-DataAnalysis-SingleTraitSingleTime.R
# Date: 12 May 2025
# Author: Ky L Mathews
# Collaborator: Katie O'Connor
# Description: Single Trait Single Time analysis
#  With Total genetic and pedigree. 
# - Output data to .RData for analysis


## Setup ----
rm(list = ls())
require(tidyverse)
# require(pedicure)
require(asreml)
# require(dwreml)
# require(dwrPlus)
# require(AGHmatrix)
# source("G:/Delivery/R&DDel/HortForestSc/Horticulture/STRWBERY/BREED/Katie/Strawberry Breeding Program/KMDR/ASBP_KMD_login.R")
# source("../../Analyses/00-functions/Pedigree4-Ngenerations.R")

## Import data ----
load(paste0("3-RData/", Sys.Date(), "-Data4Analysis.RData"))

unique(d2$TYH)

## Data summary ----
d.sum <- d2 %>% group_by(TYH) %>%
  summarise(mean = mean(Value, na.rm = TRUE),
            nUniqueValue = length(unique(Value[is.na(Value)==FALSE])),
            nGeno = length(unique(Genotype[is.na(Value)==FALSE])),
            nGKeep = length(unique(GKeep)),
            nGKeepNA = length(unique(GKeep[is.na(Value)==TRUE])),
            nPlot = length(unique(Plot[is.na(Value)==FALSE])),
            pcNAplot = round(n_distinct(Plot[is.na(Value)==TRUE]) / n_distinct(Plot) * 100, 1),
#             # Genotypes where all plots have NA:
             pcNAgeno = round(sum(sapply(split(Value, factor(GKeep)), function(x) all(is.na(x)))) / n_distinct(GKeep) * 100, 1))
d.sum

t.lt50NA <- d.sum$TYH[d.sum$pcNAplot < 50]

d3 <- droplevels(subset(d2, TYH %in% t.lt50NA))

## Loop over ainv.ls, TYH level -----
asreml.options(ai.sing = TRUE)
out.ls <- list()

for(i in 1:length(ainv.ls)){# i=1 
  a.tmp <- ainv.ls[[i]]
  a.tmp <- a.tmp[names(a.tmp) %in% c("ped", "a.ainv", "p.nrm", "agh.2", "agh.8")]  #this drops out dgh.2 cos model is different
  
  type.ls <- list()
  for(j in 2:length(a.tmp)){# j=2 #
    ped.cy <- a.tmp[[1]]
    ainv.tmp <- a.tmp[[j]]
    
    results <- data.frame()
    for(k in 1:nlevels(d3$TYH)){ # k=1  
      tyh.tmp <- levels(d3$TYH)[k]
      d.tmp <- data.frame(droplevels(subset(d3, TYH == tyh.tmp)))
      s.abar <- mean(attr(a.tmp$a.ainv,'inbreeding')[levels(d.tmp$GKeep)])

      asr1 <- asreml(fixed = Value ~ GDrop, 
                     random =~ GKeep + Block,
                     residual =~ units,
                     data = d.tmp,
                     na.action = na.method(x = "include", y = "include"))
      asr1 <- update(asr1)
      
      asr2 <- asreml(fixed = Value ~ GDrop, 
                     random =~ vm(GKeep, ainv.tmp) + ide(GKeep, ainv.tmp) + Block,
                     residual =~ units,
                     data = d.tmp,
                     na.action = na.method(x = "include", y = "include"))
      
      asr2 <- update(asr2)
      
      # Extract required info
      vc1 <- cbind.data.frame(term = rownames(summary(asr1)$varcomp), summary(asr1)$varcomp)
      vc2 <- cbind.data.frame(term = rownames(summary(asr2)$varcomp), summary(asr2)$varcomp)
      
      #Check if any terms are singular
      vc1Singular <- FALSE 
      if(length(grep("S", vc1$bound))>0) vc1Singular <- TRUE
      
      vc2Singular <- FALSE 
      if(length(grep("S", vc2$bound))>0) vc2Singular <- TRUE
      
      coef1 <- cbind(asr1$coeff$random, asr1$vcoeff$random) %>%
        as.data.frame() %>%
        filter(str_detect(row.names(.), "GKeep")) %>%
        rename("effect" = 1,
               "vcoeffrandom" = 2) %>% 
        mutate(se = sqrt(vcoeffrandom),
               Genotype = factor(str_remove(row.names(.), "GKeep_")))
      
      # Calculate % replication
      tt <- table(table(d.tmp$GKeep[is.na(d.tmp$Value)==FALSE]))
      pRep <- (sum(tt)-sum(tt[names(tt) %in% c(0, 1)]))/sum(tt)*100
      nFamily <- n_distinct(ped.cy$Family[ped.cy$Genotype %in% unique(d.tmp$GKeep[is.na(d.tmp$Value)==FALSE & 
                                                                                    is.na(d.tmp$GKeep) == FALSE])])
      # Store calculations
      nBlock = nBlock = nlevels(d.tmp$Block)
      genvar = vc1$component[grep("GKeep", vc1$term)]
      blockvar1 = vc1$component[grep("Block", vc1$term, fixed = TRUE)]
      resvar1 = vc1$component[grep("units", vc1$term)]
      tgenvar = s.abar*vc2$component[grep("vm(GKeep", vc2$term, fixed = TRUE)] + vc2$component[grep("ide(GKeep", vc2$term, fixed = TRUE)]
      addvar0 = vc2$component[grep("vm(GKeep", vc2$term, fixed = TRUE)]
      addvar = s.abar*vc2$component[grep("vm(GKeep", vc2$term, fixed = TRUE)]
      nonaddvar = vc2$component[grep("ide(GKeep", vc2$term, fixed = TRUE)]
      blockvar2 = vc2$component[grep("Block", vc2$term, fixed = TRUE)]
      resvar2 = vc2$component[grep("units", vc2$term)]
      pcAdd = (addvar / (addvar + nonaddvar)) * 100
      nrep = d.tmp %>% filter(!is.na(Value)) %>% count(Genotype) %>% filter(n > 0) %>%  
        summarise(mean_freq = mean(n)) %>% pull(mean_freq)  
      h2 = addvar/(addvar + nonaddvar + blockvar2/nBlock + resvar2/nrep)      # calculated as per Oakey 2006, and Falconer and Mackay 1996
      
      #KM0250502: I would simplify this to be for using model 1, and also, you've already pulled out most of the required terms - just use them.
      H2 = genvar /sum(genvar, blockvar1/nBlock + resvar1/nrep)      # calculated as per Oakey 2006, and Falconer and Mackay 1996 #broad sense heritability
      r2 = mean(1 - coef1$se^2 / genvar)   #KM: this is sometimes negative - eeek, which it will be when GKeep < ResidualMS
      
      # Summarise results for this iteration
      iter_res <- data.frame(
        Depth = names(ainv.ls)[i],
        Type = names(a.tmp)[j],
        TYH = tyh.tmp,
        pcRep = pRep,
        nFamily = nFamily,
        logl1 = asr1$loglik,
        singular1 = vc1Singular,
        logl2 = asr2$loglik,
        singular2 = vc2Singular,
        genvar = genvar,
        resvar1 = resvar1,
        tgenvar = tgenvar,
        addvar0 = addvar0,
        addvar = addvar, 
        nonaddvar = nonaddvar,
        resvar2 = resvar2,
        pcAdd = pcAdd,
        h2 = h2,
        H2 = H2,
        r2 = r2
      )
      
      # Combine results
      results <- rbind(results, iter_res)
      
    }
    type.ls[[j]] <- results
  }
  out.ls[[i]] <- type.ls %>% bind_rows()
}

out.df <- out.ls %>% bind_rows()



# modesl from SimulatingMissingData.R



