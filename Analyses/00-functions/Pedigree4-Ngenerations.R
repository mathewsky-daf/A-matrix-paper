
#####################
# Title: Pedigree4-Ngenerations.R
# Author: Ky L Mathews
# Date: 18 March 2025, based on Pedigree4-Ngenerations.R written for powdery mildew 2024 analyses
# Description: A function to subset a pedigree file going back N generations.
#####################

# 
# ngen <- 6
# ped.g0 <- unique(g2[, c("Genotype", "FemaleParentName", "MaleParentName")]) #1753 unique entries
# # names(ped.g0) <- c("Me", "Mum", "Dad")
# names(ped.g0)[1] <- "GenotypeName"



ped.Ngen <- function(ped.full, gen2select, me = "GenotypeName", mum = "FemaleParentName",
                     dad = "MaleParentName", ngen = 2){
             ped <- ped.full
             names(ped)[1:3] <- c(me, mum, dad)
             ped.g0 <- subset(ped, GenotypeName %in% gen2select)
             ped.ls <- list()
             ped.ls[[1]] <- ped.g0
             for(i in 1:ngen){# i <- 1
                 tmp <- ped.ls[[i]]
                 gen.tmp <- unique(c(tmp[, mum], tmp[,dad]))
                 ped.ls[[i+1]] <- subset(ped, GenotypeName %in% gen.tmp, select = c(me, mum, dad))
             }
             ped0 <- unique(ped.ls %>% bind_rows())
             return(ped0)
}
