#############################################
# Title: 2-DataAnalysis-MantelTests.R
# Date: 16 June 2025
# Author: Katie O'Connor
# Collaborator: Ky L Mathews 
# Description: Mantel tests between A-matrices
# - Output data to .RData for analysis

require(vegan)

## Import data ----
load(paste0("3-RData/", "2025-05-23", "-Data4Analysis.RData"))

# Create an empty data frame to store the results
mantel_results <- data.frame(
  Level = character(),
  Matrix1 = character(),
  Matrix2 = character(),
  Mantel_r = numeric(),
  Significance = numeric(),
  stringsAsFactors = FALSE
)

# Loop through all pedigree depths
for (X in 2:7) {
  
  # Change from sparse to full and extract the matrices for the current level
  a.ainv.f <- sp2Matrix(ainv.ls[[paste0("A", X)]]$a.ainv) # .f = full matrix, not sparse
  p.nrm.f <- sp2Matrix(ainv.ls[[paste0("A", X)]]$p.nrm)
  agh.2.f <- sp2Matrix(ainv.ls[[paste0("A", X)]]$agh.2)
  agh.8.f <- sp2Matrix(ainv.ls[[paste0("A", X)]]$agh.8)
  
  # List of matrices to compare
  matrices <- list(
    a.ainv = a.ainv.f,
    p.nrm = p.nrm.f,
    agh.2 = agh.2.f,
    agh.8 = agh.8.f
  )
  
  # Perform Mantel tests for every combination of matrices within a pedigree depth
  matrix_names <- names(matrices)
  
  for (i in 1:(length(matrices) - 1)) {
    
    for (j in (i + 1):length(matrices)) {
      
      # Perform the Mantel test
      mantel_result <- mantel(matrices[[i]], matrices[[j]], method = "pearson", permutations = 999)
      
      # Extract the Mantel statistic and significance
      mantel_r <- mantel_result$statistic
      p_value <- mantel_result$signif
      
      # Add the results to the data frame
      mantel_results <- rbind(
        mantel_results,
        data.frame(
          Level = paste0("A", X),
          Matrix1 = matrix_names[i],
          Matrix2 = matrix_names[j],
          Mantel_r = mantel_r,
          Significance = p_value,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

mantel_results

# Which pair/s of matrices has/have the lowest correlation/s?
mantel_results %>%
  filter(Mantel_r == min(Mantel_r, na.rm = TRUE))

save(list = c("mantel_results"), file = paste0("3-RData/", today(), "-MantelTests.RData"))
