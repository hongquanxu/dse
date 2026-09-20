## Tables 3-6 and Tables S2-S3
##
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

library(ggplot2)
library(dplyr)
library(tidyr)
library(UniDOE)
library(MaxPro)
source("dse_functions.R")
source("generate_soa.R")
source("dse_tables.R")

## GSOAs used in Tables 4-6 and S2-S3
gsoa32x9 <- get.soa.new(N = 32, k = 9) %/% 2 # GSOA(32, 9, 2^2, 3)
gsoa64x16 <- get.soa.new(N = 64, k = 16) %/% 2 # GSOA(64, 16, 2^2, 3)
gsoa81x9 <- read.csv("soa81x9.csv")[,-1] # GSOA(81,9, 3^2, 3)


######################
### Table 3 (a)(b) ### 
######################
options(pillar.sigfig = 4) 
tab_3a <- generate_table(s=2,p=2) # Table 3(a)
tab_3a[,c(6,1:5)] %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

tab_3b <- generate_table(s=3,p=2) # Table 3(b)
tab_3b[,c(6,1:5)] %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

###############
### Table 4 ### 
###############
phi_bar <- function(x, s, p, crit = "CD2") {
  n <- nrow(x); m <- ncol(x)
  d <- 1/s^p
  
  if (crit == "CD2") {
    g0 <- (13/12)^m
    g1 <- -2
    g2 <- 1
  } else if (crit == "MaxPro") {
    g0 <- - 1/( (n-1)*(d)^(2*m) )
    g1 <- 0
    g2 <- n/(n-1)
  }
  
  c1 <- 0
  for (a in 0:(s^p-1)) {
    for (b in 0:(s^p-1)) {
      c1 <- c1 + F_FUNC(c(a,b), type = crit, s, p)
    }
  }
  c1 <- g2*(c1 / (s^p)^2 )^m
  
  c2 <- 0
  for (a in 0:(s^p-1)) {
    c2 <- c2 + G_FUNC(a, type = crit, s, p)
  }
  c2 <- g0 + g1*(c2 / s^p)^m
  
  final_res <- generate_table(s,p)
  if (crit == "CD2") {
    if (p == 2) {
      y <- as.numeric(final_res[1,4])
      z <- as.numeric(final_res[1,5])
    } else if (p == 1) {
      y <- as.numeric(final_res[1,3])
      z <- 1
    }
  } else if (crit == "MaxPro") {
    if (p == 2) {
      y <- as.numeric(final_res[5,4])
      z <- as.numeric(final_res[5,5])
    } else if (p == 1) {
      y <- as.numeric(final_res[5,3])
      z <- 1
    }
  }
  
  if (crit == "CD2") {
    return(c1*dse(x, s, p, y, z) + c2)
  } else if (crit == "MaxPro") 
    return((c1*dse(x, s, p, y, z) + c2)^(1/m))
}

design <- gsoa32x9 
k <- 6
comb <- combinat::combn(9,k)
table_4 <- data.frame(design = 1:4,
                     CD_sp = NA,
                     CD_all = NA,
                     MP = NA,
                     MP_all = NA)

index <- c(34,49,11,84)

for (i in 1:4) {
  j <- index[i]
  design_new <- design[,comb[,j]]
  Pij <- fastPij(design_new, 2, 2, imax=6)
  Aj <- fastSP2(design_new, 4, 1)
  table_4[i,2] <- phi_bar(design_new,s=2,p=2,crit = "CD2")
  table_4[i,3] <- phi_bar(design_new,s=4,p=1,crit = "CD2")
  table_4[i,4] <- phi_bar(design_new,s=2,p=2,crit = "MaxPro")
  table_4[i,5] <- phi_bar(design_new,s=4,p=1,crit = "MaxPro")
  cat("D", i, "  ", sep="")
  cat(comb[,j], sep=","); cat("  ");
  cat(Pij[3, 4:6], sep=","); cat("  ");
  cat(Aj[3:4], sep=","); cat("  ");
  cat(round(as.matrix(table_4[i,2:5]),5), "\n" )
}


#################################
## Table 5 and 6 and S2 and S3 ##
#################################
generate_comparison <- function(design, s, p, crit = "CD2") {
  final_res <- generate_table(s, p)
  
  y <- as.numeric(final_res[which(final_res$crit == crit),4])
  z <- as.numeric(final_res[which(final_res$crit == crit),5])
  df2 <- data.frame(n = nrow(design),
  					m = 3:ncol(design),
                    S2LP = NA,
                    ALP = NA,
                    Algm = NA)	# Algm is GenUD or MaxProQQ
  
  for (k in 3:ncol(design)) {
    comb <- combinat::combn(ncol(design),k)
    
    if (k != ncol(design)) {
      values <- parallel::mclapply(
        seq_len(ncol(comb)),
        function(i) {
          idx <- comb[, i]
          dse(design[, idx], s, p, y, z)
        },
        mc.cores=parallel::detectCores()
      )
      values <- unlist(values)
      min_index <- which.min(values)
      design_new <- design[, comb[, min_index]]
    } else design_new <- design
    
   cat("k=", k, ", ", sep="") # show progress
    
    res <- main_p2(design_new, s, crit = crit)
    df2$S2LP[k-2] <- min(res[,1]) #S2LP
    df2$ALP[k-2] <- min(res[,3]) #ALP
    df2$Algm[k-2] <- min(res[,2]) #ALG
  }
  
  cat("\n")
  
  df2
}


## Note: The results may vary slightly due to the randomness of the procedure 

## Table 5
crit = "CD2"
generate_comparison(gsoa32x9, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa64x16, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa81x9, s = 3, p = 2, crit = crit) 

## Table 6
crit = "MaxPro"
generate_comparison(gsoa32x9, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa64x16, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa81x9, s = 3, p = 2, crit = crit) 

## Table S2
crit = "WD2"
generate_comparison(gsoa32x9, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa64x16, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa81x9, s = 3, p = 2, crit = crit) 

## Table S3
crit = "MD2"
generate_comparison(gsoa32x9, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa64x16, s = 2, p = 2, crit = crit) 
generate_comparison(gsoa81x9, s = 3, p = 2, crit = crit) 


