## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

## key functions used by other files
##      main_p2
##		generate_table
### based on Theorem 6: Linking the average criterion value with the double stratification enumerator

## 9/18/26


library(LHD)
library(tidyverse)
library(MaxPro)
library(UniDOE) # Archive at https://cran.r-project.org/src/contrib/Archive/UniDOE/


### Define f(x,y)
nrt_distance <- function(x, y, s = 2, p = 3) {
  i <- 1
  match <- FALSE
  distance <- 0
  while (i <= p & match == FALSE) {
    fx <- floor(x / s^(p-i)) %% s
    fy <- floor(y / s^(p-i)) %% s
    if (fx != fy) {
      distance <- p + 1 - i
      match <- TRUE
    }
    i <- i + 1
  }
  
  return(distance)
}

delta <- function(t, z, u, s) {
  sum(floor(s^u * t) == floor(s^u * z))
}

### Following functions give the f and g given different phi
cd_f <- function(x1, s = 2, p = 2) 
{
  s1 <- s^p
  x1 <- (2*x1 + 1 - s1)/(2*s1)
  x <- x1[1]
  y <- x1[2]
  1 + (abs(x) + abs(y) - abs(x-y))/2
}

cd_g <- function(x1, s = 2, p = 2) 
{
  x1 <- (2*x1 + 1 - s^p)/(2*s^p)
  1 + abs(x1)/2 - x1^2/2
}

wd_f <- function(x1, s = 2, p = 2) 
{
  s1 <- s^p
  x1 <- (2*x1 + 1 - s1)/(2*s1)
  x <- x1[1]
  y <- x1[2]
  1.5 - abs(x-y) + abs(x-y)^2
}

wd_g <- function(x1, s = 2, p = 2)
{
  x1 <- (2*x1 + 1 - s^p)/(2*s^p)
  5/3 - abs(x1)/4 - x1^2/4
}

md_f <- function(x1, s = 2, p = 2) 
{
  s1 <- s^p
  x1 <- (2*x1 + 1 - s1)/(2*s1)
  x <- x1[1]
  y <- x1[2]
  15/8 - abs(x)/4 - abs(y)/4 - 3*abs(x-y)/4 + abs(x-y)^2/2
}

sd_f <- function(x1, s = 2, p = 2)
{
  x1 <- x1/s^p
  x <- x1[1]
  y <- x1[2]
  res <- 0
  for (i in 0:p) res <- res + s^(-i)*delta(x,y,i,s)
  res
}

maxpro_f <- function(x1, s = 2, p = 2) 
{
  s1 <- s^p
  sigma <- 1/s1
  x1 <- x1/(s1-1)
  (abs(x1[1]-x1[2]) + sigma)^(-2)
}

others_g <- function(x1, s = 2, p = 2)
{
  0
}

F_FUNC <- function(x1, type = "CD2", s = 2, p = 2) { ## 9/17 changed the names
  type1 <- switch(type,
         CD2 = cd_f,
         WD2 = wd_f,
         MD2 = md_f,
         SD = sd_f,
         MaxPro = maxpro_f)
  
  type1(x1, s, p)
}

G_FUNC <- function(x1, type = "CD2", s = 2, p = 2) { ## 9/17 changed the names
  type1 <- switch(type,
                  CD2 = cd_g,
                  WD2 = wd_g,
                  MD2 = others_g,
                  SD = others_g,
                  MaxPro = others_g)
  
  type1(x1, s, p)
}

# Given a vector x, implement level permutation
permuteLevels1 <- function(x)
{
  y <- x
  lev <- unique(x)
  newlev <- sample(lev)
  for (i in 1:length(lev)) y[ x == lev[i] ] <- newlev[i]
  y
}

permuteLevels <- function(x)
{
  apply(x,2,permuteLevels1)
}

permuteLevels1_p2 <- function(x, s = 2) 
{
  x <- x - min(x)
  A <- x %/% s
  B <- x %% s
  
  A1 <- permuteLevels1(A)
  B1 <- B
  for (i in 0:(s-1)) B1[A == i] <- permuteLevels1(B[A == i])
  
  A1*s + B1
}

permuteLevels_p2 <- function(x, s = 2)
{
  apply(x,2,permuteLevels1_p2, s)
}

### Following functions give the criterion
UD_crit=function(x, crit="CD2")
{
  x = x - min(x) + 1	# levels 1, 2, ...
  UniDOE::DesignEval(x, crit=crit)
}

cd2=function(x) UD_crit(x, crit="CD2")
md2=function(x) UD_crit(x, crit="MD2")
wd2=function(x) UD_crit(x, crit="WD2")

MaxProMeasureX=function(x)
{ # scale x to [0,1] for any design
  x = x - min(x)
  x = x / max(x)
  MaxPro::MaxProMeasure(x)	
}

MaxProX = function(n, m, q)
{
  #	x0 = (x0 %% q)/(q-1) # needs to be in [0,1] not (0,1)
  
  x0 = matrix(rep(seq(from=0,to=1,length=q),each=n/q), nrow=n, ncol=m)
  x0 = rand_design_part2=apply(x0, 2, sample)
  a=MaxPro::MaxProQQ(x0); a$time; a$measure
  a$xbest = as.matrix( (q-1)*a$Design );	# *(q-1) instead of q	
  a
}

papply <- function(expr, rep = 10000) {
  seeds = sample(1e8, rep) # set rep seeds
  fn=function(i) {set.seed(seeds[i]); expr} 	# use distinct seeds for diffferent processes
  res=parallel::mclapply(1:rep, fn, mc.cores=parallel::detectCores()); 
  simplify2array(res) 
}

### This is the main function that takes the design and return SOA, Alg, OA
main_p2 <- function(x, s = 2, rep = 10000, rep2 = 100, iter = 1000, crit = "CD2")
{
  n <- nrow(x); m <- ncol(x); q <- s^2;
  crit1 <- switch(crit,
                  CD2 = cd2,
                  WD2 = wd2,
                  MD2 = md2,
                  MaxPro = MaxProMeasureX)
  
  res_oa <- papply(crit1(permuteLevels(x)), rep)
  res_soa <- papply(crit1(permuteLevels_p2(x, s)), rep)
  
  if (crit == "CD2") {
    res_alg <- papply(GenUD(n,m,q,crit = "CD2",maxiter = iter)$criterion_value, rep2)
  } else if (crit == "WD2") {
    res_alg <- papply(GenUD(n,m,q,crit = "WD2",maxiter = iter)$criterion_value, rep2)
  } else if (crit == "MD2") {
    res_alg <- papply(GenUD(n,m,q,crit = "MD2",maxiter = iter)$criterion_value, rep2)
  } else if (crit == "MaxPro") {
    res_alg <- papply(MaxProX(n,m,q)$measure, rep2)
  } 
  
  res <- cbind(res_soa, res_alg, res_oa)	
#  print(summary(res))
  res
}

#######################################
### Table 3-6, S2-S3  ### 
#######################################
generate_table <- function(s, p) { ## 9/17 changed the names and deleted a few columns
  
  final_res <- c()
  s_p <- 0:(s^p-1)
  
  for (crit in c("CD2", "WD2", "MD2", "SD", "MaxPro")) { 
    nrt <- data.frame(i = 0,
                      j = 0,
                      NRT = 0,
                      f = 0)
    
    k <- 1
    for (i in 1:length(s_p)) {
      for (j in 1:length(s_p)) {
          nrt[k,] <- c(s_p[i], s_p[j], nrt_distance(s_p[i],s_p[j],s, p), F_FUNC(c(s_p[i],s_p[j]), type = crit, s, p))
          k <- k + 1
      }
    }
    
    temp <- nrt %>%
      group_by(NRT) %>%
      summarise(Term = mean(f)) %>%
      mutate(NRT_new = paste0("NRT_", NRT))
    
    if (p == 2) {
      res <- temp %>%
        select(NRT_new, Term) %>%
        pivot_wider(names_from = NRT_new, values_from = Term) %>%
        mutate(y = (NRT_0 - NRT_1) / (NRT_0 + (s-1) * NRT_1 - s* NRT_2),
          z = (NRT_0 + (s-1) * NRT_1 - s* NRT_2)^2 / ((NRT_0 - NRT_1) * (NRT_0 + (s-1) * NRT_1 + s*(s-1)*NRT_2)),
          crit = crit)
    } else if (p == 1) {
      res <- temp %>%
        select(NRT_new, Term) %>%
        pivot_wider(names_from = NRT_new, values_from = Term) %>%
        mutate(yz = (NRT_0 - NRT_1) / (NRT_0 + (s-1)*NRT_1), 
          crit = crit)
    }
    
    final_res <- rbind(final_res, res)
  }
  
  final_res
}
