## Figure_S1.R : simulations for surrogate modeling 
##  clean vesion using Rcpp and dse_criteria.cpp 
## 
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

## 9/19/2026

library(DiceKriging)
require(gtools)
require(parallel) # mclapply
require(LHD)

library(Rcpp)
sourceCpp("dse_criteria.cpp") # Compile and load the external C++ file 
source("generate_soa.R") # get.soa.new() 
source("test_functions.R")	# test functions for surrogate modeling

## criteria using C functions in dse-criteria.cpp for speed
SD2=function(x, s=2) wsd2_yz(x, s=s, y=1, z=1, adjust=0) # similar to y=1/s^2, z=1, adjust=1
CD2=function(x)	GD2(x, crt="CD2")
WD2=function(x)	GD2(x, crt="WD2")
MD2=function(x)	GD2(x, crt="MD2")
Maxpro2=function(x,beta=0) maxpro(x, beta=beta)
DSE1=function(x, s=2) wsd2_yz(x, s=s, y=1/(s*ncol(x)+s), z=1, adjust=1) # DSE=SD2(w)
DSE2=function(x, s=2, y=1/(s*ncol(x)+s), z=1) 1+wsd2_yz(x, s=s, y=y, z=z, adjust=1)*(1-y)^ncol(x) # DSE=SD2(w)+1 
DSE=function(x, s=2, y=1/(s*ncol(x)+s), z=1) dse_yz(x, s=s, y=y, z=z) # DSE=DSE2=SD2(w)+1

## 
scale01=function(X, mid=F, rand=F)
{ # scale x to [0,1], modified on 7/15/22
	X=X-min(X)
	q = max(X)+1 # number of levels
	if(mid==T){
		if(rand==T){
			u=runif(nrow(X)*ncol(X),0, 1)
			(X+u)/q  # scale to a random point within each interval
		}
		else (X+0.5)/q  # scale to the middle point within each interval
	} 
	else X/(q-1) # scale to [0,1]
}
randomLHD=function(n,k)
{ # generate a random LHD
  x=matrix(0,n,k)
  for(j in 1:k) x[,j]= sample(1:n, n)
  x
}

gen.response1 = function(D, fun=sim.fun1)
{
   apply(D ,1, fun)  # apply fun to each row
}

gen.response = function(D, fun=sim.fun1)
{  
	if(fun$dim > ncol(D)){ # append 0.5 if needed 10/6/25
		Da = matrix(0.5, nrow=nrow(D), ncol=fun$dim-ncol(D))
		D = cbind(D, Da)
	}	
   apply(D ,1, fun$fun)  # apply fun to each row
}

randomPermuteSign=function(x)
{ # randomly permute levels
	if(max(x)>1) cat("error: max(x)>1")
  	n=ncol(x)
	sgn=sample(c(0,1), n, rep=T)
  	for(i in 1:n) if(sgn[i]==1) x[,i]=1-x[,i]  # assume x is within [0,1]
  	x
}
randomPermuteCol=function(x)
{ # randomly permute columns
  n=ncol(x)
  a=sample(1:n, n)
  xa=x[,a]
  dimnames(xa)[[2]]=dimnames(x)[[2]]  # change variable names
  xa
}


test.mse=function(fit0, X, Y)
{ # X is the test data
  pre0=predict(fit0,data.frame(X),"UK")$mean;
  mean((pre0-Y)^2)
}

fit.mse=function(x0, nugget=1e-8, xTest, yTest, sim.fun1=sim.fun1, covtype="matern5_2", permuteCol=TRUE)
{
	 # return normalized mse as in Chen et al. (2016)
	if(min(x0)<0 || max(x0)>1) D0=scale01(x0) # scale to 0-1
	else D0=x0 # x0 is within [0,1] already, do not rescale
  	if(permuteCol){
	  	D0 = randomPermuteCol(D0)
  	  	D0 = randomPermuteSign(D0)
  	  	}
	dimnames(D0)[[2]]=paste("V", 1:ncol(x0), sep="")
	y0 = gen.response(D0, fun=sim.fun1)	
#	mu0 = mean(yTest)
	mu0 = mean(y0)		# sample mean from the training data
	
	mse = NA
	for(i in 1:3){		# try 3 times if failed or mse > 0.9
 		fit0=km(~1, design=data.frame(D0), response=data.frame(y0), nugget=nugget, covtype=covtype, control = list(pop.size = 100, trace = FALSE))  # fit with D0 and y0
 
  		if(inherits(fit0, "try-error"))	next		# fitting error, try again
  		
   		mse= test.mse(fit0, xTest, yTest)/mean((yTest-mu0)^2)
   		if(mse < 0.9)	return(mse)
	}
  	mse
}  

fit.mse1=function(dummy=1, x0, nugget=1e-8, xTest, yTest, sim.fun1=sim.fun1, covtype="matern5_2", permuteCol=TRUE)
{ # used by mclapply
	fit.mse(x0, nugget=nugget, xTest=xTest, yTest=yTest, sim.fun1=sim.fun1, covtype= covtype, permuteCol= permuteCol)
}


search.subdesigns=function(x0, k, crt, maxK=10000, ...)
{ # search the best subdesign according to a criterion from all possible or maxK subdesigns
	
	m = ncol(x0)
	if(choose(m, k) > maxK){
		subs=matrix(0, maxK, k)
		for(i in 1:maxK) subs[i,]= sample(m, k)  # avoid to generate huge subsets
	}
	else subs = gtools::combinations(m, k) # all subsets
	
	val = apply(subs, 1, \(s) crt(x0[,s], ...))
	imin=which.min(val)
	isub = subs[imin,]
	list(xbest=x0[, isub], isub=isub, opt=val[imin], val=val )
}


get_design=function(N=64, k=15, type, readfile=F)
{	# start with a soa(N,k) + a LHD(N,k), then search.subdesigns by crt 8/11/26
	# levels: 1, ..., N	
	
	if(N>=64 && k<=15)	x1=get.soa.new(N,15)+1
	else x1=get.soa.new(N,k)+1		# may cause error if such an soa does not exist
	x1=LHD::OA2LHD(x1)
	x2=randomLHD(N,ncol(x1))
	x0=cbind(x1, x2)	# 
	
	if(N%%9 == 0)	s=3 else s=2
	
	x = switch(type, 
		# 8/11/26 search subdesigns
		soa = search.subdesigns(x1, k, crt=DSE, s=s)$x, # ma-type
	    DSE = search.subdesigns(x0, k, crt=DSE, s=s)$x, #
	    SD = search.subdesigns(x0, k, crt=SD2, s=s)$x, # 
	    CD = search.subdesigns(x0, k, crt=CD2)$x,
	    MD = search.subdesigns(x0, k, crt=MD2)$x,
	    WD = search.subdesigns(x0, k, crt=WD2)$x,
	    MaxPro = search.subdesigns(x0, k, crt=Maxpro2)$x, # =dmaxpro2
		lhd=randomLHD(N,k), # randomLHD 
		LHD=randomLHD(N,k) # randomLHD 
		)	
	x = x - min(x) + 1 # change levels from 1 to N 8/2/26
	return(x)	# x is an LHD
}



gen.TestX=function(k = pFac, nTest=10000)
{ 
	X=randomLHD(nTest, k)
	X=as.data.frame(X)
 	dimnames(X)[[2]]=paste("V", 1:k, sep="")
	X
}


sim1.mse=function(dummy=1, N, k=pFac, nugget=NULL, X.test=X.test, Y.test=Y.test, sim.fun1=sim.fun1, covtype="matern5_2")
{  # dummy is used for mclapply
# N is number of runs 	
## setting for Ying and Xu (2026+)

	# use soa or lhd as reference
	types=types_dse=c("CD",  "DSE", "lhd", "MaxPro", "SD", "WD") # 
	q=length(types)
	mse=rep(0,q)
	 for(i in 1:q){ 
	 	D=get_design(N, k, types[i])		
	  	mse[i] = fit.mse(D, nugget, xTest=X.test, yTest=Y.test, sim.fun1=sim.fun1, covtype=covtype, permuteCol=TRUE)
		}
	names(mse)=types
 	mse 	
}
Run.sim=function(id=1,N=64, k0=pFac, nRep=100, covtype="matern5_2", sim.mse=sim1.mse, fun.lib=AllFun.lib,  plot.mean=TRUE)
{
	nugget= 1e-8 # NULL # 1e-8; use a small nugget to avoid singular issue
  	mc.cores <- if(.Platform$OS.type == 'windows') 1 else parallel::detectCores()
	pFac=fun.lib[[id]]$dim;  
	fun.name=fun.lib[[id]]$name; sim.f1 = fun.lib[[id]]$fun;
	sim.fun1=fun.lib[[id]]; 	# 10/6/25 for isp1
	cat("===> Simulating", pFac, "-dim", fun.name, "function with", k0, "dimensions\n")
	nTest=10000 # test sample size
	X=gen.TestX(k=k0, nTest) # use randomLHD
	X.test=scale01(X) # scale to [0,1]
	Y.test=gen.response(X.test, fun=sim.fun1) # X is required in [0,1]

	print(date()) # the first argument of mcapply is X
	res <- mclapply(X=1:nRep, sim.mse, N=N, k=k0, nugget=nugget, X.test=X.test, Y.test=Y.test,  sim.fun1=sim.fun1, covtype=covtype, mc.cores = mc.cores)
	print(date())
	
	out = matrix(0, length(res), length(res[[1]])) 
	for(i in 1:length(res)) out[i,]=res[[i]]
	colnames(out)=names(res[[1]]); 
	rmse = sqrt(out) # normalized rmse
#	title = paste(fun.name, ":",  "N=", N, "k=", k0, "Rep=", nRep)
	if(plot.mean) print(plot_means_with_se(rmse, xlab="", title=fun.name)) # 8/29/26
	rmse
}



plot_means_with_se <- function(df, xlab="Design", ylab="Mean RMSE (±2 SE)" , title = "Comparison of Means with ±2 SE") 
{
	require(ggplot2)
	require(dplyr)
	require(tidyr)

 	df = as.data.frame(df) 
  # 1. Reshape and calculate mean and 2*SE for all numeric columns
  summary_df <- df %>%
    select(where(is.numeric)) %>% 
    pivot_longer(
      cols = everything(), 
      names_to = "Variable", 
      values_to = "Value"
    ) %>%
    group_by(Variable) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SE = sd(Value, na.rm = TRUE) / sqrt(sum(!is.na(Value))),
      Error = 2 * SE,
      .groups = "drop"
    )
  
  # 2. Plot using ggplot2
  p <- ggplot(summary_df, aes(x = Variable, y = Mean)) +
    geom_point(size = 3, color = "#2b5c8f") +
    geom_errorbar(   aes(ymin = Mean - Error, ymax = Mean + Error), 
      width = 0.2,    linewidth = 0.8,     color = "#333333"  ) +
    theme_minimal(base_size = 12) + labs( title = title, x = xlab, y = ylab) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      # Draws a border around the plot area:
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    )
  
  return(p)
}


Run=function()
{
#	source("Figure_S1.R") # source this file
	system.time(a<-get_design(64, 15, "DSE"))
	system.time(a<-get_design(64, 15, "SD"))
	system.time(a<-get_design(64, 15, "CD"))
	system.time(a<-get_design(64, 15, "MaxPro"))

	N=64;  nRep=200;   #  8/17/26 using soa64x15+lhd
	pdf(paste0("DSE_SOA64x15+lhd-R", nRep, ".pdf"), w=6, h=4)
	for(id in c(1:3)) { a=Run.sim(id=id, N=N, nRep=nRep, plot.mean=T);  print(apply(a,2,mean, na.rm=T)) }
	dev.off()
}

## Figure S1(a) OTL Circuit function
 a=Run.sim(id=1, N=64, nRep=200, plot.mean=T)	# about 2.5 minutes with 10 cores

## Figure S1(b) OTL Circuit function
 a=Run.sim(id=2, N=64, nRep=200, plot.mean=T)	# about 3.5 minutes with 10 cores

## Figure S1(c) OTL Circuit function
 a=Run.sim(id=3, N=64, nRep=200, plot.mean=T)	# about 6.5 minutes with 10 cores
