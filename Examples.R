## Examples 4, S1 and S2; Tables 2 and S1
##
## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

library(SOAs) 
source("dse_functions.R")
source("generate_soa.R")
	
#### Example 4/Table 2
####  compute partial Stratification Pattern for GSOA(81,9,2^3,3)

	x = read.csv("soa81x9.csv")[,-1]; dim(x) 
	s=3; p=2; n=nrow(x); m=ncol(x); 	s^(m*p)/n-1
	sd2.yz(x,s,p,adj=2) # = s^(m*p)/n-1
	system.time(pij <- fastPijFull(x, 3, 2)); # 2.4s for entire SP 
	sum(pij) # sum of Pij = 3^18/81-1
	3^18/81-1 # whole SP
	tab_2 = pij[1:6, 1:12]; tab_2

	# compare with SOAs::Spattern 
	system.time(b <- Spattern(x, 3, maxwt=6, maxdim=3) ); dim_wt_tab(b)	# about 1.7s
	system.time(b <- Spattern(x, 3, maxwt=8, maxdim=4) ); dim_wt_tab(b)	# about 21.6s
	system.time(b <- Spattern(x, 3, maxwt=10, maxdim=5) ); dim_wt_tab(b)	# about 182s
	system.time(b <- Spattern(x, 3, maxwt=12, maxdim=6) ); dim_wt_tab(b)	# about 1043s or 17m38s
#	date(); system.time(b <- Spattern(x, 3, maxwt=14, maxdim=7) ); dim_wt_tab(b)	# about 4065s or 68m

## Example S1/Table S1
##  compute partial Stratification Pattern for GSOA(81,9,2^3,3)

	x = read.csv("soa81x9.csv")[,-1]; dim(x) 
	system.time(tab_S1<-fastPij(x,3,2, u=0.1, v=0.1, K=8,L=4, imax=8, jmax=4,  digits=3)) # 0.5s for partial SP
	tab_S1

## Example S2: compute partial Stratification Pattern for SOA(128,31,8,3)
##
 	x=get.soa.new(128,31); dim(x)
	system.time(pij <- fastPijFull(x, 2, 3, imax=5)); pij # 338s, with wrong Pij<0
	
	system.time(pij <- fastPij(x, 2, 3, u=0.1, v=0.1, imax=5, K=10, L=6)); pij	 # 8s
	m=31; choose(m,3)/4 - choose(m,2)/12 # match P44=1085 Theorem 6 of Shi and Xu (2024)

	system.time(Spat <- (Spattern(x, s = 2, maxwt=4))); dim_wt_tab(Spat)	# 14s
	system.time(Spat <- (Spattern(x, s = 2, maxwt=5))); dim_wt_tab(Spat)	# 114s


## Tables 3-6, S2, and S3
## see Tables.R
