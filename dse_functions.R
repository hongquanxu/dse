## Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 
##
##  key functions: 
##		dse() -  compute the double stratification enumerator via Lemma 2
##		dse.yz() and DSE() - compute the double stratification enumerator via Theorem 2
##		fastPijFull - compute the stratification patten based on Theorem 5
##		fastPij - compute the stratification patten based on Theorem S1 
##
## Date: 9/18/2026

### functions 
###
## generate a full s^k factorial design,
full.fd=function(s,k)
{ # generate s^k full factorial design, 4/26/20
	if(k==1) return(matrix(0:(s-1),ncol=1))
	x0=full.fd(s,k-1)
	xk=cbind(0,x0)
	for(i in 1:(s-1)){
		xk=rbind(xk, cbind(i,x0))
	}
	dimnames(xk)=list(1:s^k, 1:k)
	xk
}

## kernel R() as Lemma 2 of Ying and Xu (2026+)
r.kernel <- function(s, p, y, z) {
	syk <- function(y, k) if(y != 1) (1-y^k)/(1-y) else k	# avoid dividing by 0
  Fd = full.fd(s, p)
  Ker = matrix(0, s^p, s^p)
  for(u in 1:s^p) for(v in 1:u){
    if (u == v) Ker[u,v] = Ker[v,u] = 1 - z + z * ( (1-y)*syk(s*y, p+1) + s^p*y^(p+1)) 
    else {
      k <- p + 1 - which(abs(Fd[u,] - Fd[v,]) != 0)[1]
      Ker[u,v] = Ker[v,u] = 1 - z + z * (1-y)*syk(s*y, p-k+1)
    }
  }
  Ker
}

## Definition 1
dse <- function(x, s, p, y, z) {
# x is a matrix with levels from 0,1,...,s^p
   	x = x - min(x)	# ensure starting level at 0
	if(missing(p))	p=floor(log(max(x)+1+1e-8, s)) # default s^p as number of levels
	if(s^p != max(x)+1)	stop("s^p != max(x)+1")	
 	if(missing(y))	y=1/(s*ncol(x)+s) # set default y as Theorem 4
  	if(missing(z))	z=1 # set default z=1

 	N=nrow(x); n=ncol(x); 
   	R = r.kernel(s, p, y, z)

  res = 0
  for(a in 1:N) for(b in a:N){
    rk = 1
    for(j in 1:n) rk = rk * R[x[a,j]+1, x[b,j]+1] 
    if(b>a) res = res + 2* rk  # (b,a)
    else res = res + rk  # a==b
  }
  
  res/N^2
}

### functions from SD2.R (Tian and Xu 2026, JRSSB)
###

### functions for the stratified L2-discrepancy 
sd.kernel=function(q, k, w=rep(1,k), w0=1)
{ # w is a vector of k weights, 3/7/21
# return a q^k * q^k matrix K(x,y)
	Fd = full.fd(q, k) # q^k full factorial
	Ker = matrix(0, q^k, q^k)
	for(x in 1:q^k) for(y in 1:x){
		k1 = w0
		for(j in 1:k){
			if(Fd[x,j] != Fd[y,j]) break  # no need to continue
			k1 = k1 + w[j]/q^j
		}  
		Ker[x,y] = Ker[y,x] = k1
	}
	Ker   # include the first term 1
}

sd2 = sd2.w=function(x, q=2, k=floor(log(nrow(x)+1e-8, q)), w=rep(1,k), w0=1, linear=F)
{ # stratified L2-discrepancy with weights, Theorem 1
# x: N*n design, within [0,1) or scaled to [0,1)
# q, k: stratification parameters. Each dimension is divided into q^k intervals.
# w: a vector of weights with length k

	N=nrow(x); n=ncol(x); 
	if(max(x)-min(x)>=1){ # 
		s = max(x) - min(x) + 1	 # number of levels
		x = x - min(x)  # coded as 0:(s-1)
	# change levels (0:(s-1)) to (0,1)
		x = (x+0.5)/s
	}
	
	if(max(x)-min(x)<1){ # 0<=x<1
		x = floor(q^k * x) + 1 # convert to [0,q^k-1] +1
	}	
	Kw = sd.kernel(q, k, w, w0)
	
	res = 0
	if(linear){		# linear 3/20/25
		for(a in 1:N){
			pk = 1;	 b=1	# set b=1 for linear
			for(j in 1:n) pk = pk * Kw[x[a,j], x[b,j]] 
			res = res + pk  
			}
		res = res * N
	}
	else{	# nonlinear
		for(a in 1:N) for(b in a:N){
			pk = 1
			for(j in 1:n) pk = pk * Kw[x[a,j], x[b,j]] 
			if(b>a) res = res + 2* pk  # (b,a)
			else res = res + pk  # a==b
		}
	}
	res1= w0
	for(j in 1:k) res1 = res1 + w[j]/q^(2*j)
	res/N^2 -(res1)^n
}

sd2.y=function(x, q=2, k=floor(log(nrow(x)+1e-8, q)), y=.01, adjust=F, linear=F)
{ # stratified L2-discrepancy with exponenital weighting scheme (7) or special weights in Theorem 4 
# x: N*n design, y can be a complex number
# q, k: stratification parameters. Each dimension is divided into q^k intervals.
# 0 <y <1

	if(adjust == T){	# adjust weights as in Theorem 4
		w=(q^2*y)^(1:k); 
		w[k]=(q^2*y)^k/(1-y)
	}
	else 	w=(y)^(1:k); 	# exponenital weighting scheme (7)
	sd2.w(x,q,k, w=w, w0=1, linear=linear)
}

## ----------------------------
## new functions for Ying and Xu (2026+)
##
## dse.yz and sd2.yz: extension of sd2.y with y and z 
## fastPij : fast algorithm for computing Pij, 
##

sd2.yz  <- function(x, q=2, k=floor(log(nrow(x)+1e-8, q)), y=1, z=1, adjust=0, linear=F)
{ # stratified L2-discrepancy with exponenital weighting scheme or special weights in Theorem 2 of Ying and Xu (2026+) 
# x: N*n design, within [0,1) or scaled to [0,1)
# q, k: stratification parameters. Each dimension is divided into q^k intervals.
#  y for volume, z for dimension
#  match wsd2_yz in dse-criteria.cpp 

	if(adjust==0) { w0 = 1; w=(y)^(1:k) * z; }	# exponenital weighting scheme; 
	else if(adjust != 0){	# adjust weights as in Theorem 2 of Ying and Xu (2026+) 
		w0=1-y*z
		w=(q^2*y)^(1:k) *(1-y)*z; 
		w[k]=(q^2*y)^k * z 
	}
	
	if(adjust == 1){ # normalize w0=1, 9/1/26
		w = w/w0
		w0 = 1 # 
	}
	sd2.w(x,q,k, w=w, w0=w0, linear=linear)
}

## dse.yz uses adjust=2 as in Theorem 2 of Ying and Xu (2026+)
dse.yz <- function(x, s, p, y=1, z=1,linear=F) 1 + sd2.yz(x, s, p, y=y, z=z, adjust=2, linear=linear) # use different default y as dse 
DSE <- function(x, s, p, y=1/(s*ncol(x)+s), z=1) 1+sd2.yz(x, s, p, y=y, z=z, adjust=2)	# match dse  with same s and p (for LHDs)

## -------------------
## fastSP2 using sd2.yz, 
##
fastSP2=function(x, q, k=floor(log(nrow(x)+1e-8, q)), y0=0.1, limit=5, K, linear=F, digits=3)
{ # same as fastSP.K, which uses EDz() while fastSP2 uses sd2.yz(), 3/12/25
 # see SOAs::fastSP, which is a variant of fastSP.K
	# compute the first few components of SPattern using EDz and complex numbers 
	# EDz() uses Rd.kernel() as Lemma 1 of Tian and Xu (2024)
	# x is a GSOA with q^k levels
	# q is the base and k is the length
	# limit is the number of components of SPattern computed
	# K is the number of components as in Theorem 4 of Tian and Xu (2024)
	# when k=1, it returns tha usual GWLP as in Xu and Wu (2001)

	 x= as.matrix(x) # treat a vector as a matrix
	 x = x - min(x)	# set lowest level to 0
	 n=ncol(x);  
	 if(missing(K)) K=ceiling((log(0.01)-log(q^(k*ncol(x))/nrow(x)-1))/log(y0)) # default K as in Tian and Xu (2024)
	 K = max(K, limit)	# K must be >= limit
#	 omega = as.complex(exp(1i*2*pi/K))  # omega is a Kth root of 1, omega^K=1
	 z= bz = Ez=complex(K)
	 for(j in 1:K) z[j] = as.complex(exp(1i*2*j*pi/K))  # z[j]=omega^j 
#	 for(j in 1:K) Ez[j] = EDz(x, q, k, y=z[j]*y0) -1 # y[j]=omega^j*y0
	 for(j in 1:K) Ez[j] = dse.yz(x, q, k, y=z[j]*y0, z=1, linear=linear) - 1   # z=1
	 for(i in 1:limit){
	 		ij=(i*(1:K)) %% K
	 		bz[i] = sum(z[K-ij] * Ez)/K/exp(i*log(y0)) # bz[i]=sum_{j=1}^K (omega^(-i*j)*Ez[j])/K/y0^i
	 } 
	 if(max(abs(Im(bz)))>10^-5) cat("max(abs(Im(bz))) =", max(abs(Im(bz))), '\n')
	 round(Re(bz)[1:limit], digits=digits)  # Im(bz) should be zero
}


## fast algorithm based on Theorem S1 of Ying and Xu (2026+)
fastPij=function(x, q, k=floor(log(nrow(x)+1e-8, q)), u=0.1, v=0.1, imax=5, jmax, K, L, linear=F, digits=3)
{ # compute Pij with sd2.yz via complex numbers
	# q is the base and k is the length
	# when k=1, it returns tha usual GWLP as in Xu and Wu (2001)
	# 3/22/25 add v and L to deal with designs with large number of columns
	 x= as.matrix(x) # treat a vector as a matrix
	 x = x - min(x)	# set lowest level to 0
	 m=ncol(x);  # dim
	 if(missing(jmax))	jmax=min(imax, m)	
	 if(missing(K)) K=min(m*k, 20) # 
	 if(missing(L)) L=min(m, 10) # 
	 K = max(K, imax)
	 L = max(L, jmax)
	 y = complex(K);  
	 z = complex(L)
	 Pij = Eyz = matrix(0, L, K)
	 for(a in 1:K) y[a] = as.complex(exp(1i*2*a*pi/K))  #
	 for(b in 1:L) z[b] = as.complex(exp(1i*2*b*pi/L))  # dim
	 for(a in 1:K) for(b in 1:L) Eyz[b, a] = dse.yz(x, q, k, y=y[a]*u, z=z[b]*v, linear=linear) - 1 
#	 for(a in 1:K) for(b in 1:L) Eyz[b, a] = dse(x, q, k, y=y[a]*u, z=z[b]*v) - 1 
	 for(j in 1:jmax){
	 	for(i in j:min(k*j,imax)){		# 1/5/26
	 		ia=(i*(1:K)) %% K;	yia = y[K-ia]
	 		jb=(j*(1:L)) %% L; zjb = z[L-jb]
	 		Pij[j, i] = sum( zjb %*% Eyz %*% yia)/(K * L * u^i * v^j)
	 		}
	 } 
	 if(max(abs(Im(Pij)))>1e-5) cat("Im(Pij) =", max(abs(Im(Pij))), '\n')
	 round(Re(Pij)[1:jmax, 1:imax], digits=digits)  # Im(Pij) should be zero
}


fastPij2=function(x, q, k=floor(log(nrow(x)+1e-8, q)), imax=5, u=0.1, v=0.1, ... )
 fastPij(x, q, k, u=u, v=v, imax=imax, ...) # 


## fast algorithm based on Theorem 5 of Ying and Xu (2026+)
fastPijFull=function(x, q, k=floor(log(nrow(x)+1e-8, q)), imax=k*ncol(x), digit=4)
 fastPij(x, q, k, u=1, v=1, imax=imax, K=k*ncol(x), L=ncol(x), digit=digit) # fix u=1 and v=1


## ----------------

