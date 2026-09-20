## R code from Shi and Xu (2024, JASA), "A Projection Space-Filling Criterion and Related Optimality Results"
## main function: get.soa.new() to generate SOAs according to Shi and Tang (2020) and Shi and Xu (2024)
## Shi and Xu, 3/11/2023

## functions to generate SOAs
ffd2gen = function(gen)
{ # generate 2-level ffds with generator matrix 
    if(is.vector(gen) || nrow(gen) == 1)     rbind(0, gen) %% 2
    else{
        g <- gen[nrow(gen),] # last row
        x <- ffd2gen(gen[-nrow(gen),])  
        x <- t(x)
        x<- cbind(x, x+g ) %% 2
        x<- t(x)
        dimnames(x) <- list(1:nrow(x), 1:ncol(x))
        x
        }
}

ffd.s = function(n, s=2, rev=T)
{ # a simple function to generate 2-level ffd
	if(s !=2) stop("s !=2: not implemented")
	Identity = function(n) diag(rep(1,n))
    gen <- t(ffd2gen(Identity(n)))[,-1]	# remove col of 0's
    if(rev) x <- ffd2gen(gen[n:1,])  # reverse the order so that it is easy to read
    else x <- ffd2gen(gen)  
	list(code=x)
}

get.xyz=function(k)
{ # create columns of X, Y, Z using Yates orders
	if(k==2) return(cbind(x=c(1,2,3), y=c(2,3,1), z=c(3,1,2)))
	else if(k==3) return(cbind(x=c(1:7), y=c(7,5,2,1,6,4,3), z=c(6,7,1,5,3,2,4)))
	xyz=get.xyz(k-2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	x=c(0, xyz[,1]); x=c(x, k1+x, k2+x, k3+x)[-1]
	y=c(0, xyz[,2]); y=c(y, k2+y, k3+y, k1+y)[-1]
	z=c(0, xyz[,3]); z=c(z, k3+z, k1+z, k2+z)[-1]
	cbind(x=x,y=y,z=z)
}
soa.ST5=function(k, iC=k1, printC=F)
{ # construct SOA(n=2^k,n/4-1,8, 3) as in Shi and Tang (2020, Theorem 5)
	# iC=k1=n/4 is the unique choice.
	xyz = get.xyz(k-2)
	a=ffd.s(k,2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	A=a$code[,c(k1+xyz[,1])]
	B=a$code[,c(k2+xyz[,2])]
	B1=a$code[,c(k3+xyz[,3] )] # B1 =A+B
	setAB=c(k1+xyz[,1], k2+xyz[,2], k3+xyz[,3], xyz[,1]) # A, B, B1, A2
	setC=(1:(2^k-1))[-setAB] # all possible C not in (A, A2, B, B1)
	if(printC) print(setC)
#	cc=rep(setC, ncol(A))[1:ncol(A)]
	if(missing(iC)) iC=setC[1]
	cc=rep(iC, ncol(A))[1:ncol(A)]
	C=a$code[, cc ] # avoid C in 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

soa.ST5x=function(k, iC=1, printC=F)
{ # construct SOA(n=2^k,n/4,8, 3) as in Shi and Tang (2020, corollary 1)
	# iC can be from 1:(k1-1)
	xyz = get.xyz(k-2)
	a=ffd.s(k,2)
	k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	A=a$code[,c(k1, k1+xyz[,1])]
	B=a$code[,c(k2, k2+xyz[,2])]
	B1=a$code[,c(k3, k3+xyz[,3] )] # B1 =A+B
	setAB=c(k1+xyz[,1], k2+xyz[,2], k3+xyz[,3], xyz[,1]) # A, B, B1, A2
#	setC=(1:(2^k-1))[-setAB] # all possible C not in (A, A2, B, B1)
#	if(printC) print(setC)
#	cc=rep(setC, ncol(A))[1:ncol(A)]
#	if(missing(iC)) iC=setC[1]
	cc=rep(iC, ncol(A))[1:ncol(A)]
	C=a$code[, cc ] # avoid C in 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

get.abc=function(k)
{ # initial A, B, C as in Section 4.1 for constructing SOA(n=2^k,5n/16,8, 3)	 or SOA(32,9,8,3)		
	if(k==4) return(cbind(x=c(1,2,4,8,15), y=c(12,9,3,6,5), z=c(11,7,9,5,12))) 
	else if(k==5) return(cbind(x=c(1,2,4,8,16,15,19,21,22), y=c(24,28,3,5,10,6,12,14,11), z=c(6,5,10,20,7,12,7,18,17)))
	else if(k==7){
	# create columns of A, B, C using Yates orders as Theorem 3 for N=128
A=c(1,2,4,8,15,17,18,20,24,31,33,34,36,40,47,49,50,52,56,63,65,66,68,72,79,81,82,84,88,95,97,98,100,104,111,113,114,116,120,127);
					B=c(42,37,25,3,117,74,41,10,14,102,92,69,23,6,83,90,73,71,21,86,54,28,7,5,57,61,44,26,19,53,60,12,9,13,58,55,62,35,27,38);
					C=c(60,27,90,69,54,61,119,71,59,108,75,105,77,85,53,76,78,125,59,106,61,101,122,107,53,42,119,109,46,102,74,9,3,6,73,10,69,42,5,21); # change C[40]=22
	 return(cbind(x=A, y=B, z=C))
	}

	if(k<4) stop("k<4")
	xyz=get.abc(k-2)
	k0=0; k1=2^(k-2); k2=2^(k-1); k3=k1+k2
	x=xyz[,1]; x=c(x, k1+x, k2+x, k3+x)
	y=xyz[,2]; y=c(y, k2+y, k3+y, k1+y) # original as in Shi and Tang (2020)
	z=xyz[,3]; z=c(z, k2+z, k3+z, k1+z)
	cbind(x=x,y=y,z=z)
}

soa.ST3=function(k)
{ # constructing SOA(n=2^k,5n/16,8, 3) as Theorem 3 in Section 4.1
# except for k=5, soa(32,9,8,3)
	xyz = get.abc(k)
	a=ffd.s(k,2)
	A=a$code[,(xyz[,1])]
	B=a$code[,(xyz[,2])]
	C=a$code[,(xyz[,3])] # 
	# check if C=A+B for any column
	x=(A+B+C)%%2; if(sum(apply(x,2,sum)==0)) print("A+B=C for some columns\n")
	x=4*A+2*B+C
	x
}

get.soa.new=function(N=64, k=15, rand=F)
{ 
	p=log(N,2);
	if(N != 2^p) stop("N !=2^p")
	if(k>N/4) x=soa.ST3(p) # Theorem 3, m=5N/16
	else if(k==N/4) x=soa.ST5x(p) # Shi and Tang (2020, corollary 1)
	else x=soa.ST5(p) # Theorem 6, m=N/4-1; Shi and Tang (2020, Theorem 5)
	if(rand) x=x[,sample(ncol(x), k)]
	x[,1:k]	
}


