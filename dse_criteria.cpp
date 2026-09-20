/*  Rcpp file called by Figure_S1.R for computing double stratification enumerator (DSE), stratified L2-discrepancy, and other criteria

Ying and Xu (2026+). Efficient Representation and Construction of Space-Filling Designs via a Double Stratification Enumerator. 

 Date: 9/1/2026
*/

#include <Rcpp.h> 

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix nrt_kernel(int s, int p)
{  // return s^p * s^p kernel matrix; see nrt.kernel() in gwlp.R
	int Sp = (int) pow(s, p);
	
	// first generate the nrt.kernel()
	NumericMatrix Kw(Sp, Sp);	
	for (int i=0; i<Sp; i++) 
		for(int j=0; j< Sp; j++)	
			Kw(i, j) = p; // initial value p
  
    for (int i=p; i>0; i--){
    		int si = (int) pow(s, i-1);		// block size
          for(int j=0; j< (Sp/si); j++){
          	int block=si*j;		// jth block
            for (int k=0; k<si; k++)
            	for (int l=0; l<si; l++)
            		Kw(block+k, block+l) -= 1 ;
            }
       }
    return Kw;	// Sp*Sp matrix
}

// [[Rcpp::export()]]
NumericMatrix wsd2_kernel(int s, int p, NumericVector w)
{  // return s^p * s^p kernel matrix; see sd.kernel() in gwlp.R
// w is a vector of length p+1, including w0
 // w[i] =pow(y, i)
	int Sp = (int) pow(s, p);
	
	// first generate the nrt.kernel()
	NumericMatrix Kw = nrt_kernel(s, p);	
 
    // Kw is nrt.kernel, now change it to sd.kernel using Thereom 3; see wsd2.g
    for (int i=0; i<Sp; i++){
          for(int j=0; j< Sp; j++){
          	int d = Kw(i,j);
          	Kw(i,j) = w[0]; // w0
          	if(d < p)  for (int k=1; k<= p-d; k++)	Kw(i,j) += w[k]/pow(s, k) ;  // w[k]=y^k * z
            }
       } 
   	return Kw;	// Sp*Sp matrix
}

// [[Rcpp::export()]]
NumericVector InitWeight_spyz(int s, int p, double y=1, double z=1, int adjust=0)
{ // to match wsd2.spyz
	NumericVector w(p+1, 1.0); // = rep(1.0, p+1); // w[0]-w[p]
	w[0] = 1.0; 
	for(int i=1; i<=p; i++) w[i] = pow(y, i) * z; // set exponential weights
	
	if(adjust != 0){	// set weights as in Theorem 2 (Ying and Xu 2026+)
		w[0] = 1 -y*z;
		for(int i=1; i<p; i++) w[i] = pow(s*s*y, i) * (1-y)*z;
		w[p] = pow(s*s*y, p)*z;
	}
	if(adjust == 1){	// normalize w[0]=1, 9/1/26
		w = w / w[0]; 		
	}
	return (w);
}

//// [[Rcpp::export()]]
double wsd2_g(NumericMatrix x, int s, int p=0, NumericVector w={1,1})
{  // should match wsd2.g  8/11/26
	int 	n = x.nrow(), m= x.ncol();
	if(p<=0) p = (int) floor(log(n+pow(1.0,-8))/log(s)); 
	if(w.size()<p+1) w = rep(1.0, p+1);
	NumericMatrix Kw = wsd2_kernel(s, p, w);		
	
	double res0 = 0;  // 
	for(int j=0; j<=p; j++)	res0 += w[j]/pow(s, 2*j);	
	
	if(min(x)==0)	x = x+ 1;
	double	q = max(x) - min(x) + 1;	// number of levels
	if(q < 2) Rcout << "Error in wsd_type: q < 2!";

	double res = 0.0;
		for(int i = 0; i < n; i++ ){
			for(int j = i; j < n; j++ ){
				double pk = 1;	// 			
				for(int k = 0; k < m; k++ ){
					// x() are within 1 and q, convert them to [0,1]
					double z1=(x(i,k) -0.5)/q, z2=(x(j,k) -0.5)/q;
					// convert z1 and z2 to index of Kw
					int x1=(int)(Kw.nrow() * z1), x2=(int)(Kw.nrow() * z2);
					pk *= Kw(x1, x2);	// K(x_{ik}, x_{jk})
				}	// for k
				res += pk;
				if(j > i)	res += pk;	// 
			} // for j
		}	// for i
		
		res /= pow(n, 2)  ;  // may differ by a constant, -pow(res0, m)
		res -= pow(res0,m);
		return (res);	
}

// [[Rcpp::export()]]
double wsd2_yz(NumericMatrix x, int s=2, int p=0, double y=1, double z=1, int adjust=0)
{  // should match sd2.yz, wsd2.yz and wsd2.spyz 8/11/26
	int 	n = x.nrow();  // m= x.ncol();
	if(p<=0) p = (int) floor(log(n+exp(-8))/log(s)); 

	NumericVector w = InitWeight_spyz(s, p, y, z, adjust);
	return wsd2_g(x, s, p,  w);
}

// [[Rcpp::export()]]
double dse_yz(NumericMatrix x, int s=2, int p=0, double y=1, double z=1)
{  // use adjust=2 as in Theorem 2 of Ying and Xu (2026+), 9/1/26
	return 1 + wsd2_yz(x, s, p, y, z, 2);	// fix adjust=2
}


// other criteria: maxpro and cd2, wd2, md2
// 

// [[Rcpp::export()]]
double maxpro(NumericMatrix D, double beta=0, int p_norm=2, bool scale=true)
{ // same as bid_log 
  int n=D.rows(), m=D.cols();
  int q=max(D)-min(D)+1;		// number of levels
	if(beta==0 && q<n)	return(log(0.0));		// not properly defined
	
  if(scale)	 D = D/q; 		// normalized to [0,1]

  double obj = 0;
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double prod = 1.0;
      for(int k=0; k<m; k++){
      	double dijk = (double) fabs(D(i,k) - D(j,k)); 
  		if(beta==0 && dijk==0){	 // happens if q < n
  			continue;
  			prod *= 1;	// 
  		}
        prod *= pow(beta + pow( dijk , p_norm), -2/p_norm);
      }
      obj += prod;
    }
  }
  return	pow( obj / (n*(n-1)/2), 1.0/m );
}


// Kernel functions for CD2/WD2/MD2
double Kg(double z, int type)
{
	double res = 0;
	z = fabs(z -0.5); // center
	if(type==-15)	res=1 + z*(1-z)/2; //cd2
	else if(type==-16)	res = 0; // wd2
	else if(type==-17)	res = 5.0/3 -z*(1+z)/4; // md2
	return(res);
}

double Kf(double z1, double z2, int type)
{
	double res = 0;
	double diff = fabs(z1-z2);
	z1 = fabs(z1 -0.5); z2 = fabs(z2 -0.5); // center
	if(type==-15)	res = 1+(z1+z2-diff)/2; //cd2
	else if(type==-16)	res = 1.5 - diff*(1-diff); // wd2, a typo in Tian and Xu (2026)
	else if(type==-17)	res = 15.0/8 -(z1+z2)/4 - 3*diff/4 +diff*diff/2; // md2
	return(res);
}

// [[Rcpp::export()]]
double GD2(NumericMatrix D, String crt="CD2")
{ // crt= "CD2", "WD2", "MD2"
 // implement generalized L2-discrepancy, include cd2, wd2, md2 
  int n=D.rows(), m=D.cols();
  D = D - min(D);	// starting level at 0
  int q=max(D)+1;		// number of levels
  D = (D+0.5)/q; 		// normalized to (0,1)

	double res0=0;
	int type = 0;
	if(crt=="CD2")	{type=-15; res0 = pow(13.0/12, m); } //cd2
	else if(crt=="WD2")	{type=-16; res0 = -pow(4.0/3, m); } // wd2
	else if(crt=="MD2")	{type=-17; res0 = pow(19.0/12, m); } // md2	

  double res1 = 0;
  for(int i=0; i<n; i++){
  	double pk1 = 1.0;
  	for(int k=0; k<m; k++)	pk1 *= Kg(D(i,k), type); 
  	res1 += pk1;
  	}
  		
  double res2=0;
  for(int i=0; i<n; i++){
	for(int j=i; j<n; j++){
      double pk2 = 1.0;
      for(int k=0; k<m; k++){
  			double z1=D(i,k), z2=D(j,k);
			// z1 and z2 in (0, 1)
			pk2 *= Kf(z1, z2, type);	 // K(x_{ik}, x_{jk})
			}	// for k
		res2 += pk2;
		if(j > i)	res2 += pk2;		// 
		} // for j
	}	// for i
		
	return ( res0 -2.0*res1/n + res2/pow(n, 2) ) ;  
}

