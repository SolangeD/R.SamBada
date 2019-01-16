require(Rcpp)

cppFunction("
            NumericVector exfun(NumericVector x, int i){
            x = x*i;
            return x;
            }")

cppFunction("
NumericVector MGF(NumericMatrix xx, double maxMGFAllowed){
int xrow = xx.nrow() ;
int xcol = xx.ncol();
int aa;
int Aa;
int AA;
int sum_max;
int mgf;
NumericVector yy(xcol);
NumericVector yybool(xcol);
NumericVector yysnpid(xcol);
int k=0;
for(int i = 0; i < xcol; i++) {
	aa=0;
	Aa=0;
	AA=0;
	for(int j = 0; j < xrow; j++){
		if(xx(j,i)==0){
			aa++;
		}
		else if(xx(j,i)==1){
			Aa++;
		}
		else if(xx(j,i)==2){
			AA++;
		}
	}
  if(aa>=Aa){
    sum_max=aa;
    if(AA>aa){
    sum_max=AA;
    }
  }
  else{
    if(AA>=Aa){
      sum_max=AA;
    }
    else{
      sum_max=Aa;
    }
  }
  if(aa+AA+Aa>0){
    mgf=(sum_max*100)/(aa+Aa+AA);
  }
  else{
    mgf=0;
  }
  yy(i)=mgf;
  if(mgf>maxMGFAllowed*100){
    yybool(i)=0;
    
  }
  else{
    yybool(i)=1;
    yysnpid(k)=i+1;
    k++;
  }

}

//NumericVector yysnpid2(k);
//yysnpid2(yysnpid.begin() , yysnpid.begin() + k);
if(maxMGFAllowed>=0){
  return yysnpid;
}
else{
  return yy;
}
}")

