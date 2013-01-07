errdist = function(n,errpar,cn=0){

  #Generating three-level nested data data in simulation
  #Not used in rlme
  #when no cn, it's case 1


	eps=.20
#The persentage of contamination is 20%
	sigmac=5	
#4  makes the ratio of the contaminated variance to the uncontaminated is set at 16.
	ind=rbinom(n,1,eps)		#eps contaminated percentage, sigmac is sd ratio 
	
	if(cn==0){ errdist = rnorm(n,errpar[1],errpar[2])}
	if(cn==1){ x=rnorm(n,errpar[1],errpar[2]) ; errdist =x*(1-ind)+sigmac*x*ind }
 
errdist
}

