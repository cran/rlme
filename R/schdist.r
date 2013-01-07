schdist = function(n,schpar,cn=0){
  #Generating three-level nested data data in simulation
  #Not used in rlme
  
	eps=.20
	sigmac=5
	ind=rbinom(n,1,eps)	#eps contaminated percentage, sigmac is sd ratio 
	
	if(cn==0){ schdist = rnorm(n,schpar[1],schpar[2])}
	if(cn==1){x=rnorm(n,schpar[1],schpar[2]);schdist=x*(1-ind)+sigmac*x*ind}
schdist
}
    
