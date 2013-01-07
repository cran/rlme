wilstep = function(I,sec,mat,init=F,y,x,sigmaa2=1,sigmaw2=1,sigmae2=1,thetaold=c(0),eps = .0001,iflag2=0, rprpair = 'hl-disp'){

#wwest using wil scores, yields fixed estimates(theta) abd sigmas
#or, use rfit.
  
    location = scale = 2
    if(rprpair == 'med-mad') {
      location = scale = 1
    }
  
     n = length(y) 

     if(init == T){
         fitw = wwest(x,y,print.tbl=F) #Actually, we can use rfit if faster!!
         theta = fitw$tmp1$coef		#beta estimates are stored in theta.
         ehat = fitw$tmp1$residuals	#ehat=a+w+epl
#		tauhat=fitw$tmp2$tau	#wilc scale est
#		taus=fitw$tmp2$tau1		#scale est in intercept, L1
         	fitvc = rprmeddis(I,sec,mat,ehat,location,scale, rprpair = rprpair)
         	sigmaa2 = fitvc$siga2
         	sigmaw2 = fitvc$sigw2
 	 	sigmae2 = fitvc$sigmae2
         	rea = fitvc$frei
         	rew = fitvc$frew
      	ree = fitvc$free	#epsilon error (not raw)
	
     } else {
         sigmay = sigymake(I,sec,mat,sigmaa2,sigmaw2,sigmae2)
         sigma12inv = matrix(sigmay$sigy12i,ncol=n)
         ystar = sigma12inv%*%y         
         xstar = sigma12inv%*%cbind(rep(1,n),x)
         fitw = wwest(xstar,ystar,print.tbl=F)         
         yhat = ystar - fitw$tmp1$residuals
         fitcsp = projcsp(xstar,yhat)
         theta = fitcsp$betahat
         ehat = y - cbind(rep(1,n),x)%*%theta
#		tauhat=fitw$tmp2$tau	#wilc scale est
#		taus=fitw$tmp2$tau1		#scale est in intercept, L1
         
         fitvc = rprmeddis(I,sec,mat,ehat,location,scale, rprpair = rprpair)
         sigmaa2 = fitvc$siga2
         sigmaw2 = fitvc$sigw2
         sigmae2 = fitvc$sigmae2
         rea = fitvc$frei	#a effect
         rew = fitvc$frew	#w effect
         ree = fitvc$free	#epsilon error (not raw)

         chk =  sum((theta - thetaold)^2)/sum(thetaold^2)
         if(chk < eps){iflag2 = 1}
     }
         
     list(theta=theta,ehat=ehat,sigmaa2=sigmaa2,sigmaw2=sigmaw2,sigmae2=sigmae2,rea=rea,rew=rew,ree=ree,iflag2=iflag2)
}
