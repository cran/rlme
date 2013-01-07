dispvar = function(x,score=1){
#this is scale estimate like sd depending on scores
#1-wilc score, 2-sign score (look pg.201-203, McKean&Hettm-2011)
   n = length(x)
   rx = rank(x,ties.method=c("random"))

if(score==1){
   sc = sqrt(12)*((rx/(n+1)) - .5)
   dispvar = sqrt(pi/3)*sum(x*sc)/n	#when wilc scores used. correct sqrt(n/(n-2))
   }

if(score==2){
sc=sign((rx-(n+1)/2))	#for L1
dispvar = sqrt(pi/2)*sum(x*sc)/n	#when sign scores used
   }

dispvar
}

