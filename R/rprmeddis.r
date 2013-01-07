rpr = function(ehat, school, section, rprpair = "hl-disp") {
  I = length(unique(factor(section)))   #number of uppest clusters
  sec = as.vector(sec_vec(school, section))  #number of sections in each school  
  
  mat = mat_vec(school, section)  #sample size for each section in school

  rprpair = tolower(rprpair)
  location = scale = 2
  if(rprpair == "med-mad") {
    location = scale = 1
  }
  
  return(rprmeddis(I, sec, mat, ehat, location, scale))
}

rprmeddis = function(I,sec,mat,ehat,location,scale, rprpair = "hl-disp"){

#this fn gets location and scale est 
#location: 	1-median, 	2-wilc or hl est, 3-huber, 	4-mean
#scale:	1-mad,  	2-disp,  		3-tau,  	4-sd, 5-huber est,
#In the rlme pack, only options 1 and 2 are used in rprpair="med-mad" and "hl-disp"
  
rprpair = tolower(rprpair)
if(rprpair == 'hl-disp' ) {
  location = scale = 2
}
if(rprpair == 'med-mad') {
  location = scale = 1
}
  
####location
if(location==1){  
     matre = matrix(c(0,0),ncol=2)

     uhatij = 0*mat
     ni = apply(mat,1,sum)
     reij = 0*mat
     rei = rep(0,I)
     rew = rep(0,sum(sec))
     frei = rep(0,I)
     uhati = rep(0,I)
     ehat2 = ehat
     ehat3 = 0*ehat
     ic = 0
     id = 0
     ia = 0
     for(i in 1:I){
         for(j in 1:sec[i]){
                ehattmp = ehat[(ic+1):(ic+mat[i,j])]
                ic = ic+mat[i,j]
                uhatij[i,j] = median(ehattmp)
#		    uhatij[i,j] = onesampwil(ehattmp,maktable=F,plotb=F)$est   #hl location 

         }
         uhati[i] = median(uhatij[i,1:sec[i]])
# 	   uhati[i] = onesampwil(uhatij[i,1:sec[i]],maktable=F,plotb=F)$est
         for(j in 1:sec[i]){
              reij[i,j] = uhatij[i,j] - uhati[i]
         }
         for(j in 1:sec[i]){
                ehat2[(id+1):(id+mat[i,j])] = ehat[(id+1):(id+mat[i,j])] - reij[i,j]
                id = id+mat[i,j]
         }
        rei[i] = median(ehat2[(ia+1):(ia+ni[i])]) 
#   	 rei[i]= onesampwil(ehat2[(ia+1):(ia+ni[i])],maktable=F,plotb=F)$est    
         ia = ia + ni[i]
    }


    frei = rei - median(rei)
#     frei = rei - onesampwil(rei,maktable=F,plotb=F)$est 
    ib = 0
    for(i in 1:I){
        for(j in 1:sec[i]){
            ib = ib + 1
            rew[ib] = reij[i,j]
        }
    }
    frew = rew - median(rew) 
#     frew = rew - onesampwil(rew,maktable=F,plotb=F)$est 
}

if(location==2){
     matre = matrix(c(0,0),ncol=2)

     uhatij = 0*mat
     ni = apply(mat,1,sum)
     reij = 0*mat
     rei = rep(0,I)
     rew = rep(0,sum(sec))
     frei = rep(0,I)
     uhati = rep(0,I)
     ehat2 = ehat
     ehat3 = 0*ehat
     ic = 0
     id = 0
     ia = 0
     for(i in 1:I){
         for(j in 1:sec[i]){
                ehattmp = ehat[(ic+1):(ic+mat[i,j])]
                ic = ic+mat[i,j]
#               uhatij[i,j] = median(ehattmp)

			#when 1 obs in section in school, need the following
			 if(length(ehattmp)==1){
				uhatij[i,j] =  ehattmp   #
				} else {
	  		 uhatij[i,j] = onesampwil(ehattmp,maktable=F,plotb=F)$est   #hl location 
	  		 #uhatij[i,j] = as.vector(wilcox.test(ehattmp,conf.int=T)$estimate)
				}
         }
#         uhati[i] = median(uhatij[i,1:sec[i]])
			 if(length(uhatij[i,1:sec[i]])==1){
				uhati[i] =  uhatij[i,1:sec[i]] # 
				} else {
		 	  uhati[i] = onesampwil(uhatij[i,1:sec[i]],maktable=F,plotb=F)$est
#				uhati[i] = as.vector(wilcox.test(uhatij[i,1:sec[i]],conf.int=T)$estimate)
			 }


         for(j in 1:sec[i]){
              reij[i,j] = uhatij[i,j] - uhati[i]
         }
         for(j in 1:sec[i]){
                ehat2[(id+1):(id+mat[i,j])] = ehat[(id+1):(id+mat[i,j])] - reij[i,j]
                id = id+mat[i,j]
         }
#        rei[i] = median(ehat2[(ia+1):(ia+ni[i])]) 

			 if(length(ehat2[(ia+1):(ia+ni[i])])==1){
				rei[i] =  ehat2[(ia+1):(ia+ni[i])] # 
				} else {
		    	rei[i]= onesampwil(ehat2[(ia+1):(ia+ni[i])],maktable=F,plotb=F)$est  
#				  rei[i]= as.vector(wilcox.test(ehat2[(ia+1):(ia+ni[i])],conf.int=T)$estimate)
			 }

         ia = ia + ni[i]
    }
#?Error in xpairs[, 1] : incorrect number of dimensions

#    frei = rei - median(rei)
    frei = rei - onesampwil(rei,maktable=F,plotb=F)$est 
#     frei = rei - as.vector(wilcox.test(rei,conf.int=T)$estimate)
    ib = 0
    for(i in 1:I){
        for(j in 1:sec[i]){
            ib = ib + 1
            rew[ib] = reij[i,j]
        }
    }
#    frew = rew - median(rew) 

			 if(length(rew)==1){
			     frew = rew - rew 
				} else {
			    frew = rew - onesampwil(rew,maktable=F,plotb=F)$est 
#				  frew = rew - as.vector(wilcox.test(rew,conf.int=T)$estimate)
			 }


}

if(location==3){  
     matre = matrix(c(0,0),ncol=2)

     uhatij = 0*mat
     ni = apply(mat,1,sum)
     reij = 0*mat
     rei = rep(0,I)
     rew = rep(0,sum(sec))
     frei = rep(0,I)
     uhati = rep(0,I)
     ehat2 = ehat
     ehat3 = 0*ehat
     ic = 0
     id = 0
     ia = 0
     for(i in 1:I){
         for(j in 1:sec[i]){
                ehattmp = ehat[(ic+1):(ic+mat[i,j])]
                ic = ic+mat[i,j]
               uhatij[i,j] = huber(ehattmp)$mu
#		   uhatij[i,j] = onesampwil(ehattmp,maktable=F,plotb=F)$est   #hl location 

         }
         uhati[i] = huber(uhatij[i,1:sec[i]])$mu
# 	   uhati[i] = onesampwil(uhatij[i,1:sec[i]],maktable=F,plotb=F)$est
         for(j in 1:sec[i]){
              reij[i,j] = uhatij[i,j] - uhati[i]
         }
         for(j in 1:sec[i]){
                ehat2[(id+1):(id+mat[i,j])] = ehat[(id+1):(id+mat[i,j])] - reij[i,j]
                id = id+mat[i,j]
         }
        rei[i] = huber(ehat2[(ia+1):(ia+ni[i])])$mu 
#   	  rei[i]= onesampwil(ehat2[(ia+1):(ia+ni[i])],maktable=F,plotb=F)$est    
         ia = ia + ni[i]
    }


      frei = rei - huber(rei)$mu
#     frei = rei - onesampwil(rei,maktable=F,plotb=F)$est 
    ib = 0
    for(i in 1:I){
        for(j in 1:sec[i]){
            ib = ib + 1
            rew[ib] = reij[i,j]
        }
    }
    frew = rew - huber(rew)$mu 
#   frew = rew - onesampwil(rew,maktable=F,plotb=F)$est 
}

if(location==4){  
     matre = matrix(c(0,0),ncol=2)

     uhatij = 0*mat
     ni = apply(mat,1,sum)
     reij = 0*mat
     rei = rep(0,I)
     rew = rep(0,sum(sec))
     frei = rep(0,I)
     uhati = rep(0,I)
     ehat2 = ehat
     ehat3 = 0*ehat
     ic = 0
     id = 0
     ia = 0
     for(i in 1:I){
         for(j in 1:sec[i]){
                ehattmp = ehat[(ic+1):(ic+mat[i,j])]
                ic = ic+mat[i,j]
               uhatij[i,j] = mean(ehattmp)
#		uhatij[i,j] = onesampwil(ehattmp,maktable=F,plotb=F)$est   #hl location 

         }
         uhati[i] = mean(uhatij[i,1:sec[i]])
# 	  uhati[i] = onesampwil(uhatij[i,1:sec[i]],maktable=F,plotb=F)$est
         for(j in 1:sec[i]){
              reij[i,j] = uhatij[i,j] - uhati[i]
         }
         for(j in 1:sec[i]){
                ehat2[(id+1):(id+mat[i,j])] = ehat[(id+1):(id+mat[i,j])] - reij[i,j]
                id = id+mat[i,j]
         }
        rei[i] = mean(ehat2[(ia+1):(ia+ni[i])]) 
#   	 rei[i]= onesampwil(ehat2[(ia+1):(ia+ni[i])],maktable=F,plotb=F)$est    
         ia = ia + ni[i]
    }


    frei = rei - mean(rei)
#     frei = rei - onesampwil(rei,maktable=F,plotb=F)$est 
    ib = 0
    for(i in 1:I){
        for(j in 1:sec[i]){
            ib = ib + 1
            rew[ib] = reij[i,j]
        }
    }
    frew = rew - mean(rew) 
#     frew = rew - onesampwil(rew,maktable=F,plotb=F)$est 
}


####now, scale est
if(scale==1){
     siga2 = mad(frei)^2	#mad
#    siga2 =onesampwil(frei,maktable=F,plotb=F)$tau^2	#tau
#    siga2 = dispvar(frei)^2	#disp
#    siga2 = sd(frei)^2		#var

     sigw2 = mad(frew)^2	#mad
#    sigw2 =onesampwil(frew,maktable=F,plotb=F)$tau^2	#tau
#    sigw2 = dispvar(frew)^2	#disp
#    sigw2 = sd(frew)^2		#var

    ic = 0
    id = 0
    for(i in 1:I){
        re1 = frei[i]
        for(j in 1:sec[i]){
           id = id + 1
           re2 = frew[id]
           for(k in 1:mat[i,j]){
              ic = ic + 1
              ehat3[ic] = ehat[ic] - re1 - re2
           }
        }
    }
     sigmae2 = mad(ehat3)^2 #mad
#    sigmae2 = onesampwil(ehat3,maktable=F,plotb=F)$tau^2	#tau
#    sigmae2 = dispvar(ehat3)^2	#disp
#    sigmae2 = sd(ehat3)^2  #var
}

if(scale==2){

     siga2 = dispvar(frei)^2	#.88*disp for r.eff school. coeffs are correction
     sigw2 = dispvar(frew)^2	#1.12*disp for r.eff sect in sch 

    ic = 0
    id = 0
    for(i in 1:I){
        re1 = frei[i]
        for(j in 1:sec[i]){
           id = id + 1
           re2 = frew[id]
           for(k in 1:mat[i,j]){
              ic = ic + 1
              ehat3[ic] = ehat[ic] - re1 - re2
           }
        }
    }

     sigmae2 = dispvar(ehat3)^2	#1.12*disp 
}

if(scale==3){

    siga2 =onesampwil(frei,maktable=F,plotb=F)$tau^2	#tau
    sigw2 =onesampwil(frew,maktable=F,plotb=F)$tau^2	#tau

    ic = 0
    id = 0
    for(i in 1:I){
        re1 = frei[i]
        for(j in 1:sec[i]){
           id = id + 1
           re2 = frew[id]
           for(k in 1:mat[i,j]){
              ic = ic + 1
              ehat3[ic] = ehat[ic] - re1 - re2
           }
        }
    }
    sigmae2 = onesampwil(ehat3,maktable=F,plotb=F)$tau^2	#tau

}


if(scale==4){

    siga2 = sd(frei)^2		#var
    sigw2 = sd(frew)^2		#var

    ic = 0
    id = 0
    for(i in 1:I){
        re1 = frei[i]
        for(j in 1:sec[i]){
           id = id + 1
           re2 = frew[id]
           for(k in 1:mat[i,j]){
              ic = ic + 1
              ehat3[ic] = ehat[ic] - re1 - re2
           }
        }
    }

    sigmae2 = sd(ehat3)^2  #var
}

if(scale==5){
     siga2 = huber(frei)$s^2	#huber
     sigw2 = huber(frew)$s^2	#huber

    ic = 0
    id = 0
    for(i in 1:I){
        re1 = frei[i]
        for(j in 1:sec[i]){
           id = id + 1
           re2 = frew[id]
           for(k in 1:mat[i,j]){
              ic = ic + 1
              ehat3[ic] = ehat[ic] - re1 - re2
           }
        }
    }
     sigmae2 = huber(ehat3)$s^2 #huber
}


sigmae2=sigmae2

list(frei=frei,frew=frew,free=ehat3,siga2=siga2,sigw2=sigw2,sigmae2=sigmae2)

}
