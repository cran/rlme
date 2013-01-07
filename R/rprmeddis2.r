rprmeddis2 = function(I,sec,mat,ehat,location,scale, rprpair = ""){

#this fn gets location and scale est 
#location: 	1-median (Med), 	2-wilc or hl est (HL)
#scale:	1-mad (MAD),  	2-Disp (Disp)

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
     sch_size=as.vector(apply(mat,1,sum))
     
     ahati = rep(0,I) #this is a vector in 2level
     epshatij=rep(0,sum(sch_size)) #this is epsilon in 2level so ei=ai+epsij
  
#      ehat2 = ehat
#      ehat3 = 0*ehat

     ic = 0
     id = 0
    
#      ia = 0

     for(i in 1:I){
                ehattmp = ehat[(ic+1):(ic+sch_size[i])]
                ic = ic+sch_size[i]
                ahati[i] = median(ehattmp)
         }
         amed = median(ahati)
         ahati=ahati-amed
# 	   uhati[i] = onesampwil(uhatij[i,1:sec[i]],maktable=F,plotb=F)$est
     for(j in 1:I){
        for(k in 1:sch_size[j]){
            epshatij[(id+k)] = ehat[(id+k)] - ahati[j]
                               }
            id = id+sch_size[j]
                  }
     epshatij=epshatij-median(epshatij)
}

if(location==2){  
  
  matre = matrix(c(0,0),ncol=2)
  sch_size=as.vector(apply(mat,1,sum))
  
  ahati = rep(0,I) #this is a vector in 2level
  epshatij=rep(0,sum(sch_size)) #this is epsilon in 2level so ei=ai+epsij
  
  #      ehat2 = ehat
  #      ehat3 = 0*ehat
  
  ic = 0
  id = 0
  
  #      ia = 0
  
  for(i in 1:I){
    ehattmp = ehat[(ic+1):(ic+sch_size[i])]
    ic = ic+sch_size[i]
#   ahati[i] = median(ehattmp)
    
    if(length(ehattmp)==1){
      ahati[i] =  ehattmp   #
      } else {
     ahati[i] = onesampwil(ehattmp,maktable=F,plotb=F)$est   #hl location 
#    ahati[i] =as.vector(wilcox.test(ehattmp,conf.int=T)$estimate)
      }
    }
#  amed = median(ahati)
  amed =  onesampwil(ahati,maktable=F,plotb=F)$est
#  amed =as.vector(wilcox.test(ahati,conf.int=T)$estimate)
  ahati=ahati-amed
  for(j in 1:I){
    for(k in 1:sch_size[j]){
      epshatij[(id+k)] = ehat[(id+k)] - ahati[j]
    }
    id = id+sch_size[j]
  }
# epshatij=epshatij-median(epshatij)
  epshatij=epshatij-onesampwil(epshatij,maktable=F,plotb=F)$est
#  epshatij=epshatij-as.vector(wilcox.test(epshatij,conf.int=T)$estimate)

}


####now, scale est
if(scale==1){
     siga2 = mad(ahati)^2	#mad
     sigmae2 = mad(epshatij)^2 #mad
#    sigmae2 = onesampwil(ehat3,maktable=F,plotb=F)$tau^2	#tau
#    sigmae2 = dispvar(ehat3)^2	#disp
#    sigmae2 = sd(ehat3)^2  #var
}

if(scale==2){

     siga2 = dispvar(ahati)^2	#.88*disp for r.eff school. coeffs are correction
     sigmae2 = dispvar(epshatij)^2	#1.12*disp for r.eff sect in sch 
  
	}

list(frei=ahati,free=epshatij,siga2=siga2,sigmae2=sigmae2)

}

# effect_sch = scale.fit$frei  #school random eff
# effect_sec = scale.fit$frew	#sec random eff
# effect_err = scale.fit$free
