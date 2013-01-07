sigymake = function(I,sec,mat,siga2,sigw2,sige2){
  #Generating three-level nested data data in simulation
  #Used in rlme?
  
     ss = apply(mat,1,sum)
     n = sum(ss)
     iflag = 0
     
     for(i in 1:I){
         for(j in 1:sec[i]){
                mattmp2 = jmake(mat[i,j])
                if(j == 1){
                    secpart = mattmp2
                } else {
                    secpart = bdmake(secpart,mattmp2)
                }
          }
          mattmp = siga2*jmake(ss[i]) + sigw2*secpart + sige2*diag(rep(1,ss[i]))
          see = eigen(mattmp,symmetric=T)$values
          if(min(see) < 0){iflag = 1}
          mattmp3 = sigma12i(mattmp)
          if(i == 1){
               schpart = mattmp
               schpart12i = mattmp3
          } else {
               schpart = bdmake(schpart,mattmp)
               schpart12i = bdmake(schpart12i,mattmp3)
          }
     }
     list(sigy2=schpart,sigy12i=schpart12i,iflag=iflag)
}
