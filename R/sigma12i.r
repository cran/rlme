sigma12i<-function(sigma){
temp<-eigen(sigma,symmetric=T)
sigma12i<-temp$vectors%*%diag(1/temp$values^.5)%*%t(temp$vectors)
sigma12i
}
