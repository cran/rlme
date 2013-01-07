sigma12<-function(sigma){
temp<-eigen(sigma)
sigma12<-temp$vectors%*%diag(temp$values^.5)%*%t(temp$vectors)
sigma12
}
