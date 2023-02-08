rm(list=ls())
library(rethinking)
data(Howell1)

d<-Howell1
d2 <- d[ d$age >= 18 , ]

dens(d2$height)

tmp<-function(i,data,postvals){
  tmp2<-dnorm( data$height , postvals$mu[i] , postvals$sigma[i] , log=TRUE )
  sum(tmp2) 
}



mu.list <- seq( from=150, to=160 , length.out=100 ) 
sigma.list <- seq( from=7 , to=9 , length.out=100 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )

post$LL <- sapply(1:nrow(post),tmp,d2,post)

post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )