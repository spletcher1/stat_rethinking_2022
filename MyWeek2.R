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


curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )

sample_mu <- rnorm( 1e4 , 178 , 20 )
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h)

flist<-alist(
  height~dnorm(mu,sigma),
  mu~dnorm(178,20),
  sigma~dunif(0,15)
)

answer<-quap(flist,data=d2)


### Now to the linear model

## Priors 
## Maybe assume normal for alpha and beta with mean near the linear model value
## alpha Normal(114,20)
## beta  log Normal(0,1)
## sigma uniform(0,50)

mean.weight<-mean(d2$weight)

flist2<-alist(
  height~dnorm(mu,sigma),
  mu~alpha+beta*(weight-mean.weight),
  alpha~dnorm(114,20),
  beta~dlnorm(0,1),
  sigma~dunif(0,50)
)

answer2<-quap(flist2,data=d2)

plot( height ~ weight , data=d2 , col=rangi2 )
post <- extract.samples( answer2)
a_map <- mean(post$alpha)
b_map <- mean(post$beta)
for(i in 1:nrow(post)){
  a_map<-post$alpha[i]
  b_map<-post$beta[i]
  curve( a_map + b_map*(x - mean.weight) , add=TRUE )
}

sim.height <- sim(answer2 , data=list(weight=c(140-mean.weight,160-mean.weight,175-mean.weight) ))
str(sim.height)

precis(data.frame(sim.height))
