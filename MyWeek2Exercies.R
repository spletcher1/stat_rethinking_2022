data<-read.csv("Flies.csv")
data$Genotype<-factor(data$Genotype)
data$Exposure<-factor(data$Exposure)

control<-subset(data,data$Genotype=="Dsk81_W")

deaths<-control[control$Event==1,]

deaths.expanded<-deaths[rep(1:nrow(deaths),deaths$Nobs),]
deaths.expanded$Nobs<-1

summary(aov(Age~Exposure,data=deaths.expanded))

write.csv(deaths.expanded,file="Data.csv")


###########################
## 
library(rethinking)

data<-read.csv("Data.csv")
data$Exposure<-factor(data$Exposure)

precis(data)
## Assume Age is normally distributed with mean mu and sd=sigma
## Mean prior is normal.  SD prior is uniform
data.w<-subset(data,data$Exposure=="W_Tra")
data.ok<-subset(data,data$Exposure=="OK_Tra")

##################################################S
## ## W_Tra first using 90% of the data.
mu.list<-seq(from=986, to=1824, length.out=100)
sigma.list<-seq(from=100, to=500, length.out=100)
post<-expand.grid(mu=mu.list,sigma=sigma.list)

ll.function<-function(i,post,data){
  sum(dnorm(data$Age , post$mu[i] , post$sigma[i] , log=TRUE ) )
}

apply.prior.mu<-function(mu.vector,mean,sd){
  dnorm(mu.vector,mean,sd,log=TRUE)
}

apply.prior.sigma<-function(sigma.vector,uniStart,uniEnd){
  dunif(sigma.vector,uniStart,uniEnd,log=TRUE)
}

post$LL<-sapply(1:nrow(post),ll.function,post,data.w)

post$prod<-post$LL + apply.prior.mu(post$mu,1000,200) + apply.prior.sigma(post$sigma,100,500)

post$prob<-exp(post$prod-max(post$prod))

##Visualize the posterior
contour_xyz(post$mu,post$sigma,post$prob,xlab="Mu",ylab="Sigma")
image_xyz(post$mu,post$sigma,post$prob,xlab="Mu",ylab="Sigma")


sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
samples <- data.frame(post$mu[ sample.rows ],post$sigma[ sample.rows ])

precis(samples)
summary(data.w$Age)

##################################################S
## ## OK_Tra first using 90% of the data.
mu.list<-seq(from=986, to=1824, length.out=100)
sigma.list<-seq(from=100, to=500, length.out=100)
post<-expand.grid(mu=mu.list,sigma=sigma.list)

post$LL<-sapply(1:nrow(post),ll.function,post,data.ok)
post$prod<-post$LL + apply.prior.mu(post$mu,1000,200) + apply.prior.sigma(post$sigma,100,500)
post$prob<-exp(post$prod-max(post$prod))

##Visualize the posterior
contour_xyz(post$mu,post$sigma,post$prob)
image_xyz(post$mu,post$sigma,post$prob)


sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
samples <- data.frame(post$mu[ sample.rows ],post$sigma[ sample.rows ])

precis(samples)
summary(data.ok$Age)


## How much does the estimate of the mean change with different priors
prior.means<-seq(from=5,to=5000,length.out=100)
results<-rep(NA,100)
for(i in 1:length(prior.means)){
  post$LL<-sapply(1:nrow(post),ll.function,post,data.ok)
  post$prod<-post$LL + apply.prior.mu(post$mu,prior.means[i],200) + apply.prior.sigma(post$sigma,100,500)
  post$prob<-exp(post$prod-max(post$prod))
  sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                         prob=post$prob )
  samples <- data.frame(post$mu[ sample.rows ],post$sigma[ sample.rows ])
  results[i]<-mean(samples[,1])
}

plot(prior.means,results)

post$LL<-sapply(1:nrow(post),ll.function,post,data.ok)
post$prod<-post$LL + apply.prior.mu(post$mu,1000,200) + apply.prior.sigma(post$sigma,100,500)
post$prob<-exp(post$prod-max(post$prod))
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
samples <- data.frame(post$mu[ sample.rows ],post$sigma[ sample.rows ])
dens(samples[,1])
tmp<-samples[,1]


post$LL<-sapply(1:nrow(post),ll.function,post,data.ok)
post$prod<-post$LL + apply.prior.mu(post$mu,5000,200) + apply.prior.sigma(post$sigma,100,500)
post$prob<-exp(post$prod-max(post$prod))
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
samples <- data.frame(post$mu[ sample.rows ],post$sigma[ sample.rows ])
dens(samples[,1],add=TRUE)

### Now lets try a baysian approach to asking whether there are differences between
### the means of the two cohors.

#death is Normal(mu, sigma)
#mu is alpha + beta*exposure (OK_tra = 0, w_Tra=1)
#alpha prior is normal (1500,300)
# beta prior is normal (0,500)
#signam prior is uni(200,500)

## First try grid approach
tmp<-rep(1,nrow(data))
tmp[data$Exposure=="OK_Tra"]<-0
new.data<-data.frame(data,tmp)
names(new.data)<-c(names(data),"Indicator")

alpha.list<-seq(from=986, to=1824, length.out=100)
beta.list<-seq(from=-400,to=400,length.out=100)
sigma.list<-seq(from=100, to=500, length.out=100)
post<-expand.grid(alpha=alpha.list,beta=beta.list,sigma=sigma.list)

ll2.function<-function(i,post,data){
  sum(dnorm(data$Age , post$alpha[i]+(post$beta[i]*data$Indicator) , post$sigma[i] , log=TRUE ) )
}

post$LL<-sapply(1:nrow(post),ll2.function,post,new.data)
post$prod<-post$LL + apply.prior.mu(post$alpha,1500,300) + apply.prior.mu(post$beta,0,500) + apply.prior.sigma(post$sigma,200,500)
post$prob<-exp(post$prod-max(post$prod))

sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       prob=post$prob )
samples <- data.frame(post$alpha[ sample.rows ],post$beta[ sample.rows ],post$sigma[ sample.rows ])
names(samples)<-c("Alpha","Beta","Sigma")
dens(samples$Beta)

PI(samples$Beta)
summary(samples$Beta)
## Probably the difference is greater than 0
sum(samples$Beta>0)/length(samples$Beta)









data<-read.csv("Data.csv")
data$Exposure<-factor(data$Exposure)

tmp<-rep(1,nrow(data))
tmp[data$Exposure=="OK_Tra"]<-0
new.data<-data.frame(data,tmp)
names(new.data)<-c(names(data),"Indicator")

## Let's try fitting the quadratic approximation to the prior.
answer <- quap(
  alist(
    Age ~ dnorm( mu , sigma ) ,
    mu <- alpha+beta*Indicator,
    alpha ~ dnorm( 1500 , 10 ) ,
    beta ~ dnorm( 100 , 100 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=new.data )

precis(answer)
round(vcov(answer),3)

post<-extract.samples(answer,10000)
hist(post$beta)
dens(post$beta)

## Probability that beta>0
sum(post$beta>0)/length(post$beta)
