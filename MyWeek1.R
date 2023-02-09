library(rethinking)

######################
# Question #1
######################

probs<-seq(0,1,length.out=500)

likelihoods<-dbinom(4,15,probs)

plot(probs,likelihoods)

prior<-rep(1/500,500)

posterior<-(likelihoods*prior)/(sum(likelihoods*prior))

plot(probs,posterior,xlab="Water Fraction of Earth")

sum(posterior)


######################
# Question #2
######################

probs<-seq(0,1,length.out=50000)

likelihoods<-dbinom(4,6,probs)

plot(probs,likelihoods)

tmp<-rep(0,25000)
prior<-rep(1,25000)
prior<-c(tmp,prior)

posterior<-(likelihoods*prior)/(sum(likelihoods*prior))

plot(probs,posterior,xlab="Water Fraction of Earth")

sum(posterior)

######################
# Question #3
######################

samples<-sample(probs,10000, prob=posterior,replace=TRUE)
PI(samples)
HPDI(samples)


######################
# Optional
######################

sample.size<-1000
num.repeats<-100
prop.water<-0.7

raw.water.results<-rbinom(num.repeats,sample.size,prop.water)
recorded.water.results<-rep(NA,length(raw.water.results))
for(i in 1:length(raw.water.results)){
  recorded.water.results[i]<-rbinom(1,raw.water.results[i],0.8)  
}

estimated.perc.water<-recorded.water.results/sample.size

raw.perc.water<-raw.water.results/sample.size

dens(estimated.perc.water,xlab="Estimated",xlim=c(0,1))
dens(raw.perc.water,xlab="Raw",xlim=c(0,1))
  

mean(estimated.perc.water)
  


