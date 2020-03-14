######################### SCRIPT 0: SUPPORTING FUNCTIONS #########################

runModel = function(indexu,indexk,NITERATION)
{
  dataset=list(
    N=dim(data0)[1],
    Kk=length(indexk),
    nk=data0[,indexk],
    Ku=length(indexu),
    nu=data0[,indexu],
    Sk=knownpop$real_data[indexk],
    Su=rep(NA,length(indexu)))
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model1),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=NITERATION,progress.bar="text",thin=10)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=10)
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  return(results)
}


runModel_te = function(indexu,indexk,NITERATION,x)
{
  dataset=list(
    N=dim(data0)[1],
    Kk=length(indexk),
    nk=data0[,indexk],
    Ku=length(indexu),
    nu=data0[,indexu],
    Sk=knownpop$real_data[indexk],
    Su=rep(NA,length(indexu)),
    x=x)
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model2),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda",'beta','sigma',"tau","Su"),n.iter=NITERATION,progress.bar="text",thin=10)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=10)
  
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  
  return(results)
}


runModel_te_be = function(indexu,indexk,NITERATION,x)
{
  dataset=list(
    N=dim(data0)[1],
    Kk=length(indexk),
    nk=data0[,indexk],
    Ku=length(indexu),
    nu=data0[,indexu],
    Sk=knownpop$real_data[indexk],
    Su=rep(NA,length(indexu)),
    x=x,
    age=demo$age-40.5,
    sex=demo$sex-0.49,
    malay=demo$malay-0.15,
    indian=demo$indian-0.07)
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(model3),data=dataset,inits=initialisation,n.chains=2)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda",'beta','sigma','sigmab',"tau","Su","b1u","b2u","b3u",'b4u',"b1k","b2k",'b3k','b4k'),n.iter=NITERATION,progress.bar="text",thin=10)
  dicsamples = dic.samples(jagmod,type = "pD",n.iter=20000,thin=10)
  
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior,dicsamples=dicsamples)
  return(results)
}


### Models

model1 = 'model {

for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}

for(k in 1:Ku)

{
  
  Su[k]~dunif(0,2500000)
  
}

for(k in 1:Kk)

{
  
  Sk[k]~dunif(0,2500000)
  
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

}
'

model2 ='model {

for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*exp(beta[k]*x[i,k])*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}

for(k in 1:Ku)

{
  
  Su[k]~dunif(0,2500000)
  
}

for(k in 1:Kk)

{
  
  Sk[k]~dunif(0,2500000)
  
}

for(k in 1:Ku)

{
  
  beta[k]~dnorm(0,(1/sigma)^2)
  
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

sigma ~ dgamma(1,0.01)
}
'

model3 = 'model {
  
for(i in 1:N)

{
  
  for(k in 1:Ku)
  
  {
  
  nu[i,k] ~ dpois(lambda*alpha[i]*exp(beta[k]*x[i,k])*exp(b1u[k]*age[i])*exp(b2u[k]*sex[i])*exp(b3u[k]*malay[i])*exp(b4u[k]*indian[i])*Su[k])
  
  }
  
  for(k in 1:Kk)
  
  {
  
  nk[i,k] ~ dpois(lambda*alpha[i]*exp(b1k[k]*age[i])*exp(b2k[k]*sex[i])*exp(b3k[k]*malay[i])*exp(b4k[k]*indian[i])*Sk[k])
  
  }
  
  alpha[i]~dlnorm(0,tau)
  
}

for(k in 1:Ku)

{
  
  Su[k]~dunif(0,2500000)
  
}

for(k in 1:Kk)

{
  
  Sk[k]~dunif(0,2500000)
  b1k[k]~dnorm(0,(1/sigmab)^2)
  b2k[k]~dnorm(0,(1/sigmab)^2)
  b3k[k]~dnorm(0,(1/sigmab)^2)
  b4k[k]~dnorm(0,(1/sigmab)^2)
}

for(k in 1:Ku)

{
  
  beta[k]~dnorm(0,(1/sigma)^2)
  b1u[k]~dnorm(0,(1/sigmab)^2)
  b2u[k]~dnorm(0,(1/sigmab)^2)
  b3u[k]~dnorm(0,(1/sigmab)^2)
  b4u[k]~dnorm(0,(1/sigmab)^2)
}

lambda ~ dunif(0,10)

tau ~ dunif(0,10)

sigma ~ dgamma(1,0.01)
sigmab ~ dgamma(1,0.01)
}
'