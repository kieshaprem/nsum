######################### SCRIPT 1: KNOWN POPUALTION VALIDATION #########################

# dependencies: JAGS, R packages (coda, rjags)

# Load required packages
library(rjags)

# read in data
knownpop = read.csv('data/knownpop.csv',as.is = TRUE)
nsum = read.csv('data/nsum_sg.csv',as.is = TRUE)

names_pop = c('Medical doctors','Primary school teachers','Full-time NS men',
            'Licensed taxi drivers','Hawker stalls owners',
            'Licensed property agents','Male clients of sex workers',
            'Women who had a baby in 2016','NUS graduates 2016',
            'Stroke 2016','Heart attack 2016','Bought an HDB in 2016',
            'Driving license 2016','PSLE 2016','NDP 2016','MSM',
            'O-Levels 2016','Couples who got married in 2016',
            'Dengue 2016','Female sex workers','Men >70 years old',
            'Women >70 years old','Single men >50 years old','IVDU')
hiddenpop = c(0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,0,1,0,0,0,1)
popindex = grep(pattern = '_total',x = names(nsum))


modelstring='model {\n
for(i in 1:N)\n
{\n
for(k in 1:Ku)\n
{\n
nu[i,k] ~ dpois(lambda*alpha[i]*Su[k])\n
}\n
for(k in 1:Kk)\n
{\n
nk[i,k] ~ dpois(lambda*alpha[i]*Sk[k])\n
}\n
alpha[i]~dlnorm(0,tau)\n
}\n
for(k in 1:Ku)\n
{\n
Su[k]~dunif(0,2500000)\n
}\n
for(k in 1:Kk)\n
{\n
Sk[k]~dunif(0,2500000)\n
}\n
lambda ~ dunif(0,10)\n
tau ~ dunif(0,10)\n
}\n'



data0 = nsum[,popindex]

runModel = function(indexu)
{
  indexk = which(!(1:24 %in% indexu))
  dataset=list(
    N=dim(data0)[1],
    Kk=dim(data0)[2]-length(indexu),
    nk=data0[,indexk],
    Ku=length(indexu),
    nu=data0[,indexu],
    Sk=knownpop$real_data[indexk],
    Su=rep(NA,length(indexu)))
  
  initialisation=list(lambda=0.1)
  jagmod=jags.model(textConnection(modelstring),data=dataset,inits=initialisation,n.chains=1)
  update(jagmod, n.iter=5000, progress.bar="text")
  posterior = coda.samples(jagmod, c("alpha","lambda","tau","Su"),n.iter=10000,progress.bar="text",thin=10)
  
  results = list(indexk=indexk,indexu=indexu,dataset = dataset,posterior=posterior)
  return(results)
}


indexu = c(3,7,16,20,24)

rs = runModel(indexu = c(3,7,16,20,24))
rs1 = runModel(indexu = c(1,3,7,16,20,24))
rs2 = runModel(indexu = c(2,3,7,16,20,24))
rs4 = runModel(indexu = c(3,4,7,16,20,24))
rs5 = runModel(indexu = c(3,5,7,16,20,24))
rs6 = runModel(indexu = c(3,6,7,16,20,24))
rs8 = runModel(indexu = c(3,7,8,16,20,24))
rs9 = runModel(indexu = c(3,7,9,16,20,24))
rs10 = runModel(indexu = c(3,7,10,16,20,24))
rs11 = runModel(indexu = c(3,7,11,16,20,24))
rs12 = runModel(indexu = c(3,7,12,16,20,24))
rs13 = runModel(indexu = c(3,7,13,16,20,24))
rs14 = runModel(indexu = c(3,7,14,16,20,24))
rs15 = runModel(indexu = c(3,7,15,16,20,24))
rs17 = runModel(indexu = c(3,7,16,17,20,24))
rs18 = runModel(indexu = c(3,7,16,18,20,24))
rs19 = runModel(indexu = c(3,7,16,19,20,24))
rs21 = runModel(indexu = c(3,7,16,20,21,24))
rs22 = runModel(indexu = c(3,7,16,20,22,24))
rs23 = runModel(indexu = c(3,7,16,20,23,24))

rslist = list(rs1=rs1,rs2=rs2,rs4=rs4,rs4=rs4,rs5=rs5,rs6=rs6,rs8=rs8,
              rs9=rs9,rs10=rs10,rs11=rs11,rs12=rs12,rs13=rs13,rs14=rs14,
              rs15=rs15,rs17=rs17,rs18=rs18,rs19=rs19,rs21=rs21,rs22=rs22,rs23=rs23)

getModelResults = function(RESULTS,VALIDATE)
{
  postmat = as.matrix(RESULTS$posterior) 
  su = postmat[,paste0('Su[',VALIDATE,']')]
  return(su)
}

quantile(getModelResults(RESULTS = rs1,VALIDATE = 1),probs = c(0.025,0.5,0.975))
knownpopindex = which(!(is.na(knownpop$real_data)))
knownpop$known_population[!(is.na(knownpop$real_data))]
validation = data.frame(index = knownpopindex,known_population = knownpop$known_population[!(is.na(knownpop$real_data))],
                        real = knownpop$real_data[!(is.na(knownpop$real_data))], 
                        median=NA,cril=NA,criu=NA)

for(i in 1:length(knownpopindex))
{
  LI = knownpopindex[i]
  validation$median[validation$index %in% knownpopindex[i]] = quantile(getModelResults(RESULTS = rslist[[paste0('rs',LI)]],VALIDATE = which(rslist[[paste0('rs',LI)]]$indexu %in% knownpopindex[i])),probs = c(0.5))
  validation$cril[validation$index %in% knownpopindex[i]] = quantile(getModelResults(RESULTS = rslist[[paste0('rs',LI)]],VALIDATE = which(rslist[[paste0('rs',LI)]]$indexu %in% knownpopindex[i])),probs = c(0.025))
  validation$criu[validation$index %in% knownpopindex[i]] = quantile(getModelResults(RESULTS = rslist[[paste0('rs',LI)]],VALIDATE = which(rslist[[paste0('rs',LI)]]$indexu %in% knownpopindex[i])),probs = c(0.975))
}



discrepancy = (log(validation$median/validation$real))
discrepancy_adjusted = discrepancy
discrepancy_adjusted[validation$known_population %in% "Men >70 years old"] = NA
discrepancy_rank = rank(abs(discrepancy_adjusted),na.last = TRUE)
validation$ratio = (validation$median/validation$real)
validation$discrepancy = log(validation$median/validation$real)
validation$rank = rank(abs(discrepancy_adjusted))
validation

write.csv(validation,file = 'data/validation.csv',row.names = FALSE)




