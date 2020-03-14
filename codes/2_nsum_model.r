######################### SCRIPT 2: MODEL-BUILDING #########################

# This script runs the three main NSUM models
# 1. Basic model
# 2. Transmission error model
# 3. Barrier effect (+ transmission error) model

# dependencies: JAGS, R packages (coda, rjags)

# Load required packages
library(rjags)

# load the functions from the supporting_function.r script 
source('codes/supporting_functions.r')

# read in data
knownpop = read.csv('data/knownpop.csv',as.is = TRUE)
nsum = read.csv('data/nsum_sg.csv',as.is = TRUE)
validation = read.csv('data/validation.csv',as.is = TRUE)

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


data0 = nsum[,popindex]
te = nsum[,c('socialaccept_mcfsw','socialaccept_msm','socialaccept_fsw','socialaccept_pwid')]
demo = nsum[,c('age','sex','ethnicity_x')]
demo$age[is.na(demo$age)] = mean(demo$age,na.rm = TRUE)
demo$malay = 1*(demo$ethnicity_x %in% "Malay")
demo$indian = 1*(demo$ethnicity_x %in% "Indian")

validation$known_population[validation$rank<11]
indexk = validation$index[validation$rank<11]
indexu = c(7,16,20,24)

# The models were fit using a Markov chain Monte Carlo algorithm with 50 000 iterations 
# with a burn-in of 5 000, storing 1 out of 10 iterations. 
# Convergence was assessed visually with trace plots. The data analyses and visualisations 
# were performed in R18 and the model building was done in JAGS. 
# Deviance Information Criterion (DIC) was used to compare the models. 


set.seed(666)
results = runModel(indexu = c(7,16,20,24),indexk = validation$index[validation$rank<11],NITERATION=50000)
set.seed(666)
results_te = runModel_te(indexu = c(7,16,20,24),indexk = validation$index[validation$rank<11],NITERATION=50000,x = te)
set.seed(666)
results_te_be = runModel_te_be(indexu = c(7,16,20,24),indexk = validation$index[validation$rank<11],NITERATION=500,x = te)


rs = results
rs = results_te
rs = results_te_be
G = geweke.diag(rs$posterior) # Note: you may want to run these diagnostics on the other results too: this just does it for the final model
H = heidel.diag(rs$posterior) 
H
if(FALSE)
{
  # Note: these lines create plots to assess convergence
  #       However, if you are using RStudio and your plots panel is too small, the GUI may stop the code here
  #       So you should manually run this if you want trace plots, by switching the if statement to TRUE, or 
  #       ...just by running these lines manually
  hist(G[[1]][[1]])
  
  plot(rs$posterior[[1]][,'lambda'])
  
  plot(rs$posterior[[1]][,c('beta[1]','beta[2]','beta[3]','beta[4]')])
  
  plot(rs$posterior[[1]][,c('Su[1]','Su[2]','Su[3]','Su[4]')])
  
  plot(rs$posterior[[1]][,c('bu1','bu2','bu3','bu4')])
}


posterior = list()

posterior$basic = as.matrix(results$posterior) 
posterior$transmission = as.matrix(results_te$posterior) 
posterior$barrier = as.matrix(results_te_be$posterior) 


size_posterior = list()

size_posterior$basic = posterior$basic[,grep('Su',colnames(posterior$basic))]
size_posterior$transmission = posterior$transmission[,grep('Su',colnames(posterior$transmission))]
size_posterior$barrier = posterior$barrier[,grep('Su',colnames(posterior$barrier))]

alpha = list()

alpha$basic = posterior$basic[,grep('alpha',colnames(posterior$basic))]
alpha$transmission = posterior$transmission[,grep('alpha',colnames(posterior$transmission))]
alpha$barrier = posterior$barrier[,grep('alpha',colnames(posterior$barrier))]

lambda = list()

lambda$basic = posterior$basic[,grep('lambda',colnames(posterior$basic))]
lambda$transmission = posterior$transmission[,grep('lambda',colnames(posterior$transmission))]
lambda$barrier = posterior$barrier[,grep('lambda',colnames(posterior$barrier))]

size_estimates = list()

size_estimates$basic = (signif(apply(posterior$basic[,grep('Su',colnames(posterior$basic))],2,function(x)quantile(x,probs = c(0.5,0.025,0.975))),2))
size_estimates$transmission = (signif(apply(posterior$transmission[,grep('Su',colnames(posterior$transmission))],2,function(x)quantile(x,probs = c(0.5,0.025,0.975))),2))
size_estimates$barrier = (signif(apply(posterior$barrier[,grep('Su',colnames(posterior$barrier))],2,function(x)quantile(x,probs = c(0.5,0.025,0.975))),2))

colnames(size_estimates[[1]]) = c('MCFSW','MSM','FSW','IVDU')
colnames(size_estimates[[2]]) = c('MCFSW','MSM','FSW','IVDU')
colnames(size_estimates[[3]]) = c('MCFSW','MSM','FSW','IVDU')

sizeestimate_summary = rbind(cbind(model = 1,size_estimates[[1]]),
                             cbind(model = 2,size_estimates[[2]]),
                             cbind(model = 3,size_estimates[[3]]))

write.csv(x = sizeestimate_summary,file = 'results/sizeestimates.csv')


postmat_adjusted = as.matrix(results_te$posterior) 
sizehiddenmat_adjusted = postmat_adjusted[,grep('Su',colnames(postmat_adjusted))]
lambdas_adjusted = postmat_adjusted[,grep('lambda',colnames(postmat_adjusted))]
alphas_adjusted = postmat_adjusted[,grep('alpha',colnames(postmat_adjusted))]
betas = postmat_adjusted[,grep('beta',colnames(postmat_adjusted))]
#transmission error
betasummary = t(signif(apply(betas,2,function(x)quantile(x,probs = c(0.5,0.025,0.975))),2))
rownames(betasummary) = c('MCFSW','MSM','FSW','IVDU')
write.csv(x = betasummary,file = 'results/transmissionerror.csv')
#size estimates
sizeestimates = t(signif(apply(sizehiddenmat_adjusted,2,function(x)quantile(x,probs = c(0.5,0.025,0.975))),2))
rownames(sizeestimates) = c('MCFSW','MSM','FSW','IVDU')
write.csv(x = sizeestimates,file = 'results/sizeestimates.csv')

# network size
mod = 3
alphas_mean = apply(alpha[[mod]],2,function(x)quantile(x,probs = c(0.5)))
alphas_cil = apply(alpha[[mod]],2,function(x)quantile(x,probs = c(0.025)))
alphas_ciu = apply(alpha[[mod]],2,function(x)quantile(x,probs = c(0.975)))
networksize_crude = mean(lambda[[mod]])*5E6
networksize_crude
networksize = mean(lambda[[mod]])*mean(alphas_mean)*5E6
networksize_cil = mean(lambda[[mod]])*mean(alphas_cil)*5E6
networksize_ciu = mean(lambda[[mod]])*mean(alphas_ciu)*5E6
network  = data.frame(networksize=networksize,networksize_cil=networksize_cil,networksize_ciu=networksize_ciu)
network
write.csv(round(network,1),'results/network.csv',row.names = FALSE)

