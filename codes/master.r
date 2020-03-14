######################### MASTER SCRIPT #########################
# This R project consists of 4 R scripts:

# 1)	master.r 
# Prior to running the following scripts, 
# please ensure that you have the necessary dependencies: JAGS, R packages (coda, rjags). 

# 2)	supporting_functions.r : contains additional functions, and models. 
source('codes/supporting_functions.r')

# 3)	1_knownpop_validation.r : first validation when selecting the 10 known population (computationally intensive)
source('codes/1_knownpop_validation.r')

# 4)	2_nsum_model.r : 3 models---basic model, tranmission error model, barrier effect model (computationally intensive)
source('codes/2_nsum_model.r')
