# Set the path to the main directory
mainDir <- "path to folder/Data_and_Software"

# Load user functions (make sure to have packages coda and TeachingDemos installed)
source(file.path(mainDir,"CardiopulmNet_functions.r"))

# Load metadata (ID, full name, medical semantic, name of categories)
load(file.path(mainDir,"CardiopulmNet_info.RData"))

# Load data
CardiopulmNet_data <- read.table(file.path(mainDir,"CardiopulmNet_data.csv"),sep=";",
  dec=".",na.strings="",header=T)

# Read JAGS code as text
pnet.code <- readLines(file.path(mainDir,"CardiopulmNet_mcmcModel.r"))

# Create the data list readable by JAGS through the user function 'jagsData'.
# Arguments of the function:
#  - 'code': a JAGS model loaded on R using readLines()
#  - 'data': a dataset consistent with the JAGS model
#  - 'tol': the numerical tolerance for 0 probabiity values
# A list with two components is returned: 'data' (containing data) and 'inits'
# (containing initial values for non-stochastic relationships)
pnet.data <- jagsData(code=pnet.code,data=CardiopulmNet_data,tol=1e-7)

# Load rjags package (make sure to have JAGS installed on your system)
library(rjags)
load.module("glm")

# Initialize the jags model (Mersenne-Twister random number generator with seed 123)
pnet.jags <- jags.model(file=file.path(mainDir,"CardiopulmNet_mcmcModel.r"),
  data=pnet.data$data,inits=c(pnet.data$inits,'.RNG.name'="base::Mersenne-Twister",'.RNG.seed'=123),
  n.chains=1,n.adapt=5000)

# Perform MCMC simulation (batches of 10000 iterations - each one takes approximately 13 hours with JAGS 3.4.0)
pnet.mcmc1 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc1,file=file.path(mainDir,"save/sim/pnet.mcmc1.RData"))
pnet.mcmc2 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc2,file=file.path(mainDir,"save/sim/pnet.mcmc2.RData"))
pnet.mcmc3 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc3,file=file.path(mainDir,"save/sim/pnet.mcmc3.RData"))
pnet.mcmc4 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc4,file=file.path(mainDir,"save/sim/pnet.mcmc4.RData"))
pnet.mcmc5 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc5,file=file.path(mainDir,"save/sim/pnet.mcmc5.RData"))
pnet.mcmc6 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc6,file=file.path(mainDir,"save/sim/pnet.mcmc6.RData"))
pnet.mcmc7 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc7,file=file.path(mainDir,"save/sim/pnet.mcmc7.RData"))
pnet.mcmc8 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc8,file=file.path(mainDir,"save/sim/pnet.mcmc8.RData"))
pnet.mcmc9 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc9,file=file.path(mainDir,"save/sim/pnet.mcmc9.RData"))
pnet.mcmc10 <- coda.samples(pnet.jags,pnet.data$monitor,n.iter=10000)
save(pnet.mcmc10,file=file.path(mainDir,"save/sim/pnet.mcmc10.RData"))

# Assemble batches of MCMC output + thinning
mysim <- c()
for(i in 5:10) {
  load(file.path(mainDir,paste("save/sim/pnet.mcmc",i,".RData",sep="")))
  mystr <- paste("mysim <- rbind(mysim,pnet.mcmc[[1]][seq(1,10000,by=3),])",sep="")
  eval(parse(text=mystr))
  }
pnet.mcmc <- mcmc.list(mcmc(mysim))
save(pnet.mcmc,file=file.path(mainDir,"save/pnet.mcmc.RData"))


### Run the following line to load the output of MCMC simulation
#   load(file.path(mainDir,"save/pnet.mcmc.RData"))


# Produce trace plots through the user function 'mcmcTrace'.
# Arguments of the function:
#  - 'sim': the output of a MCMC simulation
#  - 'plotDir': the path to a directory where plots are saved
mcmcTrace(pnet.mcmc,plotDir=file.path(mainDir,"save/trace"))

# Perform convergence diagnostic tests (shown in Table 4 of the manuscript)
# through the user function 'mcmcDiag'. Arguments of the function:
#  - 'sim': the output of a MCMC simulation
pnet.diag <- mcmcDiag(pnet.mcmc)
save(pnet.diag,file=file.path(mainDir,"save/pnet.diag.RData"))
pnet.diag$summary_n
pnet.diag$summary_p

# Produce summaries and density plots through the user function 'mcmcSumm'.
# Arguments of the function:
#  - 'code': a JAGS model loaded on R using readLines()
#  - 'sim': the output of a MCMC simulation
#  - 'plotDir': the path to a directory where to save plots
# A matrix is returned with prior and posterior summaries and the D-statistic for each stochastic parameter.
pnet.summ <- mcmcSumm(pnet.code,pnet.mcmc,plotDir=file.path(mainDir,"save/density"))
save(pnet.summ,file=file.path(mainDir,"save/pnet.summ.RData"))

# Produce the summary of D-statistic (shown in Table 5 of the manuscript)
pnet.div <- discretize(pnet.summ[,"D-stat"],c(0.0001,0.1,0.5,0.95))
table(pnet.div)
round(prop.table(table(pnet.div)),4)

# Implement the probabilistic network in GeNIe through the user function 'genieExport'.
# Arguments of the function:
#  - 'code': a JAGS model loaded on R using readLines()
#  - 'sim': the output of a MCMC simulation
#  - 'nclass': the number of categories for each range (lp, n, hp) of continuous variables
#  - 'node.labels': a named vector containing the full name of each variable
#  - 'state.labels': a named list containing the name of categories for each variable
#  - 'file.name': the name of the GeNIe file
#  - 'save.dir': the path to a directory where to save the GeNIe file
pnet.cpts <- genieExport(pnet.code,pnet.mcmc,nclass=3,node.labels=CardiopulmNet_info$titles,
  state.labels=CardiopulmNet_info$samspace,file.name="CardiopulmNet_BNet",save.dir=file.path(mainDir,"save"))
# save the conditional probability tables
save(pnet.cpts,file=file.path(mainDir,"save/pnet.cpts.RData"))
