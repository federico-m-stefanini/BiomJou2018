
## splitting of the original object
for(aux in seq(1000,20000,1000)){
  tmp<- pnet.mcmc[[1]][(aux-999):aux,]
  save(tmp,
       file=paste("pnet.mcmc.p",aux,".RData",sep="")
  )}


## merging of component files

ppnet.mcmc<- list()
load(file=paste("pnet.mcmc.p1000.RData",sep=""))
ppnet.mcmc[[1]] <- tmp
for(aux in seq(2000,20000,1000)){
  filename <- paste("pnet.mcmc.p",aux,".RData",sep="")
  load(file=filename)
  ppnet.mcmc[[1]] <- rbind(ppnet.mcmc[[1]], tmp)
  }


# some checking
sum(pnet.mcmc[[1]][1:20000,] != ppnet.mcmc[[1]][,])

pnet.mcmc[[1]][100,20:30] 
ppnet.mcmc[[1]][100,20:30]

pnet.mcmc[[1]][1100,20:30] 
ppnet.mcmc[[1]][1100,20:30]


## ppnet.mcmc is the reassembled-merged  matrix
