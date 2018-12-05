require(coda)
require(TeachingDemos)

cscales <- list(
  Age=c(0,35,70,105),
  O2=c(200,200,80,15),
  SaO2=c(100,100,98,25),
  FiO2=c(20,20,22,100.01),
  Emoglob=c(20,20,12,3),
  Latt=c(0,0,1.8,5),
  CO2=c(0,35,45,100),
  FC=c(0,60,80,200),
  PA=c(0,75,100,200),
  pH=c(6.8,7.35,7.45,7.8),
  FT3=c(0,2.3,4.5,13),
  FT4=c(0,0.7,1.8,7),
  TSH=c(0,0.3,5,20)
  )

# discretization
discretize <- function(x,cutpoints,begin=1,extremes=F) {
  if((is.vector(x)==F || is.numeric(x)==F) |
    (identical(cutpoints,sort(cutpoints,T))==F &
     identical(cutpoints,sort(cutpoints,F))==F)) stop("Invalid arguments")
  if(extremes==T) {
    x[which(x<min(cutpoints) | x>max(cutpoints))] <- NA
    }
  dx <- x
  if(identical(cutpoints,sort(cutpoints,decreasing=T))) {
    decr <- 1
    if(extremes==F) {
      cutaux <- c(Inf,cutpoints,-Inf)
      } else if(extremes==T) {
      cutaux <- c(Inf,cutpoints[-c(1,length(cutpoints))],-Inf)
      }
    } else {
    decr <- 0
    if(extremes==F) {
      cutaux <- c(-Inf,cutpoints,Inf)
      } else if(extremes==T) {
      cutaux <- c(-Inf,cutpoints[-c(1,length(cutpoints))],Inf)
      }
    }
  for(i in 1:length(x)) {
    for(j in 1:(length(cutaux)-1)) {
      if(!is.na(x[i])) {
        if(decr==1) {
          if(x[i]<cutaux[j] & x[i]>=cutaux[j+1]) dx[i] <- j-1+begin
          } else if(decr==0) {
          if(x[i]>=cutaux[j] & x[i]<cutaux[j+1]) dx[i] <- j-1+begin
          }
        }
      }
    }
  names(dx) <- names(x)
  dx
  }

# rescale a continuous variable within the interval (0,1)
jagsRescale <- function(v,vscale) {
  if(identical(sort(vscale,T),vscale)) {
    vscale <- -vscale
    v <- -v
    }
  vtil <- c()                                         
  for(i in 1:length(v)) {
    if(!is.na(v[i])) {
      if(v[i]<=vscale[1]) vtil[i] <- 0
      if(v[i]>vscale[1]&v[i]<vscale[2]) vtil[i] <- 1/3*(v[i]-vscale[1])/(vscale[2]-vscale[1])
      if(v[i]>=vscale[2]&v[i]<=vscale[3]) vtil[i] <- 1/3+1/3*(v[i]-vscale[2])/(vscale[3]-vscale[2])
      if(v[i]>=vscale[3]&v[i]<vscale[4]) vtil[i] <- 2/3+1/3*(v[i]-vscale[3])/(vscale[4]-vscale[3])
      if(v[i]>=vscale[4]) vtil[i] <- 1
      } else {
      vtil[i] <- NA
      }
    }
  aux0 <- which(vtil<=0)
  if(length(aux0)>0) {
    x0 <- vtil[setdiff(1:length(vtil),aux0)]
    if(length(x0)>0) {
      vtil[aux0] <- min(0.0001,min(x0,na.rm=T))
      } else {
      vtil[aux0] <- 0.0001
      }
    }
  aux1 <- which(vtil>=1)
  if(length(aux1)>0) {
    x1 <- vtil[setdiff(1:length(vtil),aux1)]
    if(length(x1)>0) {    
      vtil[aux1] <- max(0.9999,max(x1,na.rm=T))  
      } else {
      vtil[aux1] <- 0.9999
      }
    }
  vtil
  }

# inverted rescaling
jagsRescale.inv <- function(v,vscale) {
  vtil <- c()
  for(i in 1:length(v)) {
    if(!is.na(v[i])) {
      if(v[i]< 1/3) vtil[i] <- vscale[1]+v[i]*3*(vscale[2]-vscale[1])
      if(v[i]>= 1/3&v[i]<=2/3) vtil[i] <- vscale[2]+3*(v[i]-1/3)*(vscale[3]-vscale[2])
      if(v[i]>2/3) vtil[i] <- vscale[3]+3*(v[i]-2/3)*(vscale[4]-vscale[3])
      } else {
      vtil[i] <- NA
      }
    }
  vtil
  }

# clean code
cleanCode <- function(code) {
  for(i in 1:length(code)) {
    if(length(grep("#",code[i]))>0) {
      code[i] <- strsplit(code[i],"#")[[1]][1]
      }
    }
  res <- code
  for(i in 1:length(code)) {
    if(grepl(";",code[i])) {
      if(i==1) {
        res <- c(strsplit(code[i],";")[[1]],res[(i+1):length(res)])
        } else if(i==length(res)) {
        res <- c(res[1:(i-1)],strsplit(code[i],";")[[1]])
        } else {
        res <- c(res[1:(i-1)],strsplit(code[i],";")[[1]],res[(i+1):length(res)])
        }
      }
    }
  res <- gsub(" in ","\\*\\_in\\_\\*",res)
  res <- gsub(" ","",res)
  res <- gsub("\\*\\_in\\_\\*"," in ",res)
  auxind <- which(substr(res,nchar(res),nchar(res)) %in% c("+","-","*","/","(","^"))
  auxblank <- which(grepl("^\\s*$",res))
  auxdel <- c()
  if(length(grep("model\\{",res))>0) {
    auxdel <- c(auxdel,which(substr(res,1,6)=="model{"),max(which(grepl("\\}",res))))
    }
  if(length(auxind)>0) {
    auxL <- split(auxind,cumsum(c(1,diff(auxind)!= 1)))
    for(i in 1:length(auxL)) {
      iaux <- c(auxL[[i]],max(auxL[[i]])+1)
      res[iaux[1]] <- paste(res[iaux],collapse="")
      auxdel <- c(auxdel,iaux[-1])
      }
    }
  if(length(auxblank)>0) {
    auxB <- split(auxblank,cumsum(c(1,diff(auxblank)!= 1)))
    for(i in 1:length(auxB)) {
      auxdel <- c(auxdel,auxB[[i]][-1])
      }
    }
  res[setdiff(1:length(res),auxdel)]
  }

# get sub-string between two characters
strBW <- function(x,ch1,ch2) {
  strsplit(strsplit(x,ch1)[[1]][-1],ch2)[[1]][1]
  }

# find variable names by typology in jags code
findVar <- function(code) {
  code <- cleanCode(code)
  xcont <- c()
  auxcont <- code[grep("tilde\\..*dbeta",code)]
  for(i in 1:length(auxcont)) {
    xcont <- c(xcont,strBW(auxcont[i],"tilde\\.","\\["))
    }
  xbin <- c()
  auxbin <- code[grep(".*dbern",code)]
  for(i in 1:length(auxbin)) {
    xbin <- c(xbin,strsplit(auxbin[i],"\\[")[[1]][1])
    }
  xcat <- c()
  auxcat <- code[grep("z\\..*dcat",code)]
  for(i in 1:length(auxcat)) {
    xcat <- c(xcat,strBW(auxcat[i],"z\\.","\\["))
    }
  list(continuous=sort(unique(xcont)),binary=sort(unique(xbin)),multivalued=sort(unique(xcat)))
  }

# find parameters in jags code
findParams <- function(code) {
  code <- cleanCode(code)
  auxpar <- code[grep("^p\\..*~",code)]
  auxK <- code[grep("^p\\..*<-",code)]
  parnam <- parnum <- Knam <- Knum <- c()
  for(i in 1:length(auxpar)) {
    parnam <- c(parnam,paste("p.",strBW(auxpar[i],"p\\.","\\["),sep=""))
    parnum <- c(parnum,strBW(auxpar[i],"\\[","\\]"))
    }
  for(i in 1:length(auxK)) {
    Knam <- c(Knam,paste("p.",strBW(auxK[i],"p\\.","\\["),sep=""))
    Knum <- strBW(auxK[i],"\\[","\\]")
    }
  auxeta <- code[grep("^eta\\..*~",code)]
  auxetK <- code[grep("^eta\\..*<-",code)]
  for(i in 1:length(auxeta)) {
    parnam <- c(parnam,paste("eta.",strBW(auxeta[i],"eta\\.","\\["),sep=""))
    parnum <- c(parnum,strBW(auxeta[i],"\\[","\\]"))
    }
  for(i in 1:length(auxetK)) {
    Knam <- c(Knam,paste("eta.",strBW(auxetK[i],"eta\\.","\\["),sep=""))
    Knum <- strBW(auxetK[i],"\\[","\\]")
    }
  allpnam <- sort(unique(c(parnam,Knam)))
  res <- list()
  for(i in 1:length(allpnam)) {
    auxind1 <- which(parnam==allpnam[i])
    auxstr1 <- paste("c(",paste(parnum[auxind1],collapse=","),")",sep="")
    iPar <- eval(parse(text=auxstr1))
    auxind2 <- which(Knam==allpnam[i])
    auxstr2 <- paste("c(",paste(Knum[auxind2],collapse=","),")",sep="")
    iK <- eval(parse(text=auxstr2))
    iN <- length(c(iPar,iK))
    ires <- rep(F,iN)
    ires[iPar] <- T
    res[[i]] <- ires
    }
  names(res) <- allpnam
  res
  }

# find prior on parameters by typology in jags code
findPrior <- function(code) {
  code <- cleanCode(code)
  auxtheta <- c(code[grep("^p\\..*~dbeta",code)],code[grep("^eta\\..*~dnorm",code)],code[grep("^p\\..*~dgamma",code)],code[grep(".*~ddirch",code)])
  betanam <- dirnam <- gaunam <- gamnam <- betamat <- gaumat <- gammat <- c()
  dirmat <- list()
  for(i in 1:length(auxtheta)) {
    if(grepl("dbeta",auxtheta[i])) {
      betanam <- c(betanam,strsplit(auxtheta[i],"~")[[1]][1])
      ihp <- strsplit(strBW(auxtheta[i],"\\(","\\)"),",")[[1]]
      betamat <- rbind(betamat,sapply(ihp,function(x){eval(parse(text=x))}))
      } else if(grepl("ddirch",auxtheta[i])) {
      if(grepl("p\\.",auxtheta[i])) {                                            
        inam <- strBW(auxtheta[i],"p\\.","\\[")
        ialph <- strBW(auxtheta[i],"\\(","\\)")       
        inum <- eval(parse(text=strBW(auxtheta[i],"\\[","\\]")))                        
        dirnam <- c(dirnam,paste("p.",inam,"[",paste(sort(inum),collapse=","),"]",sep=""))
        ihp <- c()
        for(j in 1:length(inum)) {
          ijdir <- code[grep(paste(strsplit(ialph,"\\[")[[1]][1],"\\[",j,"\\]",sep=""),code)]   
          ijal <- strsplit(ijdir,"<-")[[1]][2]                   
          ihp[j] <- eval(parse(text=paste("tol <- 0; ",ijal,sep="")))
          }
        names(ihp) <- paste("p.",inam,"[",inum,"]",sep="")
        } else {
        iauxnam <- strsplit(auxtheta[i],"~")[[1]][1]
        inam <- strBW(auxtheta[i],"\\.","~")
        iauxhp <- strBW(auxtheta[i],"\\(","\\)")                           
        icod_hp <- setdiff(code[grep(strsplit(iauxhp,"\\[")[[1]][1],code)],auxtheta[i])
        ihpstr <- sapply(strsplit(icod_hp,"<-"),function(z){z[2]})         
        ihp <- eval(parse(text=paste("c(",paste(ihpstr,collapse=","),")",sep="")))
        icod_num <- setdiff(code[grep(strsplit(iauxnam,"\\[")[[1]][1],code)],auxtheta[i])
        inumstr <- sapply(strsplit(icod_num,"<-"),function(z){strBW(z[1],"\\[","\\]")})         
        inum <- eval(parse(text=paste("c(",paste(inumstr,collapse=","),")",sep="")))
        dirnam <- c(dirnam,paste("p.",inam,"[",paste(sort(inum),collapse=","),"]",sep=""))
        #dirnam <- c(dirnam,iauxnam)
        names(ihp) <- paste("p.",inam,"[",inum,"]",sep="")                           
        #names(ihp) <- paste(iauxnam,"[1:",length(ihp),"]",sep="")
        }
      dirmat <- c(dirmat,list(ihp))
      } else if(grepl("dnorm",auxtheta[i])) {
      gaunam <- c(gaunam,strsplit(auxtheta[i],"~")[[1]][1])
      ihp <- strsplit(strBW(auxtheta[i],"\\(","\\)"),",")[[1]]
      gaumat <- rbind(gaumat,sapply(ihp,function(x){eval(parse(text=x))}))
      } else if(grepl("dgamma",auxtheta[i])) {
      gamnam <- c(gamnam,strsplit(auxtheta[i],"~")[[1]][1])
      ihp <- strsplit(strBW(auxtheta[i],"\\(","\\)"),",")[[1]]
      gammat <- rbind(gammat,sapply(ihp,function(x){eval(parse(text=x))}))
      }
    }
  if(!is.null(betamat)) {
    rownames(betamat) <- betanam
    colnames(betamat) <- c("hp1","hp2")
    }
  if(!is.null(gaumat)) {
    rownames(gaumat) <- gaunam
    colnames(gaumat) <- c("hp1","hp2")
    }
  if(!is.null(gammat)) {
    rownames(gammat) <- gamnam
    colnames(gammat) <- c("hp1","hp2")
    }
  if(!is.null(dirmat)) names(dirmat) <- dirnam
  list(beta=betamat,dirichlet=dirmat,gaussian=gaumat,gamma=gammat)
  }

# initial values for deterministic relationships
varInits <- function(x) {
  n <- max(sapply(x,length))
  res <- list()
  res[["Aort_diss"]] <- x$Flap
  res[["Aneur"]] <- x$Dil_aort
  res[["Arit_sopra_acut"]] <- rep(NA,n)
  res[["Arit_sopra_acut"]][which(is.na(x$Arit_sopra_acut)&x$Arit_sopra_o==0)] <- 0
  res[["Arit_sopra_acut"]][which(is.na(x$Arit_sopra_acut)&x$Arit_sopra_o==1)] <- 1
  res[["Arit_sopra_cron"]] <- rep(NA,n)
  res[["Arit_sopra_cron"]][which(x$Arit_sopra_o==0)] <- 0
  res[["Arit_sopra_cron"]][which(x$Arit_sopra_o==1)] <- 1
  res[["z.Opac"]] <- rep(NA,n)
  res[["z.Opac"]][which(x$intersPol==0|x$z.OpacRX==1)] <- 1
  res[["z.Opac"]][which(x$z.OpacRX==2)] <- 2
  res[["z.Opac"]][which(x$z.OpacRX==3)] <- 3
  res[["z.Cardiomio"]] <- rep(NA,n)
  res[["z.Cardiomio"]][which(x$IpertSx==0)] <- 1
  res[["z.ConvPth"]] <- rep(NA,n)
  res[["z.ConvPth"]][which(x$Conv==0|x$PerdCosc==0)] <- 1
  res[["z.ConvPth"]][which(x$Conv==1)] <- 3
  res[["z.BPCOa"]] <- rep(NA,n)
  res[["z.BPCOa"]][which(is.na(x$z.BPCOa)&(x$BPCOa_J==0|x$z.EO_bsp==1))] <- 1
  res[["z.BPCOa"]][which(is.na(x$z.BPCOa)&(x$BPCOa_J==1|x$z.EO_bsp>1))] <- 2
  res[["Card_dil"]] <- rep(NA,n)
  res[["Card_dil"]][which(is.na(x$Card_dil)&x$DilatVsx==0)] <- 0
  res[["z.Consol"]] <- rep(NA,n)
  res[["z.Consol"]][which(is.na(x$z.Consol)&x$bronAer==0)] <- 1
  res[["z.Consol"]][which(is.na(x$z.Consol)&x$bronAer==1)] <- 2
  res[["InP"]] <- rep(NA,n)
  res[["InP"]][which(is.na(x$InP)&x$z.Consol==3)] <- 1
  res[["z.Polm"]] <- rep(NA,n)
  res[["z.Polm"]][which(is.na(x$z.Polm)&(x$Polm_J==0|res[["z.Opac"]]==1|res[["z.Consol"]]<3))] <- 1          
  res[["z.Polm"]][which(is.na(x$z.Polm)&x$Polm_J==1)] <- 2
  res[["z.EP"]] <- rep(NA,n)
  res[["z.EP"]][which(x$EP_J==0)] <- 1
  res[["z.EP"]][which(x$EP_J==1)] <- 2
  res[["z.EP"]][which((is.na(x$EP_J)|x$EP_J==1)&res[["InP"]]==1)] <- 2
  res[["z.EP"]][which((is.na(x$EP_J)|x$EP_J==1)&x$z.EO_bsp>1&x$z.Toss==1)] <- 2
  res[["z.Pnx"]] <- rep(NA,n)
  res[["z.Pnx"]][which(is.na(x$z.Pnx)&x$z.Iper_Diaf<3)] <- 1
  res[["z.Pnx"]][which(is.na(x$z.Pnx)&x$z.Iper_Diaf==3)] <- 2
  res[["z.Enf"]] <- rep(NA,n)
  res[["z.Enf"]][which(x$z.Iper_Diaf==1|res[["z.Opac"]]==1)] <- 1
  res[["z.Enf"]][which(x$z.Iper_Diaf==2|res[["z.BPCOa"]]==2|x$z.BPCOa==2|x$z.EO_bsp>1)] <- 2
  res[["z.Enf"]][which(x$z.BPCOa==3)] <- 3
  res[["Asma"]] <- rep(NA,n)
  res[["Asma"]][which(is.na(x$Asma)&(res[["z.BPCOa"]]>1&res[["z.Enf"]]==1))] <- 1
  res[["Insuf_aorAcut"]] <- rep(NA,n)
  res[["Insuf_aorAcut"]][which(is.na(x$Insuf_aorAcut)&x$Insuf_aor==0)] <- 0
  res[["Endocard"]] <- rep(NA,n)
  res[["Endocard"]][which(is.na(x$Endocard)&x$Veget_end==0)] <- 0
  res[["z.Insuf_aorCron"]] <- rep(NA,n)
  res[["z.Insuf_aorCron"]][which(x$Insuf_aor==0)] <- 1
  res[["z.Insuf_aorCron"]][which(x$Insuf_aor==1)] <- 2
  res[["InsuffMitrAcut"]] <- rep(NA,n)
  res[["InsuffMitrAcut"]][which(is.na(x$InsuffMitrAcut)&x$InsuffMitr==0)] <- 0
  res[["z.InsuffMitrCron"]] <- rep(NA,n)
  res[["z.InsuffMitrCron"]][which(x$InsuffMitr==0)] <- 1
  res[["z.InsuffMitrCron"]][which(x$InsuffMitr==1)] <- 2
  res[["ProlMitralAcut"]] <- rep(NA,n)
  res[["ProlMitralAcut"]][which(x$ProlMitral==0)] <- 0
  res[["ProlMitralCron"]] <- rep(NA,n)
  res[["ProlMitralCron"]][which(x$ProlMitral==0)] <- 0
  res[["ProlMitralCron"]][which(x$ProlMitral==1)] <- 1
  res[["z.IMA"]] <- rep(NA,n)
  res[["z.IMA"]][which(x$Trop_o==0|x$IMA_J==0)] <- 1
  res[["z.IMA"]][which(x$IMA_J==1)] <- 3
  res[["z.IMA"]][which((is.na(x$IMA_J)|x$IMA_J==0)&x$CK_MB==1)] <- 2
  res[["z.IMA"]][which(x$IMA_J==1&x$CK_MB==1)] <- 3
  res[["Miocard"]] <- rep(NA,n)
  res[["Miocard"]][which(is.na(x$Miocard)&res[["z.IMA"]]>1)] <- 1
  res[["Interst"]] <- rep(NA,n)
  res[["Interst"]][which(is.na(x$Interst)&res[["z.Opac"]]==1)] <- 0
  res[["tilde.Shu"]] <- rep(NA,n)
  res[["tilde.Shu"]][which(x$Shu_J==0|x$Dispnea==0)] <- 0.5
  res[["tilde.Shu"]][which(x$Shu_J==1)] <- 0.9
  res[["tilde.pulmPerf"]] <- rep(NA,n)
  res[["tilde.pulmPerf"]][which(x$pulmPerf_J==1)] <- 0.9
  res[["z.Ed_polm"]] <- rep(NA,n)
  res[["z.Ed_polm"]][which(x$Ed_polm_J==0)] <- 1
  res[["z.Ed_polm"]][which(x$Ed_polm_J==1|(is.na(x$Ed_polm_J)&res[["z.Opac"]]==1))] <- 2
  res[["tilde.LHI"]] <- rep(NA,n)
  res[["tilde.LHI"]][which(x$LHI_J==0|res[["z.Ed_polm"]]==1)] <- 0.5
  res[["tilde.LHI"]][which(x$LHI_J==1)] <- 0.9
  res[["ARDS"]] <- rep(NA,n)
  res[["ARDS"]][which(is.na(x$ARDS)&res[["z.Ed_polm"]]==1)] <- 0
  res[["ARDS"]][which(is.na(x$ARDS)&res[["z.Ed_polm"]]>1)] <- 1
  res[["tilde.leftPump"]] <- rep(NA,n)
  res[["tilde.leftPump"]][which(x$leftPump_J==0)] <- 0.5
  res[["tilde.leftPump"]][which(x$leftPump_J==1|x$z.BPCOa==4)] <- 0.9
  res[["tilde.RHO"]] <- rep(NA,n)
  res[["tilde.RHO"]][which(x$RHO_J==0)] <- 0.5
  res[["tilde.RHO"]][which(x$RHO_J==1)] <- 0.9
  res[["tilde.Latt"]] <- rep(NA,n)
  res[["tilde.Latt"]][which(x$Shock==0)] <- 0.5
  res[["tilde.Latt"]][which(x$Shock==1)] <- 0.9
  res[["tilde.Disidrat"]] <- rep(NA,n)
  res[["tilde.Disidrat"]][which(x$Disidrat_J==0)] <- 0.5
  res[["tilde.Disidrat"]][which(x$Disidrat_J==1)] <- 0.9
  res[["tilde.SaO2"]] <- rep(NA,n)
  res[["tilde.SaO2"]][which(is.na(x$tilde.SaO2)&x$SaO2_J==0)] <- 0.5
  res[["tilde.SaO2"]][which(is.na(x$tilde.SaO2)&(x$Cia==1|x$SaO2_J==1))] <- 0.9
  res[["tilde.CO2"]] <- rep(NA,n)
  res[["tilde.CO2"]][which(x$CO2_J==0)] <- 0.5
  res[["tilde.CO2"]][which(x$CO2_J==1)] <- 0.9
  res[["tilde.minVent"]] <- rep(NA,n)
  res[["tilde.minVent"]][which(x$Tcp==1)] <- 0.9
  res[["tilde.RHI"]] <- rep(NA,n)
  res[["tilde.RHI"]][which(x$ECG_dx==1)] <- 0.9
  res[["Scint_bias"]] <- rep(NA,n)
  res[["Scint_bias"]][which(x$S_per_r==0)] <- 0
  res[["HemoPeric"]] <- rep(NA,n)
  res[["HemoPeric"]][which(x$Vers_peric==0)] <- 0
  res[["HemoPeric"]][which(x$z.tamp>1)] <- 1
  res[["z.tamp"]] <- rep(NA,n)
  res[["z.tamp"]][which(x$DilatVdx==1)] <- 1
  res[["URI"]] <- rep(NA,n)
  res[["URI"]][which(is.na(x$URI)&x$z.EO_bsp>1)] <- 1
  res[["Fertil"]] <- rep(NA,n)
  res[["Fertil"]][which(x$Estrog==1|x$z.Grav>1)] <- 1
  res[["Infiam_b"]] <- rep(NA,n)
  res[["Infiam_b"]][which(x$FANS==1|x$Fever==1)] <- 1
  res[["Infiam_nb"]] <- rep(NA,n)
  res[["Infiam_nb"]][which(x$FANS==1|x$Fever==1)] <- 1
  res[["D_a_pu"]] <- rep(NA,n)
  res[["D_a_pu"]][which(x$z.Dol_type==1)] <- 0
  res[["D_cen"]] <- rep(NA,n)
  res[["D_cen"]][which(x$z.Dol_type<=2)] <- 0
  res[["D_pl"]] <- rep(NA,n)
  res[["D_pl"]][which(x$z.Dol_type<=3)] <- 0
  res[["D_quadr_sup"]] <- rep(NA,n)
  res[["D_quadr_sup"]][which(x$z.Dol_type<=4)] <- 0
  res[["D_par"]] <- rep(NA,n)
  res[["D_par"]][which(x$z.Dol_type<=5)] <- 0
  res[["Hemo"]] <- rep(NA,n)
  res[["Hemo"]][which(is.na(x$Hemo)&x$acuteAnem==0)] <- 0
  res[["z.Dispeps"]] <- rep(NA,n)
  res[["z.Dispeps"]][which(x$Dispeps_J==0)] <- 1
  res[["z.Dispeps"]][which(x$Dispeps_J==1)] <- 2
  res[["BBSx"]] <- rep(NA,n)
  res[["BBSx"]][which(is.na(x$BBSx)&(x$z.Q_inT>0|x$z.nonSTE>0))] <- 0
  res[["z.ArtIntraVasCoa"]] <- rep(NA,n)
  res[["z.ArtIntraVasCoa"]][which(is.na(x$z.ArtIntraVasCoa)&x$Coronar==0)] <- 1
  res[["EmbArte"]] <- rep(NA,n)
  res[["EmbArte"]][which(x$IschCer==0)] <- 0
  res[["z.Agit"]] <- rep(NA,n)
  res[["z.Agit"]][which(x$Agit_o==0)] <- 1
  res[["z.Agit"]][which(x$Agit_o==1)] <- 2
  res[["Linfoc_o"]] <- rep(NA,n)
  res[["Linfoc_o"]][which(is.na(x$Linfoc_o)&x$z.Leucoc_o<3)] <- 0
  res[["Linfoc_o"]][which(is.na(x$Linfoc_o)&x$z.Leucoc_o>=3)] <- 1
  res[["Neutrof_o"]] <- rep(NA,n)
  res[["Neutrof_o"]][which(x$z.Leucoc_o==2&(x$Linfoc_o==0|res[["Linfoc_o"]]==0))] <- 0
  res[["Neutrof_o"]][which(x$z.Leucoc_o<=3&(x$Linfoc_o==1|res[["Linfoc_o"]]==1))] <- 0
  res[["Neutrof_o"]][which(x$z.Leucoc_o==3&(x$Linfoc_o==0|res[["Linfoc_o"]]==0))] <- 1
  res[["Neutrof_o"]][which(x$z.Leucoc_o==4)] <- 1
  res[["tilde.OT"]] <- rep(NA,n)
  res[["tilde.OT"]][which(x$Tireotos==1)] <- 0.9
  res[["tilde.driveCard"]] <- rep(NA,n)
  res[["tilde.driveCard"]][which(x$z.BradiTachi==2)] <- 0.9  
  res[["tilde.pulmPerf"]][which((is.na(x$pulmPerf_J)|x$pulmPerf_J==0)&x$S_per_r==1)] <- 0.6
  res[["tilde.pulmPerf"]][which(x$pulmPerf_J==1&x$S_per_r==1)] <- 0.9
  res
  }
  
# create data for jags
jagsData <- function(code,data,tol=1e-7) {
  if(!is.numeric(tol)||tol>=0.01||tol<=0) stop("Argument 'tol' must be a positive number less than 0.01",call.=F)
  auxV <- findVar(code)
  netvar <- c()
  if(length(auxV$continuous)>0) netvar <- c(netvar,paste("tilde.",auxV$continuous,sep=""))
  if(length(auxV$multivalued)>0) netvar <- c(netvar,paste("z.",auxV$multivalued,sep=""))
  if(length(auxV$binary)>0) netvar <- c(netvar,auxV$binary)
  xcat <- auxV$multivalued
  out <- list()
  varnam <- colnames(data)
  outnam <- c()
  for(i in 1:length(varnam)) {
    if(varnam[i] %in% xcat) {
      out[[i]] <- data[,varnam[i]]+1
      outnam[i] <- paste("z.",varnam[i],sep="")
      } else if(varnam[i] %in% names(cscales)) {
      out[[i]] <- jagsRescale(data[,varnam[i]],cscales[[varnam[i]]])
      outnam[i] <- paste("tilde.",varnam[i],sep="")
      } else {
      out[[i]] <- data[,varnam[i]]
      outnam[i] <- varnam[i]
      }
    }
  names(out) <- outnam
  namOK <- intersect(outnam,netvar)  
  mydat <- c(out[namOK],list(nrow(data),tol))
  names(mydat) <- c(namOK,"n","tol")
  xpar <- findParams(code)
  myini <- varInits(out)
  list(data=mydat,monitor=names(xpar),inits=myini)
  }

# sort parameters by name
sortParams <- function(x) {
  auxnam <- auxnum <- c()
  for(i in 1:length(x)) {
    if(grepl("\\[",x[i])==F) {
      ix <- paste(x[i],"[1]",sep="")
      } else {
      ix <- x[i]
      }
    auxnam[i] <- strsplit(ix,"\\[")[[1]][1]    
    auxnum[i] <- as.numeric(strBW(ix,"\\[","\\]"))
    }
  allnam <- sort(unique(auxnam))
  res <- c()
  for(i in 1:length(allnam)) {
    auxind <- which(auxnam==allnam[i])
    if(length(auxind)==1) {
      ipar <- x[auxind]
      } else {
      ipar <- x[auxind][order(auxnum[auxind])]
      }
    res <- c(res,ipar)
    }
  res
  }

# mcmc traces
mcmcTrace <- function(sim,plotDir) {
  if(class(sim)!="mcmc.list") stop("argument 'sim' must be an object of class 'mcmc.list'")
  nomiOK <- colnames(sim[[1]])[which(apply(sim[[1]],2,var)>0)]
  pdf(file.path(plotDir,"trace%03d.pdf"),onefile=F)
  plot(sim[,sortParams(nomiOK)],density=F)
  dev.off()
  }

# mcmc diagnostic
mcmcDiag <- function(sim) {
  if(class(sim)!="mcmc.list") stop("argument 'sim' must be an object of class 'mcmc.list'")
  nomi <- sortParams(colnames(sim[[1]])[which(apply(sim[[1]],2,var)>0)])
  res <- summ <- summp <- list()
  for(i in 1:length(sim)) {
    imat <- matrix(nrow=length(nomi),ncol=3)
    rownames(imat) <- nomi
    colnames(imat) <- c("Geweke","Heidel","Raftery")
    for(j in 1:length(nomi)) {         
      ijdat <- sim[[i]][,nomi[j]]
      imat[j,1] <- 1*(abs(geweke.diag(ijdat)$z)>1.96)
      if(var(ijdat)>1e-6) {
        imat[j,2] <- 1*(heidel.diag(ijdat)[,"htest"]==0)
        imat[j,3] <- 1*(raftery.diag(ijdat)$resmatrix[,"I"]>5)
        } else {
        imat[j,2] <- imat[j,3] <- 0
        }
      }
    res[[i]] <- imat
    ifail <- apply(imat,1,sum)
    itab <- table(3-ifail)
    icum <- c()
    icum[1] <- itab[1]
    icum[2] <- sum(itab[2:4])
    icum[3] <- sum(itab[3:4])
    icum[4] <- itab[4]
    names(icum) <- c("0",">=1",">=2","3")
    summ[[i]] <- round(icum,4)
    summp[[i]] <- round(icum/length(ifail),4)
    }
  list(tests=res,summary_n=summ,summary_p=summp)
  }

# mcmc summary
mcmcSumm <- function(code,sim,p=0.95,digits=6,plotDir=NULL) {
  if(class(sim)!="mcmc.list") stop("argument 'sim' must be an object of class 'mcmc.list'")
  if(!is.numeric(digits)||digits<=0||round(digits)!=digits) stop("Argument 'digits' must be a positive integer number",call.=F)
  if(!is.numeric(p)||p<=0||p>=1) stop("Argument 'p' must be a positive number less than 1",call.=F)  
  mydat <- c()
  for(i in 1:length(sim)) {
    mydat <- rbind(mydat,sim[[i]])
    }
  pMat <- findPrior(code)
  auxnam <- colnames(mydat)
  for(i in 1:length(auxnam)) {
    if(grepl("\\[",auxnam[i])==F) auxnam[i] <- paste(auxnam[i],"[1]",sep="")
    }
  colnames(mydat) <- auxnam
  dirnam <- lapply(pMat$dirichlet,names)
  dirOne <- sapply(pMat$dirichlet,function(z){names(z)[1]})
  nomi <- sortParams(setdiff(auxnam,dirOne))
  res <- matrix(nrow=length(nomi),ncol=9)
  rownames(res) <- nomi
  colnames(res) <- c("Prior mean","Prior std dev.","Lower prior HPD","Upper prior HPD",
    "Post. mean","Post. std dev.","Lower post. HPD","Upper post. HPD","D-stat")
  if(identical(plotDir,NULL)==F) {
    pdf(paste(plotDir,"/density%03d.pdf",sep=""),width=20,height=15,onefile=F)
    par(mfrow=c(3,3))
    }
  auxdel <- c()
  for(i in 1:length(nomi)) {
    idat <- mydat[,nomi[i]]
    if(nomi[i] %in% rownames(pMat$beta)) {
      a <- pMat$beta[nomi[i],1]
      b <- pMat$beta[nomi[i],2]
      options(warn=-1)
      ihpd <- hpd(qbeta,shape1=a,shape2=b,conf=p)
      options(warn=0)
      res[i,1] <- a/(a+b)
      res[i,2] <- sqrt(res[i,1]*(1-res[i,1])/(1+a+b))
      res[i,3:4] <- ihpd      
      if(!is.null(plotDir)) {
        z <- density(idat,from=0,to=1)
        p1 <- dbeta(z$x,a,b)
        p2 <- z$y
        options(warn=-1)
        limx <- c(min(min(idat),qbeta(0.005,a,b)),max(max(idat),qbeta(0.995,a,b)))
        options(warn=0)
        limy <- c(0,max(ifelse(a<=1||b<=1,p1[which(p1<Inf)],dbeta((a-1)/(a+b-2),a,b)),p2))
        plot(z,cex.axis=1.2,cex.lab=1.1,xlim=limx,ylim=limy,main="")
        title(main=nomi[i],cex.main=1.6)
        curve(dbeta(x,a,b),add=T,lty=6)
        }
      } else if(nomi[i] %in% rownames(pMat$gamma)) {
      a <- pMat$gamma[nomi[i],1]
      b <- pMat$gamma[nomi[i],2]
      ihpd <- hpd(qgamma,shape=a,rate=b,conf=p)
      res[i,1] <- a/b
      res[i,2] <- sqrt(a)/b
      res[i,3:4] <- ihpd  
      if(!is.null(plotDir)) {
        z <- density(idat,from=0)
        p1 <- dgamma(z$x,a,b)
        p2 <- z$y
        limx <- c(min(min(idat),qgamma(0.005,a,b)),max(max(idat),qgamma(0.995,a,b)))
        limy <- c(0,max(ifelse(a<=1,p1[which(p1<Inf)],dgamma((a-1)/b,a,b)),p2))
        plot(z,cex.axis=1.2,cex.lab=1.1,xlim=limx,ylim=limy,main="")
        title(main=nomi[i],cex.main=1.6)
        curve(dgamma(x,a,b),add=T,lty=6)
        }
      } else if(nomi[i] %in% rownames(pMat$gaussian)) {
      a <- pMat$gaussian[nomi[i],1]
      b <- 1/pMat$gaussian[nomi[i],2]^2
      ihpd <- hpd(qnorm,mean=a,sd=b,conf=p)
      res[i,1] <- a
      res[i,2] <- b
      res[i,3:4] <- ihpd        
      if(!is.null(plotDir)) {
        z <- density(idat)
        p1 <- dnorm(z$x,a,b)
        p2 <- z$y
        limx <- c(min(min(idat),qnorm(0.005,a,b)),max(max(idat),qnorm(0.995,a,b)))
        limy <- c(0,max(dnorm(a,a,b),p2))
        plot(z,cex.axis=1.2,cex.lab=1.1,xlim=limx,ylim=limy,main="")
        title(main=nomi[i],cex.main=1.6)
        curve(dnorm(x,a,b),add=T,lty=6)
        }
      } else if(nomi[i] %in% unlist(dirnam)) {
      for(j in 1:length(pMat$dirichlet)) {
        if(nomi[[i]] %in% names(pMat$dirichlet[[j]])) {
          alp <- pMat$dirichlet[[j]]
          auxind <- which(names(alp)==nomi[i])
          }
        }
      if(auxind==1) auxdel <- c(auxdel,nomi[i])
      a <- alp[auxind]
      b <- sum(alp[-auxind])
      options(warn=-1)
      ihpd <- hpd(qbeta,shape1=a,shape2=b,conf=p)
      options(warn=0)
      res[i,1] <- a/(a+b)
      res[i,2] <- sqrt(res[i,1]*(1-res[i,1])/(1+a+b))
      res[i,3:4] <- ihpd
      if(!is.null(plotDir)) {
        z <- density(idat,from=0,to=1)
        p1 <- dbeta(z$x,a,b)
        p2 <- z$y
        options(warn=-1)
        limx <- c(min(min(idat),qbeta(0.005,a,b)),max(max(idat),qbeta(0.995,a,b)))
        options(warn=0)
        limy <- c(0,max(ifelse(a<=1||b<=1,p1[which(p1<Inf)],dbeta((a-1)/(a+b-2),a,b)),p2))
        plot(z,cex.axis=1.2,cex.lab=1.1,xlim=limx,ylim=limy,main="")
        title(main=nomi[i],cex.main=1.6)
        curve(dbeta(x,a,b),add=T,lty=6)
        }
      } else {
      auxdel <- c(auxdel,nomi[i])
      }
    res[i,5] <- mean(idat)
    res[i,6] <- sd(idat)
    if(nomi[i] %in% auxdel) {
      res[i,7:9] <- NA
      } else {
      res[i,7:8] <- emp.hpd(idat,conf=p)
      res[i,9] <- length(which(idat>=ihpd[1]&idat<=ihpd[2]))/nrow(mydat)
      }
    }
  if(!is.null(plotDir)) {
    par(mfrow=c(1,1))
    dev.off()
    }
  round(res[setdiff(nomi,auxdel),],digits)
  }

# compute posterior means
mcmcMean <- function(sim) {
  if(class(sim)!="mcmc.list") stop("argument 'sim' must be an object of class 'mcmc.list'")
  simdat <- c()
  for(i in 1:length(sim)) {
    simdat <- rbind(simdat,sim[[i]])
    }
  nomi <- colnames(sim[[1]])
  auxnam <- auxnum <- c()
  for(i in 1:length(nomi)) {
    if(grepl("\\[",nomi[i])) {
      if(substr(nomi[i],1,2)=="p.") {
        auxnam[i] <- paste("p.",strBW(nomi[i],"p\\.","\\["),sep="")
        } else {
        auxnam[i] <- paste("eta.",strBW(nomi[i],"eta\\.","\\["),sep="")
        }
      auxnum[i] <- strBW(nomi[i],"\\[","\\]")
      } else {
      auxnam[i] <- nomi[i]
      auxnum[i] <- 1
      }
    }
  allnam <- sort(unique(auxnam))
  res <- list()
  for(i in 1:length(allnam)) {             
    auxind <- which(auxnam==allnam[i])
    iN <- as.numeric(auxnum[auxind])
    if(length(iN)==1) {
      iP <- allnam[i]
      } else {
      iP <- paste(allnam[i],"[",sort(iN),"]",sep="")
      }
    auxdat <- matrix(simdat[,iP],ncol=length(iP))
    res[[i]] <- colMeans(auxdat)
    }
  names(res) <- allnam
  res
  }

# find the dimension of non-stochastic parameters in jags code
findConst <- function(code) {
  code <- cleanCode(code)
  auxFor <- grep("for\\(",code)
  fornam <- fornum <- c()
  for(i in 1:length(auxFor)) {
    fornam[i] <- strBW(code[auxFor[i]],"for\\("," in ")
    fornum[i] <- strBW(code[auxFor[i]],":","\\)")
    }
  auxK <- grep("<-",code)
  parnam <- parnum <- c()
  for(i in 1:length(auxK)) {
    istr <- strsplit(code[auxK[i]],"<-")[[1]][1]
    if(grepl("\\[",istr)) {
      parnam[i] <- strsplit(istr,"\\[")[[1]][1]
      parnum[i] <- strBW(istr,"\\[","\\]")
      } else {
      parnam[i] <- istr
      parnum[i] <- 1
      }
    if(grepl("\\(",istr)) parnam[i] <- strBW(parnam[i],"\\(","\\)")
    }
  indOK <- which(substr(parnam,1,2)!="p." & substr(parnam,1,5)!="alpha" & substr(parnam,1,4)!="eta.")
  allpnam <- sort(unique(parnam[indOK]))
  res <- list()
  for(i in 1:length(allpnam)) {
    auxind <- which(parnam==allpnam[i])
    ires <- c()
    for(j in 1:length(auxind)) {
      ijnum <- strsplit(parnum[auxind[j]],",")[[1]]
      for(w in 1:length(ijnum)) {
        if(suppressWarnings(is.na(as.numeric(ijnum[w])))) {
          ijind <- which(fornam==ijnum[w])
          if(length(ijind)==1) {
            ijnum[w] <- fornum[ijind]
            } else {
            ijnum[w] <- fornum[max(which(auxFor[ijind]<=auxK[auxind]))+1]
            }
          }
        }
      ires <- rbind(ires,unname(ijnum))
      }
    if(nrow(ires)>1) {
      res[[i]] <- c(ires[1,1],max(as.numeric(ires[,2])))
      } else {
      res[[i]] <- c(ires)
      }
    }
  names(res) <- allpnam
  res
  }

# find the number of values for each variable in jags code
findOmega <- function(code) {
  auxvar <- findVar(code)
  res <- auxnam <- c()
  if(length(auxvar$continuous)>0) {
    res <- c(res,rep(0,length(auxvar$continuous)))
    auxnam <- c(auxnam,paste("tilde.",auxvar$continuous,sep=""))
    }
  if(length(auxvar$binary)>0) {
    res <- c(res,rep(2,length(auxvar$binary)))
    auxnam <- c(auxnam,auxvar$binary)
    }
  if(length(auxvar$multivalued)>0) {
    dcatnam <- auxvar$multivalued
    auxnam <- c(auxnam,paste("z.",dcatnam,sep=""))
    auxconst <- findConst(code)
    for(i in 1:length(dcatnam)) {
      idim <- auxconst[[paste("theta.",dcatnam[i],sep="")]]
      res <- c(res,as.numeric(idim[2]))
      }
    }
  names(res) <- auxnam
  res
  }

# find parent sets
findParsets <- function(code) {
  auxV <- findVar(code)
  netvar <- c()
  if(length(auxV$continuous)>0) netvar <- c(netvar,paste("tilde.",auxV$continuous,sep=""))
  if(length(auxV$multivalued)>0) netvar <- c(netvar,paste("z.",auxV$multivalued,sep=""))
  if(length(auxV$binary)>0) netvar <- c(netvar,auxV$binary)
  codcl <- cleanCode(code)
  res <- list()
  for(i in 1:length(netvar)) res[[i]] <- character(0)
  names(res) <- netvar
  auxind <- c(grep("^delta\\.",codcl),grep("^wu\\.",codcl),grep("^p0\\.",codcl),
    grep("^pStar.\\.",codcl),grep("^logit\\(relmu\\.",codcl),grep("^logit\\(theta\\.",codcl),
    grep("^theta\\.",codcl),grep("^odds\\.",codcl))
  for(i in 1:length(auxind)) {
    if(grepl("^delta\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"delta\\.")[[1]][2]
      } else if(grepl("^wu\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"wu\\.")[[1]][2]
      } else if(grepl("^p0\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"p0\\.")[[1]][2]
      } else if(grepl("^pStar.\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"pStar")[[1]][2],"\\.")[[1]][2]
      } else if(grepl("^logit\\(relmu\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"logit\\(relmu\\.")[[1]][2]
      } else if(grepl("^logit\\(theta\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"logit\\(theta\\.")[[1]][2]
      } else if(grepl("^theta\\.",codcl[auxind[i]])) {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"theta\\.")[[1]][2]
      } else {
      inam <- strsplit(strsplit(codcl[auxind[i]],"\\[")[[1]][1],"odds\\.")[[1]][2]
      }
    if(inam %in% auxV$continuous) inam <- paste("tilde.",inam,sep="")
    if(inam %in% auxV$multivalued) inam <- paste("z.",inam,sep="")
    ipstr <- strsplit(gsub("[^A-Za-z0-9._ ]"," ",codcl[auxind[i]])," ")[[1]]
    res[[inam]] <- sort(unique(c(res[[inam]],intersect(ipstr,netvar))))
    }
  res[sort(netvar)]
  }

# find a topological order from the list of parent sets
topolog <- function(parSet) {
  nomi <- names(parSet)
  L <- c()
  S <- nomi[which(sapply(parSet,length)==0)]
  while(length(S)>0) {
    xaux <- S[1]
    S <- setdiff(S,xaux)
    L <- c(L,xaux)
    sch <- c()
    for(j in 1:length(parSet)) {
      if(xaux %in% parSet[[j]]) sch <- c(sch,nomi[j])
      }
    if(length(sch)>0) {
      for(j in 1:length(sch)) {
        parSet[[sch[j]]] <- setdiff(parSet[[sch[j]]],xaux)
        if(length(parSet[[sch[j]]])==0) S <- c(S,sch[j])  
        }
      }
    }
  if(sum(sapply(parSet,length))==0) L else NULL
  }
  
# export to GeNIe
genieExport <- function(code,sim,nclass=3,digits=2,node.labels=NULL,state.labels=NULL,file.name=NULL,save.dir=NULL) {
  param.val <- mcmcMean(sim)
  if(is.null(file.name)) file.name <- deparse(substitute(code))
  if(is.null(save.dir)) save.dir <- getwd()
  auxvar <- findVar(code)
  varnames <- unlist(auxvar)
  auxPS <- findParsets(code)
  parsets <- lapply(auxPS,function(z){gsub("z\\.","",gsub("tilde\\.","",z))})
  names(parsets) <- gsub("z\\.","",gsub("tilde\\.","",names(auxPS)))
  sortnames <- topolog(parsets)
  if(identical(sort(sortnames),sort(varnames))==F)
  varnames <- sortnames
  if(is.null(node.labels)) {
    node.labels <- varnames
    names(node.labels) <- varnames
    } else {
    if(is.null(names(node.labels))) names(node.labels) <- names(parsets)
    }
  xcut <- seq(0,1,length=nclass*3+1)
  xval <- c()
  for(i in 2:length(xcut)) {
    xval <- c(xval,c(xcut[i-1]+xcut[i])/2)
    }
  #####
  auxIni <- findConst(code)
  res <- cleanCode(code)
  auxbl1 <- which(grepl("\\{",res))
  auxbl2 <- which(grepl("\\}",res))
  auxOp <- c()
  for(i in 1:length(auxbl2)) {
    iaux <- setdiff(auxbl1,auxOp)
    auxOp[i] <- max(iaux[which(iaux<auxbl2[i])])
   }
  blMat <- cbind(auxOp,auxbl2)            
  auxmain <- which(diff(blMat[,1])<0)+1
  indOK <- blMat[auxmain[1],]
  blMat <- blMat[setdiff(1:nrow(blMat),auxmain),]
  auxdel <- c()
  if(nrow(blMat)>0) {
    for(i in 1:nrow(blMat)) {
      if(blMat[i,2]-blMat[i,1]>0) {
        res[blMat[i,1]] <- paste(res[blMat[i,1]:blMat[i,2]],collapse=";")
        auxdel <- c(auxdel,(blMat[i,1]+1):blMat[i,2])
        }
      }
    }
  mystr <- res[setdiff(c(indOK[1],(indOK[2]-1):(indOK[1]+1),indOK[2]),auxdel)]
  mystr <- gsub("\\|\\|","\\|",mystr)
  mystr <- gsub("\\&\\&","\\&",mystr)
  auxlog <- grep("logit\\(",mystr)
  if(length(auxlog)>0) {
    for(i in 1:length(auxlog)) {  
      inam <- strBW(mystr[auxlog[i]],"logit\\(","\\)")
      istr <- strsplit(mystr[auxlog[i]],"<-")[[1]][2]
      mystr[auxlog[i]] <- paste(inam,"<-plogis(",istr,",0,1)",sep="")
      }
    }
  getBlocks <- function(z) {
    auxind <- which(grepl("^\\s*$",z))
    out <- list()
    for(i in 1:(length(auxind)-1)) {
      out[[i]] <- z[(auxind[i]):(auxind[i+1]-1)]
      }
    out
    } 
  codebl <- getBlocks(mystr) 
  nambl <- c()
  for(i in 1:length(codebl)) {
    iind <- grep("~",codebl[[i]])             
    if(length(iind)==1) {
      inam <- strsplit(gsub(" ","",codebl[[i]][iind]),"~")[[1]][1]
      nambl[i] <- strsplit(inam,"\\[")[[1]][1]
      }
    }
  names(codebl) <- nambl
  #####
  auxnames <- auxdel <- c()
  for(i in 1:length(varnames)) {
    if(varnames[i] %in% auxvar$continuous) {
      auxnames[i] <- paste("tilde.",varnames[i],sep="")
      } else if(varnames[i] %in% auxvar$multivalued) {
      auxnames[i] <- paste("z.",varnames[i],sep="")
      } else {
      auxnames[i] <- varnames[i]      
      }
    if((auxnames[i] %in% nambl)==F) auxdel <- c(auxdel,i)
    }
  names(auxnames) <- varnames
  auxnames <- auxnames[setdiff(1:length(auxnames),auxdel)]
  varnames <- varnames[setdiff(1:length(varnames),auxdel)]
  #####
  strIni <- c()
  if(length(auxIni)>0) {
    for(i in 1:length(auxIni)) {
      if(length(auxIni[[i]])==1) {
        strIni <- c(strIni,paste(names(auxIni)[[i]]," <- array(dim=",auxIni[[i]],")",sep=""))
        } else {
        strIni <- c(strIni,paste(names(auxIni)[[i]]," <- array(dim=c(",paste(auxIni[[i]],collapse=","),"))",sep=""))
        }
      }
    }
  rstr <- paste("n<-1;",paste(strIni,collapse=";"),";",paste(paste(names(param.val),"<-",param.val,sep=""),collapse=";"),sep="")
  eval(parse(text=rstr))
  forind <- strBW(mystr[1],"for\\("," in ")
  maxnamN <- max(sapply(varnames,nchar))
  auxlev <- findOmega(code)
  if(is.null(state.labels)|!is.list(state.labels)) {
    state.labels <- list()
    for(i in 1:length(varnames)) {
      state.labels[[i]] <- 0:(auxlev[auxnames[varnames[i]]]-1)
      }
    names(state.labels) <- varnames
    } else {
    if(is.null(names(node.labels))) names(node.labels) <- names(parsets)
    if(length(auxvar$continuous)>0) {
      for(i in 1:length(auxvar$continuous)) {
        state.labels[[auxvar$continuous[i]]] <- c(paste("lp",1:nclass,sep=""),paste("n",1:nclass,sep=""),paste("hp",1:nclass,sep=""))
        }
      }
    }
  geniecod <- paste('<?xml version="1.0" encoding="ISO-8859-1"?>',"\n",
    '<smile version="1.0" id="',file.name,'" numsamples="100000" discsamples="10000">',"\n",
    "<nodes>","\n",sep="")
  betatab <- function(cutp,relmu,tau) {
    p <- c()
    p[1] <- pbeta(cutp[1],relmu*tau,(1-relmu)*tau)
    for(i in 2:length(cutp)) {
      p[i] <- pbeta(cutp[i],relmu*tau,(1-relmu)*tau)-pbeta(cutp[i-1],relmu*tau,(1-relmu)*tau)
      }
    p[length(cutp)+1] <- 1-pbeta(rev(cutp)[1],relmu*tau,(1-relmu)*tau)
    p
    }
  CPT <- list()
  for(i in 1:length(varnames)) {
    geniecod <- paste(geniecod,'<cpt id="',varnames[i],'">',"\n",sep="")
    ircode <- codebl[[auxnames[varnames[i]]]]
    iauxtil <- grep("~",ircode)
    ircode <- ircode[setdiff(1:length(ircode),iauxtil)]
    if(auxlev[auxnames[varnames[i]]]==0) {
      ilev <- xval
      } else if(auxlev[auxnames[varnames[i]]]==2) {
      ilev <- c(0,1)
      } else {
      ilev <- 1:auxlev[auxnames[varnames[i]]]
      }
    if(varnames[i] %in% auxvar$continuous) {
      if(varnames[i] %in% names(cscales)) {
        iocut <- jagsRescale.inv(xcut,cscales[[varnames[i]]])
        iauxdup <- duplicated(iocut)
        if(sum(iauxdup)>0) {
          idx <- iocut[which(iauxdup==T)[1]]
          idw <- which(iocut==idx)
          if(identical(iocut,sort(iocut))) {
            iocut[idw] <- idx+(iocut[max(idw)+1]-idx)/length(idw)*(0:(length(idw)-1))
            } else {
            iocut[idw] <- iocut[max(idw)+1]+(idx-iocut[max(idw)+1])/length(idw)*(length(idw):1)
            }
          }
        iocut <- round(iocut,digits)
        if(identical(iocut,sort(iocut))) {
          iclab <- paste("less_than_",iocut[2],sep="")
          } else {
          iclab <- paste("more_than_",iocut[2],sep="")
          }
        for(j in 3:(length(iocut)-1)) {
          iclab <- c(iclab,paste("from_",iocut[j-1],"_to_",iocut[j],sep=""))
          }
        if(identical(iocut,sort(iocut))) {
          iclab <- c(iclab,paste("more_than_",iocut[length(iocut)-1],sep=""))
          } else {
          iclab <- c(iclab,paste("less_than_",iocut[length(iocut)-1],sep=""))
          }
        iclab <- gsub("\\.","\\_",iclab)
        } else {
        iclab <- c(paste("lp",1:nclass,sep=""),paste("n",1:nclass,sep=""),paste("hp",1:nclass,sep=""))
        }
      } else {
      iclab <- state.labels[[varnames[i]]]
      }
    for(j in 1:length(iclab)) {    
      geniecod <- paste(geniecod,'<state id="',iclab[j],'" />',"\n",sep="")
      }
    ipar <- parsets[[varnames[i]]]
    if(length(ipar)>0) {
      geniecod <- paste(geniecod,'<parents>',paste(ipar,collapse=" "),'</parents>',"\n",sep="")
      auxgrid <- list()
      for(j in 1:length(ipar)) {
        ijlev <- auxlev[auxnames[ipar[j]]]
        if(ijlev==0) {
          auxgrid[[j]] <- xval
          } else if(ijlev==2) {
          auxgrid[[j]] <- c(0,1)
          } else {
          auxgrid[[j]] <- 1:ijlev
          }
        }
      names(auxgrid) <- auxnames[ipar]
      icnfg <- as.matrix(expand.grid(auxgrid[length(auxgrid):1])[,length(auxgrid):1])
      colnames(icnfg) <- names(auxgrid)
      icpt <- matrix(nrow=nrow(icnfg),ncol=length(ilev))
      colnames(icpt) <- iclab
      for(j in 1:nrow(icnfg)) {
        rstr <- paste(paste(paste(colnames(icnfg),"<-c()",sep=""),collapse=";"),";",mystr[1],
          paste(paste(colnames(icnfg),"[",forind,"]<-",icnfg[j,],sep=""),collapse=";"),";",paste(ircode,collapse=";"),"}",sep="")
        eval(parse(text=rstr))                           
        if(varnames[i] %in% auxvar$continuous) {                      
          rstr <- paste("ijmu <- relmu.",varnames[i],sep="")
          eval(parse(text=rstr))
          ijtau <- rev(param.val[[paste("p.",varnames[i],sep="")]])[1]
          icpt[j,] <- betatab(xcut[2:(length(xcut)-1)],ijmu,ijtau)
          } else if(varnames[i] %in% auxvar$multivalued) {
          rstr <- paste("icpt[",j,",] <- theta.",varnames[i],sep="")
          eval(parse(text=rstr))
          } else {
          rstr <- paste("ijpr <- theta.",varnames[i],sep="")
          eval(parse(text=rstr))
          icpt[j,] <- c(1-ijpr,ijpr)
          }
        }
      idimnam <- state.labels[c(varnames[i],ipar)]
      CPT[[i]] <- array(c(t(icpt)),dim=sapply(idimnam,length))
      dimnames(CPT[[i]]) <- idimnam
      } else {
      ipval <- param.val[[paste("p.",varnames[i],sep="")]]
      if(varnames[i] %in% auxvar$continuous) {
        icpt <- betatab(xcut[2:(length(xcut)-1)],ipval[1],rev(ipval)[1])
        } else if(varnames[i] %in% auxvar$multivalued) {        
        icpt <- ipval
        } else {
        icpt <- c(1-ipval[1],ipval[1])
        }
      idimnam <- state.labels[varnames[i]]
      CPT[[i]] <- array(icpt,dim=sapply(idimnam,length))
      dimnames(CPT[[i]]) <- idimnam
      }
    if(varnames[i] %in% auxvar$continuous) {
      if(varnames[i] %in% names(cscales)) {
        geniecod <- paste(geniecod,'<property id="cutpoints">',paste(iocut,collapse=";"),'</property>',"\n",sep="")
        } else {
        geniecod <- paste(geniecod,'<property id="cutpoints">',paste(xcut,collapse=";"),'</property>',"\n",sep="")
        }
      } else {
      geniecod <- paste(geniecod,'<property id="cutpoints"></property>',"\n",sep="")
      }
    geniecod <- paste(geniecod,'<probabilities>',paste(c(t(icpt)),collapse=" "),'</probabilities>',"\n","</cpt>","\n",sep="")
    }
  geniecod <- paste(geniecod,"</nodes>","\n","<extensions>","\n",
	  '<genie version="1.0" app="GeNIe 2.0.4843.0" name="',file.name,'" faultnameformat="nodestate">',"\n",sep="")  
  auxpos <- 1
  for(i in 2:length(varnames)) {
    if(length(intersect(varnames[which(auxpos==auxpos[i-1])],parsets[[varnames[i]]]))>0) {     
      auxpos[i] <- auxpos[i-1]+1
      } else {
      auxpos[i] <- auxpos[i-1]
      }
    }
  laydim <- as.numeric(table(auxpos))                                       
  laypos <- c()
  for(i in 1:length(laydim)) {
    laypos <- c(laypos,(1:laydim[i])*max(laydim)/(laydim[i]+1))  
    }
  names(laypos) <- varnames
  for(i in 1:length(varnames)) {
    w1 <- 60
    w2 <- 26+2*max(auxlev)
    ix <- laypos[varnames[i]]
    if(i==1) {
      iy <- 1
      } else {
      if(length(intersect(varnames[which(auxpos==auxpos[i-1])],parsets[[varnames[i]]]))>0) iy <- iy+1
      }
    ipos <- round(c(5*w1*ix+w2/2,5*w2*iy+w1/2,5*w1*ix-w2/2,5*w2*iy-w1/2))
    geniecod <- paste(geniecod,'<node id="',varnames[i],'">',"\n",
    '<name>',node.labels[varnames[i]],'</name>',"\n",
    '<interior color="e5f6f7" />',"\n",
    '<outline color="0000bb" />',"\n",
    '<font color="000000" name="Arial" size="8" />',"\n",
    '<position>',paste(ipos,collapse=" "),'</position>',"\n",
    '<barchart active="true" />',"\n",'</node>',"\n",sep="") 
    }
  geniecod <- paste(geniecod,"</genie>","\n","</extensions>","\n","</smile>","\n")
  write(geniecod,file=file.path(save.dir,paste(file.name,".xdsl",sep="")))
  names(CPT) <- varnames
  CPT
  }
