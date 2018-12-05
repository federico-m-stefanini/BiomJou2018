model {
                        
  for(ind in 1:n) {

    # Confusion (binary)
    Confus[ind] ~ dbern(theta.Confus[ind])
    theta.Confus[ind] <- (delta.Confus[ind]<1)*thetaStar.Confus[ind]+(delta.Confus[ind]>=1)*p.Confus[1]
    logit(thetaStar.Confus[ind]) <- eta.Confus[1]*delta.Confus[ind]+
      qlogis(p.Confus[2],0,1)*(z.ConvPth[ind]==2)+
      qlogis(p.Confus[3],0,1)*(z.ConvPth[ind]==3)+
      qlogis(p.Confus[4],0,1)*(z.GCS[ind]==2)+
      qlogis(p.Confus[5],0,1)*(z.GCS[ind]==3)+
      qlogis(p.Confus[6],0,1)*(z.GCS[ind]==4) 
    delta.Confus[ind] <- 1-(z.ConvPth[ind]>1)-(z.GCS[ind]>1)

    # Fall-down (binary)
    Caduta[ind] ~ dbern(theta.Caduta[ind])
    theta.Caduta[ind] <- (delta.Caduta[ind]<1)*thetaStar.Caduta[ind]+(delta.Caduta[ind]>=1)*p.Caduta[1]
    logit(thetaStar.Caduta[ind]) <- eta.Caduta[1]*delta.Caduta[ind]+
      qlogis(p.Caduta[2],0,1)*PerdCosc[ind]+
      qlogis(p.Caduta[3],0,1)*IschCer[ind]
    delta.Caduta[ind] <- 1-PerdCosc[ind]-IschCer[ind]
    
    # auxiliar node
    SaO2_J[ind] ~ dbern(theta.SaO2_J[ind])
    theta.SaO2_J[ind] <- 1*(tilde.SaO2[ind]>=2/3)
    
    # paO2 (continuous, hyper-restricted)
    tilde.O2[ind] ~ dbeta(relmu.O2[ind]*p.O2[3],(1-relmu.O2[ind])*p.O2[3]) T(tol,1-tol)
    logit(relmu.O2[ind]) <- qlogis(p.O2[1],0,1)*delta.O2[ind]+
      qlogis(p.O2[2],0,1)*(tilde.SaO2[ind]>=1/3)*(-1.5+3*tilde.SaO2[ind])
    delta.O2[ind] <- 1-(tilde.SaO2[ind]>=1/3)*(-1.5+3*tilde.SaO2[ind])

    # Glaskow coma score (ordinal, hyper-restricted)
    z.GCS[ind] ~ dcat(theta.GCS[ind,1:4])      
    theta.GCS[ind,1] <- pbeta(2/3,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6])
    theta.GCS[ind,2] <- pbeta(7/9,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6])-pbeta(2/3,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6])
    theta.GCS[ind,3] <- pbeta(8/9,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6])-pbeta(7/9,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6])
    theta.GCS[ind,4] <- 1-pbeta(8/9,relmu.GCS[ind]*p.GCS[6],(1-relmu.GCS[ind])*p.GCS[6]) 
    logit(relmu.GCS[ind]) <- qlogis(p.GCS[1],0,1)*delta.GCS[ind]+
      qlogis(p.GCS[2],0,1)*Ipoglic[ind]+
      qlogis(p.GCS[3],0,1)*massCer[ind]+
      qlogis(p.GCS[4],0,1)*(z.IpossCerPth[ind]==3)+
      qlogis(p.GCS[5],0,1)*(tilde.SaO2[ind]>=1/3)*(-1.5+3*tilde.SaO2[ind])
    delta.GCS[ind] <- 1-Ipoglic[ind]-massCer[ind]-(z.IpossCerPth[ind]==3)-(tilde.SaO2[ind]>=1/3)*(-1.5+3*tilde.SaO2[ind])

    # Central cyanosis (binary)
    Cia[ind] ~ dbern(theta.Cia[ind])
    theta.Cia[ind] <- (delta.Cia[ind]<1)*thetaStar.Cia[ind]+(delta.Cia[ind]>=1)*p.Cia[1]
    logit(thetaStar.Cia[ind]) <- eta.Cia[1]*delta.Cia[ind]+
      qlogis(p.Cia[2],0,1)*(tilde.SaO2[ind]>=2/3)
    delta.Cia[ind] <- 1-(tilde.SaO2[ind]>=2/3)
        
    # pH (continuous)
    tilde.pH[ind] ~ dbeta(relmu.pH[ind]*p.pH[7],(1-relmu.pH[ind])*p.pH[7]) T(tol,1-tol)
    logit(relmu.pH[ind]) <- qlogis(p.pH[1],0,1)*delta.pH[ind]+
      qlogis(p.pH[2],0,1)*metabAlk[ind]+
      qlogis(p.pH[3],0,1)*(metabAlk[ind]==0)*(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])+
      qlogis(p.pH[4],0,1)*(metabAlk[ind]==1)*(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])+
      qlogis(p.pH[5],0,1)*(metabAlk[ind]==0)*(-1.5+3*tilde.CO2[ind])+
      qlogis(p.pH[6],0,1)*(metabAlk[ind]==1)*(-1.5+3*tilde.CO2[ind])
    delta.pH[ind] <- 1-metabAlk[ind]-(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])-(-1.5+3*tilde.CO2[ind])
    
    # auxiliar node
    CO2_J[ind] ~ dbern(theta.CO2_J[ind])
    theta.CO2_J[ind] <- 1*(tilde.CO2[ind]<1/3||tilde.CO2[ind]>=2/3)

    # Syncope (binary)
    PerdCosc[ind] ~ dbern(theta.PerdCosc[ind])
    theta.PerdCosc[ind] <- 1*(z.ConvPth[ind]==2)+(z.ConvPth[ind]!=2)*((delta.PerdCosc[ind]<1)*thetaStar.PerdCosc[ind]+(delta.PerdCosc[ind]>=1)*p.PerdCosc[1])
    logit(thetaStar.PerdCosc[ind]) <- eta.PerdCosc[1]*delta.PerdCosc[ind]+
      qlogis(p.PerdCosc[2],0,1)*Ipoglic[ind]+
      qlogis(p.PerdCosc[3],0,1)*(z.IpossCerPth[ind]==2)
    delta.PerdCosc[ind] <- 1-(z.ConvPth[ind]==2)-Ipoglic[ind]-(z.IpossCerPth[ind]==2)
    
    # Tongue bite (binary)
    Mors[ind] ~ dbern(theta.Mors[ind])
    theta.Mors[ind] <- (delta.Mors[ind]<1)*thetaStar.Mors[ind]+(delta.Mors[ind]>=1)*p.Mors[1]
    logit(thetaStar.Mors[ind]) <- eta.Mors[1]*delta.Mors[ind]+
      qlogis(p.Mors[2],0,1)*(z.ConvPth[ind]>1)
    delta.Mors[ind] <- 1-(z.ConvPth[ind]>1)

    # Sphincter incontinence (binary)
    Inc_sfint[ind] ~ dbern(theta.Inc_sfint[ind])
    theta.Inc_sfint[ind] <- (delta.Inc_sfint[ind]<1)*thetaStar.Inc_sfint[ind]+(delta.Inc_sfint[ind]>=1)*p.Inc_sfint[1]
    logit(thetaStar.Inc_sfint[ind]) <- eta.Inc_sfint[1]*delta.Inc_sfint[ind]+
      qlogis(p.Inc_sfint[2],0,1)*(z.ConvPth[ind]>1) 
    delta.Inc_sfint[ind] <- 1-(z.ConvPth[ind]>1)

    # Generalized epileptic seizure (binary, non-stochastic)
    Conv[ind] ~ dbern(theta.Conv[ind])
    theta.Conv[ind] <- 1*(z.ConvPth[ind]==3)
    
    # auxiliar node
    Shu_J[ind] ~ dbern(theta.Shu_J[ind])
    theta.Shu_J[ind] <- 1*(tilde.Shu[ind]>=2/3)
    
    # Oxygen saturation (continuous, hyper-restricted)
    tilde.SaO2[ind] ~ dbeta(relmu.SaO2[ind]*p.SaO2[7],(1-relmu.SaO2[ind])*p.SaO2[7]) T(tol,1-tol)
    logit(relmu.SaO2[ind]) <- qlogis(p.SaO2[1],0,1)*delta.SaO2[ind]+
      qlogis(p.SaO2[2],0,1)*(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])+
      qlogis(p.SaO2[3],0,1)*(tilde.Shu[ind]>=1/3)*(-1.5+3*tilde.Shu[ind])+
      qlogis(p.SaO2[4],0,1)*(tilde.minVent[ind]<0.5)*(-1.5+3*tilde.minVent[ind])+
      qlogis(p.SaO2[5],0,1)*(tilde.minVent[ind]>0.5)*(-1.5+3*tilde.minVent[ind])+
      qlogis(p.SaO2[6],0,1)*(tilde.FiO2[ind]>=1/3)*(-1.5+3*tilde.FiO2[ind])
    delta.SaO2[ind] <- 1-(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])-(tilde.Shu[ind]>=1/3)*(-1.5+3*tilde.Shu[ind])-
      (-1.5+3*tilde.minVent[ind])-(tilde.FiO2[ind]>=1/3)*(-1.5+3*tilde.FiO2[ind])
    
    # Dyspnea (binary)
    Dispnea[ind] ~ dbern(theta.Dispnea[ind])
    theta.Dispnea[ind] <- 1*(tilde.Shu[ind]>=2/3)+(tilde.Shu[ind]<2/3)*((delta.Dispnea[ind]<1)*thetaStar.Dispnea[ind]+(delta.Dispnea[ind]>=1)*p.Dispnea[1])
    logit(thetaStar.Dispnea[ind]) <- eta.Dispnea[1]*delta.Dispnea[ind]+
      qlogis(p.Dispnea[2],0,1)*(tilde.Latt[ind]>=2/3)+
      qlogis(p.Dispnea[3],0,1)*(tilde.pulmPerf[ind]>=2/3)+
      qlogis(p.Dispnea[4],0,1)*acuteAnem[ind]+
      qlogis(p.Dispnea[5],0,1)*(z.Agit[ind]==3)
    delta.Dispnea[ind] <- 1-(tilde.Latt[ind]>=2/3)-(tilde.pulmPerf[ind]>=2/3)-(tilde.Shu[ind]>=2/3)-acuteAnem[ind]-(z.Agit[ind]==3)

    # paCO2 (continuous)
    tilde.CO2[ind] ~ dbeta(relmu.CO2[ind]*p.CO2[5],(1-relmu.CO2[ind])*p.CO2[5]) T(tol,1-tol)
    logit(relmu.CO2[ind]) <- qlogis(p.CO2[1],0,1)*delta.CO2[ind]+
      qlogis(p.CO2[2],0,1)*(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])+
      qlogis(p.CO2[3],0,1)*(tilde.Shu[ind]>=1/3)*(-1.5+3*tilde.Shu[ind])+
      qlogis(p.CO2[4],0,1)*(-1.5+3*tilde.minVent[ind])
    delta.CO2[ind] <- 1-(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])-
      (tilde.Shu[ind]>=1/3)*(-1.5+3*tilde.Shu[ind])-(-1.5+3*tilde.minVent[ind])

    # Lung perfusion scintigraphy (binary)
    S_per_r[ind] ~ dbern(theta.S_per_r[ind])
    theta.S_per_r[ind] <- 1*(Scint_bias[ind]==1)+(Scint_bias[ind]==0)*((delta.S_per_r[ind]<1)*thetaStar.S_per_r[ind]+(delta.S_per_r[ind]>=1)*p.S_per_r[1])
    logit(thetaStar.S_per_r[ind]) <- eta.S_per_r[1]*delta.S_per_r[ind]+
      qlogis(p.S_per_r[2],0,1)*(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])
    delta.S_per_r[ind] <- 1-Scint_bias[ind]-(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])

    # Vomit (binary)
    vomit[ind] ~ dbern(theta.vomit[ind])
    theta.vomit[ind] <- (delta.vomit[ind]<1)*thetaStar.vomit[ind]+(delta.vomit[ind]>=1)*p.vomit[1]
    logit(thetaStar.vomit[ind]) <- eta.vomit[1]*delta.vomit[ind]+
      qlogis(p.vomit[2],0,1)*(z.Dispeps[ind]>1)+
      qlogis(p.vomit[3],0,1)*massCer[ind]   
    delta.vomit[ind] <- 1-(z.Dispeps[ind]>1)-massCer[ind]

    # Previous transient seizure (3 categories)
    z.ConvPth[ind] ~ dcat(theta.ConvPth[ind,1:3]) 
    for(aux in 1:3) {
      theta.ConvPth[ind,aux] <- (delta.ConvPth[ind]<1)*odds.ConvPth[ind,aux]/sum(odds.ConvPth[ind,1:3])+(delta.ConvPth[ind]>=1)*p.ConvPth[aux]
      }
    odds.ConvPth[ind,1] <- 1           
    odds.ConvPth[ind,2] <- exp(eta.ConvPth[1]*delta.ConvPth[ind]+
      log(p.ConvPth[5]/p.ConvPth[4])*Ipoglic[ind]+
      log(p.ConvPth[8]/p.ConvPth[7])*(z.IpossCerPth[ind]==2)+
      log(p.ConvPth[11]/p.ConvPth[10])*(z.IpossCerPth[ind]==3)+
      log(p.ConvPth[14]/p.ConvPth[13])*massCer[ind])      
    odds.ConvPth[ind,3] <- exp(eta.ConvPth[2]*delta.ConvPth[ind]+
      log(p.ConvPth[6]/p.ConvPth[4])*Ipoglic[ind]+
      log(p.ConvPth[9]/p.ConvPth[7])*(z.IpossCerPth[ind]==2)+
      log(p.ConvPth[12]/p.ConvPth[10])*(z.IpossCerPth[ind]==3)+
      log(p.ConvPth[15]/p.ConvPth[13])*massCer[ind])  
    delta.ConvPth[ind] <- 1-Ipoglic[ind]-(z.IpossCerPth[ind]>1)-massCer[ind]

    # Tachypnea (binary)
    Tcp[ind] ~ dbern(theta.Tcp[ind])
    theta.Tcp[ind] <- (delta.Tcp[ind]<1)*thetaStar.Tcp[ind]+(delta.Tcp[ind]>=1)*p.Tcp[1]
    logit(thetaStar.Tcp[ind]) <- eta.Tcp[1]*delta.Tcp[ind]+
      qlogis(p.Tcp[2],0,1)*(tilde.minVent[ind]<1/3)+
      qlogis(p.Tcp[3],0,1)*(tilde.minVent[ind]>=2/3)
    delta.Tcp[ind] <- 1-(tilde.minVent[ind]<1/3||tilde.minVent[ind]>=2/3)
    
    # Pulmonary shunt (continuous, hyper-restricted)
    tilde.Shu[ind] ~ dbeta(relmu.Shu[ind]*p.Shu[8],(1-relmu.Shu[ind])*p.Shu[8]) T(tol,1-tol)
    logit(relmu.Shu[ind]) <- qlogis(p.Shu[1],0,1)*delta.Shu[ind]+
      qlogis(p.Shu[2],0,1)*(z.Ed_polm[ind]==2)+
      qlogis(p.Shu[3],0,1)*(z.Ed_polm[ind]==3)+
      qlogis(p.Shu[4],0,1)*(z.BPCOa[ind]==2)+
      qlogis(p.Shu[5],0,1)*(z.BPCOa[ind]==4)+
      qlogis(p.Shu[6],0,1)*(z.Pnx[ind]==3)+
      qlogis(p.Shu[7],0,1)*Atel[ind]
    delta.Shu[ind] <- 1-(z.Ed_polm[ind]>1)-(z.BPCOa[ind]==2||z.BPCOa[ind]==4)-(z.Pnx[ind]==3)-Atel[ind]

    # Bias of perfusion scintigraphy (binary)
    Scint_bias[ind] ~ dbern(theta.Scint_bias[ind])
    theta.Scint_bias[ind] <- (delta.Scint_bias[ind]<1)*thetaStar.Scint_bias[ind]+(delta.Scint_bias[ind]>=1)*p.Scint_bias[1]
    logit(thetaStar.Scint_bias[ind]) <- eta.Scint_bias[1]*delta.Scint_bias[ind]+
      qlogis(p.Scint_bias[2],0,1)*(z.Opac[ind]>1)+
      qlogis(p.Scint_bias[3],0,1)*(z.Consol[ind]>1)+
      qlogis(p.Scint_bias[4],0,1)*Atel[ind]+
      qlogis(p.Scint_bias[5],0,1)*(z.V_pl[ind]>1)
    delta.Scint_bias[ind] <- 1-(z.Opac[ind]>1)-(z.Consol[ind]>1)-Atel[ind]-(z.V_pl[ind]>1)

    # Elevated hemidiaphragm (binary)
    E_dia[ind] ~ dbern(theta.E_dia[ind])
    theta.E_dia[ind] <- (delta.E_dia[ind]<1)*thetaStar.E_dia[ind]+(delta.E_dia[ind]>=1)*p.E_dia[1]
    logit(thetaStar.E_dia[ind]) <- eta.E_dia[1]*delta.E_dia[ind]+
      qlogis(p.E_dia[2],0,1)*InP[ind]+
      qlogis(p.E_dia[3],0,1)*Atel[ind]  
    delta.E_dia[ind] <- 1-InP[ind]-Atel[ind]
    
    # Blood pressure (continuous)
    tilde.PA[ind] ~ dbeta(relmu.PA[ind]*p.PA[4],(1-relmu.PA[ind])*p.PA[4]) T(tol,1-tol)
    logit(relmu.PA[ind]) <- qlogis(p.PA[1],0,1)*delta.PA[ind]+
      qlogis(p.PA[2],0,1)*(-1.5+3*tilde.LHO[ind])+
      qlogis(p.PA[3],0,1)*(-1.5+3*tilde.resVas[ind])
    delta.PA[ind] <- 1-(-1.5+3*tilde.LHO[ind])-(-1.5+3*tilde.resVas[ind])

    # Focal neurological signs (binary)
    Neur_foc[ind] ~ dbern(theta.Neur_foc[ind])
    theta.Neur_foc[ind] <- (delta.Neur_foc[ind]<1)*thetaStar.Neur_foc[ind]+(delta.Neur_foc[ind]>=1)*p.Neur_foc[1]
    logit(thetaStar.Neur_foc[ind]) <- eta.Neur_foc[1]*delta.Neur_foc[ind]+
      qlogis(p.Neur_foc[2],0,1)*IschCer[ind]+
      qlogis(p.Neur_foc[3],0,1)*IschCerCron[ind]
    delta.Neur_foc[ind] <- 1-IschCer[ind]-IschCerCron[ind]

    # Cerebral mass (binary)
    massCer[ind] ~ dbern(theta.massCer[ind])
    theta.massCer[ind] <- (delta.massCer[ind]<1)*thetaStar.massCer[ind]+(delta.massCer[ind]>=1)*p.massCer[1]
    logit(thetaStar.massCer[ind]) <- eta.massCer[1]*delta.massCer[ind]+
      qlogis(p.massCer[2],0,1)*IschCer[ind]+
      qlogis(p.massCer[3],0,1)*Neo_o[ind]
    delta.massCer[ind] <- 1-IschCer[ind]-Neo_o[ind]

    # LDH (binary)
    LDH_o[ind] ~ dbern(theta.LDH_o[ind])
    theta.LDH_o[ind] <- (delta.LDH_o[ind]<1)*thetaStar.LDH_o[ind]+(delta.LDH_o[ind]>=1)*p.LDH_o[1]
    logit(thetaStar.LDH_o[ind]) <- eta.LDH_o[1]*delta.LDH_o[ind]+
      qlogis(p.LDH_o[2],0,1)*(z.IMA[ind]>1)+
      qlogis(p.LDH_o[3],0,1)*Neo_o[ind]+
      qlogis(p.LDH_o[4],0,1)*InP[ind]+
      qlogis(p.LDH_o[5],0,1)*Infiam_b[ind]+
      qlogis(p.LDH_o[6],0,1)*IschCer[ind]+
      qlogis(p.LDH_o[7],0,1)*(z.Tiroid[ind]==2||z.Tiroid[ind]==3)     
    delta.LDH_o[ind] <- 1-(z.IMA[ind]>1)-Neo_o[ind]-InP[ind]-Infiam_b[ind]-IschCer[ind]-(z.Tiroid[ind]==2||z.Tiroid[ind]==3)
 
    # Cerebral hypoxia (3 categories)
    z.IpossCerPth[ind] ~ dcat(theta.IpossCerPth[ind,1:3])  
    for(aux in 1:3) {
      theta.IpossCerPth[ind,aux] <- (delta.IpossCerPth[ind]<1)*odds.IpossCerPth[ind,aux]/sum(odds.IpossCerPth[ind,1:3])+(delta.IpossCerPth[ind]>=1)*p.IpossCerPth[aux]
      }
    odds.IpossCerPth[ind,1] <- 1
    odds.IpossCerPth[ind,2] <- exp(eta.IpossCerPth[1]*delta.IpossCerPth[ind]+
      log(p.IpossCerPth[5]/p.IpossCerPth[4])*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])+
      log(p.IpossCerPth[8]/p.IpossCerPth[7])*Ostr_circ_sist[ind]+
      log(p.IpossCerPth[11]/p.IpossCerPth[10])*IschCer[ind]+
      log(p.IpossCerPth[14]/p.IpossCerPth[13])*MNodSA[ind]+
      log(p.IpossCerPth[17]/p.IpossCerPth[16])*CVV[ind]+
      log(p.IpossCerPth[20]/p.IpossCerPth[19])*Arit_sopra_acut[ind]+
      log(p.IpossCerPth[23]/p.IpossCerPth[22])*(1.5-3*tilde.LHO[ind])+
      log(p.IpossCerPth[26]/p.IpossCerPth[25])*acuteAnem[ind]+
      log(p.IpossCerPth[29]/p.IpossCerPth[28])*SuspHeaDrivTmp[ind])      
    odds.IpossCerPth[ind,3] <- exp(eta.IpossCerPth[2]*delta.IpossCerPth[ind]+
      log(p.IpossCerPth[6]/p.IpossCerPth[4])*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])+
      log(p.IpossCerPth[9]/p.IpossCerPth[7])*Ostr_circ_sist[ind]+
      log(p.IpossCerPth[12]/p.IpossCerPth[10])*IschCer[ind]+
      log(p.IpossCerPth[15]/p.IpossCerPth[13])*MNodSA[ind]+
      log(p.IpossCerPth[18]/p.IpossCerPth[16])*CVV[ind]+
      log(p.IpossCerPth[21]/p.IpossCerPth[19])*Arit_sopra_acut[ind]+
      log(p.IpossCerPth[24]/p.IpossCerPth[22])*(1.5-3*tilde.LHO[ind])+
      log(p.IpossCerPth[27]/p.IpossCerPth[25])*acuteAnem[ind]+
      log(p.IpossCerPth[30]/p.IpossCerPth[28])*SuspHeaDrivTmp[ind]) 
    delta.IpossCerPth[ind] <- 1-(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])-Ostr_circ_sist[ind]-IschCer[ind]-MNodSA[ind]-
      CVV[ind]-Arit_sopra_acut[ind]-(1.5-3*tilde.LHO[ind])-acuteAnem[ind]+SuspHeaDrivTmp[ind]
 
    # Palpitations (binary)
    cpl[ind] ~ dbern(theta.cpl[ind])
    theta.cpl[ind] <- (delta.cpl[ind]<1)*thetaStar.cpl[ind]+(delta.cpl[ind]>=1)*p.cpl[1]
    logit(thetaStar.cpl[ind]) <- eta.cpl[1]*delta.cpl[ind]+
      qlogis(p.cpl[2],0,1)*extraSist[ind]+
      qlogis(p.cpl[3],0,1)*(-1.5+3*tilde.FC[ind])
    delta.cpl[ind] <- 1-extraSist[ind]-(-1.5+3*tilde.FC[ind])

    # Shock (binary)
    Shock[ind] ~ dbern(theta.Shock[ind])
    theta.Shock[ind] <- (delta.Shock[ind]<1)*thetaStar.Shock[ind]+(delta.Shock[ind]>=1)*p.Shock[1]
    logit(thetaStar.Shock[ind]) <- eta.Shock[1]*delta.Shock[ind]+
      qlogis(p.Shock[2],0,1)*(tilde.Latt[ind]>=0.75)
    delta.Shock[ind] <- 1-(tilde.Latt[ind]>=0.75)
    
    # Minute ventilation (continuous)
    tilde.minVent[ind] ~ dbeta(relmu.minVent[ind]*p.minVent[12],(1-relmu.minVent[ind])*p.minVent[12]) T(tol,1-tol)
    logit(relmu.minVent[ind]) <- qlogis(p.minVent[1],0,1)*delta.minVent[ind]+
      qlogis(p.minVent[2],0,1)*acuteAnem[ind]+
      qlogis(p.minVent[3],0,1)*(z.Agit[ind]==2)+
      qlogis(p.minVent[4],0,1)*(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])+
      qlogis(p.minVent[5],0,1)*(z.ventTrigg[ind]==2)+
      qlogis(p.minVent[6],0,1)*(z.ventTrigg[ind]==3)+
      qlogis(p.minVent[7],0,1)*Sepsi[ind]+
      qlogis(p.minVent[8],0,1)*(z.BPCOa[ind]==2)+
      qlogis(p.minVent[9],0,1)*(z.BPCOa[ind]==3)+
      qlogis(p.minVent[10],0,1)*(z.BPCOa[ind]==4)+
      qlogis(p.minVent[11],0,1)*Neuromusc[ind]
    delta.minVent[ind] <- 1-acuteAnem[ind]-(z.Agit[ind]==2)-(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])-
      (z.ventTrigg[ind]>1)-Sepsi[ind]-(z.BPCOa[ind]>1)-Neuromusc[ind]
    
    # Augmented lactates (according to clinical judgement) (binary)
    Latt_o[ind] ~ dbern(theta.Latt_o[ind])
    theta.Latt_o[ind] <- (delta.Latt_o[ind]<1)*thetaStar.Latt_o[ind]+(delta.Latt_o[ind]>=1)*p.Latt_o[1]
    logit(thetaStar.Latt_o[ind]) <- eta.Latt_o[1]*delta.Latt_o[ind]+
      qlogis(p.Latt_o[2],0,1)*(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])
    delta.Latt_o[ind] <- 1-(tilde.Latt[ind]>=1/3)*(-1.5+3*tilde.Latt[ind])

    # Opacity to chest X-rays (3 categories)
    z.OpacRX[ind] ~ dcat(theta.OpacRX[ind,1:3]) 
    for(aux in 1:3) {
      theta.OpacRX[ind,aux] <- (delta.OpacRX[ind]<1)*odds.OpacRX[ind,aux]/sum(odds.OpacRX[ind,1:3])+(delta.OpacRX[ind]>=1)*p.OpacRX[aux]
      }
    odds.OpacRX[ind,1] <- 1
    odds.OpacRX[ind,2] <- exp(eta.OpacRX[1]*delta.OpacRX[ind]+
      log(p.OpacRX[5]/p.OpacRX[4])*(z.Opac[ind]==2)+
      log(p.OpacRX[8]/p.OpacRX[7])*(z.Opac[ind]==3))
    odds.OpacRX[ind,3] <- exp(eta.OpacRX[2]*delta.OpacRX[ind]+
      log(p.OpacRX[6]/p.OpacRX[4])*(z.Opac[ind]==2)+
      log(p.OpacRX[9]/p.OpacRX[7])*(z.Opac[ind]==3))
    delta.OpacRX[ind] <- 1-(z.Opac[ind]>1)

    # Pulmonary interstitium (binary)
    intersPol[ind] ~ dbern(theta.intersPol[ind])    
    theta.intersPol[ind] <- 1*(z.Opac[ind]==2)+(z.Opac[ind]!=2)*((delta.intersPol[ind]<1)*thetaStar.intersPol[ind]+(delta.intersPol[ind]>=1)*p.intersPol[1])
    logit(thetaStar.intersPol[ind]) <- eta.intersPol[1]*delta.intersPol[ind]+
      qlogis(p.intersPol[2],0,1)*(z.Opac[ind]==3)  
    delta.intersPol[ind] <- 1-(z.Opac[ind]>1)

    # Air bronchogram (binary)
    bronAer[ind] ~ dbern(theta.bronAer[ind])
    theta.bronAer[ind] <- (delta.bronAer[ind]<1)*thetaStar.bronAer[ind]+(delta.bronAer[ind]>=1)*p.bronAer[1]
    logit(thetaStar.bronAer[ind]) <- eta.bronAer[1]*delta.bronAer[ind]+
      qlogis(p.bronAer[2],0,1)*(z.Consol[ind]>1)
    delta.bronAer[ind] <- 1-(z.Consol[ind]>1)
    
    # Atelactasis (binary)
    Atel[ind] ~ dbern(theta.Atel[ind])
    theta.Atel[ind] <- (delta.Atel[ind]<1)*thetaStar.Atel[ind]+(delta.Atel[ind]>=1)*p.Atel[1]
    logit(thetaStar.Atel[ind]) <- eta.Atel[1]*delta.Atel[ind]+
      qlogis(p.Atel[2],0,1)*ARDS[ind]+
      qlogis(p.Atel[3],0,1)*(z.Consol[ind]>1)+
      qlogis(p.Atel[4],0,1)*Neo_p[ind]+
      qlogis(p.Atel[5],0,1)*(z.Pnx[ind]==2)+
      qlogis(p.Atel[6],0,1)*(z.Pnx[ind]==3)  
    delta.Atel[ind] <- 1-ARDS[ind]-(z.Consol[ind]>1)-Neo_p[ind]-(z.Pnx[ind]>1)
    
    # Orthostatic hypotension (binary)
    Ipo_orto[ind] ~ dbern(theta.Ipo_orto[ind])
    theta.Ipo_orto[ind] <- (delta.Ipo_orto[ind]<1)*thetaStar.Ipo_orto[ind]+(delta.Ipo_orto[ind]>=1)*p.Ipo_orto[1]
    logit(thetaStar.Ipo_orto[ind]) <- eta.Ipo_orto[1]*delta.Ipo_orto[ind]+
      qlogis(p.Ipo_orto[2],0,1)*Psico[ind]+
      qlogis(p.Ipo_orto[3],0,1)*LDopa[ind]+
      qlogis(p.Ipo_orto[4],0,1)*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])+
      qlogis(p.Ipo_orto[5],0,1)*(z.EP[ind]==2)+
      qlogis(p.Ipo_orto[6],0,1)*(z.EP[ind]==3)
    delta.Ipo_orto[ind] <- 1-Psico[ind]-LDopa[ind]-(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])-(z.EP[ind]>1)
    
    # Arterial vascular resistance (continuous)
    tilde.resVas[ind] ~ dbeta(relmu.resVas[ind]*p.resVas[7],(1-relmu.resVas[ind])*p.resVas[7]) T(tol,1-tol)
    logit(relmu.resVas[ind]) <- qlogis(p.resVas[1],0,1)*delta.resVas[ind]+
      qlogis(p.resVas[2],0,1)*(z.autNerSys[ind]==3)+
      qlogis(p.resVas[3],0,1)*(z.autNerSys[ind]==4)+
      qlogis(p.resVas[4],0,1)*(z.autNerSys[ind]==6)+
      qlogis(p.resVas[5],0,1)*cronIpert[ind]+
      qlogis(p.resVas[6],0,1)*(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])
    delta.resVas[ind] <- 1-(z.autNerSys[ind]==3||z.autNerSys[ind]==4||z.autNerSys[ind]==6)-cronIpert[ind]-
      (tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])

    # Acute cerebro-vascular disease (binary)
    IschCer[ind] ~ dbern(theta.IschCer[ind])
    theta.IschCer[ind] <- 1*(EmbArte[ind]==1||z.ArtIntraVasCoa[ind]==3||z.ArtIntraVasCoa[ind]==5)+
      (EmbArte[ind]==0&&z.ArtIntraVasCoa[ind]!=3&&z.ArtIntraVasCoa[ind]!=5)*((delta.IschCer[ind]<1)*thetaStar.IschCer[ind]+(delta.IschCer[ind]>=1)*p.IschCer[1])
    logit(thetaStar.IschCer[ind]) <- eta.IschCer[1]*delta.IschCer[ind]+
      qlogis(p.IschCer[2],0,1)*IschCerCron[ind]+
      qlogis(p.IschCer[3],0,1)*(z.Prof[ind]==3)+
      qlogis(p.IschCer[4],0,1)*(z.autNerSys[ind]==4)+
      qlogis(p.IschCer[5],0,1)*(IschCerCron[ind]==1&&z.autNerSys[ind]==4)+
      qlogis(p.IschCer[6],0,1)*(z.autNerSys[ind]==5)+
      qlogis(p.IschCer[7],0,1)*(z.autNerSys[ind]==6)
    delta.IschCer[ind] <- 1-IschCerCron[ind]-(z.Prof[ind]==3)-(IschCerCron[ind]==1&&z.autNerSys[ind]==4)-
      (z.autNerSys[ind]==4||z.autNerSys[ind]==5||z.autNerSys[ind]==6)-(EmbArte[ind]==1||z.ArtIntraVasCoa[ind]==3||z.ArtIntraVasCoa[ind]==5)
    
    # Heart rate (continuous)
    tilde.FC[ind] ~ dbeta(relmu.FC[ind]*p.FC[9],(1-relmu.FC[ind])*p.FC[9]) T(tol,1-tol)
    logit(relmu.FC[ind]) <- qlogis(p.FC[1],0,1)*delta.FC[ind]+
      qlogis(p.FC[2],0,1)*(z.autNerSys[ind]==2)+
      qlogis(p.FC[3],0,1)*(z.autNerSys[ind]==3)+
      qlogis(p.FC[4],0,1)*(z.autNerSys[ind]==5)+
      qlogis(p.FC[5],0,1)*(z.autNerSys[ind]==6)+
      qlogis(p.FC[6],0,1)*(z.BradiTachi[ind]==2)+
      qlogis(p.FC[7],0,1)*(z.BradiTachi[ind]==3)+
      qlogis(p.FC[8],0,1)*(z.BradiTachi[ind]==4)
    delta.FC[ind] <- 1-(z.autNerSys[ind]==2||z.autNerSys[ind]==3||z.autNerSys[ind]==5||z.autNerSys[ind]==6)-(z.BradiTachi[ind]>1)

    # Cardiomegaly (binary)
    carMeg[ind] ~ dbern(theta.carMeg[ind])
    theta.carMeg[ind] <- (delta.carMeg[ind]<1)*thetaStar.carMeg[ind]+(delta.carMeg[ind]>=1)*p.carMeg[1]
    logit(thetaStar.carMeg[ind]) <- eta.carMeg[1]*delta.carMeg[ind]+
      qlogis(p.carMeg[2],0,1)*(z.Cardiomio[ind]==3)+
      qlogis(p.carMeg[3],0,1)*(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])+
      qlogis(p.carMeg[4],0,1)*Vers_peric[ind]
    delta.carMeg[ind] <- 1-(z.Cardiomio[ind]==3)-(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])-Vers_peric[ind]

    # Troponin I (binary)
    Trop_o[ind] ~ dbern(theta.Trop_o[ind])
    theta.Trop_o[ind] <- 1*(z.IMA[ind]>1)+(z.IMA[ind]==1)*((delta.Trop_o[ind]<1)*thetaStar.Trop_o[ind]+(delta.Trop_o[ind]>=1)*p.Trop_o[1])
    logit(thetaStar.Trop_o[ind]) <- eta.Trop_o[1]*delta.Trop_o[ind]+
      qlogis(p.Trop_o[2],0,1)*MioStretch[ind]
    delta.Trop_o[ind] <- 1-(z.IMA[ind]>1)-MioStretch[ind]

    # Oliguria/anuria (binary)
    Olig[ind] ~ dbern(theta.Olig[ind])
    theta.Olig[ind] <- (delta.Olig[ind]<1)*thetaStar.Olig[ind]+(delta.Olig[ind]>=1)*p.Olig[1]
    logit(thetaStar.Olig[ind]) <- eta.Olig[1]*delta.Olig[ind]+
      qlogis(p.Olig[2],0,1)*(1.5-3*tilde.LHO[ind])
    delta.Olig[ind] <- 1-(1.5-3*tilde.LHO[ind])
       
    # Lactates (continuous, hyper-restricted)
    tilde.Latt[ind] ~ dbeta(relmu.Latt[ind]*p.Latt[3],(1-relmu.Latt[ind])*p.Latt[3]) T(tol,1-tol)
    logit(relmu.Latt[ind]) <- qlogis(p.Latt[1],0,1)*delta.Latt[ind]+
      qlogis(p.Latt[2],0,1)*(-1.5+3*tilde.LHO[ind])
    delta.Latt[ind] <- 1-(-1.5+3*tilde.LHO[ind])

    # Abnormal ventilation trigger (3 categories)
    z.ventTrigg[ind] ~ dcat(theta.ventTrigg[ind,1:3])  
    for(aux in 1:3) {
      theta.ventTrigg[ind,aux] <- (delta.ventTrigg[ind]<1)*odds.ventTrigg[ind,aux]/sum(odds.ventTrigg[ind,1:3])+(delta.ventTrigg[ind]>=1)*p.ventTrigg[aux]
      }
    odds.ventTrigg[ind,1] <- 1           
    odds.ventTrigg[ind,2] <- exp(eta.ventTrigg[1]*delta.ventTrigg[ind]+
      log(p.ventTrigg[5]/p.ventTrigg[4])*(z.EP[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[8]/p.ventTrigg[7])*(z.Ed_polm[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[11]/p.ventTrigg[10])*(z.Polm[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[14]/p.ventTrigg[13])*(z.EP[ind]>1&&z.Enf[ind]>1)+
      log(p.ventTrigg[17]/p.ventTrigg[16])*(z.Ed_polm[ind]>1&&z.Enf[ind]>1)+
      log(p.ventTrigg[20]/p.ventTrigg[19])*(z.Polm[ind]>1&&z.Enf[ind]>1))      
    odds.ventTrigg[ind,3] <- exp(eta.ventTrigg[2]*delta.ventTrigg[ind]+
      log(p.ventTrigg[6]/p.ventTrigg[4])*(z.EP[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[9]/p.ventTrigg[7])*(z.Ed_polm[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[12]/p.ventTrigg[10])*(z.Polm[ind]>1&&z.Enf[ind]==1)+
      log(p.ventTrigg[15]/p.ventTrigg[13])*(z.EP[ind]>1&&z.Enf[ind]>1)+
      log(p.ventTrigg[18]/p.ventTrigg[16])*(z.Ed_polm[ind]>1&&z.Enf[ind]>1)+
      log(p.ventTrigg[21]/p.ventTrigg[19])*(z.Polm[ind]>1&&z.Enf[ind]>1))  
    delta.ventTrigg[ind] <- 1-(z.EP[ind]>1)-(z.Ed_polm[ind]>1)-(z.Polm[ind]>1)
 
    # Cough (3 categories)
    z.Toss[ind] ~ dcat(theta.Toss[ind,1:3]) 
    for(aux in 1:3) {
      theta.Toss[ind,aux] <- (delta.Toss[ind]<1)*odds.Toss[ind,aux]/sum(odds.Toss[ind,1:3])+(delta.Toss[ind]>=1)*p.Toss[aux]
      }
    odds.Toss[ind,1] <- 1           
    odds.Toss[ind,2] <- exp(eta.Toss[1]*delta.Toss[ind]+
      log(p.Toss[5]/p.Toss[4])*(z.BPCOa[ind]==2)+
      log(p.Toss[8]/p.Toss[7])*(z.BPCOa[ind]==3)+
      log(p.Toss[11]/p.Toss[10])*(z.BPCOa[ind]==4)+
      log(p.Toss[14]/p.Toss[13])*(z.Ed_polm[ind]==2)+
      log(p.Toss[17]/p.Toss[16])*(z.Ed_polm[ind]==3)+
      log(p.Toss[20]/p.Toss[19])*Ple[ind]+
      log(p.Toss[23]/p.Toss[22])*(z.Polm[ind]>1)+
      log(p.Toss[26]/p.Toss[25])*URI[ind]+
      log(p.Toss[29]/p.Toss[28])*(z.Enf[ind]>1))      
    odds.Toss[ind,3] <- exp(eta.Toss[2]*delta.Toss[ind]+
      log(p.Toss[6]/p.Toss[4])*(z.BPCOa[ind]==2)+
      log(p.Toss[9]/p.Toss[7])*(z.BPCOa[ind]==3)+
      log(p.Toss[12]/p.Toss[10])*(z.BPCOa[ind]==4)+
      log(p.Toss[15]/p.Toss[13])*(z.Ed_polm[ind]==2)+
      log(p.Toss[18]/p.Toss[16])*(z.Ed_polm[ind]==3)+
      log(p.Toss[21]/p.Toss[19])*Ple[ind]+
      log(p.Toss[24]/p.Toss[22])*(z.Polm[ind]>1)+
      log(p.Toss[27]/p.Toss[25])*URI[ind]+
      log(p.Toss[30]/p.Toss[28])*(z.Enf[ind]>1))  
    delta.Toss[ind] <- 1-(z.BPCOa[ind]>1)-(z.Ed_polm[ind]>1)-Ple[ind]-(z.Polm[ind]>1)-URI[ind]-(z.Enf[ind]>1)

    # Orthopnea (binary)
    Ort[ind] ~ dbern(theta.Ort[ind])
    theta.Ort[ind] <- (delta.Ort[ind]<1)*thetaStar.Ort[ind]+(delta.Ort[ind]>=1)*p.Ort[1]
    logit(thetaStar.Ort[ind]) <- eta.Ort[1]*delta.Ort[ind]+
      qlogis(p.Ort[2],0,1)*(z.Ed_polm[ind]==2)+
      qlogis(p.Ort[3],0,1)*(z.Ed_polm[ind]==3)+
      qlogis(p.Ort[4],0,1)*(z.Enf[ind]==2)+
      qlogis(p.Ort[5],0,1)*(z.Enf[ind]==3)+
      qlogis(p.Ort[6],0,1)*(z.Cardiomio[ind]==3)  
    delta.Ort[ind] <- 1-(z.Ed_polm[ind]>1)-(z.Enf[ind]>1)-(z.Cardiomio[ind]==3)

    # Pulmonary opacity (3 categories with age)
    z.Opac[ind] ~ dcat(theta.Opac[ind,1:3]) 
    for(aux in 1:3) {
      theta.Opac[ind,aux] <- (delta.Opac[ind]<1)*odds.Opac[ind,aux]/sum(odds.Opac[ind,1:3])+(delta.Opac[ind]>=1)*p0.Opac[ind,aux]
      }
    odds.Opac[ind,1] <- 1
    odds.Opac[ind,2] <- exp(eta0.Opac[ind,1]*delta.Opac[ind]+
      log(p.Opac[5]/p.Opac[4])*(z.Ed_polm[ind]==2)+
      log(p.Opac[8]/p.Opac[7])*(z.Ed_polm[ind]==3)+
      log(p.Opac[11]/p.Opac[10])*(z.Polm[ind]==2)+
      log(p.Opac[14]/p.Opac[13])*(z.Polm[ind]==3)+
      log(p.Opac[17]/p.Opac[16])*Neo_p[ind]+
      log(p.Opac[20]/p.Opac[19])*(z.Enf[ind]>1)+
      log(p.Opac[23]/p.Opac[22])*Interst[ind])      
    odds.Opac[ind,3] <- exp(eta0.Opac[ind,2]*delta.Opac[ind]+
      log(p.Opac[6]/p.Opac[4])*(z.Ed_polm[ind]==2)+
      log(p.Opac[9]/p.Opac[7])*(z.Ed_polm[ind]==3)+
      log(p.Opac[12]/p.Opac[10])*(z.Polm[ind]==2)+
      log(p.Opac[15]/p.Opac[13])*(z.Polm[ind]==3)+
      log(p.Opac[18]/p.Opac[16])*Neo_p[ind]+
      log(p.Opac[21]/p.Opac[19])*(z.Enf[ind]>1)+
      log(p.Opac[24]/p.Opac[22])*Interst[ind])     
    eta0.Opac[ind,1] <- log(wu.Opac[ind]*exp(eta.Opac[1])/(1+(1-wu.Opac[ind])*(exp(eta.Opac[1])+exp(eta.Opac[2]))))
    eta0.Opac[ind,2] <- log(wu.Opac[ind]*exp(eta.Opac[2])/(1+(1-wu.Opac[ind])*(exp(eta.Opac[1])+exp(eta.Opac[2]))))
    p0.Opac[ind,1] <- 1-wu.Opac[ind]*(1-p.Opac[1])
    p0.Opac[ind,2] <- wu.Opac[ind]*p.Opac[2]
    p0.Opac[ind,3] <- wu.Opac[ind]*p.Opac[3]
    wu.Opac[ind] <- (tilde.Age[ind]<0.86)*tilde.Age[ind]/0.86+(tilde.Age[ind]>=0.86)*(1-tilde.Age[ind])/(1-0.86)
    delta.Opac[ind] <- 1-(z.Ed_polm[ind]>1)-(z.Polm[ind]>1)-Neo_p[ind]-(z.Enf[ind]>1)-Interst[ind]

    # Ground glass (binary)
    grouGl[ind] ~ dbern(theta.grouGl[ind])
    theta.grouGl[ind] <- (delta.grouGl[ind]<1)*thetaStar.grouGl[ind]+(delta.grouGl[ind]>=1)*p.grouGl[1]
    logit(thetaStar.grouGl[ind]) <- eta.grouGl[1]*delta.grouGl[ind]+
      qlogis(p.grouGl[2],0,1)*(z.Ed_polm[ind]==2)+
      qlogis(p.grouGl[3],0,1)*(z.Ed_polm[ind]==3)+
      qlogis(p.grouGl[4],0,1)*(z.Enf[ind]==2)+
      qlogis(p.grouGl[5],0,1)*(z.Enf[ind]==3)
    delta.grouGl[ind] <- 1-(z.Ed_polm[ind]>1)-(z.Enf[ind]>1)
    
    # auxiliar node
    Ed_polm_J[ind] ~ dbern(theta.Ed_polm_J[ind])
    theta.Ed_polm_J[ind] <- 1*(z.Ed_polm[ind]>1)

    # Pulmonary consolidation (3 categories)
    z.Consol[ind] ~ dcat(theta.Consol[ind,1:3])  
    for(aux in 1:3) {
      theta.Consol[ind,aux] <- (delta.Consol[ind]<1)*odds.Consol[ind,aux]/sum(odds.Consol[ind,1:3])+(delta.Consol[ind]>=1)*p.Consol[aux]
      }
    odds.Consol[ind,1] <- 1           
    odds.Consol[ind,2] <- exp(eta.Consol[1]*delta.Consol[ind]+
      log(p.Consol[5]/p.Consol[4])*InP[ind]+
      log(p.Consol[8]/p.Consol[7])*(z.Polm[ind]==3)+
      log(p.Consol[11]/p.Consol[10])*Neo_p[ind]+
      log(p.Consol[14]/p.Consol[13])*(z.Ed_polm[ind]==2)+
      log(p.Consol[17]/p.Consol[16])*(z.Ed_polm[ind]==3))      
    odds.Consol[ind,3] <- exp(eta.Consol[2]*delta.Consol[ind]+
      log(p.Consol[6]/p.Consol[4])*InP[ind]+
      log(p.Consol[9]/p.Consol[7])*(z.Polm[ind]==3)+
      log(p.Consol[12]/p.Consol[10])*Neo_p[ind]+
      log(p.Consol[15]/p.Consol[13])*(z.Ed_polm[ind]==2)+
      log(p.Consol[18]/p.Consol[16])*(z.Ed_polm[ind]==3))    
    delta.Consol[ind] <- 1-InP[ind]-(z.Polm[ind]==3)-Neo_p[ind]-(z.Ed_polm[ind]>1)

    # Psychiatric medication (binary)
    Psico[ind] ~ dbern(theta.Psico[ind])  
    theta.Psico[ind] <- (delta.Psico[ind]<1)*thetaStar.Psico[ind]+(delta.Psico[ind]>=1)*p.Psico[1]
    logit(thetaStar.Psico[ind]) <- eta.Psico[1]*delta.Psico[ind]+
      qlogis(p.Psico[2],0,1)*(z.Agit[ind]>1)            
    delta.Psico[ind] <- 1-(z.Agit[ind]>1)

    # Autonomic nervous system status (6 categories)
    z.autNerSys[ind] ~ dcat(theta.autNerSys[ind,1:6]) 
    for(aux in 1:6) {
      theta.autNerSys[ind,aux] <- (delta.autNerSys[ind]<1)*odds.autNerSys[ind,aux]/sum(odds.autNerSys[ind,1:6])+(delta.autNerSys[ind]>=1)*p.autNerSys[aux]
      }
    odds.autNerSys[ind,1] <- 1                  
    odds.autNerSys[ind,2] <- exp(eta.autNerSys[1]*delta.autNerSys[ind]+
      log(p.autNerSys[8]/p.autNerSys[7])*(z.BPCOa[ind]>1)+
      log(p.autNerSys[14]/p.autNerSys[13])*Coronar[ind]+
      log(p.autNerSys[20]/p.autNerSys[19])*(z.Agit[ind]>1)+
      log(p.autNerSys[26]/p.autNerSys[25])*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      log(p.autNerSys[32]/p.autNerSys[31])*(z.Cardiomio[ind]==2)+
      log(p.autNerSys[38]/p.autNerSys[37])*(z.Cardiomio[ind]==3)+      
      log(p.autNerSys[44]/p.autNerSys[43])*(tilde.leftPump[ind]>=2/3)+
      log(p.autNerSys[50]/p.autNerSys[49])*(z.IMA[ind]>1)+
      log(p.autNerSys[56]/p.autNerSys[55])*Sepsi[ind])
    odds.autNerSys[ind,3] <- exp(eta.autNerSys[2]*delta.autNerSys[ind]+
      log(p.autNerSys[9]/p.autNerSys[7])*(z.BPCOa[ind]>1)+
      log(p.autNerSys[15]/p.autNerSys[13])*Coronar[ind]+
      log(p.autNerSys[21]/p.autNerSys[19])*(z.Agit[ind]>1)+
      log(p.autNerSys[27]/p.autNerSys[25])*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      log(p.autNerSys[33]/p.autNerSys[31])*(z.Cardiomio[ind]==2)+
      log(p.autNerSys[39]/p.autNerSys[37])*(z.Cardiomio[ind]==3)+      
      log(p.autNerSys[45]/p.autNerSys[43])*(tilde.leftPump[ind]>=2/3)+
      log(p.autNerSys[51]/p.autNerSys[49])*(z.IMA[ind]>1)+
      log(p.autNerSys[57]/p.autNerSys[55])*Sepsi[ind])
    odds.autNerSys[ind,4] <- exp(eta.autNerSys[3]*delta.autNerSys[ind]+
      log(p.autNerSys[10]/p.autNerSys[7])*(z.BPCOa[ind]>1)+
      log(p.autNerSys[16]/p.autNerSys[13])*Coronar[ind]+
      log(p.autNerSys[22]/p.autNerSys[19])*(z.Agit[ind]>1)+
      log(p.autNerSys[28]/p.autNerSys[25])*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      log(p.autNerSys[34]/p.autNerSys[31])*(z.Cardiomio[ind]==2)+
      log(p.autNerSys[40]/p.autNerSys[37])*(z.Cardiomio[ind]==3)+      
      log(p.autNerSys[46]/p.autNerSys[43])*(tilde.leftPump[ind]>=2/3)+
      log(p.autNerSys[52]/p.autNerSys[49])*(z.IMA[ind]>1)+
      log(p.autNerSys[58]/p.autNerSys[55])*Sepsi[ind])
    odds.autNerSys[ind,5] <- exp(eta.autNerSys[4]*delta.autNerSys[ind]+
      log(p.autNerSys[11]/p.autNerSys[7])*(z.BPCOa[ind]>1)+
      log(p.autNerSys[17]/p.autNerSys[13])*Coronar[ind]+
      log(p.autNerSys[23]/p.autNerSys[19])*(z.Agit[ind]>1)+
      log(p.autNerSys[29]/p.autNerSys[25])*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      log(p.autNerSys[35]/p.autNerSys[31])*(z.Cardiomio[ind]==2)+
      log(p.autNerSys[41]/p.autNerSys[37])*(z.Cardiomio[ind]==3)+      
      log(p.autNerSys[47]/p.autNerSys[43])*(tilde.leftPump[ind]>=2/3)+
      log(p.autNerSys[53]/p.autNerSys[49])*(z.IMA[ind]>1)+
      log(p.autNerSys[59]/p.autNerSys[55])*Sepsi[ind])
    odds.autNerSys[ind,6] <- exp(eta.autNerSys[5]*delta.autNerSys[ind]+
      log(p.autNerSys[12]/p.autNerSys[7])*(z.BPCOa[ind]>1)+
      log(p.autNerSys[18]/p.autNerSys[13])*Coronar[ind]+
      log(p.autNerSys[24]/p.autNerSys[19])*(z.Agit[ind]>1)+
      log(p.autNerSys[30]/p.autNerSys[25])*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      log(p.autNerSys[36]/p.autNerSys[31])*(z.Cardiomio[ind]==2)+
      log(p.autNerSys[42]/p.autNerSys[37])*(z.Cardiomio[ind]==3)+      
      log(p.autNerSys[48]/p.autNerSys[43])*(tilde.leftPump[ind]>=2/3)+
      log(p.autNerSys[54]/p.autNerSys[49])*(z.IMA[ind]>1)+
      log(p.autNerSys[60]/p.autNerSys[55])*Sepsi[ind])
    delta.autNerSys[ind] <- 1-(z.BPCOa[ind]>1)-Coronar[ind]-(z.Agit[ind]>1)-(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])-
      (z.Cardiomio[ind]>1)-(tilde.leftPump[ind]>=2/3)-(z.IMA[ind]>1)-Sepsi[ind]

    # Anxiety/agitation (according to clinical judgement) (binary)
    Agit_o[ind] ~ dbern(theta.Agit_o[ind])
    theta.Agit_o[ind] <- 1*(z.Agit[ind]==3)+(z.Agit[ind]!=3)*((delta.Agit_o[ind]<1)*thetaStar.Agit_o[ind]+(delta.Agit_o[ind]>=1)*p.Agit_o[1])
    logit(thetaStar.Agit_o[ind]) <- eta.Agit_o[1]*delta.Agit_o[ind]+
      qlogis(p.Agit_o[2],0,1)*(z.Agit[ind]==2)
    delta.Agit_o[ind] <- 1-(z.Agit[ind]>1)

    # Pericardial effusion (binary)
    Vers_peric[ind] ~ dbern(theta.Vers_peric[ind])
    theta.Vers_peric[ind] <- 1*(HemoPeric[ind]==1)+(HemoPeric[ind]==0)*((delta.Vers_peric[ind]<1)*thetaStar.Vers_peric[ind]+(delta.Vers_peric[ind]>=1)*p.Vers_peric[1])
    logit(thetaStar.Vers_peric[ind]) <- eta.Vers_peric[1]*delta.Vers_peric[ind]+
      qlogis(p.Vers_peric[2],0,1)*Peric[ind]+
      qlogis(p.Vers_peric[3],0,1)*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])
    delta.Vers_peric[ind] <- 1-Peric[ind]-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])-HemoPeric[ind]

    # Pleural effusion (3 categories)
    z.V_pl[ind] ~ dcat(theta.V_pl[ind,1:3])
    for(aux in 1:3) {
      theta.V_pl[ind,aux] <- (delta.V_pl[ind]<1)*odds.V_pl[ind,aux]/sum(odds.V_pl[ind,1:3])+(delta.V_pl[ind]>=1)*p.V_pl[aux]
      }
    odds.V_pl[ind,1] <- 1               
    odds.V_pl[ind,2] <- exp(eta.V_pl[1]*delta.V_pl[ind]+
      log(p.V_pl[5]/p.V_pl[4])*InP[ind]+
      log(p.V_pl[8]/p.V_pl[7])*Ple[ind]+
      log(p.V_pl[11]/p.V_pl[10])*Neo_p[ind]+
      log(p.V_pl[14]/p.V_pl[13])*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind]))
    odds.V_pl[ind,3] <- exp(eta.V_pl[2]*delta.V_pl[ind]+
      log(p.V_pl[6]/p.V_pl[4])*InP[ind]+
      log(p.V_pl[9]/p.V_pl[7])*Ple[ind]+
      log(p.V_pl[12]/p.V_pl[10])*Neo_p[ind]+
      log(p.V_pl[15]/p.V_pl[13])*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind]))
    delta.V_pl[ind] <- 1-InP[ind]-Ple[ind]-Neo_p[ind]-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Jugular venous distention (binary)
    Turg_giu[ind] ~ dbern(theta.Turg_giu[ind])
    theta.Turg_giu[ind] <- (delta.Turg_giu[ind]<1)*thetaStar.Turg_giu[ind]+(delta.Turg_giu[ind]>=1)*p.Turg_giu[1]
    logit(thetaStar.Turg_giu[ind]) <- eta.Turg_giu[1]*delta.Turg_giu[ind]+
      qlogis(p.Turg_giu[2],0,1)*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])
    delta.Turg_giu[ind] <- 1-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Reflux of contrast medium into the hepatic veins (binary)
    reflVenSEp[ind] ~ dbern(theta.reflVenSEp[ind])
    theta.reflVenSEp[ind] <- (delta.reflVenSEp[ind]<1)*thetaStar.reflVenSEp[ind]+(delta.reflVenSEp[ind]>=1)*p.reflVenSEp[1]
    logit(thetaStar.reflVenSEp[ind]) <- eta.reflVenSEp[1]*delta.reflVenSEp[ind]+
      qlogis(p.reflVenSEp[2],0,1)*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])
    delta.reflVenSEp[ind] <- 1-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Right ventricular hypokinesis (binary)
    IpoCinVdx[ind] ~ dbern(theta.IpoCinVdx[ind])
    theta.IpoCinVdx[ind] <- (delta.IpoCinVdx[ind]<1)*thetaStar.IpoCinVdx[ind]+(delta.IpoCinVdx[ind]>=1)*p.IpoCinVdx[1]
    logit(thetaStar.IpoCinVdx[ind]) <- eta.IpoCinVdx[1]*delta.IpoCinVdx[ind]+
      qlogis(p.IpoCinVdx[2],0,1)*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])
    delta.IpoCinVdx[ind] <- 1-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Hepatomegaly (binary)
    Epatom[ind] ~ dbern(theta.Epatom[ind])
    theta.Epatom[ind] <- (delta.Epatom[ind]<1)*thetaStar.Epatom[ind]+(delta.Epatom[ind]>=1)*p.Epatom[1]
    logit(thetaStar.Epatom[ind]) <- eta.Epatom[1]*delta.Epatom[ind]+
      qlogis(p.Epatom[2],0,1)*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])
    delta.Epatom[ind] <- 1-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Peripheral edema (3 categories)
    z.Edem_g[ind] ~ dcat(theta.Edem_g[ind,1:3])
    for(aux in 1:3) {
      theta.Edem_g[ind,aux] <- (delta.Edem_g[ind]<1)*odds.Edem_g[ind,aux]/sum(odds.Edem_g[ind,1:3])+(delta.Edem_g[ind]>=1)*p.Edem_g[aux]
      }
    odds.Edem_g[ind,1] <- 1           
    odds.Edem_g[ind,2] <- exp(eta.Edem_g[1]*delta.Edem_g[ind]+
      log(p.Edem_g[5]/p.Edem_g[4])*TVPi[ind]+
      log(p.Edem_g[8]/p.Edem_g[7])*Fr[ind]+
      log(p.Edem_g[11]/p.Edem_g[10])*(z.Cardiomio[ind]==3)+
      log(p.Edem_g[14]/p.Edem_g[13])*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind]))
    odds.Edem_g[ind,3] <- exp(eta.Edem_g[2]*delta.Edem_g[ind]+
      log(p.Edem_g[6]/p.Edem_g[4])*TVPi[ind]+
      log(p.Edem_g[9]/p.Edem_g[7])*Fr[ind]+
      log(p.Edem_g[12]/p.Edem_g[10])*(z.Cardiomio[ind]==3)+
      log(p.Edem_g[15]/p.Edem_g[13])*(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind]))
    delta.Edem_g[ind] <- 1-TVPi[ind]-Fr[ind]-(z.Cardiomio[ind]==3)-(tilde.RHF[ind]>=1/3)*(-1.5+3*tilde.RHF[ind])

    # Myocardial stretching (binary)
    MioStretch[ind] ~ dbern(theta.MioStretch[ind])
    theta.MioStretch[ind] <- (delta.MioStretch[ind]<1)*thetaStar.MioStretch[ind]+(delta.MioStretch[ind]>=1)*p.MioStretch[1]
    logit(thetaStar.MioStretch[ind]) <- eta.MioStretch[1]*delta.MioStretch[ind]+
      qlogis(p.MioStretch[2],0,1)*(-1.5+3*tilde.LHI[ind])+
      qlogis(p.MioStretch[3],0,1)*(z.EP[ind]==2)+
      qlogis(p.MioStretch[4],0,1)*(z.EP[ind]==3)+
      qlogis(p.MioStretch[5],0,1)*Insuf_aorAcut[ind]
    delta.MioStretch[ind] <- 1-(-1.5+3*tilde.LHI[ind])-(z.EP[ind]>1)-Insuf_aorAcut[ind]

    # Left ventricular output (continuous)
    tilde.LHO[ind] ~ dbeta(relmu.LHO[ind]*p.LHO[4],(1-relmu.LHO[ind])*p.LHO[4]) T(tol,1-tol)
    logit(relmu.LHO[ind]) <- qlogis(p.LHO[1],0,1)*delta.LHO[ind]+
      qlogis(p.LHO[2],0,1)*(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])+
      qlogis(p.LHO[3],0,1)*(-1.5+3*tilde.LHI[ind])
    delta.LHO[ind] <- 1-(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])-(-1.5+3*tilde.LHI[ind])
    
    # auxiliar node
    LHI_J[ind] ~ dbern(theta.LHI_J[ind])
    theta.LHI_J[ind] <- 1*(tilde.LHI[ind]>=2/3)

    # Pulmonary edema (3 categories)
    z.Ed_polm[ind] ~ dcat(theta.Ed_polm[ind,1:3])
    for(aux in 1:3) {
      theta.Ed_polm[ind,aux] <- (delta.Ed_polm[ind]<1)*odds.Ed_polm[ind,aux]/sum(odds.Ed_polm[ind,1:3])+(delta.Ed_polm[ind]>=1)*p.Ed_polm[aux]
      }
    odds.Ed_polm[ind,1] <- 1    
    odds.Ed_polm[ind,2] <- exp(eta.Ed_polm[1]*delta.Ed_polm[ind]+
      log(p.Ed_polm[5]/p.Ed_polm[4])*(tilde.LHI[ind]>=2/3)+
      log(p.Ed_polm[8]/p.Ed_polm[7])*ARDS[ind])
    odds.Ed_polm[ind,3] <- exp(eta.Ed_polm[2]*delta.Ed_polm[ind]+
      log(p.Ed_polm[6]/p.Ed_polm[4])*(tilde.LHI[ind]>=2/3)+
      log(p.Ed_polm[9]/p.Ed_polm[7])*ARDS[ind])
    delta.Ed_polm[ind] <- 1-(tilde.LHI[ind]>=2/3)-ARDS[ind]
   
    # Dilated left ventricle (binary)
    DilatVsx[ind] ~ dbern(theta.DilatVsx[ind])
    theta.DilatVsx[ind] <- 1*(Card_dil[ind]==1)+(Card_dil[ind]==0)*((delta.DilatVsx[ind]<1)*thetaStar.DilatVsx[ind]+(delta.DilatVsx[ind]>=1)*p.DilatVsx[1])
    logit(thetaStar.DilatVsx[ind]) <- eta.DilatVsx[1]*delta.DilatVsx[ind]+
      qlogis(p.DilatVsx[2],0,1)*(-1.5+3*tilde.LHI[ind])+
      qlogis(p.DilatVsx[3],0,1)*Miocard[ind]
    delta.DilatVsx[ind] <- 1-Card_dil[ind]-(-1.5+3*tilde.LHI[ind])-Miocard[ind]
   
    # Bronchospasm/reduced vescicolar murmur (3 categories)
    z.EO_bsp[ind] ~ dcat(theta.EO_bsp[ind,1:3])
    for(aux in 1:3) {
      theta.EO_bsp[ind,aux] <- (delta.EO_bsp[ind]<1)*odds.EO_bsp[ind,aux]/sum(odds.EO_bsp[ind,1:3])+(delta.EO_bsp[ind]>=1)*p.EO_bsp[aux]
      }
    odds.EO_bsp[ind,1] <- 1           
    odds.EO_bsp[ind,2] <- exp(eta.EO_bsp[1]*delta.EO_bsp[ind]+
      log(p.EO_bsp[5]/p.EO_bsp[4])*(z.BPCOa[ind]==2||z.BPCOa[ind]==3)+
      log(p.EO_bsp[8]/p.EO_bsp[7])*(z.BPCOa[ind]==4)+
      log(p.EO_bsp[11]/p.EO_bsp[10])*(z.Pnx[ind]==2)+
      log(p.EO_bsp[14]/p.EO_bsp[13])*(z.Pnx[ind]==3)+
      log(p.EO_bsp[17]/p.EO_bsp[16])*(z.EP[ind]==2)+
      log(p.EO_bsp[20]/p.EO_bsp[19])*(z.EP[ind]==3)+
      log(p.EO_bsp[23]/p.EO_bsp[22])*URI[ind]+
      log(p.EO_bsp[26]/p.EO_bsp[25])*(z.Enf[ind]>1))      
    odds.EO_bsp[ind,3] <- exp(eta.EO_bsp[2]*delta.EO_bsp[ind]+
      log(p.EO_bsp[6]/p.EO_bsp[4])*(z.BPCOa[ind]==2||z.BPCOa[ind]==3)+
      log(p.EO_bsp[9]/p.EO_bsp[7])*(z.BPCOa[ind]==4)+
      log(p.EO_bsp[12]/p.EO_bsp[10])*(z.Pnx[ind]==2)+
      log(p.EO_bsp[15]/p.EO_bsp[13])*(z.Pnx[ind]==3)+
      log(p.EO_bsp[18]/p.EO_bsp[16])*(z.EP[ind]==2)+
      log(p.EO_bsp[21]/p.EO_bsp[19])*(z.EP[ind]==3)+
      log(p.EO_bsp[24]/p.EO_bsp[22])*URI[ind]+
      log(p.EO_bsp[27]/p.EO_bsp[25])*(z.Enf[ind]>1))  
    delta.EO_bsp[ind] <- 1-(z.BPCOa[ind]>1)-(z.Pnx[ind]>1)-(z.EP[ind]>1)-URI[ind]-(z.Enf[ind]>1)
    
    # auxiliar node
    BPCOa_J[ind] ~ dbern(theta.BPCOa_J[ind])
    theta.BPCOa_J[ind] <- 1*(z.BPCOa[ind]>1)

    # Anxiety/agitation (3 categories)
    z.Agit[ind] ~ dcat(theta.Agit[ind,1:3])
    for(aux in 1:3) {
      theta.Agit[ind,aux] <- (delta.Agit[ind]<1)*odds.Agit[ind,aux]/sum(odds.Agit[ind,1:3])+(delta.Agit[ind]>=1)*p.Agit[aux]
      }
    odds.Agit[ind,1] <- 1
    odds.Agit[ind,2] <- exp(eta.Agit[1]*delta.Agit[ind]+
      log(p.Agit[5]/p.Agit[4])*ChestPain[ind]+
      log(p.Agit[8]/p.Agit[7])*Ipoglic[ind])
    odds.Agit[ind,3] <- exp(eta.Agit[2]*delta.Agit[ind]+
      log(p.Agit[6]/p.Agit[4])*ChestPain[ind]+
      log(p.Agit[9]/p.Agit[7])*Ipoglic[ind])
    delta.Agit[ind] <- 1-ChestPain[ind]-Ipoglic[ind]

    # Right heart failure (continuous, hyper-restricted)
    tilde.RHF[ind] ~ dbeta(relmu.RHF[ind]*p.RHF[5],(1-relmu.RHF[ind])*p.RHF[5]) T(tol,1-tol)
    logit(relmu.RHF[ind]) <- qlogis(p.RHF[1],0,1)*delta.RHF[ind]+
      qlogis(p.RHF[2],0,1)*(-1.5+3*tilde.leftPump[ind])+
      qlogis(p.RHF[3],0,1)*(-1.5+3*tilde.RHI[ind])+
      qlogis(p.RHF[4],0,1)*(z.tamp[ind]>1)
    delta.RHF[ind] <- 1-(-1.5+3*tilde.leftPump[ind])-(-1.5+3*tilde.RHI[ind])-(z.tamp[ind]>1)

    # Left ventricular pre-load (continuous)
    tilde.LHI[ind] ~ dbeta(relmu.LHI[ind]*p.LHI[7],(1-relmu.LHI[ind])*p.LHI[7]) T(tol,1-tol)
    logit(relmu.LHI[ind]) <- qlogis(p.LHI[1],0,1)*delta.LHI[ind]+
      qlogis(p.LHI[2],0,1)*(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])+
      qlogis(p.LHI[3],0,1)*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      qlogis(p.LHI[4],0,1)*(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])+
      qlogis(p.LHI[5],0,1)*(z.tamp[ind]==2)+
      qlogis(p.LHI[6],0,1)*(z.tamp[ind]==3)
    delta.LHI[ind] <- 1-(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])-(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])-
      (tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])-(z.tamp[ind]>1)
    
    # auxiliar node
    leftPump_J[ind] ~ dbern(theta.leftPump_J[ind])
    theta.leftPump_J[ind] <- 1*(tilde.leftPump[ind]>=2/3)

    # Acute pulmonary disease (4 categories)
    #  [URI reduces the effect of Enf]
    z.BPCOa[ind] ~ dcat(theta.BPCOa[ind,1:4])
    for(aux in 1:4) {
      theta.BPCOa[ind,aux] <- (delta.BPCOa[ind]<1)*odds.BPCOa[ind,aux]/sum(odds.BPCOa[ind,1:4])+(delta.BPCOa[ind]>=1)*p.BPCOa[aux]
      }
    odds.BPCOa[ind,1] <- 1                  
    odds.BPCOa[ind,2] <- exp(eta.BPCOa[1]*delta.BPCOa[ind]+
      log(p.BPCOa[6]/p.BPCOa[5])*(tilde.leftPump[ind]>=2/3)+
      log(p.BPCOa[10]/p.BPCOa[9])*Neo_p[ind]+
      log(p.BPCOa[14]/p.BPCOa[13])*(z.Polm[ind]==2)+
      log(p.BPCOa[18]/p.BPCOa[17])*(z.Polm[ind]==3)+
      log(p.BPCOa[22]/p.BPCOa[21])*(z.Enf[ind]==2&&URI[ind]==0)+
      log(p.BPCOa[26]/p.BPCOa[25])*(z.Enf[ind]==2&&URI[ind]==1)+      
      log(p.BPCOa[30]/p.BPCOa[29])*(z.Enf[ind]==3&&URI[ind]==0)+
      log(p.BPCOa[34]/p.BPCOa[33])*(z.Enf[ind]==3&&URI[ind]==1)+
      log(p.BPCOa[38]/p.BPCOa[37])*Asma[ind]) 
    odds.BPCOa[ind,3] <- exp(eta.BPCOa[2]*delta.BPCOa[ind]+
      log(p.BPCOa[7]/p.BPCOa[5])*(tilde.leftPump[ind]>=2/3)+
      log(p.BPCOa[11]/p.BPCOa[9])*Neo_p[ind]+
      log(p.BPCOa[15]/p.BPCOa[13])*(z.Polm[ind]==2)+
      log(p.BPCOa[19]/p.BPCOa[17])*(z.Polm[ind]==3)+
      log(p.BPCOa[23]/p.BPCOa[21])*(z.Enf[ind]==2&&URI[ind]==0)+
      log(p.BPCOa[27]/p.BPCOa[25])*(z.Enf[ind]==2&&URI[ind]==1)+      
      log(p.BPCOa[31]/p.BPCOa[29])*(z.Enf[ind]==3&&URI[ind]==0)+
      log(p.BPCOa[35]/p.BPCOa[33])*(z.Enf[ind]==3&&URI[ind]==1)+
      log(p.BPCOa[39]/p.BPCOa[37])*Asma[ind]) 
    odds.BPCOa[ind,4] <- exp(eta.BPCOa[3]*delta.BPCOa[ind]+
      log(p.BPCOa[8]/p.BPCOa[5])*(tilde.leftPump[ind]>=2/3)+
      log(p.BPCOa[12]/p.BPCOa[9])*Neo_p[ind]+
      log(p.BPCOa[16]/p.BPCOa[13])*(z.Polm[ind]==2)+
      log(p.BPCOa[20]/p.BPCOa[17])*(z.Polm[ind]==3)+
      log(p.BPCOa[24]/p.BPCOa[21])*(z.Enf[ind]==2&&URI[ind]==0)+
      log(p.BPCOa[28]/p.BPCOa[25])*(z.Enf[ind]==2&&URI[ind]==1)+      
      log(p.BPCOa[32]/p.BPCOa[29])*(z.Enf[ind]==3&&URI[ind]==0)+
      log(p.BPCOa[36]/p.BPCOa[33])*(z.Enf[ind]==3&&URI[ind]==1)+
      log(p.BPCOa[40]/p.BPCOa[37])*Asma[ind])
    delta.BPCOa[ind] <- 1-(tilde.leftPump[ind]>=2/3)-Neo_p[ind]-(z.Polm[ind]>1)-(z.Enf[ind]>1)-Asma[ind]

    # Brain natriuretic peptide (binary)
    BNP_o[ind] ~ dbern(theta.BNP_o[ind])
    theta.BNP_o[ind] <- (delta.BNP_o[ind]<1)*thetaStar.BNP_o[ind]+(delta.BNP_o[ind]>=1)*p.BNP_o[1]
    logit(thetaStar.BNP_o[ind]) <- eta.BNP_o[1]*delta.BNP_o[ind]+
      qlogis(p.BNP_o[2],0,1)*(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])
    delta.BNP_o[ind] <- 1-(tilde.leftPump[ind]>=1/3)*(-1.5+3*tilde.leftPump[ind])

    # Chest pain (binary, semi-stochastic)
    ChestPain[ind] ~ dbern(theta.ChestPain[ind])
    theta.ChestPain[ind] <- p.ChestPain[1]*(z.Dol_type[ind]==1)+1*(z.Dol_type[ind]>1)
    
    # Left heart failure (continuous, hyper-restricted)
    tilde.leftPump[ind] ~ dbeta(relmu.leftPump[ind]*p.leftPump[7],(1-relmu.leftPump[ind])*p.leftPump[7]) T(tol,1-tol)
    logit(relmu.leftPump[ind]) <- qlogis(p.leftPump[1],0,1)*delta.leftPump[ind]+
      qlogis(p.leftPump[2],0,1)*(z.IMA[ind]==3)+
      qlogis(p.leftPump[3],0,1)*(z.IMA[ind]<3&&z.Cardiomio[ind]==3)+
      qlogis(p.leftPump[4],0,1)*(z.IMA[ind]==3&&z.Cardiomio[ind]==3)+
      qlogis(p.leftPump[5],0,1)*(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])+
      qlogis(p.leftPump[6],0,1)*(tilde.postLoad[ind]>=1/3)*(-1.5+3*tilde.postLoad[ind])
    delta.leftPump[ind] <- 1-(z.IMA[ind]==3)-(z.Cardiomio[ind]==3)-(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])-
      (tilde.postLoad[ind]>=1/3)*(-1.5+3*tilde.postLoad[ind])

    # Chest pain type (6 categories)
    z.Dol_type[ind] ~ dcat(theta.Dol_type[ind,1:6]) 
    for(aux in 1:6) {
      theta.Dol_type[ind,aux] <- (delta.Dol_type[ind]<1)*odds.Dol_type[ind,aux]/sum(odds.Dol_type[ind,1:6])+(delta.Dol_type[ind]>=1)*p.Dol_type[aux]
      }
    odds.Dol_type[ind,1] <- 1                  
    odds.Dol_type[ind,2] <- exp(eta.Dol_type[1]*delta.Dol_type[ind]+
      log(p.Dol_type[8]/p.Dol_type[7])*Coronar[ind]+
      log(p.Dol_type[14]/p.Dol_type[13])*D_a_pu[ind]+
      log(p.Dol_type[20]/p.Dol_type[19])*D_cen[ind]+
      log(p.Dol_type[26]/p.Dol_type[25])*D_pl[ind]+
      log(p.Dol_type[32]/p.Dol_type[31])*D_quadr_sup[ind]+
      log(p.Dol_type[38]/p.Dol_type[37])*D_par[ind])
    odds.Dol_type[ind,3] <- exp(eta.Dol_type[2]*delta.Dol_type[ind]+
      log(p.Dol_type[9]/p.Dol_type[7])*Coronar[ind]+
      log(p.Dol_type[15]/p.Dol_type[13])*D_a_pu[ind]+
      log(p.Dol_type[21]/p.Dol_type[19])*D_cen[ind]+
      log(p.Dol_type[27]/p.Dol_type[25])*D_pl[ind]+
      log(p.Dol_type[33]/p.Dol_type[31])*D_quadr_sup[ind]+
      log(p.Dol_type[39]/p.Dol_type[37])*D_par[ind])
    odds.Dol_type[ind,4] <- exp(eta.Dol_type[3]*delta.Dol_type[ind]+
      log(p.Dol_type[10]/p.Dol_type[7])*Coronar[ind]+
      log(p.Dol_type[16]/p.Dol_type[13])*D_a_pu[ind]+
      log(p.Dol_type[22]/p.Dol_type[19])*D_cen[ind]+
      log(p.Dol_type[28]/p.Dol_type[25])*D_pl[ind]+
      log(p.Dol_type[34]/p.Dol_type[31])*D_quadr_sup[ind]+
      log(p.Dol_type[40]/p.Dol_type[37])*D_par[ind])
    odds.Dol_type[ind,5] <- exp(eta.Dol_type[4]*delta.Dol_type[ind]+
      log(p.Dol_type[11]/p.Dol_type[7])*Coronar[ind]+
      log(p.Dol_type[17]/p.Dol_type[13])*D_a_pu[ind]+
      log(p.Dol_type[23]/p.Dol_type[19])*D_cen[ind]+
      log(p.Dol_type[29]/p.Dol_type[25])*D_pl[ind]+
      log(p.Dol_type[35]/p.Dol_type[31])*D_quadr_sup[ind]+
      log(p.Dol_type[41]/p.Dol_type[37])*D_par[ind])
    odds.Dol_type[ind,6] <- exp(eta.Dol_type[5]*delta.Dol_type[ind]+
      log(p.Dol_type[12]/p.Dol_type[7])*Coronar[ind]+
      log(p.Dol_type[18]/p.Dol_type[13])*D_a_pu[ind]+
      log(p.Dol_type[24]/p.Dol_type[19])*D_cen[ind]+
      log(p.Dol_type[30]/p.Dol_type[25])*D_pl[ind]+
      log(p.Dol_type[36]/p.Dol_type[31])*D_quadr_sup[ind]+
      log(p.Dol_type[42]/p.Dol_type[37])*D_par[ind])
    delta.Dol_type[ind] <- 1-Coronar[ind]-D_a_pu[ind]-D_cen[ind]-D_pl[ind]-D_quadr_sup[ind]-D_par[ind]

    # Dilated right ventricle (absent if cardiac tamponade is present) (binary)
    DilatVdx[ind] ~ dbern(theta.DilatVdx[ind])
    theta.DilatVdx[ind] <- (z.tamp[ind]==1)*((delta.DilatVdx[ind]<1)*thetaStar.DilatVdx[ind]+(delta.DilatVdx[ind]>=1)*p.DilatVdx[1])
    logit(thetaStar.DilatVdx[ind]) <- eta.DilatVdx[1]*delta.DilatVdx[ind]+
      qlogis(p.DilatVdx[2],0,1)*(-1.5+3*tilde.RHI[ind])+
      qlogis(p.DilatVdx[3],0,1)*Card_dil[ind]
    delta.DilatVdx[ind] <- 1-(-1.5+3*tilde.RHI[ind])-Card_dil[ind]          
    
     # Heart post-load ((continuous, hyper-restricted)
    tilde.postLoad[ind] ~ dbeta(relmu.postLoad[ind]*p.postLoad[5],(1-relmu.postLoad[ind])*p.postLoad[5]) T(tol,1-tol)
    logit(relmu.postLoad[ind]) <- qlogis(p.postLoad[1],0,1)*delta.postLoad[ind]+
      qlogis(p.postLoad[2],0,1)*Insuf_aorAcut[ind]+
      qlogis(p.postLoad[3],0,1)*InsuffMitrAcut[ind]+
      qlogis(p.postLoad[4],0,1)*(z.Card_obstr[ind]==3)
    delta.postLoad[ind] <- 1-Insuf_aorAcut[ind]-InsuffMitrAcut[ind]-(z.Card_obstr[ind]==3)

    # Mitral valve failure (generic) (binary, non-stochastic)
    InsuffMitr[ind] ~ dbern(theta.InsuffMitr[ind])
    theta.InsuffMitr[ind] <- 1*(z.InsuffMitrCron[ind]>1||InsuffMitrAcut[ind]==1)
   
    # Chest pain (pleuritic origin) (binary)
    D_pl[ind] ~ dbern(theta.D_pl[ind])
    theta.D_pl[ind] <- (delta.D_pl[ind]<1)*thetaStar.D_pl[ind]+(delta.D_pl[ind]>=1)*p.D_pl[1]
    logit(thetaStar.D_pl[ind]) <- eta.D_pl[1]*delta.D_pl[ind]+
      qlogis(p.D_pl[2],0,1)*Fr[ind]+
      qlogis(p.D_pl[3],0,1)*Ple[ind]+
      qlogis(p.D_pl[4],0,1)*Costcndr[ind]+
      qlogis(p.D_pl[5],0,1)*(z.Polm[ind]==2)+
      qlogis(p.D_pl[6],0,1)*(z.Polm[ind]==3)+
      qlogis(p.D_pl[7],0,1)*(z.Pnx[ind]>1)
    delta.D_pl[ind] <- 1-Fr[ind]-Ple[ind]-Costcndr[ind]-(z.Polm[ind]>1)-(z.Pnx[ind]>1)

    # Cardiac tamponade (3 categories)
    z.tamp[ind] ~ dcat(theta.tamp[ind,1:3])
    for(aux in 1:3) {
      theta.tamp[ind,aux] <- (delta.tamp[ind]<1)*odds.tamp[ind,aux]/sum(odds.tamp[ind,1:3])+(delta.tamp[ind]>=1)*p.tamp[aux]
      }
    odds.tamp[ind,1] <- 1
    odds.tamp[ind,2] <- exp(eta.tamp[1]*delta.tamp[ind]+
      log(p.tamp[5]/p.tamp[4])*HemoPeric[ind]+
      log(p.tamp[8]/p.tamp[7])*Peric[ind])
    odds.tamp[ind,3] <- exp(eta.tamp[2]*delta.tamp[ind]+
      log(p.tamp[6]/p.tamp[4])*HemoPeric[ind]+
      log(p.tamp[9]/p.tamp[7])*Peric[ind])
    delta.tamp[ind] <- 1-HemoPeric[ind]-Peric[ind]

    # ST segment elevation (3 categories)
    #   [absent if left branch block is present]                                                                
    z.Q_inT[ind] ~ dcat(theta.Q_inT[ind,1:3])  
    for(aux in 1:3) {
      theta.Q_inT[ind,aux] <- (delta.Q_inT[ind]<1)*odds.Q_inT[ind,aux]/sum(odds.Q_inT[ind,1:3])+(delta.Q_inT[ind]>=1)*p.Q_inT[aux]
      }
    odds.Q_inT[ind,1] <- 1
    odds.Q_inT[ind,2] <- exp(eta.Q_inT[1]*delta.Q_inT[ind]+
      log(p.Q_inT[5]/p.Q_inT[4])*(z.IMA[ind]==3)+
      log(p.Q_inT[8]/p.Q_inT[7])*Peric[ind])*(BBSx[ind]==0)
    odds.Q_inT[ind,3] <- exp(eta.Q_inT[2]*delta.Q_inT[ind]+
      log(p.Q_inT[6]/p.Q_inT[4])*(z.IMA[ind]==3)+
      log(p.Q_inT[9]/p.Q_inT[7])*Peric[ind])*(BBSx[ind]==0)
    delta.Q_inT[ind] <- 1-(z.IMA[ind]==3)-Peric[ind]

    # Chest pain (retro-sternal origin) (binary)
    D_cen[ind] ~ dbern(theta.D_cen[ind])
    theta.D_cen[ind] <- (delta.D_cen[ind]<1)*thetaStar.D_cen[ind]+(delta.D_cen[ind]>=1)*p.D_cen[1]
    logit(thetaStar.D_cen[ind]) <- eta.D_cen[1]*delta.D_cen[ind]+
      qlogis(p.D_cen[2],0,1)*Aort_diss[ind]+
      qlogis(p.D_cen[3],0,1)*(z.DilArtPolm[ind]==2)+
      qlogis(p.D_cen[4],0,1)*Peric[ind]+
      qlogis(p.D_cen[5],0,1)*(z.Polm[ind]==3)+
      qlogis(p.D_cen[6],0,1)*ProlMitralCron[ind]+
      qlogis(p.D_cen[7],0,1)*Spas_es[ind]+
      qlogis(p.D_cen[8],0,1)*(z.Card_obstr[ind]>1)
    delta.D_cen[ind] <- 1-Aort_diss[ind]-(z.DilArtPolm[ind]==2)-Peric[ind]-(z.Polm[ind]==3)-ProlMitralCron[ind]-Spas_es[ind]-(z.Card_obstr[ind]>1)

    # Extrasystoles (binary)
    extraSist[ind] ~ dbern(theta.extraSist[ind])
    theta.extraSist[ind] <- (delta.extraSist[ind]<1)*thetaStar.extraSist[ind]+(delta.extraSist[ind]>=1)*p.extraSist[1]
    logit(thetaStar.extraSist[ind]) <- eta.extraSist[1]*delta.extraSist[ind]+
      qlogis(p.extraSist[2],0,1)*Arit_sotto[ind]+
      qlogis(p.extraSist[3],0,1)*ProlMitralCron[ind]
    delta.extraSist[ind] <- 1-Arit_sotto[ind]-ProlMitralCron[ind]

    # Ventricular segmental dyssynergia (4 categories)
    z.Diss_ventr[ind] ~ dcat(theta.Diss_ventr[ind,1:4])
    for(aux in 1:4) {
      theta.Diss_ventr[ind,aux] <- (delta.Diss_ventr[ind]<1)*odds.Diss_ventr[ind,aux]/sum(odds.Diss_ventr[ind,1:4])+(delta.Diss_ventr[ind]>=1)*p.Diss_ventr[aux]
      }
    odds.Diss_ventr[ind,1] <- 1
    odds.Diss_ventr[ind,2] <- exp(eta.Diss_ventr[1]*delta.Diss_ventr[ind]+
      log(p.Diss_ventr[6]/p.Diss_ventr[5])*Arit_sotto[ind]+
      log(p.Diss_ventr[10]/p.Diss_ventr[9])*Coronar[ind]+
      log(p.Diss_ventr[14]/p.Diss_ventr[13])*Card_dil[ind]+
      log(p.Diss_ventr[18]/p.Diss_ventr[17])*(z.IMA[ind]==2)+
      log(p.Diss_ventr[22]/p.Diss_ventr[21])*(z.IMA[ind]==3))
    odds.Diss_ventr[ind,3] <- exp(eta.Diss_ventr[2]*delta.Diss_ventr[ind]+
      log(p.Diss_ventr[7]/p.Diss_ventr[5])*Arit_sotto[ind]+
      log(p.Diss_ventr[11]/p.Diss_ventr[9])*Coronar[ind]+
      log(p.Diss_ventr[15]/p.Diss_ventr[13])*Card_dil[ind]+
      log(p.Diss_ventr[19]/p.Diss_ventr[17])*(z.IMA[ind]==2)+
      log(p.Diss_ventr[23]/p.Diss_ventr[21])*(z.IMA[ind]==3))
    odds.Diss_ventr[ind,4] <- exp(eta.Diss_ventr[3]*delta.Diss_ventr[ind]+
      log(p.Diss_ventr[8]/p.Diss_ventr[5])*Arit_sotto[ind]+
      log(p.Diss_ventr[12]/p.Diss_ventr[9])*Coronar[ind]+
      log(p.Diss_ventr[16]/p.Diss_ventr[13])*Card_dil[ind]+
      log(p.Diss_ventr[20]/p.Diss_ventr[17])*(z.IMA[ind]==2)+
      log(p.Diss_ventr[24]/p.Diss_ventr[21])*(z.IMA[ind]==3))
    delta.Diss_ventr[ind] <- 1-Arit_sotto[ind]-Coronar[ind]-Card_dil[ind]-(z.IMA[ind]>1)

    # Atrial arrhythmia (generic) (binary, non-stochastic)
    Arit_sopra_o[ind] ~ dbern(theta.Arit_sopra_o[ind])
    theta.Arit_sopra_o[ind] <- 1*(Arit_sopra_acut[ind]==1||Arit_sopra_cron[ind]==1)

    # Mitral valve prolapse (generic) (binary, non-stochastic)
    ProlMitral[ind] ~ dbern(theta.ProlMitral[ind])
    theta.ProlMitral[ind] <- 1*(ProlMitralCron[ind]==1||ProlMitralAcut[ind]==1)
    
    # Acute mitral valve failure (binary)
    InsuffMitrAcut[ind] ~ dbern(theta.InsuffMitrAcut[ind])
    theta.InsuffMitrAcut[ind] <- (delta.InsuffMitrAcut[ind]<1)*thetaStar.InsuffMitrAcut[ind]+(delta.InsuffMitrAcut[ind]>=1)*p.InsuffMitrAcut[1]
    logit(thetaStar.InsuffMitrAcut[ind]) <- eta.InsuffMitrAcut[1]*delta.InsuffMitrAcut[ind]+
      qlogis(p.InsuffMitrAcut[2],0,1)*(z.IMA[ind]==3)+
      qlogis(p.InsuffMitrAcut[3],0,1)*Endocard[ind]+
      qlogis(p.InsuffMitrAcut[4],0,1)*ProlMitralAcut[ind]+
      qlogis(p.InsuffMitrAcut[5],0,1)*Card_dil[ind]
    delta.InsuffMitrAcut[ind] <- 1-(z.IMA[ind]==3)-Endocard[ind]-ProlMitralAcut[ind]-Card_dil[ind]

    # Ruptured chordae tendineae (binary)
    CordTend[ind] ~ dbern(theta.CordTend[ind])
    theta.CordTend[ind] <- (delta.CordTend[ind]<1)*thetaStar.CordTend[ind]+(delta.CordTend[ind]>=1)*p.CordTend[1]
    logit(thetaStar.CordTend[ind]) <- eta.CordTend[1]*delta.CordTend[ind]+
      qlogis(p.CordTend[2],0,1)*ProlMitralCron[ind]+
      qlogis(p.CordTend[3],0,1)*ProlMitralAcut[ind]
    delta.CordTend[ind] <- 1-ProlMitralCron[ind]-ProlMitralAcut[ind]

    # Pleurisy (binary)
    Ple[ind] ~ dbern(theta.Ple[ind])
    theta.Ple[ind] <- (delta.Ple[ind]<1)*thetaStar.Ple[ind]+(delta.Ple[ind]>=1)*p.Ple[1]
    logit(thetaStar.Ple[ind]) <- eta.Ple[1]*delta.Ple[ind]+
      qlogis(p.Ple[2],0,1)*(z.Polm[ind]>1)+
      qlogis(p.Ple[3],0,1)*InP[ind]+
      qlogis(p.Ple[4],0,1)*Neo_p[ind]+
      qlogis(p.Ple[5],0,1)*Dressl[ind]+
      qlogis(p.Ple[6],0,1)*Peric_ni[ind]
    delta.Ple[ind] <- 1-(z.Polm[ind]>1)-InP[ind]-Neo_p[ind]-Dressl[ind]-Peric_ni[ind]

    # Pericarditis (binary, semi-stochastic)
    Peric[ind] ~ dbern(theta.Peric[ind])
    theta.Peric[ind] <- p.Peric[1]*(Dressl[ind]+Peric_ni[ind]==0)+1*(Dressl[ind]+Peric_ni[ind]>0)

    # Heartburn (binary)
    Piro[ind] ~ dbern(theta.Piro[ind])
    theta.Piro[ind] <- (delta.Piro[ind]<1)*thetaStar.Piro[ind]+(delta.Piro[ind]>=1)*p.Piro[1]
    logit(thetaStar.Piro[ind]) <- eta.Piro[1]*delta.Piro[ind]+
      qlogis(p.Piro[2],0,1)*(z.Dispeps[ind]==2)+
      qlogis(p.Piro[3],0,1)*(z.Dispeps[ind]==3)  
    delta.Piro[ind] <- 1-(z.Dispeps[ind]>1)

    # Nausea (binary)
    naus[ind] ~ dbern(theta.naus[ind])
    theta.naus[ind] <- (delta.naus[ind]<1)*thetaStar.naus[ind]+(delta.naus[ind]>=1)*p.naus[1]
    logit(thetaStar.naus[ind]) <- eta.naus[1]*delta.naus[ind]+
      qlogis(p.naus[2],0,1)*(z.Dispeps[ind]>1)
    delta.naus[ind] <- 1-(z.Dispeps[ind]>1)
    
    # auxiliar node
    Dispeps_J[ind] ~ dbern(theta.Dispeps_J[ind])
    theta.Dispeps_J[ind] <- 1*(z.Dispeps[ind]>1)
    
    # Ventricular arrhythmia (binary)
    Arit_sotto[ind] ~ dbern(theta.Arit_sotto[ind])   
    theta.Arit_sotto[ind] <- (delta.Arit_sotto[ind]<1)*thetaStar.Arit_sotto[ind]+(delta.Arit_sotto[ind]>=1)*p.Arit_sotto[1]
    logit(thetaStar.Arit_sotto[ind]) <- eta.Arit_sotto[1]*delta.Arit_sotto[ind]+
      qlogis(p.Arit_sotto[2],0,1)*(z.BradiTachi[ind]==2)+
      qlogis(p.Arit_sotto[3],0,1)*(z.BradiTachi[ind]==3)+
      qlogis(p.Arit_sotto[4],0,1)*(z.BradiTachi[ind]==4)
    delta.Arit_sotto[ind] <- 1-(z.BradiTachi[ind]>1)

    # Acute atrial arrhythmia (binary)
    Arit_sopra_acut[ind] ~ dbern(theta.Arit_sopra_acut[ind])
    theta.Arit_sopra_acut[ind] <- (delta.Arit_sopra_acut[ind]<1)*thetaStar.Arit_sopra_acut[ind]+(delta.Arit_sopra_acut[ind]>=1)*p.Arit_sopra_acut[1]
    logit(thetaStar.Arit_sopra_acut[ind]) <- eta.Arit_sopra_acut[1]*delta.Arit_sopra_acut[ind]+
      qlogis(p.Arit_sopra_acut[2],0,1)*Arit_sopra_cron[ind]+
      qlogis(p.Arit_sopra_acut[3],0,1)*(z.BradiTachi[ind]==3)+
      qlogis(p.Arit_sopra_acut[4],0,1)*(z.BradiTachi[ind]==4)
    delta.Arit_sopra_acut[ind] <- 1-Arit_sopra_cron[ind]-(z.BradiTachi[ind]==3||z.BradiTachi[ind]==4)

    # Acute mitral valve prolapse (binary)
    ProlMitralAcut[ind] ~ dbern(theta.ProlMitralAcut[ind])
    theta.ProlMitralAcut[ind] <- (delta.ProlMitralAcut[ind]<1)*thetaStar.ProlMitralAcut[ind]+(delta.ProlMitralAcut[ind]>=1)*p.ProlMitralAcut[1]
    logit(thetaStar.ProlMitralAcut[ind]) <- eta.ProlMitralAcut[1]*delta.ProlMitralAcut[ind]+
      qlogis(p.ProlMitralAcut[2],0,1)*(z.IMA[ind]==3)
    delta.ProlMitralAcut[ind] <- 1-(z.IMA[ind]==3)

    # Myoglobin (binary)
    Mioglob[ind] ~ dbern(theta.Mioglob[ind])
    theta.Mioglob[ind] <- (delta.Mioglob[ind]<1)*thetaStar.Mioglob[ind]+(delta.Mioglob[ind]>=1)*p.Mioglob[1]
    logit(thetaStar.Mioglob[ind]) <- eta.Mioglob[1]*delta.Mioglob[ind]+
      qlogis(p.Mioglob[2],0,1)*(z.IMA[ind]>1) 
    delta.Mioglob[ind] <- 1-(z.IMA[ind]>1)
    
    # auxiliar node
    IMA_J[ind] ~ dbern(theta.IMA_J[ind])
    theta.IMA_J[ind] <- 1*(z.IMA[ind]>1)

    # Hemopericardium (binary)
    HemoPeric[ind] ~ dbern(theta.HemoPeric[ind])
    theta.HemoPeric[ind] <- (delta.HemoPeric[ind]<1)*thetaStar.HemoPeric[ind]+(delta.HemoPeric[ind]>=1)*p.HemoPeric[1]
    logit(thetaStar.HemoPeric[ind]) <- eta.HemoPeric[1]*delta.HemoPeric[ind]+
      qlogis(p.HemoPeric[2],0,1)*(z.IMA[ind]==3)+
      qlogis(p.HemoPeric[3],0,1)*Aort_diss[ind]
    delta.HemoPeric[ind] <- 1-(z.IMA[ind]==3)-Aort_diss[ind]

    # Non-infective pericarditis (binary)
    Dressl[ind] ~ dbern(theta.Dressl[ind])
    theta.Dressl[ind] <- (delta.Dressl[ind]<1)*thetaStar.Dressl[ind]+(delta.Dressl[ind]>=1)*p.Dressl[1]
    logit(thetaStar.Dressl[ind]) <- eta.Dressl[1]*delta.Dressl[ind]+
      qlogis(p.Dressl[2],0,1)*(z.IMA[ind]==2)+
      qlogis(p.Dressl[3],0,1)*(z.IMA[ind]==3)
    delta.Dressl[ind] <- 1-(z.IMA[ind]>1)

    # Dyspepsia (3 categories)
    z.Dispeps[ind] ~ dcat(theta.Dispeps[ind,1:3])  
    for(aux in 1:3) {
      theta.Dispeps[ind,aux] <- (delta.Dispeps[ind]<1)*odds.Dispeps[ind,aux]/sum(odds.Dispeps[ind,1:3])+(delta.Dispeps[ind]>=1)*p.Dispeps[aux]
      }
    odds.Dispeps[ind,1] <- 1           
    odds.Dispeps[ind,2] <- exp(eta.Dispeps[1]*delta.Dispeps[ind]+
      log(p.Dispeps[5]/p.Dispeps[4])*Colic[ind]+
      log(p.Dispeps[8]/p.Dispeps[7])*Refl_ge[ind]+
      log(p.Dispeps[11]/p.Dispeps[10])*Ulc_pep[ind]+
      log(p.Dispeps[14]/p.Dispeps[13])*(z.IMA[ind]==2)+
      log(p.Dispeps[17]/p.Dispeps[16])*(z.IMA[ind]==3))      
    odds.Dispeps[ind,3] <- exp(eta.Dispeps[2]*delta.Dispeps[ind]+
      log(p.Dispeps[6]/p.Dispeps[4])*Colic[ind]+
      log(p.Dispeps[9]/p.Dispeps[7])*Refl_ge[ind]+
      log(p.Dispeps[12]/p.Dispeps[10])*Ulc_pep[ind]+
      log(p.Dispeps[15]/p.Dispeps[13])*(z.IMA[ind]==2)+
      log(p.Dispeps[18]/p.Dispeps[16])*(z.IMA[ind]==3))   
    delta.Dispeps[ind] <- 1-Colic[ind]-Refl_ge[ind]-Ulc_pep[ind]-(z.IMA[ind]>1)

    # CK-MB (binary)
    CK_MB[ind] ~ dbern(theta.CK_MB[ind])
    theta.CK_MB[ind] <- (delta.CK_MB[ind]<1)*thetaStar.CK_MB[ind]+(delta.CK_MB[ind]>=1)*p.CK_MB[1]
    logit(thetaStar.CK_MB[ind]) <- eta.CK_MB[1]*delta.CK_MB[ind]+
      qlogis(p.CK_MB[2],0,1)*(z.IMA[ind]==2)+
      qlogis(p.CK_MB[3],0,1)*(z.IMA[ind]==3)
    delta.CK_MB[ind] <- 1-(z.IMA[ind]>1)
 
    # Angina without infarction (binary, non-stochastic)
    Angina[ind] ~ dbern(theta.Angina[ind])
    theta.Angina[ind] <- 1*(Coronar[ind]==1&&z.IMA[ind]==1)

    # T-wave inversion in V1-V3 (binary)
    Sovracc[ind] ~ dbern(theta.Sovracc[ind])
    theta.Sovracc[ind] <- (delta.Sovracc[ind]<1)*thetaStar.Sovracc[ind]+(delta.Sovracc[ind]>=1)*p.Sovracc[1]
    logit(thetaStar.Sovracc[ind]) <- eta.Sovracc[1]*delta.Sovracc[ind]+
      qlogis(p.Sovracc[2],0,1)*ECG_dx[ind] 
    delta.Sovracc[ind] <- 1-ECG_dx[ind]

    # Cardiac axis right deviation (binary)
    DevAsse[ind] ~ dbern(theta.DevAsse[ind])
    theta.DevAsse[ind] <- (delta.DevAsse[ind]<1)*thetaStar.DevAsse[ind]+(delta.DevAsse[ind]>=1)*p.DevAsse[1]
    logit(thetaStar.DevAsse[ind]) <- eta.DevAsse[1]*delta.DevAsse[ind]+
      qlogis(p.DevAsse[2],0,1)*ECG_dx[ind]
    delta.DevAsse[ind] <- 1-ECG_dx[ind]

    # Bradycardia/Tachycardia (4 categories)
    z.BradiTachi[ind] ~ dcat(theta.BradiTachi[ind,1:4])
    for(aux in 1:4) {
      theta.BradiTachi[ind,aux] <- (delta.BradiTachi[ind]<1)*odds.BradiTachi[ind,aux]/sum(odds.BradiTachi[ind,1:4])+(delta.BradiTachi[ind]>=1)*p.BradiTachi[aux]
      }
    odds.BradiTachi[ind,1] <- 1
    odds.BradiTachi[ind,2] <- exp(eta.BradiTachi[1]*delta.BradiTachi[ind]+
      log(p.BradiTachi[6]/p.BradiTachi[5])*(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])+
      log(p.BradiTachi[10]/p.BradiTachi[9])*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind]))
    odds.BradiTachi[ind,3] <- exp(eta.BradiTachi[2]*delta.BradiTachi[ind]+
      log(p.BradiTachi[7]/p.BradiTachi[5])*(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])+
      log(p.BradiTachi[11]/p.BradiTachi[9])*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind]))
    odds.BradiTachi[ind,4] <- exp(eta.BradiTachi[3]*delta.BradiTachi[ind]+
      log(p.BradiTachi[8]/p.BradiTachi[5])*(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])+
      log(p.BradiTachi[12]/p.BradiTachi[9])*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind]))
    delta.BradiTachi[ind] <- 1-(tilde.driveCard[ind]>=1/3)*(-1.5+3*tilde.driveCard[ind])-(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])

    # Non ST segment elevation (3 categories)
    z.nonSTE[ind] ~ dcat(theta.nonSTE[ind,1:3])
    for(aux in 1:3) {
      theta.nonSTE[ind,aux] <- (delta.nonSTE[ind]<1)*odds.nonSTE[ind,aux]/sum(odds.nonSTE[ind,1:3])+(delta.nonSTE[ind]>=1)*p.nonSTE[aux]
      }
    odds.nonSTE[ind,1] <- 1
    odds.nonSTE[ind,2] <- exp(eta.nonSTE[1]*delta.nonSTE[ind]+
      log(p.nonSTE[5]/p.nonSTE[4])*Coronar[ind]+
      log(p.nonSTE[8]/p.nonSTE[7])*Miocard[ind]+
      log(p.nonSTE[11]/p.nonSTE[10])*(z.Cardiomio[ind]==2)+
      log(p.nonSTE[14]/p.nonSTE[13])*(z.Cardiomio[ind]==3))*(BBSx[ind]==0)
    odds.nonSTE[ind,3] <- exp(eta.nonSTE[2]*delta.nonSTE[ind]+
      log(p.nonSTE[6]/p.nonSTE[4])*Coronar[ind]+
      log(p.nonSTE[9]/p.nonSTE[7])*Miocard[ind]+
      log(p.nonSTE[12]/p.nonSTE[10])*(z.Cardiomio[ind]==2)+
      log(p.nonSTE[15]/p.nonSTE[13])*(z.Cardiomio[ind]==3))*(BBSx[ind]==0)
    delta.nonSTE[ind] <- 1-Coronar[ind]-Miocard[ind]-(z.Cardiomio[ind]>1)

    # Temporary suspension of heart drive (binary)
    SuspHeaDrivTmp[ind] ~ dbern(theta.SuspHeaDrivTmp[ind])
    theta.SuspHeaDrivTmp[ind] <- (delta.SuspHeaDrivTmp[ind]<1)*thetaStar.SuspHeaDrivTmp[ind]+(delta.SuspHeaDrivTmp[ind]>=1)*p.SuspHeaDrivTmp[1]
    logit(thetaStar.SuspHeaDrivTmp[ind]) <- eta.SuspHeaDrivTmp[1]*delta.SuspHeaDrivTmp[ind]+
      qlogis(p.SuspHeaDrivTmp[2],0,1)*Coronar[ind]+
      qlogis(p.SuspHeaDrivTmp[3],0,1)*Feocrom[ind]+
      qlogis(p.SuspHeaDrivTmp[4],0,1)*Tireotos[ind]+
      qlogis(p.SuspHeaDrivTmp[5],0,1)*WPW[ind]
    delta.SuspHeaDrivTmp[ind] <- 1-Coronar[ind]-Feocrom[ind]-Tireotos[ind]-WPW[ind]

    # Acute myocardial infarction (3 categories)
    z.IMA[ind] ~ dcat(theta.IMA[ind,1:3])
    for(aux in 1:3) {
      theta.IMA[ind,aux] <- (delta.IMA[ind]<1)*odds.IMA[ind,aux]/sum(odds.IMA[ind,1:3])+(delta.IMA[ind]>=1)*p.IMA[aux]
      }
    odds.IMA[ind,1] <- 1
    odds.IMA[ind,2] <- exp(eta.IMA[1]*delta.IMA[ind]+
      log(p.IMA[5]/p.IMA[4])*Coronar[ind]+
      log(p.IMA[8]/p.IMA[7])*Miocard[ind])
    odds.IMA[ind,3] <- exp(eta.IMA[2]*delta.IMA[ind]+
      log(p.IMA[6]/p.IMA[4])*Coronar[ind]+
      log(p.IMA[9]/p.IMA[7])*Miocard[ind])
    delta.IMA[ind] <- 1-Coronar[ind]-Miocard[ind]

    # ECG right heart findings (binary)
    ECG_dx[ind] ~ dbern(theta.ECG_dx[ind])
    theta.ECG_dx[ind] <- (delta.ECG_dx[ind]<1)*thetaStar.ECG_dx[ind]+(delta.ECG_dx[ind]>=1)*p.ECG_dx[1]
    logit(thetaStar.ECG_dx[ind]) <- eta.ECG_dx[1]*delta.ECG_dx[ind]+
      qlogis(p.ECG_dx[2],0,1)*Coronar[ind]+
      qlogis(p.ECG_dx[3],0,1)*(-1.5+3*tilde.RHI[ind])
    delta.ECG_dx[ind] <- 1-Coronar[ind]-(-1.5+3*tilde.RHI[ind])

    # Heart drive (continuous, hyper-restricted)
    tilde.driveCard[ind] ~ dbeta(relmu.driveCard[ind]*p.driveCard[6],(1-relmu.driveCard[ind])*p.driveCard[6]) T(tol,1-tol)
    logit(relmu.driveCard[ind]) <- qlogis(p.driveCard[1],0,1)*delta.driveCard[ind]+
      qlogis(p.driveCard[2],0,1)*Coronar[ind]+
      qlogis(p.driveCard[3],0,1)*Feocrom[ind]+
      qlogis(p.driveCard[4],0,1)*Tireotos[ind]+
      qlogis(p.driveCard[5],0,1)*WPW[ind]
    delta.driveCard[ind] <- 1-Coronar[ind]-Feocrom[ind]-Tireotos[ind]-WPW[ind]

    # Left bundle branch block (binary)
    BBSx[ind] ~ dbern(theta.BBSx[ind])
    theta.BBSx[ind] <- (delta.BBSx[ind]<1)*thetaStar.BBSx[ind]+(delta.BBSx[ind]>=1)*p.BBSx[1]
    logit(thetaStar.BBSx[ind]) <- eta.BBSx[1]*delta.BBSx[ind]+
      qlogis(p.BBSx[2],0,1)*(z.Cardiomio[ind]==2)+
      qlogis(p.BBSx[3],0,1)*(z.Cardiomio[ind]==3)+
      qlogis(p.BBSx[4],0,1)*Coronar[ind]
    delta.BBSx[ind] <- 1-(z.Cardiomio[ind]>1)-Coronar[ind]
    
    # auxiliar node
    RHO_J[ind] ~ dbern(theta.RHO_J[ind])
    theta.RHO_J[ind] <- 1*(tilde.RHO[ind]>=2/3)
    
    # Acute coronary event (binary)
    Coronar[ind] ~ dbern(theta.Coronar[ind])
    theta.Coronar[ind] <- 1*(z.ArtIntraVasCoa[ind]==2||z.ArtIntraVasCoa[ind]==5)+
      (z.ArtIntraVasCoa[ind]!=2&&z.ArtIntraVasCoa[ind]!=5)*((delta.Coronar[ind]<1)*thetaStar.Coronar[ind]+(delta.Coronar[ind]>=1)*p.Coronar[1])
    logit(thetaStar.Coronar[ind]) <- eta.Coronar[1]*delta.Coronar[ind]+
      qlogis(p.Coronar[2],0,1)*A_coc[ind]+
      qlogis(p.Coronar[3],0,1)*(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])+
      qlogis(p.Coronar[4],0,1)*(z.Cardiomio[ind]==3)
    delta.Coronar[ind] <- 1-A_coc[ind]-(tilde.RHO[ind]>=1/3)*(-1.5+3*tilde.RHO[ind])-(z.Cardiomio[ind]==3)-
      (z.ArtIntraVasCoa[ind]==2||z.ArtIntraVasCoa[ind]==5)
        
    # Leukocytosis (4 categories, semi-stochastic)
    z.Leucoc_o[ind] ~ dcat(theta.Leucoc_o[ind,1:4])
    theta.Leucoc_o[ind,1] <- (1-p.Leucoc_o[1])*(Linfoc_o[ind]+Neutrof_o[ind]==0)+(1-p.Leucoc_o[2])*(Linfoc_o[ind]+Neutrof_o[ind]==1)
    theta.Leucoc_o[ind,2] <- p.Leucoc_o[1]*(Linfoc_o[ind]+Neutrof_o[ind]==0)
    theta.Leucoc_o[ind,3] <- p.Leucoc_o[2]*(Linfoc_o[ind]+Neutrof_o[ind]==1)
    theta.Leucoc_o[ind,4] <- 1*(Linfoc_o[ind]+Neutrof_o[ind]==2)

    # Fever (binary)
    Fever[ind] ~ dbern(theta.Fever[ind])
    theta.Fever[ind] <- (delta.Fever[ind]<1)*thetaStar.Fever[ind]+(delta.Fever[ind]>=1)*p.Fever[1]
    logit(thetaStar.Fever[ind]) <- eta.Fever[1]*delta.Fever[ind]+
      qlogis(pStar2.Fever[ind],0,1)*Infiam_b[ind]+
      qlogis(pStar3.Fever[ind],0,1)*Infiam_nb[ind]
    pStar2.Fever[ind] <- p.Fever[2]*(1-p.Fever[4]*FANS[ind])
    pStar3.Fever[ind] <- p.Fever[3]*(1-p.Fever[4]*FANS[ind])
    delta.Fever[ind] <- 1-Infiam_b[ind]-Infiam_nb[ind]    
    
    # Paradoxical interventricular septum (binary)
    Setto_par[ind] ~ dbern(theta.Setto_par[ind])
    theta.Setto_par[ind] <- (delta.Setto_par[ind]<1)*thetaStar.Setto_par[ind]+(delta.Setto_par[ind]>=1)*p.Setto_par[1]
    logit(thetaStar.Setto_par[ind]) <- eta.Setto_par[1]*delta.Setto_par[ind]+
      qlogis(p.Setto_par[2],0,1)*(-1.5+3*tilde.RHI[ind])
    delta.Setto_par[ind] <- 1-(-1.5+3*tilde.RHI[ind])
    
    # Right ventricular output (continuous, hyper-restricted)
    tilde.RHO[ind] ~ dbeta(relmu.RHO[ind]*p.RHO[5],(1-relmu.RHO[ind])*p.RHO[5]) T(tol,1-tol)
    logit(relmu.RHO[ind]) <- qlogis(p.RHO[1],0,1)*delta.RHO[ind]+
      qlogis(p.RHO[2],0,1)*(tilde.RHI[ind]<0.5)*(-1.5+3*tilde.RHI[ind])+
      qlogis(p.RHO[3],0,1)*(tilde.RHI[ind]>0.5)*(-1.5+3*tilde.RHI[ind])+
      qlogis(p.RHO[4],0,1)*RHtrigg[ind]
    delta.RHO[ind] <- 1-(-1.5+3*tilde.RHI[ind])-RHtrigg[ind]

    # Tricuspid valve insufficiency (binary)
    InsTric_o[ind] ~ dbern(theta.InsTric_o[ind])
    theta.InsTric_o[ind] <- (delta.InsTric_o[ind]<1)*thetaStar.InsTric_o[ind]+(delta.InsTric_o[ind]>=1)*p.InsTric_o[1]
    logit(thetaStar.InsTric_o[ind]) <- eta.InsTric_o[1]*delta.InsTric_o[ind]+
      qlogis(p.InsTric_o[2],0,1)*CuPC[ind]+
      qlogis(p.InsTric_o[3],0,1)*Endocard[ind]+
      qlogis(p.InsTric_o[4],0,1)*(z.InsuffMitrCron[ind]==3)+
      qlogis(p.InsTric_o[5],0,1)*(-1.5+3*tilde.RHI[ind])
    delta.InsTric_o[ind] <- 1-CuPC[ind]-Endocard[ind]-(z.InsuffMitrCron[ind]==3)-(-1.5+3*tilde.RHI[ind])

    # Right bundle branch block (binary)
    BBD[ind] ~ dbern(theta.BBD[ind])
    theta.BBD[ind] <- (delta.BBD[ind]<1)*thetaStar.BBD[ind]+(delta.BBD[ind]>=1)*p.BBD[1]
    logit(thetaStar.BBD[ind]) <- eta.BBD[1]*delta.BBD[ind]+
      qlogis(p.BBD[2],0,1)*CuPC[ind]+
      qlogis(p.BBD[3],0,1)*(-1.5+3*tilde.RHI[ind])
    delta.BBD[ind] <- 1-CuPC[ind]-(-1.5+3*tilde.RHI[ind])
  
    # Lymphocytosis (binary)
    Linfoc_o[ind] ~ dbern(theta.Linfoc_o[ind])
    theta.Linfoc_o[ind] <- (delta.Linfoc_o[ind]<1)*thetaStar.Linfoc_o[ind]+(delta.Linfoc_o[ind]>=1)*p.Linfoc_o[1]
    logit(thetaStar.Linfoc_o[ind]) <- eta.Linfoc_o[1]*delta.Linfoc_o[ind]+
      qlogis(p.Linfoc_o[2],0,1)*Infiam_nb[ind]
    delta.Linfoc_o[ind] <- 1-Infiam_nb[ind]  
            
    # Anti-inflammatory drugs recent intake (binary)
    FANS[ind] ~ dbern(theta.FANS[ind])
    theta.FANS[ind] <- (delta.FANS[ind]<1)*thetaStar.FANS[ind]+(delta.FANS[ind]>=1)*p.FANS[1]
    logit(thetaStar.FANS[ind]) <- eta.FANS[1]*delta.FANS[ind]+
      qlogis(p.FANS[2],0,1)*Infiam_b[ind]+
      qlogis(p.FANS[3],0,1)*Infiam_nb[ind]  
    delta.FANS[ind] <- 1-Infiam_b[ind]-Infiam_nb[ind]

    # Neutrophilia (binary)
    Neutrof_o[ind] ~ dbern(theta.Neutrof_o[ind])
    theta.Neutrof_o[ind] <- (delta.Neutrof_o[ind]<1)*thetaStar.Neutrof_o[ind]+(delta.Neutrof_o[ind]>=1)*p.Neutrof_o[1]
    logit(thetaStar.Neutrof_o[ind]) <- eta.Neutrof_o[1]*delta.Neutrof_o[ind]+
      qlogis(p.Neutrof_o[2],0,1)*Infiam_b[ind]
    delta.Neutrof_o[ind] <- 1-Infiam_b[ind]

    # Right ventricular pre-load (continuous)
    tilde.RHI[ind] ~ dbeta(relmu.RHI[ind]*p.RHI[5],(1-relmu.RHI[ind])*p.RHI[5]) T(tol,1-tol)
    logit(relmu.RHI[ind]) <- qlogis(p.RHI[1],0,1)*delta.RHI[ind]+
      qlogis(p.RHI[2],0,1)*(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])+
      qlogis(p.RHI[3],0,1)*CuPC[ind]+
      qlogis(p.RHI[4],0,1)*RHtrigg[ind]
    delta.RHI[ind] <- 1-(tilde.Disidrat[ind]>=1/3)*(-1.5+3*tilde.Disidrat[ind])-CuPC[ind]-RHtrigg[ind]
    
    # auxiliar node    
    pulmPerf_J[ind] ~ dbern(theta.pulmPerf_J[ind])
    theta.pulmPerf_J[ind] <- 1*(tilde.pulmPerf[ind]>=2/3)

    # Obstruction of the systemic circulation (binary)
    Ostr_circ_sist[ind] ~ dbern(theta.Ostr_circ_sist[ind])
    theta.Ostr_circ_sist[ind] <- (delta.Ostr_circ_sist[ind]<1)*thetaStar.Ostr_circ_sist[ind]+(delta.Ostr_circ_sist[ind]>=1)*p.Ostr_circ_sist[1]
    logit(thetaStar.Ostr_circ_sist[ind]) <- eta.Ostr_circ_sist[1]*delta.Ostr_circ_sist[ind]+
      qlogis(p.Ostr_circ_sist[2],0,1)*(z.Card_obstr[ind]==2)+
      qlogis(p.Ostr_circ_sist[3],0,1)*(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])
    delta.Ostr_circ_sist[ind] <- 1-(z.Card_obstr[ind]==2)-(tilde.pulmPerf[ind]>=1/3)*(-1.5+3*tilde.pulmPerf[ind])
    
    # Non-bacterial infection (binary)
    Infiam_nb[ind] ~ dbern(theta.Infiam_nb[ind])
    theta.Infiam_nb[ind] <- (delta.Infiam_nb[ind]<1)*thetaStar.Infiam_nb[ind]+(delta.Infiam_nb[ind]>=1)*p.Infiam_nb[1]
    logit(thetaStar.Infiam_nb[ind]) <- eta.Infiam_nb[1]*delta.Infiam_nb[ind]+
      qlogis(p.Infiam_nb[2],0,1)*InP[ind]+
      qlogis(p.Infiam_nb[3],0,1)*Miocard[ind]+
      qlogis(p.Infiam_nb[4],0,1)*Pancr[ind]+
      qlogis(p.Infiam_nb[5],0,1)*Peric_ni[ind]+
      qlogis(p.Infiam_nb[6],0,1)*URI[ind]
    delta.Infiam_nb[ind] <- 1-InP[ind]-Miocard[ind]-Pancr[ind]-Peric_ni[ind]-URI[ind]

    # Bacterial infection (binary)
    #   [Defed decreases all conditional probabilities by a fixed quantity] 
    Infiam_b[ind] ~ dbern(theta.Infiam_b[ind])
    theta.Infiam_b[ind] <- (delta.Infiam_b[ind]<1)*thetaStar.Infiam_b[ind]+(delta.Infiam_b[ind]>=1)*p.Infiam_b[1]
    logit(thetaStar.Infiam_b[ind]) <- eta.Infiam_b[1]*delta.Infiam_b[ind]+
      qlogis(pStar2.Infiam_b[ind],0,1)*Endocard[ind]+
      qlogis(pStar3.Infiam_b[ind],0,1)*InP[ind]+
      qlogis(pStar4.Infiam_b[ind],0,1)*Pancr[ind]+
      qlogis(pStar5.Infiam_b[ind],0,1)*Perit[ind]+
      qlogis(pStar6.Infiam_b[ind],0,1)*Miocard[ind]+
      qlogis(pStar7.Infiam_b[ind],0,1)*Peric_ni[ind]+      
      qlogis(pStar8.Infiam_b[ind],0,1)*(z.Polm[ind]>1||Colec[ind]==1||Sepsi[ind]==1)
    pStar2.Infiam_b[ind] <- p.Infiam_b[2]*(1-p.Infiam_b[9]*Defed[ind])
    pStar3.Infiam_b[ind] <- p.Infiam_b[3]*(1-p.Infiam_b[9]*Defed[ind])
    pStar4.Infiam_b[ind] <- p.Infiam_b[4]*(1-p.Infiam_b[9]*Defed[ind])
    pStar5.Infiam_b[ind] <- p.Infiam_b[5]*(1-p.Infiam_b[9]*Defed[ind])
    pStar6.Infiam_b[ind] <- p.Infiam_b[6]*(1-p.Infiam_b[9]*Defed[ind])
    pStar7.Infiam_b[ind] <- p.Infiam_b[7]*(1-p.Infiam_b[9]*Defed[ind])
    pStar8.Infiam_b[ind] <- p.Infiam_b[8]*(1-p.Infiam_b[9]*Defed[ind])
    delta.Infiam_b[ind] <- 1-Endocard[ind]-InP[ind]-Pancr[ind]-Perit[ind]-Miocard[ind]-Peric_ni[ind]-(z.Polm[ind]>1)-Colec[ind]-Sepsi[ind]

    # Hemoptysis (binary)
    Emo[ind] ~ dbern(theta.Emo[ind])
    theta.Emo[ind] <- (delta.Emo[ind]<1)*thetaStar.Emo[ind]+(delta.Emo[ind]>=1)*p.Emo[1]
    logit(thetaStar.Emo[ind]) <- eta.Emo[1]*delta.Emo[ind]+
      qlogis(p.Emo[2],0,1)*InP[ind]+
      qlogis(p.Emo[3],0,1)*Neo_p[ind]+
      qlogis(p.Emo[4],0,1)*URI[ind]+
      qlogis(p.Emo[5],0,1)*(z.Polm[ind]==2)+
      qlogis(p.Emo[6],0,1)*(z.Polm[ind]==3)
    delta.Emo[ind] <- 1-InP[ind]-Neo_p[ind]-URI[ind]-(z.Polm[ind]>1)
    
    # Acute respiratory distress syndrome (binary)
    ARDS[ind] ~ dbern(theta.ARDS[ind])
    theta.ARDS[ind] <- (delta.ARDS[ind]<1)*thetaStar.ARDS[ind]+(delta.ARDS[ind]>=1)*p.ARDS[1]
    logit(thetaStar.ARDS[ind]) <- eta.ARDS[1]*delta.ARDS[ind]+
      qlogis(p.ARDS[2],0,1)*Sepsi[ind]+
      qlogis(p.ARDS[3],0,1)*Pancr[ind]+
      qlogis(p.ARDS[4],0,1)*(z.Polm[ind]==2)+
      qlogis(p.ARDS[5],0,1)*InP[ind]+
      qlogis(p.ARDS[6],0,1)*Neo_p[ind]
    delta.ARDS[ind] <- 1-Sepsi[ind]-Pancr[ind]-(z.Polm[ind]==2)-InP[ind]-Neo_p[ind]
    
    # Pulmonary artery thrombosis (binary)
    E_tap_te[ind] ~ dbern(theta.E_tap_te[ind])
    theta.E_tap_te[ind] <- (delta.E_tap_te[ind]<1)*thetaStar.E_tap_te[ind]+(delta.E_tap_te[ind]>=1)*p.E_tap_te[1]
    logit(thetaStar.E_tap_te[ind]) <- eta.E_tap_te[1]*delta.E_tap_te[ind]+
      qlogis(p.E_tap_te[2],0,1)*(z.DilArtPolm[ind]==2)+
      qlogis(p.E_tap_te[3],0,1)*(z.DilArtPolm[ind]==3)+
      qlogis(p.E_tap_te[4],0,1)*(z.EP[ind]==2)+
      qlogis(p.E_tap_te[5],0,1)*(z.EP[ind]==3)   
    delta.E_tap_te[ind] <- 1-(z.DilArtPolm[ind]>1)-(z.EP[ind]>1)

    # Dilated pulmonary artery (binary)
    DilArtPolm_ec[ind] ~ dbern(theta.DilArtPolm_ec[ind])
    theta.DilArtPolm_ec[ind] <- (delta.DilArtPolm_ec[ind]<1)*thetaStar.DilArtPolm_ec[ind]+(delta.DilArtPolm_ec[ind]>=1)*p.DilArtPolm_ec[1]
    logit(thetaStar.DilArtPolm_ec[ind]) <- eta.DilArtPolm_ec[1]*delta.DilArtPolm_ec[ind]+
      qlogis(p.DilArtPolm_ec[2],0,1)*(z.DilArtPolm[ind]==2)+
      qlogis(p.DilArtPolm_ec[3],0,1)*(z.DilArtPolm[ind]==3)
    delta.DilArtPolm_ec[ind] <- 1-(z.DilArtPolm[ind]>1)
    
    # Pulmonary artery diameter (binary)
    DiamArtPolm[ind] ~ dbern(theta.DiamArtPolm[ind])
    theta.DiamArtPolm[ind] <- (delta.DiamArtPolm[ind]<1)*thetaStar.DiamArtPolm[ind]+(delta.DiamArtPolm[ind]>=1)*p.DiamArtPolm[1]
    logit(thetaStar.DiamArtPolm[ind]) <- eta.DiamArtPolm[1]*delta.DiamArtPolm[ind]+
      qlogis(p.DiamArtPolm[2],0,1)*(z.DilArtPolm[ind]==2)+
      qlogis(p.DiamArtPolm[3],0,1)*(z.DilArtPolm[ind]==3)
    delta.DiamArtPolm[ind] <- 1-(z.DilArtPolm[ind]>1)

    # Small pulmonary vessel diameter (binary)
    calVasPol[ind] ~ dbern(theta.calVasPol[ind])
    theta.calVasPol[ind] <- (delta.calVasPol[ind]<1)*thetaStar.calVasPol[ind]+(delta.calVasPol[ind]>=1)*p.calVasPol[1]
    logit(thetaStar.calVasPol[ind]) <- eta.calVasPol[1]*delta.calVasPol[ind]+
      qlogis(p.calVasPol[2],0,1)*(z.DilArtPolm[ind]==2)+
      qlogis(p.calVasPol[3],0,1)*(z.DilArtPolm[ind]==3)
    delta.calVasPol[ind] <- 1-(z.DilArtPolm[ind]>1)

    # Endoluminal thrombus (3 categories)
    z.TC_s_new[ind] ~ dcat(theta.TC_s_new[ind,1:3])
    for(aux in 1:3) {
      theta.TC_s_new[ind,aux] <- (delta.TC_s_new[ind]<1)*odds.TC_s_new[ind,aux]/sum(odds.TC_s_new[ind,1:3])+(delta.TC_s_new[ind]>=1)*p.TC_s_new[aux]
      }
    odds.TC_s_new[ind,1] <- 1
    odds.TC_s_new[ind,2] <- exp(eta.TC_s_new[1]*delta.TC_s_new[ind]+
      log(p.TC_s_new[5]/p.TC_s_new[4])*(z.EP[ind]==1)+
      log(p.TC_s_new[8]/p.TC_s_new[7])*(z.EP[ind]==2))
    odds.TC_s_new[ind,3] <- exp(eta.TC_s_new[2]*delta.TC_s_new[ind]+
      log(p.TC_s_new[6]/p.TC_s_new[4])*(z.EP[ind]==1)+
      log(p.TC_s_new[9]/p.TC_s_new[7])*(z.EP[ind]==2))
    delta.TC_s_new[ind] <- 1-(z.EP[ind]>1)

    # Right circolatory obstruction trigger (binary)
    RHtrigg[ind] ~ dbern(theta.RHtrigg[ind])
    theta.RHtrigg[ind] <- (delta.RHtrigg[ind]<1)*thetaStar.RHtrigg[ind]+(delta.RHtrigg[ind]>=1)*p.RHtrigg[1]
    logit(thetaStar.RHtrigg[ind]) <- eta.RHtrigg[1]*delta.RHtrigg[ind]+
      qlogis(p.RHtrigg[2],0,1)*(z.EP[ind]==2)+
      qlogis(p.RHtrigg[3],0,1)*(z.EP[ind]==3)+
      qlogis(p.RHtrigg[4],0,1)*(z.Pnx[ind]==3)
    delta.RHtrigg[ind] <- 1-(z.EP[ind]>1)-(z.Pnx[ind]==3)

    # Lung perfusion (continuous, hyper-restricted)
    tilde.pulmPerf[ind] ~ dbeta(relmu.pulmPerf[ind]*p.pulmPerf[5],(1-relmu.pulmPerf[ind])*p.pulmPerf[5]) T(tol,1-tol)
    logit(relmu.pulmPerf[ind]) <- qlogis(p.pulmPerf[1],0,1)*delta.pulmPerf[ind]+
      qlogis(p.pulmPerf[2],0,1)*(z.EP[ind]==2)+
      qlogis(p.pulmPerf[3],0,1)*(z.EP[ind]==3)+
      qlogis(p.pulmPerf[4],0,1)*Ipertpolm[ind]
    delta.pulmPerf[ind] <- 1-(z.EP[ind]>1)-Ipertpolm[ind]

    # Miller index (3 categories, non-stochastic)
    z.MI[ind] ~ dcat(theta.MI[ind,1:3])
    theta.MI[ind,1] <- 1*(z.EP[ind]==1)
    theta.MI[ind,2] <- 1*(z.EP[ind]==2)
    theta.MI[ind,3] <- 1*(z.EP[ind]==3)

    # Pulmonary infarction (binary)
    InP[ind] ~ dbern(theta.InP[ind])  
    theta.InP[ind] <- (delta.InP[ind]<1)*thetaStar.InP[ind]+(delta.InP[ind]>=1)*p.InP[1]
    logit(thetaStar.InP[ind]) <- eta.InP[1]*delta.InP[ind]+
      qlogis(p.InP[2],0,1)*(z.EP[ind]==2)+
      qlogis(p.InP[3],0,1)*(z.EP[ind]==3)
    delta.InP[ind] <- 1-(z.EP[ind]>1)
    
    # auxiliar node
    EP_J[ind] ~ dbern(theta.EP_J[ind])
    theta.EP_J[ind] <- 1*(z.EP[ind]>1)
    
    # Arterial embolism (binary)
    EmbArte[ind] ~ dbern(theta.EmbArte[ind])
    theta.EmbArte[ind] <- (delta.EmbArte[ind]<1)*thetaStar.EmbArte[ind]+(delta.EmbArte[ind]>=1)*p.EmbArte[1]
    logit(thetaStar.EmbArte[ind]) <- eta.EmbArte[1]*delta.EmbArte[ind]+
      qlogis(p.EmbArte[2],0,1)*(z.ArtIntraVasCoa[ind]==4)+
      qlogis(p.EmbArte[3],0,1)*(z.EP[ind]>1&&For_ov[ind]==1)
    delta.EmbArte[ind] <- 1-(z.ArtIntraVasCoa[ind]==4)-(z.EP[ind]>1&&For_ov[ind]==1)

    # Dilatated pulmonary artery disease (3 categories)
    z.DilArtPolm[ind] ~ dcat(theta.DilArtPolm[ind,1:3]) 
    for(aux in 1:3) {
      theta.DilArtPolm[ind,aux] <- (delta.DilArtPolm[ind]<1)*odds.DilArtPolm[ind,aux]/sum(odds.DilArtPolm[ind,1:3])+(delta.DilArtPolm[ind]>=1)*p.DilArtPolm[aux]
      }
    odds.DilArtPolm[ind,1] <- 1
    odds.DilArtPolm[ind,2] <- exp(eta.DilArtPolm[1]*delta.DilArtPolm[ind]+
      log(p.DilArtPolm[5]/p.DilArtPolm[4])*Ipertpolm[ind]+
      log(p.DilArtPolm[8]/p.DilArtPolm[7])*(z.EP[ind]==2)+
      log(p.DilArtPolm[11]/p.DilArtPolm[10])*(z.EP[ind]==3))
    odds.DilArtPolm[ind,3] <- exp(eta.DilArtPolm[2]*delta.DilArtPolm[ind]+
      log(p.DilArtPolm[6]/p.DilArtPolm[4])*Ipertpolm[ind]+
      log(p.DilArtPolm[9]/p.DilArtPolm[7])*(z.EP[ind]==2)+
      log(p.DilArtPolm[12]/p.DilArtPolm[10])*(z.EP[ind]==3))
    delta.DilArtPolm[ind] <- 1-Ipertpolm[ind]-(z.EP[ind]>1)

    # Pulmonary embolism (3 categories)
    z.EP[ind] ~ dcat(theta.EP[ind,1:3])
    for(aux in 1:3) {
      theta.EP[ind,aux] <- (delta.EP[ind]<1)*odds.EP[ind,aux]/sum(odds.EP[ind,1:3])+(delta.EP[ind]>=1)*p.EP[aux]
      }
    odds.EP[ind,1] <- 1
    odds.EP[ind,2] <- exp(eta.EP[1]*delta.EP[ind]+
      log(p.EP[5]/p.EP[4])*TrCdx[ind]+
      log(p.EP[8]/p.EP[7])*TVPi[ind]+
      log(p.EP[11]/p.EP[10])*TVPsup[ind])
    odds.EP[ind,3] <- exp(eta.EP[2]*delta.EP[ind]+
      log(p.EP[6]/p.EP[4])*TrCdx[ind]+
      log(p.EP[9]/p.EP[7])*TVPi[ind]+
      log(p.EP[12]/p.EP[10])*TVPsup[ind])
    delta.EP[ind] <- 1-TrCdx[ind]-TVPi[ind]-TVPsup[ind]
    
    # Right heart thrombus (binary)
    TrCdx[ind] ~ dbern(theta.TrCdx[ind])
    theta.TrCdx[ind] <- 1*(z.VenIntraVasCoa[ind]==4)+(z.VenIntraVasCoa[ind]!=4)*((delta.TrCdx[ind]<1)*thetaStar.TrCdx[ind]+(delta.TrCdx[ind]>=1)*p.TrCdx[1])
    logit(thetaStar.TrCdx[ind]) <- eta.TrCdx[1]*delta.TrCdx[ind]+
      qlogis(p.TrCdx[2],0,1)*TVPi[ind]+
      qlogis(p.TrCdx[3],0,1)*TVPsup[ind]
    delta.TrCdx[ind] <- 1-TVPi[ind]-TVPsup[ind]-(z.VenIntraVasCoa[ind]==4)

    # Lower limbs magnetic resonance phlebography (binary)
    FlebArtInf_R[ind] ~ dbern(theta.FlebArtInf_R[ind])
    theta.FlebArtInf_R[ind] <- (delta.FlebArtInf_R[ind]<1)*thetaStar.FlebArtInf_R[ind]+(delta.FlebArtInf_R[ind]>=1)*p.FlebArtInf_R[1]
    logit(thetaStar.FlebArtInf_R[ind]) <- eta.FlebArtInf_R[1]*delta.FlebArtInf_R[ind]+
      qlogis(p.FlebArtInf_R[2],0,1)*TVPi[ind]
    delta.FlebArtInf_R[ind] <- 1-TVPi[ind]
    
    # Iliac phlebography  (binary)
    FlebArtInf[ind] ~ dbern(theta.FlebArtInf[ind])
    theta.FlebArtInf[ind] <- (delta.FlebArtInf[ind]<1)*thetaStar.FlebArtInf[ind]+(delta.FlebArtInf[ind]>=1)*p.FlebArtInf[1]
    logit(thetaStar.FlebArtInf[ind]) <- eta.FlebArtInf[1]*delta.FlebArtInf[ind]+
      qlogis(p.FlebArtInf[2],0,1)*TVPi[ind]
    delta.FlebArtInf[ind] <- 1-TVPi[ind]

    # Lower limbs echo-color doppler (binary)
    Eco_d[ind] ~ dbern(theta.Eco_d[ind])
    theta.Eco_d[ind] <- (delta.Eco_d[ind]<1)*thetaStar.Eco_d[ind]+(delta.Eco_d[ind]>=1)*p.Eco_d[1]
    logit(thetaStar.Eco_d[ind]) <- eta.Eco_d[1]*delta.Eco_d[ind]+
      qlogis(p.Eco_d[2],0,1)*TVPi[ind]
    delta.Eco_d[ind] <- 1-TVPi[ind]
    
    # Lower limbs compression ultrasounds (binary)
    Eco_c[ind] ~ dbern(theta.Eco_c[ind])        
    theta.Eco_c[ind] <- (delta.Eco_c[ind]<1)*thetaStar.Eco_c[ind]+(delta.Eco_c[ind]>=1)*p.Eco_c[1]
    logit(thetaStar.Eco_c[ind]) <- eta.Eco_c[1]*delta.Eco_c[ind]+
      qlogis(p.Eco_c[2],0,1)*TVPi[ind]
    delta.Eco_c[ind] <- 1-TVPi[ind]

    # Lower limbs pain (binary)
    Dol_g[ind] ~ dbern(theta.Dol_g[ind])
    theta.Dol_g[ind] <- (delta.Dol_g[ind]<1)*thetaStar.Dol_g[ind]+(delta.Dol_g[ind]>=1)*p.Dol_g[1]
    logit(thetaStar.Dol_g[ind]) <- eta.Dol_g[1]*delta.Dol_g[ind]+
      qlogis(p.Dol_g[2],0,1)*TVPi[ind]
    delta.Dol_g[ind] <- 1-TVPi[ind]

    # D-dimer test (binary with age)
    #  [Prof taking the 2nd non-neutral value decreases all conditional probabilities by a fixed quantity]
    Ddimer[ind] ~ dbern(theta.Ddimer[ind])   
    theta.Ddimer[ind] <- (delta.Ddimer[ind]<1)*thetaStar.Ddimer[ind]+(delta.Ddimer[ind]>=1)*wu.Ddimer[ind]*p.Ddimer[1]
    logit(thetaStar.Ddimer[ind]) <- qlogis(plogis(eta.Ddimer[1],0,1)*wu.Ddimer[ind],0,1)*delta.Ddimer[ind]+
      qlogis(pStar2.Ddimer[ind],0,1)*(z.Grav[ind]>1)+
      qlogis(pStar3.Ddimer[ind],0,1)*Fibrin[ind]
    pStar2.Ddimer[ind] <- p.Ddimer[2]*(1-p.Ddimer[4]*(z.Prof[ind]==3))
    pStar3.Ddimer[ind] <- p.Ddimer[3]*(1-p.Ddimer[4]*(z.Prof[ind]==3))
    wu.Ddimer[ind] <- (tilde.Age[ind]<0.81)*tilde.Age[ind]/0.81+(tilde.Age[ind]>=0.81)*(1-tilde.Age[ind])/(1-0.81)
    delta.Ddimer[ind] <- 1-(z.Grav[ind]>1)-Fibrin[ind]
        
    # Upper caval circle deep vein thrombosis (binary, non-stochastic)
    TVPsup[ind] ~ dbern(theta.TVPsup[ind])
    theta.TVPsup[ind] <- 1*(z.VenIntraVasCoa[ind]==2)  

    # Lower limbs deep vein thrombosis (binary, non-stochastic)
    TVPi[ind] ~ dbern(theta.TVPi[ind])  
    theta.TVPi[ind] <- 1*(z.VenIntraVasCoa[ind]==3)

    # Fibrinolysis (binary, non-stochastic)
    Fibrin[ind] ~ dbern(theta.Fibrin[ind])
    theta.Fibrin[ind] <- 1*(z.ArtIntraVasCoa[ind]>1||z.VenIntraVasCoa[ind]>1)

    # Venous intra-vascular coagulation (4 categories)
    #  [Calze decreases all conditional probabilities by a fixed quantity]
    z.VenIntraVasCoa[ind] ~ dcat(theta.VenIntraVasCoa[ind,1:4])
    theta.VenIntraVasCoa[ind,1] <- 1-(1-p.VenIntraVasCoa[41]*Calze[ind])*(1-thetaStar.VenIntraVasCoa[ind,1])
    theta.VenIntraVasCoa[ind,2] <- thetaStar.VenIntraVasCoa[ind,2]*(1-p.VenIntraVasCoa[41]*Calze[ind])
    theta.VenIntraVasCoa[ind,3] <- thetaStar.VenIntraVasCoa[ind,3]*(1-p.VenIntraVasCoa[41]*Calze[ind])
    theta.VenIntraVasCoa[ind,4] <- thetaStar.VenIntraVasCoa[ind,4]*(1-p.VenIntraVasCoa[41]*Calze[ind])
    for(aux in 1:4) {
      thetaStar.VenIntraVasCoa[ind,aux] <- (delta.VenIntraVasCoa[ind]<1)*odds.VenIntraVasCoa[ind,aux]/sum(odds.VenIntraVasCoa[ind,1:4])+(delta.VenIntraVasCoa[ind]>=1)*p.VenIntraVasCoa[aux]
      }
    odds.VenIntraVasCoa[ind,1] <- 1
    odds.VenIntraVasCoa[ind,2] <- exp(eta.VenIntraVasCoa[1]*delta.VenIntraVasCoa[ind]+
      log(p.VenIntraVasCoa[6]/p.VenIntraVasCoa[5])*Allett[ind]+
      log(p.VenIntraVasCoa[10]/p.VenIntraVasCoa[9])*(z.Cardiomio[ind]==2)+
      log(p.VenIntraVasCoa[14]/p.VenIntraVasCoa[13])*(z.Chir[ind]==2)+
      log(p.VenIntraVasCoa[18]/p.VenIntraVasCoa[17])*(z.Chir[ind]==3)+
      log(p.VenIntraVasCoa[22]/p.VenIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.VenIntraVasCoa[26]/p.VenIntraVasCoa[25])*CatVen[ind]+
      log(p.VenIntraVasCoa[30]/p.VenIntraVasCoa[29])*(z.Chir[ind]>1)+
      log(p.VenIntraVasCoa[34]/p.VenIntraVasCoa[33])*TVPscore[ind]+
      log(p.VenIntraVasCoa[38]/p.VenIntraVasCoa[37])*preTVP[ind]) 
    odds.VenIntraVasCoa[ind,3] <- exp(eta.VenIntraVasCoa[2]*delta.VenIntraVasCoa[ind]+
      log(p.VenIntraVasCoa[7]/p.VenIntraVasCoa[5])*Allett[ind]+
      log(p.VenIntraVasCoa[11]/p.VenIntraVasCoa[9])*(z.Cardiomio[ind]==2)+
      log(p.VenIntraVasCoa[15]/p.VenIntraVasCoa[13])*(z.Chir[ind]==2)+
      log(p.VenIntraVasCoa[19]/p.VenIntraVasCoa[17])*(z.Chir[ind]==3)+
      log(p.VenIntraVasCoa[23]/p.VenIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.VenIntraVasCoa[27]/p.VenIntraVasCoa[25])*CatVen[ind]+
      log(p.VenIntraVasCoa[31]/p.VenIntraVasCoa[29])*(z.Chir[ind]>1)+
      log(p.VenIntraVasCoa[35]/p.VenIntraVasCoa[33])*TVPscore[ind]+
      log(p.VenIntraVasCoa[39]/p.VenIntraVasCoa[37])*preTVP[ind]) 
    odds.VenIntraVasCoa[ind,4] <- exp(eta.VenIntraVasCoa[3]*delta.VenIntraVasCoa[ind]+
      log(p.VenIntraVasCoa[8]/p.VenIntraVasCoa[5])*Allett[ind]+
      log(p.VenIntraVasCoa[12]/p.VenIntraVasCoa[9])*(z.Cardiomio[ind]==2)+
      log(p.VenIntraVasCoa[16]/p.VenIntraVasCoa[13])*(z.Chir[ind]==2)+
      log(p.VenIntraVasCoa[20]/p.VenIntraVasCoa[17])*(z.Chir[ind]==3)+
      log(p.VenIntraVasCoa[24]/p.VenIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.VenIntraVasCoa[28]/p.VenIntraVasCoa[25])*CatVen[ind]+
      log(p.VenIntraVasCoa[32]/p.VenIntraVasCoa[29])*(z.Chir[ind]>1)+
      log(p.VenIntraVasCoa[36]/p.VenIntraVasCoa[33])*TVPscore[ind]+
      log(p.VenIntraVasCoa[40]/p.VenIntraVasCoa[37])*preTVP[ind]) 
    delta.VenIntraVasCoa[ind] <- 1-Allett[ind]-(z.Cardiomio[ind]==3)-(z.Chir[ind]>1)-(z.IntraVasCoa[ind]==2)-
      CatVen[ind]-(z.Grav[ind]>1)-TVPscore[ind]-preTVP[ind]

    # Arterial intra-vascular coagulation (5 categories)
    z.ArtIntraVasCoa[ind] ~ dcat(theta.ArtIntraVasCoa[ind,1:5])
    for(aux in 1:5) {
      theta.ArtIntraVasCoa[ind,aux] <- (delta.ArtIntraVasCoa[ind]<1)*odds.ArtIntraVasCoa[ind,aux]/sum(odds.ArtIntraVasCoa[ind,1:5])+(delta.ArtIntraVasCoa[ind]>=1)*p.ArtIntraVasCoa[aux]
      }
    odds.ArtIntraVasCoa[ind,1] <- 1
    odds.ArtIntraVasCoa[ind,2] <- exp(eta.ArtIntraVasCoa[1]*delta.ArtIntraVasCoa[ind]+
      log(p.ArtIntraVasCoa[7]/p.ArtIntraVasCoa[6])*(Arit_sopra_cron[ind]==1&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[12]/p.ArtIntraVasCoa[11])*Aort_diss[ind]+
      log(p.ArtIntraVasCoa[17]/p.ArtIntraVasCoa[16])*(z.Cardiomio[ind]==2&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[22]/p.ArtIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.ArtIntraVasCoa[27]/p.ArtIntraVasCoa[26])*Pancr[ind]) 
    odds.ArtIntraVasCoa[ind,3] <- exp(eta.ArtIntraVasCoa[2]*delta.ArtIntraVasCoa[ind]+
      log(p.ArtIntraVasCoa[8]/p.ArtIntraVasCoa[6])*(Arit_sopra_cron[ind]==1&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[13]/p.ArtIntraVasCoa[11])*Aort_diss[ind]+
      log(p.ArtIntraVasCoa[18]/p.ArtIntraVasCoa[16])*(z.Cardiomio[ind]==2&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[23]/p.ArtIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.ArtIntraVasCoa[28]/p.ArtIntraVasCoa[26])*Pancr[ind]) 
    odds.ArtIntraVasCoa[ind,4] <- exp(eta.ArtIntraVasCoa[3]*delta.ArtIntraVasCoa[ind]+
      log(p.ArtIntraVasCoa[9]/p.ArtIntraVasCoa[6])*(Arit_sopra_cron[ind]==1&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[14]/p.ArtIntraVasCoa[11])*Aort_diss[ind]+
      log(p.ArtIntraVasCoa[19]/p.ArtIntraVasCoa[16])*(z.Cardiomio[ind]==2&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[24]/p.ArtIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.ArtIntraVasCoa[29]/p.ArtIntraVasCoa[26])*Pancr[ind]) 
    odds.ArtIntraVasCoa[ind,5] <- exp(eta.ArtIntraVasCoa[4]*delta.ArtIntraVasCoa[ind]+
      log(p.ArtIntraVasCoa[10]/p.ArtIntraVasCoa[6])*(Arit_sopra_cron[ind]==1&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[15]/p.ArtIntraVasCoa[11])*Aort_diss[ind]+
      log(p.ArtIntraVasCoa[20]/p.ArtIntraVasCoa[16])*(z.Cardiomio[ind]==2&&z.IntraVasCoa[ind]<3)+
      log(p.ArtIntraVasCoa[25]/p.ArtIntraVasCoa[21])*(z.IntraVasCoa[ind]==2)+
      log(p.ArtIntraVasCoa[30]/p.ArtIntraVasCoa[26])*Pancr[ind]) 
    delta.ArtIntraVasCoa[ind] <- 1-(Arit_sopra_cron[ind]==1&&z.IntraVasCoa[ind]<3)-Aort_diss[ind]-
      (z.Cardiomio[ind]==2&&z.IntraVasCoa[ind]<3)-(z.IntraVasCoa[ind]==3)-Pancr[ind]
    
    # Dehydration (according to clinical judgement) (binary, non-stochastic)
    Disidrat_J[ind] ~ dbern(theta.Disidrat_J[ind])
    theta.Disidrat_J[ind] <- 1*(tilde.Disidrat[ind]>=2/3)

    # Myocarditis (binary)
    Miocard[ind] ~ dbern(theta.Miocard[ind])
    theta.Miocard[ind] <- (delta.Miocard[ind]<1)*thetaStar.Miocard[ind]+(delta.Miocard[ind]>=1)*p.Miocard[1]
    logit(thetaStar.Miocard[ind]) <- eta.Miocard[1]*delta.Miocard[ind]+
      qlogis(p.Miocard[2],0,1)*Peric_ni[ind]+
      qlogis(p.Miocard[3],0,1)*Sepsi[ind]+
      qlogis(p.Miocard[4],0,1)*Endocard[ind]
    delta.Miocard[ind] <- 1-Peric_ni[ind]-Sepsi[ind]-Endocard[ind]

    # Intra-vascular coagulation (3 categories)
    z.IntraVasCoa[ind] ~ dcat(theta.IntraVasCoa[ind,1:3]) 
    for(aux in 1:3) {
      theta.IntraVasCoa[ind,aux] <- (delta.IntraVasCoa[ind]<1)*odds.IntraVasCoa[ind,aux]/sum(odds.IntraVasCoa[ind,1:3])+(delta.IntraVasCoa[ind]>=1)*p.IntraVasCoa[aux]
      }
    odds.IntraVasCoa[ind,1] <- 1
    odds.IntraVasCoa[ind,2] <- exp(eta.IntraVasCoa[1]*delta.IntraVasCoa[ind]+
      log(p.IntraVasCoa[5]/p.IntraVasCoa[4])*Estrog[ind]+
      log(p.IntraVasCoa[8]/p.IntraVasCoa[7])*Neo_o[ind]+
      log(p.IntraVasCoa[11]/p.IntraVasCoa[10])*Obes[ind]+
      log(p.IntraVasCoa[14]/p.IntraVasCoa[13])*Sepsi[ind]+
      log(p.IntraVasCoa[17]/p.IntraVasCoa[16])*Trmphil[ind]+
      log(p.IntraVasCoa[20]/p.IntraVasCoa[19])*(z.Prof[ind]==2)+      
      log(p.IntraVasCoa[23]/p.IntraVasCoa[22])*(z.Prof[ind]==3)) 
    odds.IntraVasCoa[ind,3] <- exp(eta.IntraVasCoa[2]*delta.IntraVasCoa[ind]+
      log(p.IntraVasCoa[6]/p.IntraVasCoa[4])*Estrog[ind]+
      log(p.IntraVasCoa[9]/p.IntraVasCoa[7])*Neo_o[ind]+
      log(p.IntraVasCoa[12]/p.IntraVasCoa[10])*Obes[ind]+
      log(p.IntraVasCoa[15]/p.IntraVasCoa[13])*Sepsi[ind]+
      log(p.IntraVasCoa[18]/p.IntraVasCoa[16])*Trmphil[ind]+                                                                                             
      log(p.IntraVasCoa[21]/p.IntraVasCoa[19])*(z.Prof[ind]==2)+      
      log(p.IntraVasCoa[24]/p.IntraVasCoa[22])*(z.Prof[ind]==3)) 
    delta.IntraVasCoa[ind] <- 1-Estrog[ind]-Neo_o[ind]-Obes[ind]-Sepsi[ind]-Trmphil[ind]-(z.Prof[ind]>1)

    # Dehydration (continuous, hyper-restricted)
    tilde.Disidrat[ind] ~ dbeta(relmu.Disidrat[ind]*p.Disidrat[4],(1-relmu.Disidrat[ind])*p.Disidrat[4]) T(tol,1-tol)
    logit(relmu.Disidrat[ind]) <- qlogis(p.Disidrat[1],0,1)*delta.Disidrat[ind]+
      qlogis(p.Disidrat[2],0,1)*Pancr[ind]+
      qlogis(p.Disidrat[3],0,1)*Sepsi[ind]
    delta.Disidrat[ind] <- 1-Pancr[ind]-Sepsi[ind]
    
    # Sepsis (binary)
    Sepsi[ind] ~ dbern(theta.Sepsi[ind])
    theta.Sepsi[ind] <- (delta.Sepsi[ind]<1)*thetaStar.Sepsi[ind]+(delta.Sepsi[ind]>=1)*p.Sepsi[1]
    logit(thetaStar.Sepsi[ind]) <- eta.Sepsi[1]*delta.Sepsi[ind]+
      qlogis(p.Sepsi[2],0,1)*(z.Polm[ind]>1)+
      qlogis(p.Sepsi[3],0,1)*Perit[ind]
    delta.Sepsi[ind] <- 1-(z.Polm[ind]>1)-Perit[ind]
    
    # auxiliar node    
    Polm_J[ind] ~ dbern(theta.Polm_J[ind])
    theta.Polm_J[ind] <- 1*(z.Polm[ind]>1)
   
    # Chest pain (parietal origin) (binary)
    D_par[ind] ~ dbern(theta.D_par[ind]) 
    theta.D_par[ind] <- (delta.D_par[ind]<1)*thetaStar.D_par[ind]+(delta.D_par[ind]>=1)*p.D_par[1]
    logit(thetaStar.D_par[ind]) <- eta.D_par[1]*delta.D_par[ind]+
      qlogis(p.D_par[2],0,1)*Costcndr[ind]+
      qlogis(p.D_par[3],0,1)*Fr[ind]+
      qlogis(p.D_par[4],0,1)*Fr_cost[ind]+
      qlogis(p.D_par[5],0,1)*HrepZ[ind]+
      qlogis(p.D_par[6],0,1)*Neo_p[ind]+
      qlogis(p.D_par[7],0,1)*(z.Polm[ind]>1)+
      qlogis(p.D_par[8],0,1)*(z.Pnx[ind]>1)
    delta.D_par[ind] <- 1-Costcndr[ind]-Fr[ind]-Fr_cost[ind]-HrepZ[ind]-Neo_p[ind]-(z.Polm[ind]>1)-(z.Pnx[ind]>1)

    # Cavitation/colliquation (binary)
    CavColl[ind] ~ dbern(theta.CavColl[ind])
    theta.CavColl[ind] <- (delta.CavColl[ind]<1)*thetaStar.CavColl[ind]+(delta.CavColl[ind]>=1)*p.CavColl[1]
    logit(thetaStar.CavColl[ind]) <- eta.CavColl[1]*delta.CavColl[ind]+
      qlogis(p.CavColl[2],0,1)*(z.Polm[ind]==2)+
      qlogis(p.CavColl[3],0,1)*(z.Polm[ind]==3)
    delta.CavColl[ind] <- 1-(z.Polm[ind]>1)
    
    # Hilar adenopathy (binary)
    AdenIlar[ind] ~ dbern(theta.AdenIlar[ind])
    theta.AdenIlar[ind] <- (delta.AdenIlar[ind]<1)*thetaStar.AdenIlar[ind]+(delta.AdenIlar[ind]>=1)*p.AdenIlar[1]
    logit(thetaStar.AdenIlar[ind]) <- eta.AdenIlar[1]*delta.AdenIlar[ind]+
      qlogis(p.AdenIlar[2],0,1)*Neo_p[ind]+
      qlogis(p.AdenIlar[3],0,1)*Interst[ind]+
      qlogis(p.AdenIlar[4],0,1)*(z.Polm[ind]==2)+
      qlogis(p.AdenIlar[5],0,1)*(z.Polm[ind]==3)
    delta.AdenIlar[ind] <- 1-Neo_p[ind]-Interst[ind]-(z.Polm[ind]>1)

    # Prophylaxis/anticoagulation  (3 categories)
    z.Prof[ind] ~ dcat(theta.Prof[ind,1:3]) 
    for(aux in 1:3) {
      theta.Prof[ind,aux] <- (delta.Prof[ind]<1)*odds.Prof[ind,aux]/sum(odds.Prof[ind,1:3])+(delta.Prof[ind]>=1)*p.Prof[aux]
      }
    odds.Prof[ind,1] <- 1
    odds.Prof[ind,2] <- exp(eta.Prof[1]*delta.Prof[ind]+
      log(p.Prof[5]/p.Prof[4])*Arit_sopra_cron[ind]+
      log(p.Prof[8]/p.Prof[7])*(z.Chir[ind]>1)+
      log(p.Prof[11]/p.Prof[10])*preTVP[ind])
    odds.Prof[ind,3] <- exp(eta.Prof[2]*delta.Prof[ind]+
      log(p.Prof[6]/p.Prof[4])*Arit_sopra_cron[ind]+
      log(p.Prof[9]/p.Prof[7])*(z.Chir[ind]>1)+
      log(p.Prof[12]/p.Prof[10])*preTVP[ind])
    delta.Prof[ind] <- 1-Arit_sopra_cron[ind]-(z.Chir[ind]>1)-preTVP[ind]

    # Pneumonia (3 categories with age and gender)
    z.Polm[ind] ~ dcat(theta.Polm[ind,1:3])
    for(aux in 1:3) {
      theta.Polm[ind,aux] <- (delta.Polm[ind]<1)*odds.Polm[ind,aux]/sum(odds.Polm[ind,1:3])+(delta.Polm[ind]>=1)*p0.Polm[ind,aux]
      }
    odds.Polm[ind,1] <- 1
    odds.Polm[ind,2] <- exp(eta0.Polm[ind,1]*delta.Polm[ind]+
      log(p.Polm[8]/p.Polm[7])*(z.Enf[ind]==2)+
      log(p.Polm[11]/p.Polm[10])*(z.Enf[ind]==3)+
      log(p.Polm[14]/p.Polm[13])*Neo_p[ind]+
      log(p.Polm[17]/p.Polm[16])*(z.Cardiomio[ind]==3)+
      log(p.Polm[20]/p.Polm[19])*Defed[ind]+
      log(p.Polm[23]/p.Polm[22])*Alcol[ind]+
      log(p.Polm[26]/p.Polm[25])*IschCerCron[ind])      
    odds.Polm[ind,3] <- exp(eta0.Polm[ind,2]*delta.Polm[ind]+
      log(p.Polm[9]/p.Polm[7])*(z.Enf[ind]==2)+
      log(p.Polm[12]/p.Polm[10])*(z.Enf[ind]==3)+
      log(p.Polm[15]/p.Polm[13])*Neo_p[ind]+
      log(p.Polm[18]/p.Polm[16])*(z.Cardiomio[ind]==3)+
      log(p.Polm[21]/p.Polm[19])*Defed[ind]+
      log(p.Polm[24]/p.Polm[22])*Alcol[ind]+
      log(p.Polm[27]/p.Polm[25])*IschCerCron[ind])
    p0.Polm[ind,1] <- 1-wu.Polm[ind]*(1-(Sex[ind]==0)*p.Polm[1]-(Sex[ind]==1)*p.Polm[4])
    p0.Polm[ind,2] <- ((Sex[ind]==0)*p.Polm[2]+(Sex[ind]==1)*p.Polm[5])*wu.Polm[ind]
    p0.Polm[ind,3] <- ((Sex[ind]==0)*p.Polm[3]+(Sex[ind]==1)*p.Polm[6])*wu.Polm[ind]
    eta0.Polm[ind,1] <- (Sex[ind]==0)*log(wu.Polm[ind]*exp(eta.Polm[1])/(1+(1-wu.Polm[ind])*(exp(eta.Polm[1])+exp(eta.Polm[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Polm[ind]*exp(eta.Polm[2])/(1+(1-wu.Polm[ind])*(exp(eta.Polm[2])+exp(eta.Polm[4]))))
    eta0.Polm[ind,2] <- (Sex[ind]==0)*log(wu.Polm[ind]*exp(eta.Polm[3])/(1+(1-wu.Polm[ind])*(exp(eta.Polm[1])+exp(eta.Polm[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Polm[ind]*exp(eta.Polm[4])/(1+(1-wu.Polm[ind])*(exp(eta.Polm[2])+exp(eta.Polm[4]))))
    wu.Polm[ind] <- (tilde.Age[ind]<0.86)*tilde.Age[ind]/0.86+(tilde.Age[ind]>=0.86)*(1-tilde.Age[ind])/(1-0.86)
    delta.Polm[ind] <- 1-(z.Enf[ind]>1)-Neo_p[ind]-(z.Cardiomio[ind]==3)-Defed[ind]-Alcol[ind]-IschCerCron[ind]

    # Chronic atrial arrhythmia (binary with age)
    Arit_sopra_cron[ind] ~ dbern(theta.Arit_sopra_cron[ind])  
    theta.Arit_sopra_cron[ind] <- (delta.Arit_sopra_cron[ind]<1)*thetaStar.Arit_sopra_cron[ind]+(delta.Arit_sopra_cron[ind]>=1)*wu.Arit_sopra_cron[ind]*p.Arit_sopra_cron[1]
    logit(thetaStar.Arit_sopra_cron[ind]) <- qlogis(plogis(eta.Arit_sopra_cron[1],0,1)*wu.Arit_sopra_cron[ind],0,1)*delta.Arit_sopra_cron[ind]+
      qlogis(p.Arit_sopra_cron[2],0,1)*(z.Cardiomio[ind]==2)+
      qlogis(p.Arit_sopra_cron[3],0,1)*(z.Cardiomio[ind]==3)+
      qlogis(p.Arit_sopra_cron[4],0,1)*(z.Enf[ind]>1)
    wu.Arit_sopra_cron[ind] <- (tilde.Age[ind]<0.76)*tilde.Age[ind]/0.76+(tilde.Age[ind]>=0.76)*(1-tilde.Age[ind])/(1-0.76)
    delta.Arit_sopra_cron[ind] <- 1-(z.Cardiomio[ind]>1)-(z.Enf[ind]>1)
 
    # Chronic cardiac muscle disease (3 categories with age and gender)
    z.Cardiomio[ind] ~ dcat(theta.Cardiomio[ind,1:3])
    for(aux in 1:3) {
      theta.Cardiomio[ind,aux] <- (delta.Cardiomio[ind]<1)*odds.Cardiomio[ind,aux]/sum(odds.Cardiomio[ind,1:3])+(delta.Cardiomio[ind]>=1)*p0.Cardiomio[ind,aux]
      }
    odds.Cardiomio[ind,1] <- 1
    odds.Cardiomio[ind,2] <- exp(eta0.Cardiomio[ind,1]*delta.Cardiomio[ind]+
      log(p.Cardiomio[8]/p.Cardiomio[7])*Card_dil[ind]+
      log(p.Cardiomio[11]/p.Cardiomio[10])*cronIpert[ind]+
      log(p.Cardiomio[14]/p.Cardiomio[13])*CuPC[ind]+
      log(p.Cardiomio[17]/p.Cardiomio[16])*IpertSx[ind])      
    odds.Cardiomio[ind,3] <- exp(eta0.Cardiomio[ind,2]*delta.Cardiomio[ind]+    
      log(p.Cardiomio[9]/p.Cardiomio[7])*Card_dil[ind]+
      log(p.Cardiomio[12]/p.Cardiomio[10])*cronIpert[ind]+
      log(p.Cardiomio[15]/p.Cardiomio[13])*CuPC[ind]+
      log(p.Cardiomio[18]/p.Cardiomio[16])*IpertSx[ind])
    p0.Cardiomio[ind,1] <- 1-wu.Cardiomio[ind]*(1-(Sex[ind]==0)*p.Cardiomio[1]-(Sex[ind]==1)*p.Cardiomio[4])
    p0.Cardiomio[ind,2] <- ((Sex[ind]==0)*p.Cardiomio[2]+(Sex[ind]==1)*p.Cardiomio[5])*wu.Cardiomio[ind]
    p0.Cardiomio[ind,3] <- ((Sex[ind]==0)*p.Cardiomio[3]+(Sex[ind]==1)*p.Cardiomio[6])*wu.Cardiomio[ind]
    eta0.Cardiomio[ind,1] <- (Sex[ind]==0)*log(wu.Cardiomio[ind]*exp(eta.Cardiomio[1])/(1+(1-wu.Cardiomio[ind])*(exp(eta.Cardiomio[1])+exp(eta.Cardiomio[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Cardiomio[ind]*exp(eta.Cardiomio[2])/(1+(1-wu.Cardiomio[ind])*(exp(eta.Cardiomio[2])+exp(eta.Cardiomio[4]))))
    eta0.Cardiomio[ind,2] <- (Sex[ind]==0)*log(wu.Cardiomio[ind]*exp(eta.Cardiomio[3])/(1+(1-wu.Cardiomio[ind])*(exp(eta.Cardiomio[1])+exp(eta.Cardiomio[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Cardiomio[ind]*exp(eta.Cardiomio[4])/(1+(1-wu.Cardiomio[ind])*(exp(eta.Cardiomio[2])+exp(eta.Cardiomio[4]))))
    wu.Cardiomio[ind] <- (tilde.Age[ind]<0.57)*tilde.Age[ind]/0.57+(tilde.Age[ind]>=0.57)*(1-tilde.Age[ind])/(1-0.57)
    delta.Cardiomio[ind] <- 1-Card_dil[ind]-cronIpert[ind]-CuPC[ind]-IpertSx[ind]
        
    # Obstructive cardiomyopathy (3 categories with age)
    z.Card_obstr[ind] ~ dcat(theta.Card_obstr[ind,1:3]) 
    for(aux in 1:3) {
      theta.Card_obstr[ind,aux] <- (delta.Card_obstr[ind]<1)*odds.Card_obstr[ind,aux]/sum(odds.Card_obstr[ind,1:3])+(delta.Card_obstr[ind]>=1)*p0.Card_obstr[ind,aux]          
      }
    odds.Card_obstr[ind,1] <- 1
    odds.Card_obstr[ind,2] <- exp(eta0.Card_obstr[ind,1]*delta.Card_obstr[ind]+
      log(p.Card_obstr[5]/p.Card_obstr[4])*IpertSx[ind]+
      log(p.Card_obstr[8]/p.Card_obstr[7])*Sten_aor[ind])      
    odds.Card_obstr[ind,3] <- exp(eta0.Card_obstr[ind,2]*delta.Card_obstr[ind]+
      log(p.Card_obstr[6]/p.Card_obstr[4])*IpertSx[ind]+
      log(p.Card_obstr[9]/p.Card_obstr[7])*Sten_aor[ind])
    eta0.Card_obstr[ind,1] <- log(wu.Card_obstr[ind]*exp(eta.Card_obstr[1])/(1+(1-wu.Card_obstr[ind])*(exp(eta.Card_obstr[1])+exp(eta.Card_obstr[2]))))
    eta0.Card_obstr[ind,2] <- log(wu.Card_obstr[ind]*exp(eta.Card_obstr[2])/(1+(1-wu.Card_obstr[ind])*(exp(eta.Card_obstr[1])+exp(eta.Card_obstr[2]))))
    p0.Card_obstr[ind,1] <- 1-wu.Card_obstr[ind]*(1-p.Card_obstr[1])
    p0.Card_obstr[ind,2] <- wu.Card_obstr[ind]*p.Card_obstr[2]
    p0.Card_obstr[ind,3] <- wu.Card_obstr[ind]*p.Card_obstr[3]
    wu.Card_obstr[ind] <- (tilde.Age[ind]<0.24)*tilde.Age[ind]/0.24+(tilde.Age[ind]>=0.24)*(1-tilde.Age[ind])/(1-0.24)
    delta.Card_obstr[ind] <- 1-IpertSx[ind]-Sten_aor[ind]
    
    # Hemoglobin (continuous, hyper-restricted)
    tilde.Emoglob[ind] ~ dbeta(relmu.Emoglob[ind]*p.Emoglob[4],(1-relmu.Emoglob[ind])*p.Emoglob[4]) T(tol,1-tol)
    logit(relmu.Emoglob[ind]) <- qlogis(p.Emoglob[1],0,1)*delta.Emoglob[ind]+
      qlogis(p.Emoglob[2],0,1)*acuteAnem[ind]+
      qlogis(p.Emoglob[3],0,1)*cronAnem[ind]
    delta.Emoglob[ind] <- 1-acuteAnem[ind]-cronAnem[ind]

    # Left ventricular hypertrophy (binary with age)
    IpertSx[ind] ~ dbern(theta.IpertSx[ind]) 
    theta.IpertSx[ind] <- (delta.IpertSx[ind]<1)*thetaStar.IpertSx[ind]+(delta.IpertSx[ind]>=1)*wu.IpertSx[ind]*p.IpertSx[1]    
    logit(thetaStar.IpertSx[ind]) <- qlogis(plogis(eta.IpertSx[1],0,1)*wu.IpertSx[ind],0,1)*delta.IpertSx[ind]+
      qlogis(p.IpertSx[2],0,1)*cronIpert[ind]+
      qlogis(p.IpertSx[3],0,1)*(z.Insuf_aorCron[ind]==3)+
      qlogis(p.IpertSx[4],0,1)*(z.InsuffMitrCron[ind]==3)
    wu.IpertSx[ind] <- (tilde.Age[ind]<0.86)*tilde.Age[ind]/0.86+(tilde.Age[ind]>=0.86)*(1-tilde.Age[ind])/(1-0.86)
    delta.IpertSx[ind] <- 1-cronIpert[ind]-(z.Insuf_aorCron[ind]==3)-(z.InsuffMitrCron[ind]==3)
    
    # Aortic valve failure (generic) (binary, non-stochastic)
    Insuf_aor[ind] ~ dbern(theta.Insuf_aor[ind])
    theta.Insuf_aor[ind] <- 1*(z.Insuf_aorCron[ind]>1||Insuf_aorAcut[ind]==1)

    # Ventricular moderator band thickness (binary)
    Ip_ban_mod[ind] ~ dbern(theta.Ip_ban_mod[ind])
    theta.Ip_ban_mod[ind] <- (delta.Ip_ban_mod[ind]<1)*thetaStar.Ip_ban_mod[ind]+(delta.Ip_ban_mod[ind]>=1)*p.Ip_ban_mod[1]
    logit(thetaStar.Ip_ban_mod[ind]) <- eta.Ip_ban_mod[1]*delta.Ip_ban_mod[ind]+
      qlogis(p.Ip_ban_mod[2],0,1)*CuPC[ind]
    delta.Ip_ban_mod[ind] <- 1-CuPC[ind]

    # Left ventricular thickness (binary)
    I_par_lib[ind] ~ dbern(theta.I_par_lib[ind])
    theta.I_par_lib[ind] <- (delta.I_par_lib[ind]<1)*thetaStar.I_par_lib[ind]+(delta.I_par_lib[ind]>=1)*p.I_par_lib[1]
    logit(thetaStar.I_par_lib[ind]) <- eta.I_par_lib[1]*delta.I_par_lib[ind]+
      qlogis(p.I_par_lib[2],0,1)*CuPC[ind]   
    delta.I_par_lib[ind] <- 1-CuPC[ind]
    
    # auxiliar node
    Pnx_J[ind] ~ dbern(theta.Pnx_J[ind])
    theta.Pnx_J[ind] <- 1*(z.Pnx[ind]>1)

    # Pulmonary hyperlucency (3 categories)
    z.Iper_Diaf[ind] ~ dcat(theta.Iper_Diaf[ind,1:3])  
    for(aux in 1:3) {
      theta.Iper_Diaf[ind,aux] <- (delta.Iper_Diaf[ind]<1)*odds.Iper_Diaf[ind,aux]/sum(odds.Iper_Diaf[ind,1:3])+(delta.Iper_Diaf[ind]>=1)*p.Iper_Diaf[aux]
      }
    odds.Iper_Diaf[ind,1] <- 1           
    odds.Iper_Diaf[ind,2] <- exp(eta.Iper_Diaf[1]*delta.Iper_Diaf[ind]+
      log(p.Iper_Diaf[5]/p.Iper_Diaf[4])*(z.Enf[ind]==2&&z.Pnx[ind]==1)+
      log(p.Iper_Diaf[8]/p.Iper_Diaf[7])*(z.Enf[ind]==3&&z.Pnx[ind]==1)+
      log(p.Iper_Diaf[11]/p.Iper_Diaf[10])*(z.Pnx[ind]>1))      
    odds.Iper_Diaf[ind,3] <- exp(eta.Iper_Diaf[2]*delta.Iper_Diaf[ind]+
      log(p.Iper_Diaf[6]/p.Iper_Diaf[4])*(z.Enf[ind]==2&&z.Pnx[ind]==1)+
      log(p.Iper_Diaf[9]/p.Iper_Diaf[7])*(z.Enf[ind]==3&&z.Pnx[ind]==1)+
      log(p.Iper_Diaf[12]/p.Iper_Diaf[10])*(z.Pnx[ind]>1))  
    delta.Iper_Diaf[ind] <- 1-(z.Enf[ind]==2&&z.Pnx[ind]==1)-(z.Enf[ind]==3&&z.Pnx[ind]==1)-(z.Pnx[ind]>1)

    # Chest pain (stabbing origin) (binary)
    D_a_pu[ind] ~ dbern(theta.D_a_pu[ind])
    theta.D_a_pu[ind] <- (delta.D_a_pu[ind]<1)*thetaStar.D_a_pu[ind]+(delta.D_a_pu[ind]>=1)*p.D_a_pu[1]
    logit(thetaStar.D_a_pu[ind]) <- eta.D_a_pu[1]*delta.D_a_pu[ind]+
      qlogis(p.D_a_pu[2],0,1)*Aort_diss[ind]+
      qlogis(p.D_a_pu[3],0,1)*(z.Pnx[ind]>1)
    delta.D_a_pu[ind] <- 1-Aort_diss[ind]-(z.Pnx[ind]>1)
    
    # Acute anemia (binary, semi-stochastic)
    acuteAnem[ind] ~ dbern(theta.acuteAnem[ind])
    theta.acuteAnem[ind] <- p.acuteAnem[1]*(Hemo[ind]==0)+1*(Hemo[ind]==1)
    
    # Chest pain (upper-abdominal origin) (binary)
    D_quadr_sup[ind] ~ dbern(theta.D_quadr_sup[ind]) 
    theta.D_quadr_sup[ind] <- 1*(Perit[ind]==1)+(Perit[ind]==0)*((delta.D_quadr_sup[ind]<1)*thetaStar.D_quadr_sup[ind]+(delta.D_quadr_sup[ind]>=1)*p.D_quadr_sup[1])
    logit(thetaStar.D_quadr_sup[ind]) <- eta.D_quadr_sup[1]*delta.D_quadr_sup[ind]+
      qlogis(p.D_quadr_sup[2],0,1)*Colic[ind]+
      qlogis(p.D_quadr_sup[3],0,1)*Pancr[ind]+
      qlogis(p.D_quadr_sup[4],0,1)*Spas_es[ind]+
      qlogis(p.D_quadr_sup[5],0,1)*Ulc_pep[ind]
    delta.D_quadr_sup[ind] <- 1-Colic[ind]-Pancr[ind]-Spas_es[ind]-Ulc_pep[ind]-Perit[ind]

    # Chronic mitral valve failure (3 categories with age)
    z.InsuffMitrCron[ind] ~ dcat(theta.InsuffMitrCron[ind,1:3])
    for(aux in 1:3) {
      theta.InsuffMitrCron[ind,aux] <- (delta.InsuffMitrCron[ind]<1)*odds.InsuffMitrCron[ind,aux]/sum(odds.InsuffMitrCron[ind,1:3])+(delta.InsuffMitrCron[ind]>=1)*p0.InsuffMitrCron[ind,aux]
      }
    odds.InsuffMitrCron[ind,1] <- 1
    odds.InsuffMitrCron[ind,2] <- exp(eta0.InsuffMitrCron[ind,1]*delta.InsuffMitrCron[ind]+
      log(p.InsuffMitrCron[5]/p.InsuffMitrCron[4])*Card_dil[ind]+
      log(p.InsuffMitrCron[8]/p.InsuffMitrCron[7])*(z.Insuf_aorCron[ind]==2)+
      log(p.InsuffMitrCron[11]/p.InsuffMitrCron[10])*(z.Insuf_aorCron[ind]==3)+
      log(p.InsuffMitrCron[14]/p.InsuffMitrCron[13])*ProlMitralCron[ind]+
      log(p.InsuffMitrCron[17]/p.InsuffMitrCron[16])*Sten_aor[ind])      
    odds.InsuffMitrCron[ind,3] <- exp(eta0.InsuffMitrCron[ind,2]*delta.InsuffMitrCron[ind]+
      log(p.InsuffMitrCron[6]/p.InsuffMitrCron[4])*Card_dil[ind]+
      log(p.InsuffMitrCron[9]/p.InsuffMitrCron[7])*(z.Insuf_aorCron[ind]==2)+
      log(p.InsuffMitrCron[12]/p.InsuffMitrCron[10])*(z.Insuf_aorCron[ind]==3)+
      log(p.InsuffMitrCron[15]/p.InsuffMitrCron[13])*ProlMitralCron[ind]+
      log(p.InsuffMitrCron[18]/p.InsuffMitrCron[16])*Sten_aor[ind])    
    eta0.InsuffMitrCron[ind,1] <- log(wu.InsuffMitrCron[ind]*exp(eta.InsuffMitrCron[1])/(1+(1-wu.InsuffMitrCron[ind])*(exp(eta.InsuffMitrCron[1])+exp(eta.InsuffMitrCron[2]))))
    eta0.InsuffMitrCron[ind,2] <- log(wu.InsuffMitrCron[ind]*exp(eta.InsuffMitrCron[2])/(1+(1-wu.InsuffMitrCron[ind])*(exp(eta.InsuffMitrCron[1])+exp(eta.InsuffMitrCron[2]))))
    p0.InsuffMitrCron[ind,1] <- 1-wu.InsuffMitrCron[ind]*(1-p.InsuffMitrCron[1])
    p0.InsuffMitrCron[ind,2] <- wu.InsuffMitrCron[ind]*p.InsuffMitrCron[2]
    p0.InsuffMitrCron[ind,3] <- wu.InsuffMitrCron[ind]*p.InsuffMitrCron[3]
    wu.InsuffMitrCron[ind] <- (tilde.Age[ind]<0.76)*tilde.Age[ind]/0.76+(tilde.Age[ind]>=0.76)*(1-tilde.Age[ind])/(1-0.76)
    delta.InsuffMitrCron[ind] <- 1-Card_dil[ind]-(z.Insuf_aorCron[ind]>1)-ProlMitralCron[ind]-Sten_aor[ind]

    # Acute aortic valve failure (binary)
    Insuf_aorAcut[ind] ~ dbern(theta.Insuf_aorAcut[ind])
    theta.Insuf_aorAcut[ind] <- (delta.Insuf_aorAcut[ind]<1)*thetaStar.Insuf_aorAcut[ind]+(delta.Insuf_aorAcut[ind]>=1)*p.Insuf_aorAcut[1]
    logit(thetaStar.Insuf_aorAcut[ind]) <- eta.Insuf_aorAcut[1]*delta.Insuf_aorAcut[ind]+
      qlogis(p.Insuf_aorAcut[2],0,1)*Aort_diss[ind]+
      qlogis(p.Insuf_aorAcut[3],0,1)*Endocard[ind]+
      qlogis(p.Insuf_aorAcut[4],0,1)*Card_dil[ind]
    delta.Insuf_aorAcut[ind] <- 1-Aort_diss[ind]-Endocard[ind]-Card_dil[ind]

    # Ascending aorta intimal flap (binary, non-stochastic)
    Flap[ind] ~ dbern(theta.Flap[ind])
    theta.Flap[ind] <- 1*(Aort_diss[ind]==1)

    # Aortic intramural hematoma (binary)
    EmatIntram[ind] ~ dbern(theta.EmatIntram[ind])
    theta.EmatIntram[ind] <- (delta.EmatIntram[ind]<1)*thetaStar.EmatIntram[ind]+(delta.EmatIntram[ind]>=1)*p.EmatIntram[1]
    logit(thetaStar.EmatIntram[ind]) <- eta.EmatIntram[1]*delta.EmatIntram[ind]+
      qlogis(p.EmatIntram[2],0,1)*Aort_diss[ind]
    delta.EmatIntram[ind] <- 1-Aort_diss[ind]

    # Anisosfigmia (binary)
    Anisosfig[ind] ~ dbern(theta.Anisosfig[ind])
    theta.Anisosfig[ind] <- (delta.Anisosfig[ind]<1)*thetaStar.Anisosfig[ind]+(delta.Anisosfig[ind]>=1)*p.Anisosfig[1]
    logit(thetaStar.Anisosfig[ind]) <- eta.Anisosfig[1]*delta.Anisosfig[ind]+
      qlogis(p.Anisosfig[2],0,1)*Aort_diss[ind]    
    delta.Anisosfig[ind] <- 1-Aort_diss[ind]

    # Cor pulmonale (binary)
    CuPC[ind] ~ dbern(theta.CuPC[ind])
    theta.CuPC[ind] <- (delta.CuPC[ind]<1)*thetaStar.CuPC[ind]+(delta.CuPC[ind]>=1)*p.CuPC[1]
    logit(thetaStar.CuPC[ind]) <- eta.CuPC[1]*delta.CuPC[ind]+
      qlogis(p.CuPC[2],0,1)*Ipertpolm[ind]   
    delta.CuPC[ind] <- 1-Ipertpolm[ind]

    # Thyrotoxicosis (binary)
    Tireotos[ind] ~ dbern(theta.Tireotos[ind])
    theta.Tireotos[ind] <- (delta.Tireotos[ind]<1)*thetaStar.Tireotos[ind]+(delta.Tireotos[ind]>=1)*p.Tireotos[1]
    logit(thetaStar.Tireotos[ind]) <- eta.Tireotos[1]*delta.Tireotos[ind]+
      qlogis(p.Tireotos[2],0,1)*(-1.5+3*tilde.OT[ind])
    delta.Tireotos[ind] <- 1-(-1.5+3*tilde.OT[ind])
    
    # FT4 (continuous)
    tilde.FT4[ind] ~ dbeta(relmu.FT4[ind]*p.FT4[4],(1-relmu.FT4[ind])*p.FT4[4]) T(tol,1-tol)
    logit(relmu.FT4[ind]) <- qlogis(p.FT4[1],0,1)*delta.FT4[ind]+
      qlogis(p.FT4[2],0,1)*(tilde.OT[ind]<0.5)*(-1.5+3*tilde.OT[ind])+
      qlogis(p.FT4[3],0,1)*(tilde.OT[ind]>0.5)*(-1.5+3*tilde.OT[ind])
    delta.FT4[ind] <- 1-(-1.5+3*tilde.OT[ind])
        
    # FT3 (continuous)
    tilde.FT3[ind] ~ dbeta(relmu.FT3[ind]*p.FT3[4],(1-relmu.FT3[ind])*p.FT3[4]) T(tol,1-tol)
    logit(relmu.FT3[ind]) <- qlogis(p.FT3[1],0,1)*delta.FT3[ind]+
      qlogis(p.FT3[2],0,1)*(tilde.OT[ind]<0.5)*(-1.5+3*tilde.OT[ind])+
      qlogis(p.FT3[3],0,1)*(tilde.OT[ind]>0.5)*(-1.5+3*tilde.OT[ind])
    delta.FT3[ind] <- 1-(-1.5+3*tilde.OT[ind])

    # Spontaneous pneumothorax (3 categories with age)
    z.Pnx[ind] ~ dcat(theta.Pnx[ind,1:3])   
    for(aux in 1:3) {
      theta.Pnx[ind,aux] <- (delta.Pnx[ind]<1)*odds.Pnx[ind,aux]/sum(odds.Pnx[ind,1:3])+(delta.Pnx[ind]>=1)*p0.Pnx[ind,aux]
      }
    odds.Pnx[ind,1] <- 1
    odds.Pnx[ind,2] <- exp(eta0.Pnx[ind,1]*delta.Pnx[ind]+
      log(p.Pnx[5]/p.Pnx[4])*(z.Enf[ind]>1)+
      log(p.Pnx[8]/p.Pnx[7])*Neo_p[ind])      
    odds.Pnx[ind,3] <- exp(eta0.Pnx[ind,2]*delta.Pnx[ind]+
      log(p.Pnx[6]/p.Pnx[4])*(z.Enf[ind]>1)+
      log(p.Pnx[9]/p.Pnx[7])*Neo_p[ind])
    eta0.Pnx[ind,1] <- log(wu.Pnx[ind]*exp(eta.Pnx[1])/(1+(1-wu.Pnx[ind])*(exp(eta.Pnx[1])+exp(eta.Pnx[2]))))
    eta0.Pnx[ind,2] <- log(wu.Pnx[ind]*exp(eta.Pnx[2])/(1+(1-wu.Pnx[ind])*(exp(eta.Pnx[1])+exp(eta.Pnx[2]))))
    p0.Pnx[ind,1] <- 1-wu.Pnx[ind]*(1-p.Pnx[1])
    p0.Pnx[ind,2] <- wu.Pnx[ind]*p.Pnx[2]
    p0.Pnx[ind,3] <- wu.Pnx[ind]*p.Pnx[3]
    wu.Pnx[ind] <- (tilde.Age[ind]<0.24)*tilde.Age[ind]/0.24+(tilde.Age[ind]>=0.24)*(1-tilde.Age[ind])/(1-0.24)
    delta.Pnx[ind] <- 1-(z.Enf[ind]>1)-Neo_p[ind]  

    # Nodule (binary)
    Nod[ind] ~ dbern(theta.Nod[ind])
    theta.Nod[ind] <- (delta.Nod[ind]<1)*thetaStar.Nod[ind]+(delta.Nod[ind]>=1)*p.Nod[1]
    logit(thetaStar.Nod[ind]) <- eta.Nod[1]*delta.Nod[ind]+
      qlogis(p.Nod[2],0,1)*Neo_p[ind]
    delta.Nod[ind] <- 1-Neo_p[ind]

    # Central mass (thoracic) (binary)
    masCen[ind] ~ dbern(theta.masCen[ind])
    theta.masCen[ind] <- (delta.masCen[ind]<1)*thetaStar.masCen[ind]+(delta.masCen[ind]>=1)*p.masCen[1]
    logit(thetaStar.masCen[ind]) <- eta.masCen[1]*delta.masCen[ind]+
      qlogis(p.masCen[2],0,1)*Neo_p[ind]
    delta.masCen[ind] <- 1-Neo_p[ind]
    
    # Leukemic blast brisis (binary)
    blastCris[ind] ~ dbern(theta.blastCris[ind])
    theta.blastCris[ind] <- (delta.blastCris[ind]<1)*thetaStar.blastCris[ind]+(delta.blastCris[ind]>=1)*p.blastCris[1]
    logit(thetaStar.blastCris[ind]) <- eta.blastCris[1]*delta.blastCris[ind]+
      qlogis(p.blastCris[2],0,1)*Leucem[ind]
    delta.blastCris[ind] <- 1-Leucem[ind]

    # Hemorrhage (binary)
    Hemo[ind] ~ dbern(theta.Hemo[ind])
    theta.Hemo[ind] <- (delta.Hemo[ind]<1)*thetaStar.Hemo[ind]+(delta.Hemo[ind]>=1)*p.Hemo[1]
    logit(thetaStar.Hemo[ind]) <- eta.Hemo[1]*delta.Hemo[ind]+
      qlogis(p.Hemo[2],0,1)*Pancr[ind]
    delta.Hemo[ind] <- 1-Pancr[ind]

    # Amylase (binary)
    Amil_o[ind] ~ dbern(theta.Amil_o[ind])
    theta.Amil_o[ind] <- (delta.Amil_o[ind]<1)*thetaStar.Amil_o[ind]+(delta.Amil_o[ind]>=1)*p.Amil_o[1]
    logit(thetaStar.Amil_o[ind]) <- eta.Amil_o[1]*delta.Amil_o[ind]+
      qlogis(p.Amil_o[2],0,1)*Pancr[ind]
    delta.Amil_o[ind] <- 1-Pancr[ind]

    # Peritonitis (binary)
    Perit[ind] ~ dbern(theta.Perit[ind])   
    theta.Perit[ind] <- (delta.Perit[ind]<1)*thetaStar.Perit[ind]+(delta.Perit[ind]>=1)*p.Perit[1]
    logit(thetaStar.Perit[ind]) <- eta.Perit[1]*delta.Perit[ind]+
      qlogis(p.Perit[2],0,1)*Colec[ind]+
      qlogis(p.Perit[3],0,1)*Ulc_pep[ind]    
    delta.Perit[ind] <- 1-Colec[ind]-Ulc_pep[ind]

    # Pregnancy (3 categories, semi-stochastic)
    z.Grav[ind] ~ dcat(theta.Grav[ind,1:3]) 
    theta.Grav[ind,1] <- 1-theta.Grav[ind,2]-theta.Grav[ind,3]
    theta.Grav[ind,2] <- p.Grav[2]*(Estrog[ind]==0&&Fertil[ind]==1)
    theta.Grav[ind,3] <- p.Grav[3]*(Estrog[ind]==0&&Fertil[ind]==1)

    # Chronic aortic valve failure (3 categories with age)
    z.Insuf_aorCron[ind] ~ dcat(theta.Insuf_aorCron[ind,1:3]) 
    for(aux in 1:3) {
      theta.Insuf_aorCron[ind,aux] <- (delta.Insuf_aorCron[ind]<1)*odds.Insuf_aorCron[ind,aux]/sum(odds.Insuf_aorCron[ind,1:3])+(delta.Insuf_aorCron[ind]>=1)*p0.Insuf_aorCron[ind,aux]
      }
    odds.Insuf_aorCron[ind,1] <- 1
    odds.Insuf_aorCron[ind,2] <- exp(eta0.Insuf_aorCron[ind,1]*delta.Insuf_aorCron[ind]+
      log(p.Insuf_aorCron[5]/p.Insuf_aorCron[4])*Card_dil[ind]+
      log(p.Insuf_aorCron[8]/p.Insuf_aorCron[7])*Aneur[ind])      
    odds.Insuf_aorCron[ind,3] <- exp(eta0.Insuf_aorCron[ind,2]*delta.Insuf_aorCron[ind]+
      log(p.Insuf_aorCron[6]/p.Insuf_aorCron[4])*Card_dil[ind]+
      log(p.Insuf_aorCron[9]/p.Insuf_aorCron[7])*Aneur[ind])    
    eta0.Insuf_aorCron[ind,1] <- log(wu.Insuf_aorCron[ind]*exp(eta.Insuf_aorCron[1])/(1+(1-wu.Insuf_aorCron[ind])*(exp(eta.Insuf_aorCron[1])+exp(eta.Insuf_aorCron[2]))))
    eta0.Insuf_aorCron[ind,2] <- log(wu.Insuf_aorCron[ind]*exp(eta.Insuf_aorCron[2])/(1+(1-wu.Insuf_aorCron[ind])*(exp(eta.Insuf_aorCron[1])+exp(eta.Insuf_aorCron[2]))))
    p0.Insuf_aorCron[ind,1] <- 1-wu.Insuf_aorCron[ind]*(1-p.Insuf_aorCron[1])
    p0.Insuf_aorCron[ind,2] <- wu.Insuf_aorCron[ind]*p.Insuf_aorCron[2]
    p0.Insuf_aorCron[ind,3] <- wu.Insuf_aorCron[ind]*p.Insuf_aorCron[3]
    wu.Insuf_aorCron[ind] <- (tilde.Age[ind]<0.76)*tilde.Age[ind]/0.76+(tilde.Age[ind]>=0.76)*(1-tilde.Age[ind])/(1-0.76)
    delta.Insuf_aorCron[ind] <- 1-Card_dil[ind]-Aneur[ind]

    # Dilated ascending aorta (binary, non-stochastic)
    Dil_aort[ind] ~ dbern(theta.Dil_aort[ind])
    theta.Dil_aort[ind] <- 1*(Aneur[ind]==1)

    # Aortic dissection (binary with age and gender)
    Aort_diss[ind] ~ dbern(theta.Aort_diss[ind])  
    theta.Aort_diss[ind] <- (delta.Aort_diss[ind]<1)*thetaStar.Aort_diss[ind]+(delta.Aort_diss[ind]>=1)*wu.Aort_diss[ind]*p0.Aort_diss[ind]
    logit(thetaStar.Aort_diss[ind]) <- eta0.Aort_diss[ind]*delta.Aort_diss[ind]+
      qlogis(p.Aort_diss[3],0,1)*Aneur[ind]
    p0.Aort_diss[ind] <- (Sex[ind]==0)*p.Aort_diss[1]+(Sex[ind]==1)*p.Aort_diss[2]
    eta0.Aort_diss[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.Aort_diss[1],0,1)*wu.Aort_diss[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.Aort_diss[2],0,1)*wu.Aort_diss[ind],0,1)
    wu.Aort_diss[ind] <- (tilde.Age[ind]<0.52)*tilde.Age[ind]/0.52+(tilde.Age[ind]>=0.52)*(1-tilde.Age[ind])/(1-0.52)
    delta.Aort_diss[ind] <- 1-Aneur[ind]
    
    # Chest pain (gastro-oesophageal origin) (binary)
    Spas_es[ind] ~ dbern(theta.Spas_es[ind])
    theta.Spas_es[ind] <- (delta.Spas_es[ind]<1)*thetaStar.Spas_es[ind]+(delta.Spas_es[ind]>=1)*p.Spas_es[1]
    logit(thetaStar.Spas_es[ind]) <- eta.Spas_es[1]*delta.Spas_es[ind]+
      qlogis(p.Spas_es[2],0,1)*Refl_ge[ind]+
      qlogis(p.Spas_es[3],0,1)*Ern_ia[ind]+
      qlogis(p.Spas_es[4],0,1)*sindrMW[ind]
    delta.Spas_es[ind] <- 1-Refl_ge[ind]-Ern_ia[ind]-sindrMW[ind]

    # Pulmonary hypertension (binary)
    Ipertpolm[ind] ~ dbern(theta.Ipertpolm[ind])
    theta.Ipertpolm[ind] <- (delta.Ipertpolm[ind]<1)*thetaStar.Ipertpolm[ind]+(delta.Ipertpolm[ind]>=1)*p.Ipertpolm[1]
    logit(thetaStar.Ipertpolm[ind]) <- eta.Ipertpolm[1]*delta.Ipertpolm[ind]+
      qlogis(p.Ipertpolm[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.Ipertpolm[3],0,1)*(z.Enf[ind]==3)+
      qlogis(p.Ipertpolm[4],0,1)*preTVP[ind]  
    delta.Ipertpolm[ind] <- 1-(z.Enf[ind]>1)-preTVP[ind]
        
    # Thyroid-stimulating hormone (continuous)
    tilde.TSH[ind] ~ dbeta(relmu.TSH[ind]*p.TSH[6],(1-relmu.TSH[ind])*p.TSH[6]) T(tol,1-tol)
    logit(relmu.TSH[ind]) <- qlogis(p.TSH[1],0,1)*delta.TSH[ind]+
      qlogis(p.TSH[2],0,1)*(z.Tiroid[ind]==2)+
      qlogis(p.TSH[3],0,1)*(z.Tiroid[ind]==3)+
      qlogis(p.TSH[4],0,1)*(z.Tiroid[ind]==4)+
      qlogis(p.TSH[5],0,1)*(z.Tiroid[ind]==5)
    delta.TSH[ind] <- 1-(z.Tiroid[ind]>1)
        
    # Thyroid hormones (continuous)
    tilde.OT[ind] ~ dbeta(relmu.OT[ind]*p.OT[6],(1-relmu.OT[ind])*p.OT[6]) T(tol,1-tol)
    logit(relmu.OT[ind]) <- qlogis(p.OT[1],0,1)*delta.OT[ind]+
      qlogis(p.OT[2],0,1)*(z.Tiroid[ind]==2)+
      qlogis(p.OT[3],0,1)*(z.Tiroid[ind]==3)+
      qlogis(p.OT[4],0,1)*(z.Tiroid[ind]==4)+
      qlogis(p.OT[5],0,1)*(z.Tiroid[ind]==5)
    delta.OT[ind] <- 1-(z.Tiroid[ind]>1)

    # Deep vein thrombosis risk (binary)
    TVPscore[ind] ~ dbern(theta.TVPscore[ind])
    theta.TVPscore[ind] <- (delta.TVPscore[ind]<1)*thetaStar.TVPscore[ind]+(delta.TVPscore[ind]>=1)*p.TVPscore[1]
    logit(thetaStar.TVPscore[ind]) <- eta.TVPscore[1]*delta.TVPscore[ind]+
      qlogis(p.TVPscore[2],0,1)*Fr[ind]+
      qlogis(p.TVPscore[3],0,1)*Fum[ind]+
      qlogis(p.TVPscore[4],0,1)*InsVenCr[ind]+
      qlogis(p.TVPscore[5],0,1)*Neo_o[ind]
    delta.TVPscore[ind] <- 1-Fr[ind]-Fum[ind]-InsVenCr[ind]-Neo_o[ind]

    # Lung cancer (binary)
    Neo_p[ind] ~ dbern(theta.Neo_p[ind])
    theta.Neo_p[ind] <- (delta.Neo_p[ind]<1)*thetaStar.Neo_p[ind]+(delta.Neo_p[ind]>=1)*p.Neo_p[1]
    logit(thetaStar.Neo_p[ind]) <- eta.Neo_p[1]*delta.Neo_p[ind]+
      qlogis(p.Neo_p[2],0,1)*Neo_o[ind]
    delta.Neo_p[ind] <- 1-Neo_o[ind]

    # Leukemia (binary)
    Leucem[ind] ~ dbern(theta.Leucem[ind])
    theta.Leucem[ind] <- (delta.Leucem[ind]<1)*thetaStar.Leucem[ind]+(delta.Leucem[ind]>=1)*p.Leucem[1]
    logit(thetaStar.Leucem[ind]) <- eta.Leucem[1]*delta.Leucem[ind]+
      qlogis(p.Leucem[2],0,1)*Neo_o[ind]
    delta.Leucem[ind] <- 1-Neo_o[ind]

    # Rib fracture (binary)
    Fr_cost[ind] ~ dbern(theta.Fr_cost[ind])
    theta.Fr_cost[ind] <- (delta.Fr_cost[ind]<1)*thetaStar.Fr_cost[ind]+(delta.Fr_cost[ind]>=1)*p.Fr_cost[1]
    logit(thetaStar.Fr_cost[ind]) <- eta.Fr_cost[1]*delta.Fr_cost[ind]+
      qlogis(p.Fr_cost[2],0,1)*Neo_o[ind]
    delta.Fr_cost[ind] <- 1-Neo_o[ind]

    # Immunocompromission (binary with age)
    Defed[ind] ~ dbern(theta.Defed[ind]) 
    theta.Defed[ind] <- (delta.Defed[ind]<1)*thetaStar.Defed[ind]+(delta.Defed[ind]>=1)*wu.Defed[ind]*p.Defed[1]
    logit(thetaStar.Defed[ind]) <- qlogis(plogis(eta.Defed[1],0,1)*wu.Defed[ind],0,1)*delta.Defed[ind]+
      qlogis(p.Defed[2],0,1)*Neo_o[ind]
    wu.Defed[ind] <- (tilde.Age[ind]<0.95)*tilde.Age[ind]/0.95+(tilde.Age[ind]>=0.95)*(1-tilde.Age[ind])/(1-0.95)
    delta.Defed[ind] <- 1-Neo_o[ind]

    # Pancreatitis (binary with age and gender)
    Pancr[ind] ~ dbern(theta.Pancr[ind]) 
    theta.Pancr[ind] <- (delta.Pancr[ind]<1)*thetaStar.Pancr[ind]+(delta.Pancr[ind]>=1)*wu.Pancr[ind]*p0.Pancr[ind]
    logit(thetaStar.Pancr[ind]) <- eta0.Pancr[ind]*delta.Pancr[ind]+
      qlogis(p.Pancr[3],0,1)*Litias[ind]+
      qlogis(p.Pancr[4],0,1)*Alcol[ind]
    p0.Pancr[ind] <- (Sex[ind]==0)*p.Pancr[1]+(Sex[ind]==1)*p.Pancr[2]
    eta0.Pancr[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.Pancr[1],0,1)*wu.Pancr[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.Pancr[2],0,1)*wu.Pancr[ind],0,1)
    wu.Pancr[ind] <- (tilde.Age[ind]<0.52)*tilde.Age[ind]/0.52+(tilde.Age[ind]>=0.52)*(1-tilde.Age[ind])/(1-0.52)
    delta.Pancr[ind] <- 1-Litias[ind]-Alcol[ind]  

    # Biliary colic (binary)
    Colic[ind] ~ dbern(theta.Colic[ind])
    theta.Colic[ind] <- (delta.Colic[ind]<1)*thetaStar.Colic[ind]+(delta.Colic[ind]>=1)*p.Colic[1]
    logit(thetaStar.Colic[ind]) <- eta.Colic[1]*delta.Colic[ind]+
      qlogis(p.Colic[2],0,1)*Litias[ind]
    delta.Colic[ind] <- 1-Litias[ind]

    # Cholecystitis (binary)
    Colec[ind] ~ dbern(theta.Colec[ind]) 
    theta.Colec[ind] <- (delta.Colec[ind]<1)*thetaStar.Colec[ind]+(delta.Colec[ind]>=1)*p.Colec[1]
    logit(thetaStar.Colec[ind]) <- eta.Colec[1]*delta.Colec[ind]+
      qlogis(p.Colec[2],0,1)*Litias[ind]   
    delta.Colec[ind] <- 1-Litias[ind]
    
    # Extrogens use (binary)
    Estrog[ind] ~ dbern(theta.Estrog[ind])
    theta.Estrog[ind] <- (delta.Estrog[ind]<1)*thetaStar.Estrog[ind]+(delta.Estrog[ind]>=1)*p.Estrog[1]
    logit(thetaStar.Estrog[ind]) <- eta.Estrog[1]*delta.Estrog[ind]+
      qlogis(p.Estrog[2],0,1)*Fertil[ind]
    delta.Estrog[ind] <- 1-Fertil[ind]

    # Chronic anemia (binary with age)
    cronAnem[ind] ~ dbern(theta.cronAnem[ind]) 
    theta.cronAnem[ind] <- (delta.cronAnem[ind]<1)*thetaStar.cronAnem[ind]+(delta.cronAnem[ind]>=1)*wu.cronAnem[ind]*p.cronAnem[1]    
    logit(thetaStar.cronAnem[ind]) <- qlogis(plogis(eta.cronAnem[1],0,1)*wu.cronAnem[ind],0,1)*delta.cronAnem[ind]+
      qlogis(p.cronAnem[2],0,1)*Fertil[ind]  
    wu.cronAnem[ind] <- (tilde.Age[ind]<0.71)*tilde.Age[ind]/0.71+(tilde.Age[ind]>=0.71)*(1-tilde.Age[ind])/(1-0.71)
    delta.cronAnem[ind] <- 1-Fertil[ind]

    # Bronchial walls (binary)
    parBron[ind] ~ dbern(theta.parBron[ind])
    theta.parBron[ind] <- (delta.parBron[ind]<1)*thetaStar.parBron[ind]+(delta.parBron[ind]>=1)*p.parBron[1]
    logit(thetaStar.parBron[ind]) <- eta.parBron[1]*delta.parBron[ind]+
      qlogis(p.parBron[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.parBron[3],0,1)*(z.Enf[ind]==3)
    delta.parBron[ind] <- 1-(z.Enf[ind]>1)

    # Cystic areas/bullae (binary)
    CistBoll[ind] ~ dbern(theta.CistBoll[ind])
    theta.CistBoll[ind] <- (delta.CistBoll[ind]<1)*thetaStar.CistBoll[ind]+(delta.CistBoll[ind]>=1)*p.CistBoll[1]
    logit(thetaStar.CistBoll[ind]) <- eta.CistBoll[1]*delta.CistBoll[ind]+
      qlogis(p.CistBoll[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.CistBoll[3],0,1)*(z.Enf[ind]==3)
    delta.CistBoll[ind] <- 1-(z.Enf[ind]>1)
 
    # Bronchial diameter (binary)
    calBron[ind] ~ dbern(theta.calBron[ind])
    theta.calBron[ind] <- (delta.calBron[ind]<1)*thetaStar.calBron[ind]+(delta.calBron[ind]>=1)*p.calBron[1]
    logit(thetaStar.calBron[ind]) <- eta.calBron[1]*delta.calBron[ind]+
      qlogis(p.calBron[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.calBron[3],0,1)*(z.Enf[ind]==3)
    delta.calBron[ind] <- 1-(z.Enf[ind]>1)

    # Bronchiectasis (binary)
    bronect[ind] ~ dbern(theta.bronect[ind])
    theta.bronect[ind] <- (delta.bronect[ind]<1)*thetaStar.bronect[ind]+(delta.bronect[ind]>=1)*p.bronect[1]
    logit(thetaStar.bronect[ind]) <- eta.bronect[1]*delta.bronect[ind]+
      qlogis(p.bronect[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.bronect[3],0,1)*(z.Enf[ind]==3)
    delta.bronect[ind] <- 1-(z.Enf[ind]>1) 

    # Chronic obstructive pulmonary disease (3 categories)
    z.BPCO[ind] ~ dcat(theta.BPCO[ind,1:3])
    for(aux in 1:3) {
      theta.BPCO[ind,aux] <- (delta.BPCO[ind]<1)*odds.BPCO[ind,aux]/sum(odds.BPCO[ind,1:3])+(delta.BPCO[ind]>=1)*p.BPCO[aux]
      }
    odds.BPCO[ind,1] <- 1
    odds.BPCO[ind,2] <- exp(eta.BPCO[1]*delta.BPCO[ind]+
      log(p.BPCO[5]/p.BPCO[4])*(z.Enf[ind]==2)+
      log(p.BPCO[8]/p.BPCO[7])*(z.Enf[ind]==3))
    odds.BPCO[ind,3] <- exp(eta.BPCO[2]*delta.BPCO[ind]+
      log(p.BPCO[6]/p.BPCO[4])*(z.Enf[ind]==2)+
      log(p.BPCO[9]/p.BPCO[7])*(z.Enf[ind]==3))
    delta.BPCO[ind] <- 1-(z.Enf[ind]>1)

    # Air trapping (binary)
    airTrap[ind] ~ dbern(theta.airTrap[ind])
    theta.airTrap[ind] <- (delta.airTrap[ind]<1)*thetaStar.airTrap[ind]+(delta.airTrap[ind]>=1)*p.airTrap[1]
    logit(thetaStar.airTrap[ind]) <- eta.airTrap[1]*delta.airTrap[ind]+
      qlogis(p.airTrap[2],0,1)*(z.Enf[ind]==2)+
      qlogis(p.airTrap[3],0,1)*(z.Enf[ind]==3)
    delta.airTrap[ind] <- 1-(z.Enf[ind]>1)

    # Chronic cerebro-vascular disease (binary with age)
    IschCerCron[ind] ~ dbern(theta.IschCerCron[ind])  
    theta.IschCerCron[ind] <- (delta.IschCerCron[ind]<1)*thetaStar.IschCerCron[ind]+(delta.IschCerCron[ind]>=1)*wu.IschCerCron[ind]*p.IschCerCron[1]    
    logit(thetaStar.IschCerCron[ind]) <- qlogis(plogis(eta.IschCerCron[1],0,1)*wu.IschCerCron[ind],0,1)*delta.IschCerCron[ind]+
      qlogis(p.IschCerCron[2],0,1)*cronIpert[ind]
    wu.IschCerCron[ind] <- (tilde.Age[ind]<0.99)*tilde.Age[ind]/0.99+(tilde.Age[ind]>=0.99)*(1-tilde.Age[ind])/(1-0.99)
    delta.IschCerCron[ind] <- 1-cronIpert[ind]

    # Aortic aneurysm (binary with age and gender)
    Aneur[ind] ~ dbern(theta.Aneur[ind]) 
    theta.Aneur[ind] <- (delta.Aneur[ind]<1)*thetaStar.Aneur[ind]+(delta.Aneur[ind]>=1)*wu.Aneur[ind]*p0.Aneur[ind]
    logit(thetaStar.Aneur[ind]) <- eta0.Aneur[ind]*delta.Aneur[ind]+
      qlogis(p.Aneur[3],0,1)*cronIpert[ind]
    p0.Aneur[ind] <- (Sex[ind]==0)*p.Aneur[1]+(Sex[ind]==1)*p.Aneur[2]
    eta0.Aneur[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.Aneur[1],0,1)*wu.Aneur[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.Aneur[2],0,1)*wu.Aneur[ind],0,1)
    wu.Aneur[ind] <- (tilde.Age[ind]<0.76)*tilde.Age[ind]/0.76+(tilde.Age[ind]>=0.76)*(1-tilde.Age[ind])/(1-0.76)
    delta.Aneur[ind] <- 1-cronIpert[ind]

    # Compression stockings (binary)
    Calze[ind] ~ dbern(theta.Calze[ind])    
    theta.Calze[ind] <- (delta.Calze[ind]<1)*thetaStar.Calze[ind]+(delta.Calze[ind]>=1)*p.Calze[1]
    logit(thetaStar.Calze[ind]) <- eta.Calze[1]*delta.Calze[ind]+
      qlogis(p.Calze[2],0,1)*(z.Chir[ind]>1)
    delta.Calze[ind] <- 1-(z.Chir[ind]>1)
    
    # Carotid sinus massage test (binary)
    MSC[ind] ~ dbern(theta.MSC[ind])
    theta.MSC[ind] <- (delta.MSC[ind]<1)*thetaStar.MSC[ind]+(delta.MSC[ind]>=1)*p.MSC[1]
    logit(thetaStar.MSC[ind]) <- eta.MSC[1]*delta.MSC[ind]+
      qlogis(p.MSC[2],0,1)*MNodSA[ind]
    delta.MSC[ind] <- 1-MNodSA[ind]

    # Urinary catecholamines (binary)
    catecUr[ind] ~ dbern(theta.catecUr[ind])
    theta.catecUr[ind] <- (delta.catecUr[ind]<1)*thetaStar.catecUr[ind]+(delta.catecUr[ind]>=1)*p.catecUr[1]
    logit(thetaStar.catecUr[ind]) <- eta.catecUr[1]*delta.catecUr[ind]+
      qlogis(p.catecUr[2],0,1)*Feocrom[ind]
    delta.catecUr[ind] <- 1-Feocrom[ind]

    # Mallory-Weiss syndrome (binary)
    sindrMW[ind] ~ dbern(theta.sindrMW[ind])
    theta.sindrMW[ind] <- (delta.sindrMW[ind]<1)*thetaStar.sindrMW[ind]+(delta.sindrMW[ind]>=1)*p.sindrMW[1]
    logit(thetaStar.sindrMW[ind]) <- eta.sindrMW[1]*delta.sindrMW[ind]+
      qlogis(p.sindrMW[2],0,1)*Alcol[ind]
    delta.sindrMW[ind] <- 1-Alcol[ind]

    # Protein S (binary)
    protS[ind] ~ dbern(theta.protS[ind])
    theta.protS[ind] <- (delta.protS[ind]<1)*thetaStar.protS[ind]+(delta.protS[ind]>=1)*p.protS[1]
    logit(thetaStar.protS[ind]) <- eta.protS[1]*delta.protS[ind]+
      qlogis(p.protS[2],0,1)*Trmphil[ind]
    delta.protS[ind] <- 1-Trmphil[ind]

    # Protein C (binary)
    protC[ind] ~ dbern(theta.protC[ind])
    theta.protC[ind] <- (delta.protC[ind]<1)*thetaStar.protC[ind]+(delta.protC[ind]>=1)*p.protC[1]
    logit(thetaStar.protC[ind]) <- eta.protC[1]*delta.protC[ind]+
      qlogis(p.protC[2],0,1)*Trmphil[ind]
    delta.protC[ind] <- 1-Trmphil[ind]
    
    # Previous episode of deep venous thrombosis/pulmonary embolism (binary)
    preTVP[ind] ~ dbern(theta.preTVP[ind])
    theta.preTVP[ind] <- (delta.preTVP[ind]<1)*thetaStar.preTVP[ind]+(delta.preTVP[ind]>=1)*p.preTVP[1]
    logit(thetaStar.preTVP[ind]) <- eta.preTVP[1]*delta.preTVP[ind]+
      qlogis(p.preTVP[2],0,1)*Trmphil[ind]
    delta.preTVP[ind] <- 1-Trmphil[ind]
    
    # Hyperhomocysteinemia (binary)
    omocist[ind] ~ dbern(theta.omocist[ind])
    theta.omocist[ind] <- (delta.omocist[ind]<1)*thetaStar.omocist[ind]+(delta.omocist[ind]>=1)*p.omocist[1]
    logit(thetaStar.omocist[ind]) <- eta.omocist[1]*delta.omocist[ind]+
      qlogis(p.omocist[2],0,1)*Trmphil[ind]
    delta.omocist[ind] <- 1-Trmphil[ind]
    
    # Leiden factor V (binary)
    Leid[ind] ~ dbern(theta.Leid[ind])
    theta.Leid[ind] <- (delta.Leid[ind]<1)*thetaStar.Leid[ind]+(delta.Leid[ind]>=1)*p.Leid[1]
    logit(thetaStar.Leid[ind]) <- eta.Leid[1]*delta.Leid[ind]+
      qlogis(p.Leid[2],0,1)*Trmphil[ind]
    delta.Leid[ind] <- 1-Trmphil[ind]
    
    # Factor II G20210A (binary)
    fIIG[ind] ~ dbern(theta.fIIG[ind])
    theta.fIIG[ind] <- (delta.fIIG[ind]<1)*thetaStar.fIIG[ind]+(delta.fIIG[ind]>=1)*p.fIIG[1]
    logit(thetaStar.fIIG[ind]) <- eta.fIIG[1]*delta.fIIG[ind]+
      qlogis(p.fIIG[2],0,1)*Trmphil[ind]
    delta.fIIG[ind] <- 1-Trmphil[ind]

    # Factor VIII (binary)
    FattVIII[ind] ~ dbern(theta.FattVIII[ind])
    theta.FattVIII[ind] <- (delta.FattVIII[ind]<1)*thetaStar.FattVIII[ind]+(delta.FattVIII[ind]>=1)*p.FattVIII[1]
    logit(thetaStar.FattVIII[ind]) <- eta.FattVIII[1]*delta.FattVIII[ind]+
      qlogis(p.FattVIII[2],0,1)*Trmphil[ind]
    delta.FattVIII[ind] <- 1-Trmphil[ind]

    # Antithrombin III (binary)
    antiTromb[ind] ~ dbern(theta.antiTromb[ind])
    theta.antiTromb[ind] <- (delta.antiTromb[ind]<1)*thetaStar.antiTromb[ind]+(delta.antiTromb[ind]>=1)*p.antiTromb[1]
    logit(thetaStar.antiTromb[ind]) <- eta.antiTromb[1]*delta.antiTromb[ind]+
      qlogis(p.antiTromb[2],0,1)*Trmphil[ind]
    delta.antiTromb[ind] <- 1-Trmphil[ind]

    # Antiphospholipids (binary)
    abPhos[ind] ~ dbern(theta.abPhos[ind])
    theta.abPhos[ind] <- (delta.abPhos[ind]<1)*thetaStar.abPhos[ind]+(delta.abPhos[ind]>=1)*p.abPhos[1]
    logit(thetaStar.abPhos[ind]) <- eta.abPhos[1]*delta.abPhos[ind]+
      qlogis(p.abPhos[2],0,1)*Trmphil[ind]
    delta.abPhos[ind] <- 1-Trmphil[ind]
    
    # Thyroid disease (5 categories with age and gender)
    z.Tiroid[ind] ~ dcat(theta.Tiroid[ind,1:5]) 
    theta.Tiroid[ind,1] <- 1-wu.Tiroid[ind]*(1-(Sex[ind]==0)*p.Tiroid[1]-(Sex[ind]==1)*p.Tiroid[6])
    theta.Tiroid[ind,2] <- ((Sex[ind]==0)*p.Tiroid[2]+(Sex[ind]==1)*p.Tiroid[7])*wu.Tiroid[ind]
    theta.Tiroid[ind,3] <- ((Sex[ind]==0)*p.Tiroid[3]+(Sex[ind]==1)*p.Tiroid[8])*wu.Tiroid[ind]
    theta.Tiroid[ind,4] <- ((Sex[ind]==0)*p.Tiroid[4]+(Sex[ind]==1)*p.Tiroid[9])*wu.Tiroid[ind]
    theta.Tiroid[ind,5] <- ((Sex[ind]==0)*p.Tiroid[5]+(Sex[ind]==1)*p.Tiroid[10])*wu.Tiroid[ind]
    wu.Tiroid[ind] <- (tilde.Age[ind]<0.43)*tilde.Age[ind]/0.43+(tilde.Age[ind]>=0.43)*(1-tilde.Age[ind])/(1-0.43)

    # Neoplastic disease (generic) (binary with age and gender)
    Neo_o[ind] ~ dbern(theta.Neo_o[ind]) 
    theta.Neo_o[ind] <- (delta.Neo_o[ind]<1)*thetaStar.Neo_o[ind]+(delta.Neo_o[ind]>=1)*wu.Neo_o[ind]*p0.Neo_o[ind]
    logit(thetaStar.Neo_o[ind]) <- eta0.Neo_o[ind]*delta.Neo_o[ind]+
      qlogis(p.Neo_o[3],0,1)*Fum[ind]
    p0.Neo_o[ind] <- (Sex[ind]==0)*p.Neo_o[1]+(Sex[ind]==1)*p.Neo_o[2]
    eta0.Neo_o[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.Neo_o[1],0,1)*wu.Neo_o[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.Neo_o[2],0,1)*wu.Neo_o[ind],0,1)
    wu.Neo_o[ind] <- (tilde.Age[ind]<0.67)*tilde.Age[ind]/0.67+(tilde.Age[ind]>=0.67)*(1-tilde.Age[ind])/(1-0.67)
    delta.Neo_o[ind] <- 1-Fum[ind]

    # Cholelithiasis (binary with age and gender)
    Litias[ind] ~ dbern(theta.Litias[ind])  
    theta.Litias[ind] <- (delta.Litias[ind]<1)*thetaStar.Litias[ind]+(delta.Litias[ind]>=1)*wu.Litias[ind]*p0.Litias[ind]
    logit(thetaStar.Litias[ind]) <- eta0.Litias[ind]*delta.Litias[ind]+
      qlogis(p.Litias[3],0,1)*Obes[ind]
    p0.Litias[ind] <- (Sex[ind]==0)*p.Litias[1]+(Sex[ind]==1)*p.Litias[2]
    eta0.Litias[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.Litias[1],0,1)*wu.Litias[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.Litias[2],0,1)*wu.Litias[ind],0,1)
    wu.Litias[ind] <- (tilde.Age[ind]<0.38)*tilde.Age[ind]/0.38+(tilde.Age[ind]>=0.38)*(1-tilde.Age[ind])/(1-0.38)
    delta.Litias[ind] <- 1-Obes[ind]

    # Chronic venous insufficiency (binary with age and gender)
    InsVenCr[ind] ~ dbern(theta.InsVenCr[ind])
    theta.InsVenCr[ind] <- ((Sex[ind]==0)*p.InsVenCr[1]+(Sex[ind]==1)*p.InsVenCr[2])*wu.InsVenCr[ind]
    wu.InsVenCr[ind] <- (tilde.Age[ind]<0.57)*tilde.Age[ind]/0.57+(tilde.Age[ind]>=0.57)*(1-tilde.Age[ind])/(1-0.57)
        
    # Fertility (binary, semi-stochastic)
    #   [absent if male or rescaled age >0.48]
    Fertil[ind] ~ dbern(theta.Fertil[ind])
    theta.Fertil[ind] <- p.Fertil[1]*(Sex[ind]==1&&tilde.Age[ind]<0.48)*(1-1/0.48*tilde.Age[ind])
    
    # Pulmonary emphysema (3 categories with age and gender)
    z.Enf[ind] ~ dcat(theta.Enf[ind,1:3]) 
    for(aux in 1:3) {
      theta.Enf[ind,aux] <- (delta.Enf[ind]<1)*odds.Enf[ind,aux]/sum(odds.Enf[ind,1:3])+(delta.Enf[ind]>=1)*p0.Enf[ind,aux]
      }
    odds.Enf[ind,1] <- 1
    odds.Enf[ind,2] <- exp(eta0.Enf[ind,1]*delta.Enf[ind]+
      log(p.Enf[8]/p.Enf[7])*Fum[ind])      
    odds.Enf[ind,3] <- exp(eta0.Enf[ind,2]*delta.Enf[ind]+    
      log(p.Enf[9]/p.Enf[7])*Fum[ind])    
    eta0.Enf[ind,1] <- (Sex[ind]==0)*log(wu.Enf[ind]*exp(eta.Enf[1])/(1+(1-wu.Enf[ind])*(exp(eta.Enf[1])+exp(eta.Enf[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Enf[ind]*exp(eta.Enf[2])/(1+(1-wu.Enf[ind])*(exp(eta.Enf[2])+exp(eta.Enf[4]))))
    eta0.Enf[ind,2] <- (Sex[ind]==0)*log(wu.Enf[ind]*exp(eta.Enf[3])/(1+(1-wu.Enf[ind])*(exp(eta.Enf[1])+exp(eta.Enf[3]))))+
      (Sex[ind]==1)*(Sex[ind]==0)*log(wu.Enf[ind]*exp(eta.Enf[4])/(1+(1-wu.Enf[ind])*(exp(eta.Enf[2])+exp(eta.Enf[4]))))
    p0.Enf[ind,1] <- 1-wu.Enf[ind]*(1-(Sex[ind]==0)*p.Enf[1]-(Sex[ind]==1)*p.Enf[4])
    p0.Enf[ind,2] <- ((Sex[ind]==0)*p.Enf[2]+(Sex[ind]==1)*p.Enf[5])*wu.Enf[ind]
    p0.Enf[ind,3] <- ((Sex[ind]==0)*p.Enf[3]+(Sex[ind]==1)*p.Enf[6])*wu.Enf[ind]
    wu.Enf[ind] <- (tilde.Age[ind]<0.9)*tilde.Age[ind]/0.9+(tilde.Age[ind]>=0.9)*(1-tilde.Age[ind])/(1-0.9)
    delta.Enf[ind] <- 1-Fum[ind]
  
    # Chronic arterial hypertension (binary with age and gender)
    cronIpert[ind] ~ dbern(theta.cronIpert[ind]) 
    theta.cronIpert[ind] <- (delta.cronIpert[ind]<1)*thetaStar.cronIpert[ind]+(delta.cronIpert[ind]>=1)*wu.cronIpert[ind]*p0.cronIpert[ind]
    logit(thetaStar.cronIpert[ind]) <- eta0.cronIpert[ind]*delta.cronIpert[ind]+
      qlogis(p.cronIpert[3],0,1)*Fum[ind]+
      qlogis(p.cronIpert[4],0,1)*Obes[ind]
    p0.cronIpert[ind] <- (Sex[ind]==0)*p.cronIpert[1]+(Sex[ind]==1)*p.cronIpert[2]
    eta0.cronIpert[ind] <- (Sex[ind]==0)*qlogis(plogis(eta.cronIpert[1],0,1)*wu.cronIpert[ind],0,1)+(Sex[ind]==1)*qlogis(plogis(eta.cronIpert[2],0,1)*wu.cronIpert[ind],0,1)
    wu.cronIpert[ind] <- (tilde.Age[ind]<0.9)*tilde.Age[ind]/0.9+(tilde.Age[ind]>=0.9)*(1-tilde.Age[ind])/(1-0.9)
    delta.cronIpert[ind] <- 1-Fum[ind]-Obes[ind]
    
    # Dilated cardiomyopathy (binary with age and gender)
    Card_dil[ind] ~ dbern(theta.Card_dil[ind])
    theta.Card_dil[ind] <- ((Sex[ind]==0)*p.Card_dil[1]+(Sex[ind]==1)*p.Card_dil[2])*wu.Card_dil[ind]
    wu.Card_dil[ind] <- (tilde.Age[ind]<0.33)*tilde.Age[ind]/0.33+(tilde.Age[ind]>=0.33)*(1-tilde.Age[ind])/(1-0.33)
    
    # Ashtma (binary with age and gender)
    Asma[ind] ~ dbern(theta.Asma[ind])
    theta.Asma[ind] <- ((Sex[ind]==0)*p.Asma[1]+(Sex[ind]==1)*p.Asma[2])*wu.Asma[ind]
    wu.Asma[ind] <- (tilde.Age[ind]<0.19)*tilde.Age[ind]/0.19+(tilde.Age[ind]>=0.19)*(1-tilde.Age[ind])/(1-0.19)

    # Immobilisation (binary)
    Allett[ind] ~ dbern(theta.Allett[ind])
    theta.Allett[ind] <- (delta.Allett[ind]<1)*thetaStar.Allett[ind]+(delta.Allett[ind]>=1)*p.Allett[1]
    logit(thetaStar.Allett[ind]) <- eta.Allett[1]*delta.Allett[ind]+
      qlogis(p.Allett[2],0,1)*Neuromusc[ind]
    delta.Allett[ind] <- 1-Neuromusc[ind]

    # Surgery (3 categories)
    z.Chir[ind] ~ dcat(theta.Chir[ind,1:3]) 
    for(aux in 1:3) {
      theta.Chir[ind,aux] <- (delta.Chir[ind]<1)*odds.Chir[ind,aux]/sum(odds.Chir[ind,1:3])+(delta.Chir[ind]>=1)*p.Chir[aux]
      }
    odds.Chir[ind,1] <- 1
    odds.Chir[ind,2] <- exp(eta.Chir[1]*delta.Chir[ind]+
      log(p.Chir[5]/p.Chir[4])*Fr[ind])
    odds.Chir[ind,3] <- exp(eta.Chir[2]*delta.Chir[ind]+
      log(p.Chir[6]/p.Chir[4])*Fr[ind])    
    delta.Chir[ind] <- 1-Fr[ind]

    # Gastro-oesophageal reflux (binary with age)
    Refl_ge[ind] ~ dbern(theta.Refl_ge[ind])
    theta.Refl_ge[ind] <- (delta.Refl_ge[ind]<1)*thetaStar.Refl_ge[ind]+(delta.Refl_ge[ind]>=1)*wu.Refl_ge[ind]*p.Refl_ge[1]
    logit(thetaStar.Refl_ge[ind]) <- qlogis(plogis(eta.Refl_ge[1],0,1)*wu.Refl_ge[ind],0,1)*delta.Refl_ge[ind]+
      qlogis(p.Refl_ge[2],0,1)*Ern_ia[ind]    
    wu.Refl_ge[ind] <- (tilde.Age[ind]<0.62)*tilde.Age[ind]/0.62+(tilde.Age[ind]>=0.62)*(1-tilde.Age[ind])/(1-0.62)
    delta.Refl_ge[ind] <- 1-Ern_ia[ind]
    
    # Endocardial vegetations (binary, semi-stochastic)
    Veget_end[ind] ~ dbern(theta.Veget_end[ind])
    theta.Veget_end[ind] <- p.Veget_end[1]*(Endocard[ind]==0)+1*(Endocard[ind]==1)
    
    # Sick sinus syndrome (binary with age)       
    MNodSA[ind] ~ dbern(theta.MNodSA[ind])
    theta.MNodSA[ind] <- p.MNodSA[1]*wu.MNodSA[ind]
    wu.MNodSA[ind] <- (tilde.Age[ind]<0.84)*tilde.Age[ind]/0.84+(tilde.Age[ind]>=0.84)*(1-tilde.Age[ind])/(1-0.84)
        
    # Chronic metabolic alkalosis (binary with age)
    metabAlk[ind] ~ dbern(theta.metabAlk[ind])
    theta.metabAlk[ind] <- p.metabAlk[1]*wu.metabAlk[ind]
    wu.metabAlk[ind] <- (tilde.Age[ind]<0.81)*tilde.Age[ind]/0.81+(tilde.Age[ind]>=0.81)*(1-tilde.Age[ind])/(1-0.81)

    # Chronic interstitial lung disease (binary with age)
    Interst[ind] ~ dbern(theta.Interst[ind])
    theta.Interst[ind] <- p.Interst[1]*wu.Interst[ind]
    wu.Interst[ind] <- (tilde.Age[ind]<0.67)*tilde.Age[ind]/0.67+(tilde.Age[ind]>=0.67)*(1-tilde.Age[ind])/(1-0.67)
    
    # Pheochromocytoma (binary with age)
    Feocrom[ind] ~ dbern(theta.Feocrom[ind])
    theta.Feocrom[ind] <- p.Feocrom[1]*wu.Feocrom[ind]
    wu.Feocrom[ind] <- (tilde.Age[ind]<0.43)*tilde.Age[ind]/0.43+(tilde.Age[ind]>=0.43)*(1-tilde.Age[ind])/(1-0.43)

    # Alcoholism (binary with age)
    Alcol[ind] ~ dbern(theta.Alcol[ind])
    theta.Alcol[ind] <- p.Alcol[1]*wu.Alcol[ind]
    wu.Alcol[ind] <- (tilde.Age[ind]<0.19)*tilde.Age[ind]/0.19+(tilde.Age[ind]>=0.19)*(1-tilde.Age[ind])/(1-0.19)

    # Ventricular pre-excitation (binary)
    WPW[ind] ~ dbern(p.WPW[1])
    
    # Upper airways infection (binary)
    URI[ind] ~ dbern(p.URI[1])

    # Peptic ulcer (binary)
    Ulc_pep[ind] ~ dbern(p.Ulc_pep[1])

    # Thrombophilia (binary)
    Trmphil[ind] ~ dbern(p.Trmphil[1])

    # Aortic stenosis (binary)
    Sten_aor[ind] ~ dbern(p.Sten_aor[1])
    
    # Gender (binary)
    Sex[ind] ~ dbern(p.Sex[1])
    
    # Chronic mitral valve prolapse (binary)
    ProlMitralCron[ind] ~ dbern(p.ProlMitralCron[1])
   
    # Non-infarctual pericarditis (binary)
    Peric_ni[ind] ~ dbern(p.Peric_ni[1])
    
    # Obesity (binary)
    Obes[ind] ~ dbern(p.Obes[1])
    
    # Neuromuscular disease (binary)
    Neuromusc[ind] ~ dbern(p.Neuromusc[1])
    
    # L-dopa use (binary)
    LDopa[ind] ~ dbern(p.LDopa[1])  

    # Hypoglycemia (binary)
    Ipoglic[ind] ~ dbern(p.Ipoglic[1])
    
    # Herpes Zooster (binary)
    HrepZ[ind] ~ dbern(p.HrepZ[1])

    # Smoke (binary)
    Fum[ind] ~ dbern(p.Fum[1])

    # Lower limbs fractures (binary)                     
    Fr[ind] ~ dbern(p.Fr[1])
    
    # Patent foramen ovale (binary)
    For_ov[ind] ~ dbern(p.For_ov[1])

    # Inspired oxygen fraction (continuous, hyper-restricted)
    tilde.FiO2[ind] ~ dbeta(p.FiO2[1]*p.FiO2[2],(1-p.FiO2[1])*p.FiO2[2]) T(tol,1-tol)
    
    # Hiatal hernia (binary)
    Ern_ia[ind] ~ dbern(p.Ern_ia[1])
    
    # Endocarditis (binary)
    Endocard[ind] ~ dbern(p.Endocard[1])
    
    # Vasovagal syncope (binary)
    CVV[ind] ~ dbern(p.CVV[1])
    
    # Costochondritis (binary)
    Costcndr[ind] ~ dbern(p.Costcndr[1])
    
    # Central venous catheter (binary)
    CatVen[ind] ~ dbern(p.CatVen[1])

    # Age (continuous)
    tilde.Age[ind] ~ dbeta(p.Age[1]*p.Age[2],(1-p.Age[1])*p.Age[2]) T(tol,1-tol)
   
    # Cocaine/amphetamines use (binary)
    A_coc[ind] ~ dbern(p.A_coc[1])
           
    }


  #####  prior for continuous nodes
  
  p.pulmPerf[1] <- 0.5  
  p.pulmPerf[2] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.pulmPerf[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.pulmPerf[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.pulmPerf[5] ~ dgamma(89.4917,2.0304)

  p.Shu[1] <- 0.5
  p.Shu[2] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.Shu[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Shu[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Shu[5] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.Shu[6] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.Shu[7] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.Shu[8] ~ dgamma(89.4917,2.0304)

  p.driveCard[1] <- 0.5
  p.driveCard[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.driveCard[3] ~ dbeta(2/3*10,1/3*10) T(tol,1-tol)
  p.driveCard[4] ~ dbeta(2/3*10,1/3*10) T(tol,1-tol)
  p.driveCard[5] ~ dbeta(2/3*10,1/3*10) T(tol,1-tol)
  p.driveCard[6] ~ dgamma(89.4917,2.0304)

  p.minVent[1] <- 0.5
  p.minVent[2] ~ dbeta(4/5*10,1/5*10) T(tol,1-tol)
  p.minVent[3] ~ dbeta(4/5*10,1/5*10) T(tol,1-tol)
  p.minVent[4] ~ dbeta(11/15*10,4/15*10) T(tol,1-tol)
  p.minVent[5] <- 1/6
  p.minVent[6] <- 5/6
  p.minVent[7] <- 5/6
  p.minVent[8] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.minVent[9] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.minVent[10] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.minVent[11] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.minVent[12] ~ dgamma(89.4917,2.0304)

  p.leftPump[1] <- 0.5
  p.leftPump[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.leftPump[3] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.leftPump[4] ~ dbeta(13/15*10,2/15*10) T(p.leftPump[3],1-tol)
  p.leftPump[5] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.leftPump[6] <- 5/6
  p.leftPump[7] ~ dgamma(89.4917,2.0304)

  p.LHI[1] <- 0.5
  p.LHI[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.LHI[3] <- 1/6
  p.LHI[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.LHI[5] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.LHI[6] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.LHI[7] ~ dgamma(89.4917,2.0304)

  p.LHO[1] <- 0.5
  p.LHO[2] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.LHO[3] <- 5/6
  p.LHO[4] ~ dgamma(89.4917,2.0304)

  p.RHF[1] <- 0.5  
  p.RHF[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.RHF[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.RHF[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.RHF[5] ~ dgamma(89.4917,2.0304)

  p.RHI[1] <- 0.5
  p.RHI[2] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.RHI[3] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.RHI[4] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.RHI[5] ~ dgamma(89.4917,2.0304)

  p.RHO[1] <- 0.5  
  p.RHO[2] <- 1/6
  p.RHO[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.RHO[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.RHO[5] ~ dgamma(89.4917,2.0304)

  p.postLoad[1] <- 0.5  
  p.postLoad[2] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.postLoad[3] <- 5/6 
  p.postLoad[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.postLoad[5] ~ dgamma(89.4917,2.0304)

  p.resVas[1] <- 0.5
  p.resVas[2] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.resVas[3] <- 5/6
  p.resVas[4] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.resVas[5] ~ dbeta(19/30*10,11/30*10) T(tol,1-tol)
  p.resVas[6] ~ dbeta(9/10*10,1/10*10) T(tol,1-tol)
  p.resVas[7] ~ dgamma(89.4917,2.0304)

  p.Disidrat[1] ~ dbeta(0.5*30,0.5*30) T(tol,1-tol)
  p.Disidrat[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Disidrat[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Disidrat[4] ~ dgamma(89.4917,2.0304)
  
  p.Latt[1] ~ dbeta(0.5*30,0.5*30) T(tol,1-tol)
  p.Latt[2] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.Latt[3] ~ dgamma(89.4917,2.0304)  

  p.Emoglob[1] <- 0.5    
  p.Emoglob[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Emoglob[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.Emoglob[4] ~ dgamma(89.4917,2.0304)
  
  p.CO2[1] <- 0.5
  p.CO2[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.CO2[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.CO2[4] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.CO2[5] ~ dgamma(89.4917,2.0304)

  p.SaO2[1] <- 0.5  
  p.SaO2[2] <- p.CO2[2]  ### ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.SaO2[3] <- p.CO2[3]  ### ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.SaO2[4] <- p.CO2[4]  ### ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.SaO2[5] <- 7/30      ### ~ dbeta(7/30*10,23/30*10) T(tol,1-tol)
  p.SaO2[6] ~ dbeta(7/30*10,23/30*10) T(tol,1-tol)
  p.SaO2[7] ~ dgamma(89.4917,2.0304)
        
  p.O2[1] ~ dbeta(0.5*30,0.5*30) T(tol,1-tol)
  p.O2[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.O2[3] ~ dgamma(89.4917,2.0304)

  p.FC[1] <- 0.5
  p.FC[2] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.FC[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FC[4] ~ dbeta(7/30*10,23/30*10) T(tol,1-tol)
  p.FC[5] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.FC[6] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.FC[7] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.FC[8] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FC[9] ~ dgamma(89.4917,2.0304)

  p.PA[1] ~ dbeta(0.5*30,0.5*30) T(tol,1-tol)
  p.PA[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.PA[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.PA[4] ~ dgamma(89.4917,2.0304)

  p.pH[1] <- 0.5
  p.pH[2] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.pH[3] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.pH[4] ~ dbeta(11/30*10,19/30*10) T(tol,p.pH[3])
  p.pH[5] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.pH[6] ~ dbeta(11/30*10,19/30*10) T(tol,p.pH[5])
  p.pH[7] ~ dgamma(89.4917,2.0304)

  p.OT[1] <- 0.5
  p.OT[2] ~ dbeta(1/6*50,5/6*50) T(tol,1-tol)
  p.OT[3] ~ dbeta(1/6*50,5/6*50) T(tol,1-tol)
  p.OT[4] ~ dbeta(5/6*50,1/6*50) T(tol,1-tol)
  p.OT[5] ~ dbeta(5/6*50,1/6*50) T(tol,1-tol)
  p.OT[6] ~ dgamma(89.4917,2.0304)

  p.FT3[1] <- 0.5
  p.FT3[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FT3[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FT3[4] ~ dgamma(89.4917,2.0304)

  p.FT4[1] <- 0.5
  p.FT4[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FT4[3] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.FT4[4] ~ dgamma(89.4917,2.0304)

  p.TSH[1] <- 0.5
  p.TSH[2] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.TSH[3] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.TSH[4] ~ dbeta(1/6*10,5/6*10) T(tol,1-tol)
  p.TSH[5] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.TSH[6] ~ dgamma(89.4917,2.0304)

  p.GCS[1] ~ dbeta(0.5*30,0.5*30) T(tol,1-tol)
  p.GCS[2] ~ dbeta(23/30*10,7/30*10) T(tol,1-tol)
  p.GCS[3] ~ dbeta(19/30*10,11/30*10) T(tol,1-tol)
  p.GCS[4] ~ dbeta(5/6*10,1/6*10) T(tol,1-tol)
  p.GCS[5] ~ dbeta(7/10*10,3/10*10) T(tol,1-tol)
  p.GCS[6] ~ dgamma(89.4917,2.0304)

  p.FiO2[1] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.FiO2[2] ~ dgamma(0.01,0.01)

  p.Age[1] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Age[2] ~ dgamma(0.01,0.01)
  

  #####  prior for semi-stochastic variables

  p.ChestPain[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)

  p.acuteAnem[1] ~ dbeta(0.05*50,0.95*50) T(tol,1-tol)

  p.Peric[1] ~ dbeta(0.002*50,0.998*50) T(tol,1-tol)
  
  p.Veget_end[1] ~ dbeta(0.0001*30,0.9999*30) T(tol,1-tol)
  
  p.Fertil[1] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)

  alpha1.Grav[1] <- 0.97*20
  alpha1.Grav[2] <- 0.01*20
  alpha1.Grav[3] <- 0.02*20
  p.Grav[1:3] ~ ddirch(alpha1.Grav[1:3])

  p.Leucoc_o[1] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.Leucoc_o[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  
    
  #####  prior for binary nodes
  
  p.Sex[1] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  
  p.Fr[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)

  p.Fum[1] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  
  p.Trmphil[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  
  p.Obes[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)

  p.Neuromusc[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
      
  p.CatVen[1] ~ dbeta(0.005*30,0.995*30) T(tol,1-tol)
   
  p.URI[1] ~ dbeta(0.09*30,0.91*30) T(tol,1-tol)

  p.Peric_ni[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
    
  p.ProlMitralCron[1] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
 
  p.Ulc_pep[1] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
    
  p.For_ov[1] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)

  p.WPW[1] ~ dbeta(0.001*50,0.999*50) T(tol,1-tol)
  
  p.A_coc[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
             
  p.Sten_aor[1] ~ dbeta(0.05*400,0.95*400)  T(tol,1-tol)

  p.Endocard[1] ~ dbeta(0.0005*400,0.9995*400) T(tol,1-tol)
     
  p.Ipoglic[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
      
  p.CVV[1] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
    
  p.Ern_ia[1] ~ dbeta(0.075*30,0.925*30) T(tol,1-tol)
  
  p.Costcndr[1] ~ dbeta(0.005*100,0.995*100) T(tol,1-tol)
      
  p.HrepZ[1] ~ dbeta(0.008*100,0.992*100) T(tol,1-tol)
  
  p.LDopa[1] ~ dbeta(0.005*30,0.995*30) T(tol,1-tol)

  p.Allett[1] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.Allett[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Allett[1] <- 0

  p.Calze[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.Calze[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Calze[1] <- 0

  p.FlebArtInf[1] ~ dbeta(0.01*50,0.99*50) T(tol,1-tol)
  p.FlebArtInf[2] ~ dbeta(0.8*50,0.2*50) T(tol,1-tol)
  eta.FlebArtInf[1] <- 0

  p.FlebArtInf_R[1] ~ dbeta(0.01*50,0.99*50) T(tol,1-tol)
  p.FlebArtInf_R[2] ~ dbeta(0.8*50,0.2*50) T(tol,1-tol)
  eta.FlebArtInf_R[1] <- 0

  p.Dol_g[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.Dol_g[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.Dol_g[1] <- 0

  p.Eco_c[1] ~ dbeta(0.03*50,0.97*50) T(tol,1-tol)
  p.Eco_c[2] ~ dbeta(0.7*50,0.3*50) T(tol,1-tol)
  eta.Eco_c[1] <- 0

  p.Eco_d[1] ~ dbeta(0.03*50,0.97*50) T(tol,1-tol)
  p.Eco_d[2] ~ dbeta(0.7*50,0.3*50) T(tol,1-tol)
  eta.Eco_d[1] <- 0

  p.blastCris[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.blastCris[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.blastCris[1] <- 0

  p.Arit_sotto[1] ~ dbeta(0.035*30,0.965*30) T(tol,1-tol)
  p.Arit_sotto[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Arit_sotto[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.Arit_sotto[4] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Arit_sotto[1] <- 0

  p.preTVP[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.preTVP[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.preTVP[1] <- 0

  p.Linfoc_o[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Linfoc_o[2] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  eta.Linfoc_o[1] <- 0

  p.Neutrof_o[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Neutrof_o[2] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  eta.Neutrof_o[1] <- 0

  p.MSC[1] ~ dbeta(0.0021*30,0.9979*30) T(tol,1-tol)
  p.MSC[2] ~ dbeta(0.65*10,0.35*10) T(tol,1-tol)
  eta.MSC[1] <- 0

  p.Mors[1] ~ dbeta(0.00025*30,0.99975*30) T(tol,1-tol)
  p.Mors[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Mors[1] <- 0

  p.Inc_sfint[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Inc_sfint[2] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  eta.Inc_sfint[1] <- 0

  p.Neur_foc[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Neur_foc[2] ~ dbeta(0.9*30,0.1*30) T(tol,1-tol)
  p.Neur_foc[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.Neur_foc[1] <- 0

  p.DilArtPolm_ec[1] ~ dbeta(0.012*30,0.988*30) T(tol,1-tol)
  p.DilArtPolm_ec[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.DilArtPolm_ec[3] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.DilArtPolm_ec[1] <- 0

  p.Fr_cost[1] ~ dbeta(0.02*50,0.98*50) T(tol,1-tol)
  p.Fr_cost[2] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  eta.Fr_cost[1] <- 0

  p.sindrMW[1] ~ dbeta(0.0001*30,0.9999*30) T(tol,1-tol)
  p.sindrMW[2] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  eta.sindrMW[1] <- 0

  p.Colic[1] ~ dbeta(0.0001*100,0.9999*100)  T(tol,1-tol)
  p.Colic[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Colic[1] <- 0

  p.Hemo[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Hemo[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Hemo[1] <- 0

  p.antiTromb[1] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.antiTromb[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.antiTromb[1] <- 0

  p.protS[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.protS[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.protS[1] <- 0

  p.protC[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.protC[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.protC[1] <- 0

  p.Leid[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.Leid[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Leid[1] <- 0

  p.fIIG[1] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.fIIG[2] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  eta.fIIG[1] <- 0

  p.omocist[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.omocist[2] ~ dbeta(0.18*10,0.82*10) T(tol,1-tol)
  eta.omocist[1] <- 0

  p.abPhos[1] ~ dbeta(0.006*30,0.994*30) T(tol,1-tol)
  p.abPhos[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.abPhos[1] <- 0

  p.FattVIII[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.FattVIII[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.FattVIII[1] <- 0

  p.intersPol[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.intersPol[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.intersPol[1] <- 0

  p.Psico[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Psico[2] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.Psico[1] <- 0

  p.CistBoll[1] ~ dbeta(0.002*30,0.998*30) T(tol,1-tol)
  p.CistBoll[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.CistBoll[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.CistBoll[1] <- 0

  p.bronect[1] ~ dbeta(0.001*50,0.999*50) T(tol,1-tol)
  p.bronect[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.bronect[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.bronect[1] <- 0

  p.CavColl[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.CavColl[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.CavColl[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.CavColl[1] <- 0

  p.Dressl[1] ~ dbeta(0.00001*30,0.99999*30) T(tol,1-tol)
  p.Dressl[2] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.Dressl[3] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.Dressl[1] <- 0

  p.Mioglob[1] ~ dbeta(0.01*100,0.99*100) T(tol,1-tol)
  p.Mioglob[2] ~ dbeta(0.98*50,0.02*50) T(tol,1-tol)
  eta.Mioglob[1] <- 0

  p.Latt_o[1] ~ dbeta(0.05*50,0.95*50) T(tol,1-tol)
  p.Latt_o[2] ~ dbeta(0.95*50,0.05*50) T(tol,1-tol)
  eta.Latt_o[1] <- 0
  
  p.Setto_par[1] ~ dbeta(0.00005*100,0.99995*100) T(tol,1-tol)
  p.Setto_par[2] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.Setto_par[1] <- 0

  p.IpoCinVdx[1] ~ dbeta(0.0015*50,0.9985*50) T(tol,1-tol)
  p.IpoCinVdx[2] ~ dbeta(0.25*10,0.75*10) T(tol,1-tol)
  eta.IpoCinVdx[1] <- 0

  p.Epatom[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.Epatom[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Epatom[1] <- 0

  p.Turg_giu[1] ~ dbeta(0.15*20,0.85*20) T(tol,1-tol)
  p.Turg_giu[2] ~ dbeta(0.89*10,0.11*10) T(tol,1-tol)
  eta.Turg_giu[1] <- 0

  p.reflVenSEp[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.reflVenSEp[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.reflVenSEp[1] <- 0

  p.Olig[1] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.Olig[2] ~ dbeta(0.85*10,0.15*10) T(tol,1-tol)
  eta.Olig[1] <- 0

  p.BNP_o[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.BNP_o[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.BNP_o[1] <- 0

  p.DiamArtPolm[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.DiamArtPolm[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.DiamArtPolm[3] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.DiamArtPolm[1] <- 0

  p.calVasPol[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.calVasPol[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.calVasPol[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.calVasPol[1] <- 0

  p.calBron[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.calBron[2] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  p.calBron[3] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  eta.calBron[1] <- 0

  p.airTrap[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.airTrap[2] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  p.airTrap[3] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  eta.airTrap[1] <- 0

  p.parBron[1] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.parBron[2] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  p.parBron[3] ~ dbeta(0.68*10,0.32*10) T(tol,1-tol)
  eta.parBron[1] <- 0

  p.Amil_o[1] ~ dbeta(0.015*30,0.985*30) T(tol,1-tol)
  p.Amil_o[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Amil_o[1] <- 0

  p.I_par_lib[1] ~ dbeta(0.00021*30,0.99979*30) T(tol,1-tol)
  p.I_par_lib[2] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.I_par_lib[1] <- 0

  p.Ip_ban_mod[1] ~ dbeta(0.0021*30,0.9979*30) T(tol,1-tol)
  p.Ip_ban_mod[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Ip_ban_mod[1] <- 0

  p.Sovracc[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.Sovracc[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Sovracc[1] <- 0

  p.DevAsse[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.DevAsse[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.DevAsse[1] <- 0

  p.Anisosfig[1] ~ dbeta(0.002*30,0.998*30) T(tol,1-tol)
  p.Anisosfig[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.Anisosfig[1] <- 0

  p.EmatIntram[1] ~ dbeta(0.0005*30,0.9995*30) T(tol,1-tol)
  p.EmatIntram[2] ~ dbeta(0.84*10,0.16*10) T(tol,1-tol)
  eta.EmatIntram[1] <- 0

  p.Piro[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.Piro[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Piro[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Piro[1] <- 0

  p.naus[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.naus[2] ~ dbeta(0.075*30,0.925*30) T(tol,1-tol)
  eta.naus[1] <- 0

  p.catecUr[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.catecUr[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.catecUr[1] <- 0

  p.masCen[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.masCen[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.masCen[1] <- 0

  p.Nod[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.Nod[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.Nod[1] <- 0

  p.Shock[1] <- 0    
  p.Shock[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Shock[1] <- 0

  p.Estrog[1] <- 0
  p.Estrog[2] ~ dbeta(0.16*10,0.84*10) T(tol,1-tol)
  eta.Estrog[1] <- 0
  
  p.Tcp[1] <- 0
  p.Tcp[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Tcp[3] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.Tcp[1] <- 0

  p.Neo_p[1] <- 0
  p.Neo_p[2] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.Neo_p[1] <- 0
  
  p.Leucem[1] <- 0
  p.Leucem[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.Leucem[1] <- 0
  
  p.ProlMitralAcut[1] <- 0
  p.ProlMitralAcut[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  eta.ProlMitralAcut[1] <- 0
    
  p.Colec[1] <- 0
  p.Colec[2] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.Colec[1] <- 0
  
  p.CuPC[1] <- 0
  p.CuPC[2] ~ dbeta(0.85*10,0.15*10) T(tol,1-tol)
  eta.CuPC[1] <- 0
    
  p.Tireotos[1] <- 0
  p.Tireotos[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Tireotos[1] <- 0
  
  p.InP[1] <- 0
  p.InP[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.InP[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.InP[1] <- 0
    
  p.Cia[1] <- 0
  p.Cia[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Cia[1] <- 0
    
  p.Agit_o[1] <- 0
  p.Agit_o[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Agit_o[1] <- 0
  
  p.bronAer[1] <- 0
  p.bronAer[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.bronAer[1] <- 0
  
  p.CK_MB[1] <- 0
  p.CK_MB[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.CK_MB[3] ~ dbeta(0.85*10,0.15*10) T(tol,1-tol)
  eta.CK_MB[1] <- 0

  p.TVPscore[1] <- 0
  p.TVPscore[2] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.TVPscore[3] ~ dbeta(0.05*10,0.95*10) T(tol,1-tol)
  p.TVPscore[4] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.TVPscore[5] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.TVPscore[1] ~ dnorm(0,1)

  p.TrCdx[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.TrCdx[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.TrCdx[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.TrCdx[1] ~ dnorm(0,1)

  p.massCer[1] ~ dbeta(0.002*30,0.998*30) T(tol,1-tol)
  p.massCer[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.massCer[3] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  eta.massCer[1] ~ dnorm(0,1)

  p.Sepsi[1] ~ dbeta(0.007*30,0.993*30) T(tol,1-tol)
  p.Sepsi[2] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Sepsi[3] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  eta.Sepsi[1] ~ dnorm(0,1)

  p.Arit_sopra_acut[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Arit_sopra_acut[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Arit_sopra_acut[3] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Arit_sopra_acut[4] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Arit_sopra_acut[1] ~ dnorm(0,1)

  p.EmbArte[1] ~ dbeta(0.0015*50,0.9985*50) T(tol,1-tol)
  p.EmbArte[2] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.EmbArte[3] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  eta.EmbArte[1] ~ dnorm(0,1)

  p.Coronar[1] ~ dbeta(0.07*30,0.93*30) T(tol,1-tol)
  p.Coronar[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.Coronar[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Coronar[4] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Coronar[1] ~ dnorm(0,1)

  p.BBSx[1] ~ dbeta(0.07*30,0.93*30) T(tol,1-tol)
  p.BBSx[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.BBSx[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.BBSx[4] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.BBSx[1] ~ dnorm(0,1)

  p.Miocard[1] ~ dbeta(0.00008*400,0.99992*400) T(tol,1-tol)
  p.Miocard[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Miocard[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.Miocard[4] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Miocard[1] ~ dnorm(0,1)

  p.Insuf_aorAcut[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Insuf_aorAcut[2] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  p.Insuf_aorAcut[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Insuf_aorAcut[4] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.Insuf_aorAcut[1] ~ dnorm(0,1)

  p.InsuffMitrAcut[1] ~ dbeta(0.05*50,0.95*50) T(tol,1-tol)
  p.InsuffMitrAcut[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.InsuffMitrAcut[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.InsuffMitrAcut[4] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.InsuffMitrAcut[5] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.InsuffMitrAcut[1] ~ dnorm(0,1)

  p.Perit[1] ~ dbeta(0.002*100,0.998*100) T(tol,1-tol)
  p.Perit[2] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.Perit[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Perit[1] ~ dnorm(0,1)

  p.CordTend[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.CordTend[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.CordTend[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.CordTend[1] ~ dnorm(0,1)

  p.Ipertpolm[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Ipertpolm[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Ipertpolm[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Ipertpolm[4] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  eta.Ipertpolm[1] ~ dnorm(0,1)

  p.Infiam_b[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Infiam_b[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.Infiam_b[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Infiam_b[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Infiam_b[5] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Infiam_b[6] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Infiam_b[7] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.Infiam_b[8] <- 1-tol
  p.Infiam_b[9] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Infiam_b[1] ~ dnorm(0,1)

  p.Infiam_nb[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Infiam_nb[2] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Infiam_nb[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Infiam_nb[4] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Infiam_nb[5] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.Infiam_nb[6] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Infiam_nb[1] ~ dnorm(0,1)

  p.extraSist[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.extraSist[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.extraSist[3] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.extraSist[1] ~ dnorm(0,1)

  p.cpl[1] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.cpl[2] ~ dbeta(0.48*10,0.52*10) T(tol,1-tol)
  p.cpl[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.cpl[1] ~ dnorm(0,1)

  p.HemoPeric[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.HemoPeric[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.HemoPeric[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.HemoPeric[1] ~ dnorm(0,1)

  p.DilatVsx[1] ~ dbeta(0.00001*30,0.99999*30) T(tol,1-tol)
  p.DilatVsx[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.DilatVsx[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.DilatVsx[1] ~ dnorm(0,1)

  p.Vers_peric[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.Vers_peric[2] ~ dbeta(0.9*100,0.1*100) T(tol,1-tol)
  p.Vers_peric[3] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Vers_peric[1] ~ dnorm(0,1)

  p.carMeg[1] ~ dbeta(0.005*30,0.995*30) T(tol,1-tol)
  p.carMeg[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.carMeg[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.carMeg[4] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.carMeg[1] ~ dnorm(0,1)

  p.IschCer[1] ~ dbeta(0.00001*30,0.99999*30) T(tol,1-tol)
  p.IschCer[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.IschCer[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.IschCer[4] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.IschCer[5] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.IschCer[6] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.IschCer[7] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.IschCer[1] ~ dnorm(0,1)

  p.SuspHeaDrivTmp[1] ~ dbeta(0.1*50,0.9*50) T(tol,1-tol)
  p.SuspHeaDrivTmp[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.SuspHeaDrivTmp[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.SuspHeaDrivTmp[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.SuspHeaDrivTmp[5] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.SuspHeaDrivTmp[1] ~ dnorm(0,1)

  p.PerdCosc[1] ~ dbeta(0.01*50,0.99*50)  T(tol,1-tol)
  p.PerdCosc[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.PerdCosc[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.PerdCosc[1] ~ dnorm(0,1)

  p.Confus[1] ~ dbeta(0.02*30,0.98*30)  T(tol,1-tol)
  p.Confus[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Confus[3] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Confus[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Confus[5] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Confus[6] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Confus[1] ~ dnorm(0,1)

  p.LDH_o[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.LDH_o[2] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  p.LDH_o[3] ~ dbeta(0.76*10,0.24*10) T(tol,1-tol)
  p.LDH_o[4] ~ dbeta(0.94*30,0.06*30) T(tol,1-tol)
  p.LDH_o[5] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  p.LDH_o[6] ~ dbeta(0.85*10,0.15*10) T(tol,1-tol)
  p.LDH_o[7] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.LDH_o[1] ~ dnorm(0,1)

  p.ARDS[1] ~ dbeta(0.007*50,0.993*50) T(tol,1-tol)
  p.ARDS[2] ~ dbeta(0.34*10,0.66*10) T(tol,1-tol)
  p.ARDS[3] ~ dbeta(0.25*10,0.75*10) T(tol,1-tol)
  p.ARDS[4] ~ dbeta(0.05*50,0.95*50) T(tol,1-tol)
  p.ARDS[5] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.ARDS[6] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  eta.ARDS[1] ~ dnorm(0,1)

  p.Spas_es[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Spas_es[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.Spas_es[3] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Spas_es[4] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.Spas_es[1] ~ dnorm(0,1)

  p.vomit[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.vomit[2] ~ dbeta(0.45*10,0.55*10) T(tol,1-tol)
  p.vomit[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.vomit[1] ~ dnorm(0,1)

  p.Ple[1] ~ dbeta(0.0001*100,0.9999*100) T(tol,1-tol)
  p.Ple[2] ~ dbeta(0.08*30,0.92*30) T(tol,1-tol)
  p.Ple[3] ~ dbeta(0.18*10,0.82*10) T(tol,1-tol)
  p.Ple[4] ~ dbeta(0.065*30,0.935*30) T(tol,1-tol)
  p.Ple[5] ~ dbeta(0.065*30,0.935*30) T(tol,1-tol)
  p.Ple[6] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Ple[1] ~ dnorm(0,1)

  p.Ort[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.Ort[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Ort[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Ort[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Ort[5] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Ort[6] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.Ort[1] ~ dnorm(0,1)

  p.Atel[1] ~ dbeta(0.0025*30,0.9975*30) T(tol,1-tol)
  p.Atel[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Atel[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.Atel[4] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Atel[5] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Atel[6] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.Atel[1] ~ dnorm(0,1)

  p.Ostr_circ_sist[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.Ostr_circ_sist[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Ostr_circ_sist[3] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Ostr_circ_sist[1] ~ dnorm(0,1)

  p.Ipo_orto[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.Ipo_orto[2] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.Ipo_orto[3] ~ dbeta(0.07*30,0.93*30) T(tol,1-tol)
  p.Ipo_orto[4] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Ipo_orto[5] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.Ipo_orto[6] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Ipo_orto[1] ~ dnorm(0,1)

  p.Dispnea[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.Dispnea[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Dispnea[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Dispnea[4] ~ dbeta(0.9*50,0.1*50) T(tol,1-tol)
  p.Dispnea[5] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.Dispnea[1] ~ dnorm(0,1)

  p.E_dia[1] ~ dbeta(0.02*100,0.98*100) T(tol,1-tol)
  p.E_dia[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.E_dia[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.E_dia[1] ~ dnorm(0,1)

  p.Emo[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.Emo[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.Emo[3] ~ dbeta(0.35*10,0.65*10) T(tol,1-tol)
  p.Emo[4] ~ dbeta(0.0025*30,0.9975*30) T(tol,1-tol)
  p.Emo[5] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.Emo[6] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  eta.Emo[1] ~ dnorm(0,1)

  p.AdenIlar[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.AdenIlar[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.AdenIlar[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.AdenIlar[4] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  p.AdenIlar[5] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.AdenIlar[1] ~ dnorm(0,1)

  p.Trop_o[1] ~ dbeta(0.0005*30,0.9995*30) T(tol,1-tol)
  p.Trop_o[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  eta.Trop_o[1] <- 0

  p.InsTric_o[1] ~ dbeta(0.001*100,0.999*100) T(tol,1-tol)
  p.InsTric_o[2] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  p.InsTric_o[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.InsTric_o[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.InsTric_o[5] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.InsTric_o[1] ~ dnorm(0,1)

  p.DilatVdx[1] ~ dbeta(0.00002*100,0.99998*100) T(tol,1-tol)
  p.DilatVdx[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.DilatVdx[3] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.DilatVdx[1] ~ dnorm(0,1)

  p.BBD[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.BBD[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.BBD[3] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.BBD[1] ~ dnorm(0,1)

  p.Caduta[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.Caduta[2] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.Caduta[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Caduta[1] ~ dnorm(0,1)

  p.grouGl[1] ~ dbeta(0.001*30,0.999*30) T(tol,1-tol)
  p.grouGl[2] ~ dbeta(0.75*10,0.25*10) T(tol,1-tol)
  p.grouGl[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.grouGl[4] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.grouGl[5] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.grouGl[1] ~ dnorm(0,1)
    
  p.RHtrigg[1] <- 0
  p.RHtrigg[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.RHtrigg[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.RHtrigg[4] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.RHtrigg[1] ~ dnorm(0,1)
    
  p.FANS[1] <- 0
  p.FANS[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.FANS[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.FANS[1] ~ dnorm(0,1)
          
  p.Fever[1] <- 0
  p.Fever[2] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.Fever[3] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.Fever[4] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  eta.Fever[1] ~ dnorm(0,1)
  
  p.MioStretch[1] <- 0
  p.MioStretch[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.MioStretch[3] ~ dbeta(0.12*10,0.88*10) T(tol,1-tol)
  p.MioStretch[4] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.MioStretch[5] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.MioStretch[1] ~ dnorm(0,1)
  
  p.D_a_pu[1] <- 0
  p.D_a_pu[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.D_a_pu[3] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.D_a_pu[1] ~ dnorm(0,1)
  
  p.D_cen[1] <- 0
  p.D_cen[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.D_cen[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.D_cen[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.D_cen[5] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.D_cen[6] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.D_cen[7] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.D_cen[8] ~ dbeta(0.25*10,0.75*10) T(tol,1-tol)
  eta.D_cen[1] ~ dnorm(0,1)
  
  p.D_pl[1] <- 0
  p.D_pl[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.D_pl[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.D_pl[4] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  p.D_pl[5] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.D_pl[6] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.D_pl[7] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.D_pl[1] ~ dnorm(0,1)
  
  p.D_quadr_sup[1] <- 0
  p.D_quadr_sup[2] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  p.D_quadr_sup[3] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.D_quadr_sup[4] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.D_quadr_sup[5] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.D_quadr_sup[1] ~ dnorm(0,1)
  
  p.D_par[1] <- 0
  p.D_par[2] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.D_par[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.D_par[4] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.D_par[5] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  p.D_par[6] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.D_par[7] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.D_par[8] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  eta.D_par[1] ~ dnorm(0,1)
  
  p.E_tap_te[1] <- 0
  p.E_tap_te[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.E_tap_te[3] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.E_tap_te[4] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.E_tap_te[5] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  eta.E_tap_te[1] ~ dnorm(0,1)
  
  p.Scint_bias[1] <- 0
  p.Scint_bias[2] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Scint_bias[3] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Scint_bias[4] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Scint_bias[5] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  eta.Scint_bias[1] ~ dnorm(0,1)
  
  p.S_per_r[1] <- 0
  p.S_per_r[2] ~ dbeta(0.98*30,0.02*30) T(tol,1-tol)
  eta.S_per_r[1] <- 0
  
  p.ECG_dx[1] <- 0
  p.ECG_dx[2] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.ECG_dx[3] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  eta.ECG_dx[1] ~ dnorm(0,1)

  p.Alcol[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)

  p.Feocrom[1] ~ dbeta(0.0002*400,0.9998*400) T(tol,1-tol)
  
  p.Interst[1] ~ dbeta(0.01*250,0.99*250) T(tol,1-tol)

  p.MNodSA[1] ~ dbeta(0.1*30,0.9*30) T(tol,1-tol)

  p.metabAlk[1] ~ dbeta(0.3*50,0.7*50) T(tol,1-tol)

  p.InsVenCr[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.InsVenCr[2] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)

  p.Card_dil[1] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  p.Card_dil[2] ~ dbeta(0.07*30,0.93*30) T(tol,1-tol)
 
  p.Asma[1] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)
  p.Asma[2] ~ dbeta(0.025*30,0.975*30) T(tol,1-tol)

  p.IschCerCron[1] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.IschCerCron[2] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.IschCerCron[1] <- 0
  
  p.Refl_ge[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Refl_ge[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  eta.Refl_ge[1] <- 0
  
  p.cronAnem[1] ~ dbeta(0.1*50,0.9*50) T(tol,1-tol)
  p.cronAnem[2] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  eta.cronAnem[1] <- 0
    
  p.Defed[1] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Defed[2] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.Defed[1] <- 0
  
  p.Aneur[1] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  p.Aneur[2] ~ dbeta(0.008*30,0.992*30) T(tol,1-tol)
  p.Aneur[3] ~ dbeta(0.03*30,0.97*30) T(tol,1-tol)
  eta.Aneur[1] <- 0
  eta.Aneur[2] <- 0
        
  p.Aort_diss[1] ~ dbeta(0.0007*30,0.9993*30) T(tol,1-tol)
  p.Aort_diss[2] ~ dbeta(0.0003*30,0.9997*30) T(tol,1-tol)
  p.Aort_diss[3] ~ dbeta(0.14*10,0.86*10) T(tol,1-tol)
  eta.Aort_diss[1] <- 0
  eta.Aort_diss[2] <- 0
    
  p.Neo_o[1] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Neo_o[2] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.Neo_o[3] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  eta.Neo_o[1] <- 0
  eta.Neo_o[2] <- 0
    
  p.Litias[1] ~ dbeta(0.008*30,0.992*30) T(tol,1-tol)
  p.Litias[2] ~ dbeta(0.012*30,0.988*30) T(tol,1-tol)
  p.Litias[3] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  eta.Litias[1] <- 0    
  eta.Litias[2] <- 0  
    
  p.Ddimer[1] ~ dbeta(0.25*10,0.75*10) T(tol,1-tol)
  p.Ddimer[2] ~ dbeta(0.93*30,0.07*30) T(tol,1-tol)
  p.Ddimer[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Ddimer[4] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  eta.Ddimer[1] ~ dnorm(0,1)
  
  p.IpertSx[1] ~ dbeta(0.0001*30,0.9999*30) T(tol,1-tol)
  p.IpertSx[2] ~ dbeta(0.05*30,0.95*30) T(tol,1-tol)
  p.IpertSx[3] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.IpertSx[4] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  eta.IpertSx[1] ~ dnorm(0,1)
    
  p.Arit_sopra_cron[1] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.Arit_sopra_cron[2] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Arit_sopra_cron[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Arit_sopra_cron[4] ~ dbeta(0.45*10,0.55*10) T(tol,1-tol)
  eta.Arit_sopra_cron[1] ~ dnorm(0,1)

  p.Pancr[1] ~ dbeta(0.0017*400,0.9983*400) T(tol,1-tol)
  p.Pancr[2] ~ dbeta(0.0017*400,0.9983*400) T(tol,1-tol)
  p.Pancr[3] ~ dbeta(0.02*30,0.98*30) T(tol,1-tol)
  p.Pancr[4] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  eta.Pancr[1] ~ dnorm(0,1)
  eta.Pancr[2] ~ dnorm(0,1)
    
  p.cronIpert[1] ~ dbeta(0.3*20,0.7*20) T(tol,1-tol)
  p.cronIpert[2] ~ dbeta(0.3*20,0.7*20) T(tol,1-tol)
  p.cronIpert[3] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.cronIpert[4] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  eta.cronIpert[1] ~ dnorm(0,1)
  eta.cronIpert[2] ~ dnorm(0,1)
      
    
  #####  prior for nominal nodes
  
  alpha1.Chir[1] <- 0.75*20
  alpha1.Chir[2] <- 0.15*20
  alpha1.Chir[3] <- 0.10*20
  p.Chir[1:3] ~ ddirch(alpha1.Chir[1:3])  
  #
  alpha2.Chir[1] <- 0.15*10
  alpha2.Chir[2] <- 0.15*10
  alpha2.Chir[3] <- 0.70*10
  p.Chir[4:6] ~ ddirch(alpha2.Chir[1:3])
  #
  eta.Chir[1] <- 0
  eta.Chir[2] <- 0

  alpha1.TC_s_new[1] <- 0.965*20
  alpha1.TC_s_new[2] <- 0.025*20
  alpha1.TC_s_new[3] <- 0.01*20
  p.TC_s_new[1:3] ~ ddirch(alpha1.TC_s_new[1:3])
  #
  alpha2.TC_s_new[1] <- 0.025*10
  alpha2.TC_s_new[2] <- 0.95*10
  alpha2.TC_s_new[3] <- 0.025*10
  p.TC_s_new[4:6] ~ ddirch(alpha2.TC_s_new[1:3])
  #
  alpha3.TC_s_new[1] <- 0.001*10
  alpha3.TC_s_new[2] <- 0.05*10
  alpha3.TC_s_new[3] <- 0.949*10
  p.TC_s_new[7:9] ~ ddirch(alpha3.TC_s_new[1:3])
  #
  eta.TC_s_new[1] <- 0
  eta.TC_s_new[2] <- 0
  
  alpha1.BPCO[1] <- 0.99*20
  alpha1.BPCO[2] <- 0.005*20
  alpha1.BPCO[3] <- 0.005*20
  p.BPCO[1:3] ~ ddirch(alpha1.BPCO[1:3])
  #
  alpha2.BPCO[1] <- 0.01*10
  alpha2.BPCO[2] <- 0.9*10
  alpha2.BPCO[3] <- 0.09*10
  p.BPCO[4:6] ~ ddirch(alpha2.BPCO[1:3])
  #
  alpha3.BPCO[1] <- 0.005*10
  alpha3.BPCO[2] <- 0.295*10
  alpha3.BPCO[3] <- 0.7*10
  p.BPCO[7:9] ~ ddirch(alpha3.BPCO[1:3])   
  #
  eta.BPCO[1] <- 0
  eta.BPCO[2] <- 0
   
  p.OpacRX[1] <- 1
  p.OpacRX[2] <- 0
  p.OpacRX[3] <- 0
  #
  p.OpacRX[4] <- 1-p.OpacRX[5]
  p.OpacRX[5] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.OpacRX[6] <- tol
  #
  p.OpacRX[7] <- tol
  p.OpacRX[8] <- tol 
  p.OpacRX[9] <- 1-tol 
  #
  eta.OpacRX[1] <- 0
  eta.OpacRX[2] <- 0

  alpha1.Prof[1] <- 0.95*10
  alpha1.Prof[2] <- 0.03*10
  alpha1.Prof[3] <- 0.02*10
  p.Prof[1:3] ~ ddirch(alpha1.Prof[1:3])
  #
  p.Prof[4] <- 1-p.Prof[6]
  p.Prof[5] <- tol
  p.Prof[6] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  #
  p.Prof[7] <- 1-p.Prof[8]
  p.Prof[8] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Prof[9] <- tol  
  #
  alpha4.Prof[1] <- 0.80*10
  alpha4.Prof[2] <- 0.05*10
  alpha4.Prof[3] <- 0.15*10
  p.Prof[10:12] ~ ddirch(alpha4.Prof[1:3])
  #
  eta.Prof[1] ~ dnorm(0,1)
  eta.Prof[2] ~ dnorm(0,1)
  
  alpha1.Agit[1] <- 0.8*10
  alpha1.Agit[2] <- 0.14*10
  alpha1.Agit[3] <- 0.06*10
  p.Agit[1:3] ~ ddirch(alpha1.Agit[1:3])
  #
  alpha2.Agit[1] <- 0.2*10
  alpha2.Agit[2] <- 0.4*10
  alpha2.Agit[3] <- 0.4*10
  p.Agit[4:6] ~ ddirch(alpha2.Agit[1:3])
  #
  alpha3.Agit[1] <- 0.3*10
  alpha3.Agit[2] <- 0.15*10
  alpha3.Agit[3] <- 0.55*10
  p.Agit[7:9] ~ ddirch(alpha3.Agit[1:3])
  #
  eta.Agit[1] ~ dnorm(0,1)
  eta.Agit[2] ~ dnorm(0,1)
 
  p.EP[1] <- 1
  p.EP[2] <- 0
  p.EP[3] <- 0
  #
  alpha2.EP[1] <- 0.5*10
  alpha2.EP[2] <- 0.15*10
  alpha2.EP[3] <- 0.35*10
  p.EP[4:6] ~ ddirch(alpha2.EP[1:3])
  #
  alpha3.EP[1] <- 0.7*10
  alpha3.EP[2] <- 0.25*10
  alpha3.EP[3] <- 0.05*10
  p.EP[7:9] ~ ddirch(alpha3.EP[1:3])
  #
  alpha4.EP[1] <- 0.65*10
  alpha4.EP[2] <- 0.15*10
  alpha4.EP[3] <- 0.2*10
  p.EP[10:12] ~ ddirch(alpha4.EP[1:3])
  #
  eta.EP[1] ~ dnorm(0,1)
  eta.EP[2] ~ dnorm(0,1)

  p.Ed_polm[1] <- 1
  p.Ed_polm[2] <- 0
  p.Ed_polm[3] <- 0
  #
  p.Ed_polm[4] <- tol
  p.Ed_polm[5] <- 1-p.Ed_polm[6]
  p.Ed_polm[6] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  p.Ed_polm[7] <- tol
  p.Ed_polm[8] <- 1-p.Ed_polm[9]
  p.Ed_polm[9] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  #
  eta.Ed_polm[1] ~ dnorm(0,1)
  eta.Ed_polm[2] ~ dnorm(0,1)
  
  p.IMA[1] <- 1
  p.IMA[2] <- 0
  p.IMA[3] <- 0
  #
  alpha2.IMA[1] <- 0.1*10
  alpha2.IMA[2] <- 0.3*10
  alpha2.IMA[3] <- 0.6*10
  p.IMA[4:6] ~ ddirch(alpha2.IMA[1:3])
  #
  alpha3.IMA[1] <- 0.3*10
  alpha3.IMA[2] <- 0.63*10
  alpha3.IMA[3] <- 0.07*10
  p.IMA[7:9] ~ ddirch(alpha3.IMA[1:3])
  #
  eta.IMA[1] ~ dnorm(0,1)
  eta.IMA[2] ~ dnorm(0,1)

  p.tamp[1] <- 1
  p.tamp[2] <- 0
  p.tamp[3] <- 0
  #
  alpha2.tamp[1] <- 0.2*10
  alpha2.tamp[2] <- 0.56*10
  alpha2.tamp[3] <- 0.24*10
  p.tamp[4:6] ~ ddirch(alpha2.tamp[1:3])
  #
  alpha3.tamp[1] <- 0.85*10
  alpha3.tamp[2] <- 0.075*10
  alpha3.tamp[3] <- 0.075*10
  p.tamp[7:9] ~ ddirch(alpha3.tamp[1:3])
  #
  eta.tamp[1] ~ dnorm(0,1)
  eta.tamp[2] ~ dnorm(0,1)

  p.DilArtPolm[1] <- 1
  p.DilArtPolm[2] <- 0
  p.DilArtPolm[3] <- 0
  #
  p.DilArtPolm[4] <- 1-p.DilArtPolm[6]
  p.DilArtPolm[5] <- tol
  p.DilArtPolm[6] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  #
  p.DilArtPolm[7] <- 1-p.DilArtPolm[8]
  p.DilArtPolm[8] ~ dbeta(0.6*10,0.4*10) T(tol,1-tol)
  p.DilArtPolm[9] <- tol
  #
  p.DilArtPolm[10] <- 1-p.DilArtPolm[11]
  p.DilArtPolm[11] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.DilArtPolm[12] <- tol
  #
  eta.DilArtPolm[1] ~ dnorm(0,1)
  eta.DilArtPolm[2] ~ dnorm(0,1)
   
  p.Iper_Diaf[1] <- 1
  p.Iper_Diaf[2] <- 0
  p.Iper_Diaf[3] <- 0
  #
  p.Iper_Diaf[4] <- 1-p.Iper_Diaf[5]
  p.Iper_Diaf[5] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Iper_Diaf[6] <- tol   
  #
  p.Iper_Diaf[7] <- 1-p.Iper_Diaf[8]
  p.Iper_Diaf[8] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.Iper_Diaf[9] <- tol
  #
  p.Iper_Diaf[10] <- tol 
  p.Iper_Diaf[11] <- tol 
  p.Iper_Diaf[12] <- 1-tol 
  #
  eta.Iper_Diaf[1] ~ dnorm(0,1)
  eta.Iper_Diaf[2] ~ dnorm(0,1)

  p.Consol[1] <- 1-p.Consol[2]
  p.Consol[2] ~ dbeta(0.0003*100,0.9997*100) T(tol,1-tol)
  p.Consol[3] <- tol
  #
  alpha2.Consol[1] <- 0.05*10
  alpha2.Consol[2] <- 0.81*10
  alpha2.Consol[3] <- 0.14*10
  p.Consol[4:6] ~ ddirch(alpha2.Consol[1:3])
  #
  p.Consol[7] <- tol 
  p.Consol[8] <- tol 
  p.Consol[9] <- 1-tol
  #
  p.Consol[10] <- 1-p.Consol[11]
  p.Consol[11] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.Consol[12] <- tol
  #
  p.Consol[13] <- 1-p.Consol[14]
  p.Consol[14] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.Consol[15] <- tol
  #
  p.Consol[16] <- 1-p.Consol[17]
  p.Consol[17] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.Consol[18] <- tol
  #
  eta.Consol[1] ~ dnorm(0,1)
  eta.Consol[2] ~ dnorm(0,1)

  p.EO_bsp[1] <- 1
  p.EO_bsp[2] <- 0
  p.EO_bsp[3] <- 0
  #
  alpha2.EO_bsp[1] <- 0.2*10
  alpha2.EO_bsp[2] <- 0.64*10
  alpha2.EO_bsp[3] <- 0.16*10
  p.EO_bsp[4:6] ~ ddirch(alpha2.EO_bsp[1:3])
  #
  alpha3.EO_bsp[1] <- 0.1*10
  alpha3.EO_bsp[2] <- 0.55*10
  alpha3.EO_bsp[3] <- 0.35*10
  p.EO_bsp[7:9] ~ ddirch(alpha3.EO_bsp[1:3])
  #
  alpha4.EO_bsp[1] <- 0.6*10
  alpha4.EO_bsp[2] <- 0.2*10
  alpha4.EO_bsp[3] <- 0.2*10
  p.EO_bsp[10:12] ~ ddirch(alpha4.EO_bsp[1:3])
  #
  alpha5.EO_bsp[1] <- 0.4*10
  alpha5.EO_bsp[2] <- 0.3*10
  alpha5.EO_bsp[3] <- 0.3*10
  p.EO_bsp[13:15] ~ ddirch(alpha5.EO_bsp[1:3])  
  #
  alpha6.EO_bsp[1] <- 0.8*10
  alpha6.EO_bsp[2] <- 0.15*10
  alpha6.EO_bsp[3] <- 0.05*10
  p.EO_bsp[16:18] ~ ddirch(alpha6.EO_bsp[1:3])
  #
  alpha7.EO_bsp[1] <- 0.7*10
  alpha7.EO_bsp[2] <- 0.25*10
  alpha7.EO_bsp[3] <- 0.05*10
  p.EO_bsp[19:21] ~ ddirch(alpha7.EO_bsp[1:3])  
  #
  alpha8.EO_bsp[1] <- 0.6*10
  alpha8.EO_bsp[2] <- 0.32*10
  alpha8.EO_bsp[3] <- 0.08*10
  p.EO_bsp[22:24] ~ ddirch(alpha8.EO_bsp[1:3]) 
  #
  alpha9.EO_bsp[1] <- 0.575*10
  alpha9.EO_bsp[2] <- 0.4*10
  alpha9.EO_bsp[3] <- 0.025*10
  p.EO_bsp[25:27] ~ ddirch(alpha9.EO_bsp[1:3])
  #
  eta.EO_bsp[1] ~ dnorm(0,1)
  eta.EO_bsp[2] ~ dnorm(0,1)

  p.ventTrigg[1] <- 1
  p.ventTrigg[2] <- 0
  p.ventTrigg[3] <- 0
  #
  p.ventTrigg[4] <- 1-p.ventTrigg[6]
  p.ventTrigg[5] <- tol
  p.ventTrigg[6] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  p.ventTrigg[7] <- 1-p.ventTrigg[9]
  p.ventTrigg[8] <- tol
  p.ventTrigg[9] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  p.ventTrigg[10] <- 1-p.ventTrigg[12]
  p.ventTrigg[11] <- tol
  p.ventTrigg[12] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  alpha5.ventTrigg[1] <- 0.34*10
  alpha5.ventTrigg[2] <- 0.33*10
  alpha5.ventTrigg[3] <- 0.33*10
  p.ventTrigg[13:15] ~ ddirch(alpha5.ventTrigg[1:3])  
  #
  alpha6.ventTrigg[1] <- 0.34*10
  alpha6.ventTrigg[2] <- 0.33*10
  alpha6.ventTrigg[3] <- 0.33*10
  p.ventTrigg[16:18] ~ ddirch(alpha6.ventTrigg[1:3])
  #
  alpha7.ventTrigg[1] <- 0.34*10
  alpha7.ventTrigg[2] <- 0.33*10
  alpha7.ventTrigg[3] <- 0.33*10
  p.ventTrigg[19:21] ~ ddirch(alpha7.ventTrigg[1:3])  
  #
  eta.ventTrigg[1] ~ dnorm(0,1)
  eta.ventTrigg[2] ~ dnorm(0,1)
 
  alpha1.Edem_g[1] <- 0.95*10
  alpha1.Edem_g[2] <- 0.0075*10
  alpha1.Edem_g[3] <- 0.0425*10
  p.Edem_g[1:3] ~ ddirch(alpha1.Edem_g[1:3])
  #
  p.Edem_g[4] <- 1-p.Edem_g[5]
  p.Edem_g[5] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.Edem_g[6] <- tol
  #
  p.Edem_g[7] <- 1-p.Edem_g[8]
  p.Edem_g[8] ~ dbeta(0.95*30,0.05*30) T(tol,1-tol)
  p.Edem_g[9] <- tol  
  #
  p.Edem_g[10] <- 1-p.Edem_g[12]
  p.Edem_g[11] <- tol
  p.Edem_g[12] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  #
  p.Edem_g[13] <- 1-p.Edem_g[15]
  p.Edem_g[14] <- tol
  p.Edem_g[15] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  eta.Edem_g[1] ~ dnorm(0,1)
  eta.Edem_g[2] ~ dnorm(0,1)
  
  alpha1.V_pl[1] <- 0.94*10
  alpha1.V_pl[2] <- 0.03*10
  alpha1.V_pl[3] <- 0.03*10
  p.V_pl[1:3] ~ ddirch(alpha1.V_pl[1:3])
  #
  p.V_pl[4] <- 1-p.V_pl[5]
  p.V_pl[5] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.V_pl[6] <- tol
  #
  alpha3.V_pl[1] <- 0.15*10
  alpha3.V_pl[2] <- 0.68*10
  alpha3.V_pl[3] <- 0.17*10
  p.V_pl[7:9] ~ ddirch(alpha3.V_pl[1:3])
  #
  p.V_pl[10] <- 1-p.V_pl[12]
  p.V_pl[11] ~ dbeta(0.07*30,0.93*30) T(tol,1-tol)
  p.V_pl[12] <- tol
  #
  alpha5.V_pl[1] <- 0.25*10
  alpha5.V_pl[2] <- 0.65*10
  alpha5.V_pl[3] <- 0.1*10
  p.V_pl[13:15] ~ ddirch(alpha5.V_pl[1:3])
  #
  eta.V_pl[1] ~ dnorm(0,1)
  eta.V_pl[2] ~ dnorm(0,1)
  
  alpha1.nonSTE[1] <- 0.989*20
  alpha1.nonSTE[2] <- 0.01*20
  alpha1.nonSTE[3] <- 0.001*20
  p.nonSTE[1:3] ~ ddirch(alpha1.nonSTE[1:3])
  #
  alpha2.nonSTE[1] <- 0.2*10
  alpha2.nonSTE[2] <- 0.56*10
  alpha2.nonSTE[3] <- 0.24*10
  p.nonSTE[4:6] ~ ddirch(alpha2.nonSTE[1:3])
  #
  alpha3.nonSTE[1] <- 0.1*10
  alpha3.nonSTE[2] <- 0.72*10
  alpha3.nonSTE[3] <- 0.18*10
  p.nonSTE[7:9] ~ ddirch(alpha3.nonSTE[1:3])
  #
  alpha4.nonSTE[1] <- 0.3*10
  alpha4.nonSTE[2] <- 0.2*10
  alpha4.nonSTE[3] <- 0.5*10
  p.nonSTE[10:12] ~ ddirch(alpha4.nonSTE[1:3])
  #
  alpha5.nonSTE[1] <- 0.4*10
  alpha5.nonSTE[2] <- 0.2*10
  alpha5.nonSTE[3] <- 0.4*10
  p.nonSTE[13:15] ~ ddirch(alpha5.nonSTE[1:3])
  #
  eta.nonSTE[1] ~ dnorm(0,1)
  eta.nonSTE[2] ~ dnorm(0,1)
  
  alpha1.Q_inT[1] <- 0.99*20
  alpha1.Q_inT[2] <- 0.009*20
  alpha1.Q_inT[3] <- 0.001*20
  p.Q_inT[1:3] ~ ddirch(alpha1.Q_inT[1:3])
  #
  p.Q_inT[4] <- 1-p.Q_inT[5]
  p.Q_inT[5] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Q_inT[6] <- tol
  #
  p.Q_inT[7] <- 1-p.Q_inT[9]
  p.Q_inT[8] <- tol
  p.Q_inT[9] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  eta.Q_inT[1] ~ dnorm(0,1)
  eta.Q_inT[2] ~ dnorm(0,1)
      
  alpha1.Dispeps[1] <- 0.95*10
  alpha1.Dispeps[2] <- 0.03*10
  alpha1.Dispeps[3] <- 0.02*10
  p.Dispeps[1:3] ~ ddirch(alpha1.Dispeps[1:3])
  #
  alpha2.Dispeps[1] <- 0.1*10
  alpha2.Dispeps[2] <- 0.45*10
  alpha2.Dispeps[3] <- 0.45*10
  p.Dispeps[4:6] ~ ddirch(alpha2.Dispeps[1:3])
  #
  alpha3.Dispeps[1] <- 0.05*10
  alpha3.Dispeps[2] <- 0.095*10
  alpha3.Dispeps[3] <- 0.855*10
  p.Dispeps[7:9] ~ ddirch(alpha3.Dispeps[1:3])
  #
  alpha4.Dispeps[1] <- 0.1*10
  alpha4.Dispeps[2] <- 0.18*10
  alpha4.Dispeps[3] <- 0.72*10
  p.Dispeps[10:12] ~ ddirch(alpha4.Dispeps[1:3])
  #
  alpha5.Dispeps[1] <- 0.8*10
  alpha5.Dispeps[2] <- 0.1*10
  alpha5.Dispeps[3] <- 0.1*10
  p.Dispeps[13:15] ~ ddirch(alpha5.Dispeps[1:3])  
  #
  alpha6.Dispeps[1] <- 0.5*10
  alpha6.Dispeps[2] <- 0.45*10
  alpha6.Dispeps[3] <- 0.05*10
  p.Dispeps[16:18] ~ ddirch(alpha6.Dispeps[1:3])  
  #
  eta.Dispeps[1] ~ dnorm(0,1)
  eta.Dispeps[2] ~ dnorm(0,1)
    
  alpha1.Toss[1] <- 0.9*10
  alpha1.Toss[2] <- 0.05*10
  alpha1.Toss[3] <- 0.05*10
  p.Toss[1:3] ~ ddirch(alpha1.Toss[1:3])
  #
  alpha2.Toss[1] <- 0.5*10
  alpha2.Toss[2] <- 0.25*10
  alpha2.Toss[3] <- 0.25*10
  p.Toss[4:6] ~ ddirch(alpha2.Toss[1:3])
  #
  alpha3.Toss[1] <- 0.5*10
  alpha3.Toss[2] <- 0.25*10
  alpha3.Toss[3] <- 0.25*10
  p.Toss[7:9] ~ ddirch(alpha3.Toss[1:3])
  #
  alpha4.Toss[1] <- 0.5*10
  alpha4.Toss[2] <- 0.25*10
  alpha4.Toss[3] <- 0.25*10
  p.Toss[10:12] ~ ddirch(alpha4.Toss[1:3])
  #
  p.Toss[13] <- 1-p.Toss[14]
  p.Toss[14] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.Toss[15] <- tol
  #
  p.Toss[16] <- 1-p.Toss[17]
  p.Toss[17] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.Toss[18] <- tol
  #
  p.Toss[19] <- 1-p.Toss[20]
  p.Toss[20] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.Toss[21] <- tol
  #
  alpha8.Toss[1] <- 0.1*10
  alpha8.Toss[2] <- 0.63*10
  alpha8.Toss[3] <- 0.27*10
  p.Toss[22:24] ~ ddirch(alpha8.Toss[1:3]) 
  #
  alpha9.Toss[1] <- 0.01*10
  alpha9.Toss[2] <- 0.89*10
  alpha9.Toss[3] <- 0.1*10
  p.Toss[25:27] ~ ddirch(alpha9.Toss[1:3])
  #
  p.Toss[28] <- 1-p.Toss[30]
  p.Toss[29] <- tol
  p.Toss[30] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  eta.Toss[1] ~ dnorm(0,1)
  eta.Toss[2] ~ dnorm(0,1)
  
  p.IpossCerPth[1] <- 1-p.IpossCerPth[3]
  p.IpossCerPth[2] <- tol
  p.IpossCerPth[3] ~ dbeta(0.01*30,0.99*30) T(tol,1-tol)
  #
  alpha2.IpossCerPth[1] <- 0.4*10
  alpha2.IpossCerPth[2] <- 0.5*10
  alpha2.IpossCerPth[3] <- 0.1*10
  p.IpossCerPth[4:6] ~ ddirch(alpha2.IpossCerPth[1:3])
  #
  alpha3.IpossCerPth[1] <- 0.2*10
  alpha3.IpossCerPth[2] <- 0.65*10
  alpha3.IpossCerPth[3] <- 0.15*10
  p.IpossCerPth[7:9] ~ ddirch(alpha3.IpossCerPth[1:3])
  #
  alpha4.IpossCerPth[1] <- 0.6*10
  alpha4.IpossCerPth[2] <- 0.3*10
  alpha4.IpossCerPth[3] <- 0.1*10
  p.IpossCerPth[10:12] ~ ddirch(alpha4.IpossCerPth[1:3])
  #
  p.IpossCerPth[13] <- 1-p.IpossCerPth[14]
  p.IpossCerPth[14] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.IpossCerPth[15] <- tol  
  #
  p.IpossCerPth[16] <- 1-p.IpossCerPth[17]
  p.IpossCerPth[17] ~ dbeta(0.35*10,0.65*10) T(tol,1-tol)
  p.IpossCerPth[18] <- tol
  #
  p.IpossCerPth[19] <- tol
  p.IpossCerPth[20] <- tol
  p.IpossCerPth[21] <- 1-tol
  #
  p.IpossCerPth[22] <- 1-p.IpossCerPth[24]
  p.IpossCerPth[23] <- tol
  p.IpossCerPth[24] ~ dbeta(0.25*10,0.75*10) T(tol,1-tol)
  #
  p.IpossCerPth[25] <- 1-p.IpossCerPth[26]
  p.IpossCerPth[26] <- tol
  p.IpossCerPth[27] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  #
  p.IpossCerPth[28] <- tol
  p.IpossCerPth[29] <- 1-tol
  p.IpossCerPth[30] <- tol
  #
  eta.IpossCerPth[1] ~ dnorm(0,1)
  eta.IpossCerPth[2] ~ dnorm(0,1)
  
  alpha1.ConvPth[1] <- 0.975*20
  alpha1.ConvPth[2] <- 0.015*20
  alpha1.ConvPth[3] <- 0.01*20
  p.ConvPth[1:3] ~ ddirch(alpha1.ConvPth[1:3])
  #
  alpha2.ConvPth[1] <- 0.7*10
  alpha2.ConvPth[2] <- 0.05*10
  alpha2.ConvPth[3] <- 0.25*10
  p.ConvPth[4:6] ~ ddirch(alpha2.ConvPth[1:3])
  #
  p.ConvPth[7] <- 1-p.ConvPth[8]
  p.ConvPth[8] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.ConvPth[9] <- tol
  #
  p.ConvPth[10] <- 1-p.ConvPth[12]
  p.ConvPth[11] <- tol
  p.ConvPth[12] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  #
  p.ConvPth[13] <- 1-p.ConvPth[14]
  p.ConvPth[14] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.ConvPth[15] <- tol
  #
  eta.ConvPth[1] ~ dnorm(0,1)
  eta.ConvPth[2] ~ dnorm(0,1)
  
  alpha1.Diss_ventr[1] <- 0.992*100
  alpha1.Diss_ventr[2] <- 0.0015*100
  alpha1.Diss_ventr[3] <- 0.0015*100
  alpha1.Diss_ventr[4] <- 0.005*100
  p.Diss_ventr[1:4] ~ ddirch(alpha1.Diss_ventr[1:4])
  #
  alpha2.Diss_ventr[1] <- 0.58*10
  alpha2.Diss_ventr[2] <- 0.24*10
  alpha2.Diss_ventr[3] <- 0.12*10
  alpha2.Diss_ventr[4] <- 0.06*10
  p.Diss_ventr[5:8] ~ ddirch(alpha2.Diss_ventr[1:4])
  #
  alpha3.Diss_ventr[1] <- 0.5*10
  alpha3.Diss_ventr[2] <- 0.2*10
  alpha3.Diss_ventr[3] <- 0.25*10
  alpha3.Diss_ventr[4] <- 0.05*10
  p.Diss_ventr[9:12] ~ ddirch(alpha3.Diss_ventr[1:4])
  #
  alpha4.Diss_ventr[1] <- 0.2*10
  alpha4.Diss_ventr[2] <- 0.32*10
  alpha4.Diss_ventr[3] <- 0.28*10
  alpha4.Diss_ventr[4] <- 0.2*10
  p.Diss_ventr[13:16] ~ ddirch(alpha4.Diss_ventr[1:4])
  #
  alpha5.Diss_ventr[1] <- 0.2*10
  alpha5.Diss_ventr[2] <- 0.7*10
  alpha5.Diss_ventr[3] <- 0.05*10
  alpha5.Diss_ventr[4] <- 0.05*10
  p.Diss_ventr[17:20] ~ ddirch(alpha5.Diss_ventr[1:4])
  #
  alpha6.Diss_ventr[1] <- 0.05*10
  alpha6.Diss_ventr[2] <- 0.85*10
  alpha6.Diss_ventr[3] <- 0.05*10
  alpha6.Diss_ventr[4] <- 0.05*10
  p.Diss_ventr[21:24] ~ ddirch(alpha6.Diss_ventr[1:4])  
  #
  eta.Diss_ventr[1] ~ dnorm(0,1)
  eta.Diss_ventr[2] ~ dnorm(0,1)
  eta.Diss_ventr[3] ~ dnorm(0,1)

  p.IntraVasCoa[1] <- 1
  p.IntraVasCoa[2] <- 0
  p.IntraVasCoa[3] <- 0
  #
  p.IntraVasCoa[4] <- 1-p.IntraVasCoa[5]
  p.IntraVasCoa[5] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.IntraVasCoa[6] <- tol
  #
  p.IntraVasCoa[7] <- 1-p.IntraVasCoa[8]
  p.IntraVasCoa[8] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.IntraVasCoa[9] <- tol
  #
  p.IntraVasCoa[10] <- 1-p.IntraVasCoa[11]
  p.IntraVasCoa[11] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.IntraVasCoa[12] <- tol
  #
  p.IntraVasCoa[13] <- 1-p.IntraVasCoa[14]
  p.IntraVasCoa[14] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.IntraVasCoa[15] <- tol
  #
  p.IntraVasCoa[16] <- 1-p.IntraVasCoa[17]
  p.IntraVasCoa[17] ~ dbeta(0.06*30,0.94*30) T(tol,1-tol)
  p.IntraVasCoa[18] <- tol
  #
  p.IntraVasCoa[19] <- 1-p.IntraVasCoa[21]
  p.IntraVasCoa[20] <- tol
  p.IntraVasCoa[21] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  #
  p.IntraVasCoa[22] <- 1-p.IntraVasCoa[24]
  p.IntraVasCoa[23] <- tol
  p.IntraVasCoa[24] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  #
  eta.IntraVasCoa[1] ~ dnorm(0,1)
  eta.IntraVasCoa[2] ~ dnorm(0,1)
       
  alpha1.ArtIntraVasCoa[1] <- 0.89*20
  alpha1.ArtIntraVasCoa[2] <- 0.025*20
  alpha1.ArtIntraVasCoa[3] <- 0.025*20
  alpha1.ArtIntraVasCoa[4] <- 0.06*20
  p.ArtIntraVasCoa[1:4] ~ ddirch(alpha1.ArtIntraVasCoa[1:4])
  p.ArtIntraVasCoa[5] <- 0
  #
  p.ArtIntraVasCoa[6] <- 1-p.ArtIntraVasCoa[9]
  p.ArtIntraVasCoa[7] <- tol
  p.ArtIntraVasCoa[8] <- tol
  p.ArtIntraVasCoa[9] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.ArtIntraVasCoa[10] <- tol
  #
  alpha3.ArtIntraVasCoa[1] <- 0.3*10
  alpha3.ArtIntraVasCoa[2] <- 0.28*10
  alpha3.ArtIntraVasCoa[3] <- 0.21*10
  alpha3.ArtIntraVasCoa[4] <- 0.21*10
  dir3.ArtIntraVasCoa[1:4] ~ ddirch(alpha3.ArtIntraVasCoa[1:4])  #####
  p.ArtIntraVasCoa[11] <- dir3.ArtIntraVasCoa[1]
  p.ArtIntraVasCoa[12] <- dir3.ArtIntraVasCoa[2]
  p.ArtIntraVasCoa[13] <- dir3.ArtIntraVasCoa[3]
  p.ArtIntraVasCoa[14] <- tol
  p.ArtIntraVasCoa[15] <- dir3.ArtIntraVasCoa[4]
  #
  p.ArtIntraVasCoa[16] <- 1-p.ArtIntraVasCoa[17]
  p.ArtIntraVasCoa[17] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.ArtIntraVasCoa[18] <- tol
  p.ArtIntraVasCoa[19] <- tol
  p.ArtIntraVasCoa[20] <- tol
  #
  alpha5.ArtIntraVasCoa[1] <- 0.6*10
  alpha5.ArtIntraVasCoa[2] <- 0.3*10
  alpha5.ArtIntraVasCoa[3] <- 0.2*10
  p.ArtIntraVasCoa[21:23] ~ ddirch(alpha5.ArtIntraVasCoa[1:3])
  p.ArtIntraVasCoa[24] <- tol
  p.ArtIntraVasCoa[25] <- tol
  #
  p.ArtIntraVasCoa[26] <- 1-p.ArtIntraVasCoa[30]
  p.ArtIntraVasCoa[27] <- tol
  p.ArtIntraVasCoa[28] <- tol
  p.ArtIntraVasCoa[29] <- tol
  p.ArtIntraVasCoa[30] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  #
  eta.ArtIntraVasCoa[1] ~ dnorm(0,1)
  eta.ArtIntraVasCoa[2] ~ dnorm(0,1)
  eta.ArtIntraVasCoa[3] ~ dnorm(0,1)
  eta.ArtIntraVasCoa[4] ~ dnorm(0,1)

  p.BPCOa[1] <- 1
  p.BPCOa[2] <- 0
  p.BPCOa[3] <- 0
  p.BPCOa[4] <- 0
  #
  p.BPCOa[5] <- 1-p.BPCOa[8]
  p.BPCOa[6] <- tol
  p.BPCOa[7] <- tol
  p.BPCOa[8] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  #
  p.BPCOa[9] <- 1-p.BPCOa[10]
  p.BPCOa[10] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.BPCOa[11] <- tol
  p.BPCOa[12] <- tol
  #
  p.BPCOa[13] <- 1-p.BPCOa[14]
  p.BPCOa[14] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.BPCOa[15] <- tol
  p.BPCOa[16] <- tol
  #
  p.BPCOa[17] <- 1-p.BPCOa[18]
  p.BPCOa[18] ~ dbeta(0.1*10,0.9*10) T(tol,1-tol)
  p.BPCOa[19] <- tol
  p.BPCOa[20] <- tol
  #
  p.BPCOa[21] <- 1-p.BPCOa[22]
  p.BPCOa[22] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.BPCOa[23] <- tol
  p.BPCOa[24] <- tol
  #
  p.BPCOa[25] <- 1-p.BPCOa[26]
  p.BPCOa[26] ~ dbeta(0.35*10,0.65*10) T(tol,p.BPCOa[22])
  p.BPCOa[27] <- tol
  p.BPCOa[28] <- tol
  #
  p.BPCOa[29] <- 1-p.BPCOa[31]
  p.BPCOa[30] <- tol
  p.BPCOa[31] ~ dbeta(0.7*10,0.3*10) T(tol,1-tol)
  p.BPCOa[32] <- tol
  #
  p.BPCOa[33] <- 1-p.BPCOa[35]
  p.BPCOa[34] <- tol
  p.BPCOa[35] ~ dbeta(0.49*10,0.51*10) T(tol,p.BPCOa[31])
  p.BPCOa[36] <- tol
  #
  alpha10.BPCOa[1] <- 0.2*10
  alpha10.BPCOa[2] <- 0.64*10
  alpha10.BPCOa[3] <- 0.16*10
  p.BPCOa[37:39] ~ ddirch(alpha10.BPCOa[1:3])      
  p.BPCOa[40] <- tol
  #
  eta.BPCOa[1] ~ dnorm(0,1)
  eta.BPCOa[2] ~ dnorm(0,1)
  eta.BPCOa[3] ~ dnorm(0,1)

  alpha1.BradiTachi[1] <- 0.96*20
  alpha1.BradiTachi[2] <- 0.01*20
  alpha1.BradiTachi[3] <- 0.03*20
  dir1.BradiTachi[1:3] ~ ddirch(alpha1.BradiTachi[1:3]) #####
  p.BradiTachi[1] <- dir1.BradiTachi[1] 
  p.BradiTachi[2] <- 0
  p.BradiTachi[3] <- dir1.BradiTachi[2]
  p.BradiTachi[4] <- dir1.BradiTachi[3]
  #
  alpha2.BradiTachi[1] <- 0.7*10
  alpha2.BradiTachi[2] <- 0.1*10
  alpha2.BradiTachi[3] <- 0.015*10
  alpha2.BradiTachi[4] <- 0.185*10
  p.BradiTachi[5:8] ~ ddirch(alpha2.BradiTachi[1:4])
  #
  alpha3.BradiTachi[1] <- 0.5*10
  alpha3.BradiTachi[2] <- 0.05*10
  alpha3.BradiTachi[3] <- 0.35*10
  alpha3.BradiTachi[4] <- 0.1*10
  p.BradiTachi[9:12] ~ ddirch(alpha3.BradiTachi[1:4])
  #
  eta.BradiTachi[1] ~ dnorm(0,1)
  eta.BradiTachi[2] ~ dnorm(0,1)
  eta.BradiTachi[3] ~ dnorm(0,1)
      
  alpha1.VenIntraVasCoa[1] <- 0.88*20
  alpha1.VenIntraVasCoa[2] <- 0.048*20
  alpha1.VenIntraVasCoa[3] <- 0.06*20
  alpha1.VenIntraVasCoa[4] <- 0.012*20
  p.VenIntraVasCoa[1:4] ~ ddirch(alpha1.VenIntraVasCoa[1:4])
  #
  alpha2.VenIntraVasCoa[1] <- 0.9*10
  alpha2.VenIntraVasCoa[2] <- 0.01*10
  alpha2.VenIntraVasCoa[3] <- 0.08*10
  alpha2.VenIntraVasCoa[4] <- 0.01*10
  p.VenIntraVasCoa[5:8] ~ ddirch(alpha2.VenIntraVasCoa[1:4])
  #
  alpha3.VenIntraVasCoa[1] <- 0.66*10
  alpha3.VenIntraVasCoa[2] <- 0.035*10
  alpha3.VenIntraVasCoa[3] <- 0.045*10
  alpha3.VenIntraVasCoa[4] <- 0.010*10
  p.VenIntraVasCoa[9:12] ~ ddirch(alpha3.VenIntraVasCoa[1:4])
  #
  alpha4.VenIntraVasCoa[1] <- 0.3*10
  alpha4.VenIntraVasCoa[2] <- 0.1*10
  alpha4.VenIntraVasCoa[3] <- 0.5*10
  alpha4.VenIntraVasCoa[4] <- 0.1*10
  p.VenIntraVasCoa[13:16] ~ ddirch(alpha4.VenIntraVasCoa[1:4])
  #
  alpha5.VenIntraVasCoa[1] <- 0.9*10
  alpha5.VenIntraVasCoa[2] <- 0.01*10
  alpha5.VenIntraVasCoa[3] <- 0.08*10
  alpha5.VenIntraVasCoa[4] <- 0.01*10
  p.VenIntraVasCoa[17:20] ~ ddirch(alpha5.VenIntraVasCoa[1:4])
  #
  alpha6.VenIntraVasCoa[1] <- 0.2*10
  alpha6.VenIntraVasCoa[2] <- 0.1*10
  alpha6.VenIntraVasCoa[3] <- 0.7*10
  p.VenIntraVasCoa[21:23] ~ ddirch(alpha6.VenIntraVasCoa[1:3])
  p.VenIntraVasCoa[24] <- tol
  #
  p.VenIntraVasCoa[25] <- 1-p.VenIntraVasCoa[26]
  p.VenIntraVasCoa[26] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  p.VenIntraVasCoa[27] <- tol
  p.VenIntraVasCoa[28] <- tol
  #
  p.VenIntraVasCoa[29] <- 1-p.VenIntraVasCoa[31]
  p.VenIntraVasCoa[30] <- tol
  p.VenIntraVasCoa[31] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.VenIntraVasCoa[32] <- tol
  #
  p.VenIntraVasCoa[33] <- 1-p.VenIntraVasCoa[35]
  p.VenIntraVasCoa[34] <- tol
  p.VenIntraVasCoa[35] ~ dbeta(0.3*10,0.7*10) T(tol,1-tol)
  p.VenIntraVasCoa[36] <- tol
  #
  p.VenIntraVasCoa[37] <- 1-p.VenIntraVasCoa[39]
  p.VenIntraVasCoa[38] <- tol
  p.VenIntraVasCoa[39] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  p.VenIntraVasCoa[40] <- tol
  #
  p.VenIntraVasCoa[41] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  #
  eta.VenIntraVasCoa[1] ~ dnorm(0,1)
  eta.VenIntraVasCoa[2] ~ dnorm(0,1)
  eta.VenIntraVasCoa[3] ~ dnorm(0,1)

  alpha1.autNerSys[1] <- 0.85*20
  alpha1.autNerSys[2] <- 0.03*20
  alpha1.autNerSys[3] <- 0.015*20
  alpha1.autNerSys[4] <- 0.0375*20
  alpha1.autNerSys[5] <- 0.0375*20
  alpha1.autNerSys[6] <- 0.03*20
  p.autNerSys[1:6] ~ ddirch(alpha1.autNerSys[1:6])
  #
  alpha2.autNerSys[1] <- 0.2*10
  alpha2.autNerSys[2] <- 0.36*10
  alpha2.autNerSys[3] <- 0.36*10
  alpha2.autNerSys[4] <- 0.08*10
  p.autNerSys[7:10] ~ ddirch(alpha2.autNerSys[1:4])
  p.autNerSys[11] <- tol
  p.autNerSys[12] <- tol
  #
  alpha3.autNerSys[1] <- 0.32*10
  alpha3.autNerSys[2] <- 0.315*10
  alpha3.autNerSys[3] <- 0.315*10
  alpha3.autNerSys[4] <- 0.05*10
  p.autNerSys[13:16] ~ ddirch(alpha3.autNerSys[1:4])
  p.autNerSys[17] <- tol
  p.autNerSys[18] <- tol
  #
  alpha4.autNerSys[1] <- 0.1*10
  alpha4.autNerSys[2] <- 0.8*10
  alpha4.autNerSys[3] <- 0.1*10
  p.autNerSys[19:21] ~ ddirch(alpha4.autNerSys[1:3])
  p.autNerSys[22] <- tol
  p.autNerSys[23] <- tol
  p.autNerSys[24] <- tol
  #
  alpha5.autNerSys[1] <- 0.02*10
  alpha5.autNerSys[2] <- 0.7*10
  alpha5.autNerSys[3] <- 0.28*10
  p.autNerSys[25:27] ~ ddirch(alpha5.autNerSys[1:3])
  p.autNerSys[28] <- tol
  p.autNerSys[29] <- tol
  p.autNerSys[30] <- tol
  #
  p.autNerSys[31] <- 1-p.autNerSys[34]
  p.autNerSys[32] <- tol
  p.autNerSys[33] <- tol
  p.autNerSys[34] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  p.autNerSys[35] <- tol
  p.autNerSys[36] <- tol
  #
  alpha7.autNerSys[1] <- 0.6*10
  alpha7.autNerSys[2] <- 0.32*10 
  alpha7.autNerSys[3] <- 0.08*10
  dir7.autNerSys[1:3] ~ ddirch(alpha7.autNerSys[1:3])  #####
  p.autNerSys[37] <- dir7.autNerSys[1]
  p.autNerSys[38] <- tol
  p.autNerSys[39] <- dir7.autNerSys[2]
  p.autNerSys[40] <- dir7.autNerSys[3]
  p.autNerSys[41] <- tol
  p.autNerSys[42] <- tol
  #
  p.autNerSys[43] <- 1-p.autNerSys[44]
  p.autNerSys[44] ~ dbeta(0.99*30,0.01*30) T(tol,1-tol)
  p.autNerSys[45] <- tol
  p.autNerSys[46] <- tol
  p.autNerSys[47] <- tol
  p.autNerSys[48] <- tol
  #
  p.autNerSys[49] <- 1-p.autNerSys[53]
  p.autNerSys[50] <- tol
  p.autNerSys[51] <- tol
  p.autNerSys[52] <- tol
  p.autNerSys[53] ~ dbeta(0.75*10,0.25*10) T(tol,1-tol)
  p.autNerSys[54] <- tol
  #
  p.autNerSys[55] <- 1-p.autNerSys[60]
  p.autNerSys[56] <- tol
  p.autNerSys[57] <- tol
  p.autNerSys[58] <- tol
  p.autNerSys[59] <- tol
  p.autNerSys[60] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  #
  eta.autNerSys[1] ~ dnorm(0,1)
  eta.autNerSys[2] ~ dnorm(0,1)
  eta.autNerSys[3] ~ dnorm(0,1)
  eta.autNerSys[4] ~ dnorm(0,1)
  eta.autNerSys[5] ~ dnorm(0,1)

  alpha1.Dol_type[1] <- 0.9*20
  alpha1.Dol_type[2] <- 0.05*20
  alpha1.Dol_type[3] <- 0.025*20
  alpha1.Dol_type[4] <- 0.01*20
  alpha1.Dol_type[5] <- 0.01*20
  alpha1.Dol_type[6] <- 0.005*20
  p.Dol_type[1:6] ~ ddirch(alpha1.Dol_type[1:6])
  #
  alpha2.Dol_type[1] <- 0.02*10
  alpha2.Dol_type[2] <- 0.06*10
  alpha2.Dol_type[3] <- 0.26*10
  alpha2.Dol_type[4] <- 0.11*10
  alpha2.Dol_type[5] <- 0.44*10
  alpha2.Dol_type[6] <- 0.11*10
  p.Dol_type[7:12] ~ ddirch(alpha2.Dol_type[1:6])
  #
  p.Dol_type[13] <- tol
  p.Dol_type[14] <- 1-tol
  p.Dol_type[15] <- tol
  p.Dol_type[16] <- tol
  p.Dol_type[17] <- tol
  p.Dol_type[18] <- tol
  #
  p.Dol_type[19] <- tol
  p.Dol_type[20] <- tol
  p.Dol_type[21] <- 1-tol
  p.Dol_type[22] <- tol
  p.Dol_type[23] <- tol
  p.Dol_type[24] <- tol
  #
  p.Dol_type[25] <- tol
  p.Dol_type[26] <- tol
  p.Dol_type[27] <- tol
  p.Dol_type[28] <- 1-tol
  p.Dol_type[29] <- tol
  p.Dol_type[30] <- tol
  #
  p.Dol_type[31] <- tol
  p.Dol_type[32] <- tol
  p.Dol_type[33] <- tol
  p.Dol_type[34] <- tol
  p.Dol_type[35] <- 1-tol
  p.Dol_type[36] <- tol
  #
  p.Dol_type[37] <- tol
  p.Dol_type[38] <- tol
  p.Dol_type[39] <- tol
  p.Dol_type[40] <- tol
  p.Dol_type[41] <- tol
  p.Dol_type[42] <- 1-tol
  #
  eta.Dol_type[1] ~ dnorm(0,1)
  eta.Dol_type[2] ~ dnorm(0,1)
  eta.Dol_type[3] ~ dnorm(0,1)
  eta.Dol_type[4] ~ dnorm(0,1)
  eta.Dol_type[5] ~ dnorm(0,1)
    
  alpha1.Tiroid[1] <- 0.965*10
  alpha1.Tiroid[2] <- 0.0133*10
  alpha1.Tiroid[3] <- 0.0007*10
  alpha1.Tiroid[4] <- 0.02*10
  alpha1.Tiroid[5] <- 0.001*10
  p.Tiroid[1:5] ~ ddirch(alpha1.Tiroid[1:5])
  #
  alpha2.Tiroid[1] <- 0.935*10
  alpha2.Tiroid[2] <- 0.0247*10
  alpha2.Tiroid[3] <- 0.0013*10
  alpha2.Tiroid[4] <- 0.037*10
  alpha2.Tiroid[5] <- 0.002*10
  p.Tiroid[6:10] ~ ddirch(alpha2.Tiroid[1:5])
  
  alpha1.Card_obstr[1] <- 0.99*20
  alpha1.Card_obstr[2] <- 0.005*20
  alpha1.Card_obstr[3] <- 0.005*20
  p.Card_obstr[1:3] ~ ddirch(alpha1.Card_obstr[1:3])
  #
  alpha2.Card_obstr[1] <- 0.1*10
  alpha2.Card_obstr[2] <- 0.45*10
  alpha2.Card_obstr[3] <- 0.45*10
  p.Card_obstr[4:6] ~ ddirch(alpha2.Card_obstr[1:3])
  #
  alpha3.Card_obstr[1] <- 0.9*10
  alpha3.Card_obstr[2] <- 0.05*10
  alpha3.Card_obstr[3] <- 0.05*10
  p.Card_obstr[7:9] ~ ddirch(alpha3.Card_obstr[1:3])
  #
  eta.Card_obstr[1] ~ dnorm(0,1)
  eta.Card_obstr[2] ~ dnorm(0,1)

  alpha1.Pnx[1] <- 0.956*10
  alpha1.Pnx[2] <- 0.036*10
  alpha1.Pnx[3] <- 0.008*10
  p.Pnx[1:3] ~ ddirch(alpha1.Pnx[1:3])
  #
  alpha2.Pnx[1] <- 0.9*10
  alpha2.Pnx[2] <- 0.05*10
  alpha2.Pnx[3] <- 0.05*10
  p.Pnx[4:6] ~ ddirch(alpha2.Pnx[1:3])
  #
  alpha3.Pnx[1] <- 0.9*10
  alpha3.Pnx[2] <- 0.05*10
  alpha3.Pnx[3] <- 0.05*10
  p.Pnx[7:9] ~ ddirch(alpha3.Pnx[1:3])
  #
  eta.Pnx[1] ~ dnorm(0,1)
  eta.Pnx[2] ~ dnorm(0,1)

  alpha1.Insuf_aorCron[1] <- 0.88*90
  alpha1.Insuf_aorCron[2] <- 0.11*90
  alpha1.Insuf_aorCron[3] <- 0.01*90
  p.Insuf_aorCron[1:3] ~ ddirch(alpha1.Insuf_aorCron[1:3])
  #
  p.Insuf_aorCron[4] <- 1-p.Insuf_aorCron[6]
  p.Insuf_aorCron[5] <- tol
  p.Insuf_aorCron[6] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  #
  alpha3.Insuf_aorCron[1] <- 0.35*10
  alpha3.Insuf_aorCron[2] <- 0.15*10
  alpha3.Insuf_aorCron[3] <- 0.5*10
  p.Insuf_aorCron[7:9] ~ ddirch(alpha3.Insuf_aorCron[1:3])
  #
  eta.Insuf_aorCron[1] ~ dnorm(0,1)
  eta.Insuf_aorCron[2] ~ dnorm(0,1)

  alpha1.InsuffMitrCron[1] <- 0.93*100
  alpha1.InsuffMitrCron[2] <- 0.02*100
  alpha1.InsuffMitrCron[3] <- 0.05*100
  p.InsuffMitrCron[1:3] ~ ddirch(alpha1.InsuffMitrCron[1:3])
  #
  p.InsuffMitrCron[4] <- 1-p.InsuffMitrCron[6]
  p.InsuffMitrCron[5] <- tol
  p.InsuffMitrCron[6] ~ dbeta(0.4*10,0.6*10) T(tol,1-tol)
  #
  p.InsuffMitrCron[7] <- 1-p.InsuffMitrCron[8]
  p.InsuffMitrCron[8] ~ dbeta(0.15*10,0.85*10) T(tol,1-tol)
  p.InsuffMitrCron[9] <- tol
  #
  alpha4.InsuffMitrCron[1] <- 0.7*10
  alpha4.InsuffMitrCron[2] <- 0.2*10
  alpha4.InsuffMitrCron[3] <- 0.1*10
  p.InsuffMitrCron[10:12] ~ ddirch(alpha4.InsuffMitrCron[1:3])
  #
  p.InsuffMitrCron[13] <- 1-p.InsuffMitrCron[14]
  p.InsuffMitrCron[14] ~ dbeta(0.9*10,0.1*10) T(tol,1-tol)
  p.InsuffMitrCron[15] <- tol
  #
  alpha6.InsuffMitrCron[1] <- 0.3*10
  alpha6.InsuffMitrCron[2] <- 0.15*10
  alpha6.InsuffMitrCron[3] <- 0.55*10
  p.InsuffMitrCron[16:18] ~ ddirch(alpha6.InsuffMitrCron[1:3])
  #
  eta.InsuffMitrCron[1] ~ dnorm(0,1)
  eta.InsuffMitrCron[2] ~ dnorm(0,1)
    
  alpha1.Enf[1] <- 0.82*20
  alpha1.Enf[2] <- 0.12*20
  alpha1.Enf[3] <- 0.06*20
  p.Enf[1:3] ~ ddirch(alpha1.Enf[1:3])
  #
  alpha2.Enf[1] <- 0.94*10
  alpha2.Enf[2] <- 0.04*10
  alpha2.Enf[3] <- 0.02*10
  p.Enf[4:6] ~ ddirch(alpha2.Enf[1:3])
  #
  alpha3.Enf[1] <- 0.65*10
  alpha3.Enf[2] <- 0.21*10
  alpha3.Enf[3] <- 0.14*10
  p.Enf[7:9] ~ ddirch(alpha3.Enf[1:3])
  #
  eta.Enf[1] <- 0
  eta.Enf[2] <- 0
  eta.Enf[3] <- 0
  eta.Enf[4] <- 0
              
  p.Opac[1] <- 1-p.Opac[2]
  p.Opac[2] ~ dbeta(0.08*30,0.92*30) T(tol,1-tol)
  p.Opac[3] <- tol      
  #
  alpha2.Opac[1] <- 0.8*10
  alpha2.Opac[2] <- 0.15*10
  alpha2.Opac[3] <- 0.05*10
  p.Opac[4:6] ~ ddirch(alpha2.Opac[1:3])
  #
  p.Opac[7] <- 1-p.Opac[9]
  p.Opac[8] <- tol
  p.Opac[9] ~ dbeta(0.2*10,0.8*10) T(tol,1-tol)
  #
  p.Opac[10] <- tol
  p.Opac[11] <- 1-tol
  p.Opac[12] <- tol
  #
  alpha5.Opac[1] <- 0.5*10
  alpha5.Opac[2] <- 0.25*10
  alpha5.Opac[3] <- 0.25*10
  p.Opac[13:15] ~ ddirch(alpha5.Opac[1:3])
  #
  alpha6.Opac[1] <- 0.11*10
  alpha6.Opac[2] <- 0.79*10
  alpha6.Opac[3] <- 0.2*10
  p.Opac[16:18] ~ ddirch(alpha6.Opac[1:3]) 
  #
  p.Opac[19] <- 1-p.Opac[20]
  p.Opac[20] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Opac[21] <- tol
  #
  p.Opac[22] <- tol
  p.Opac[23] <- 1-tol
  p.Opac[24] <- tol
  #
  eta.Opac[1] ~ dnorm(0,1)
  eta.Opac[2] ~ dnorm(0,1)

  p.Cardiomio[1] <- 1-p.Cardiomio[2]
  p.Cardiomio[2] ~ dbeta(0.053*30,0.947*30) T(tol,1-tol)
  p.Cardiomio[3] <- tol
  #
  p.Cardiomio[4] <- 1-p.Cardiomio[5]
  p.Cardiomio[5] ~ dbeta(0.04*30,0.96*30) T(tol,1-tol)
  p.Cardiomio[6] <- tol
  #
  p.Cardiomio[7] <- 1-p.Cardiomio[9]
  p.Cardiomio[8] <- tol
  p.Cardiomio[9] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  #
  p.Cardiomio[10] <- 1-p.Cardiomio[11]
  p.Cardiomio[11] ~ dbeta(0.8*10,0.2*10) T(tol,1-tol)
  p.Cardiomio[12] <- tol
  #
  p.Cardiomio[13] <- 1-p.Cardiomio[15]
  p.Cardiomio[14] <- tol
  p.Cardiomio[15] ~ dbeta(0.5*10,0.5*10) T(tol,1-tol)
  #
  p.Cardiomio[16] <- tol
  p.Cardiomio[17] <- tol
  p.Cardiomio[18] <- 1-tol
  #
  eta.Cardiomio[1] ~ dnorm(0,1)
  eta.Cardiomio[2] ~ dnorm(0,1)
  eta.Cardiomio[3] ~ dnorm(0,1)
  eta.Cardiomio[4] ~ dnorm(0,1)
  
  alpha1.Polm[1] <- 0.985*20
  alpha1.Polm[2] <- 0.01*20
  alpha1.Polm[3] <- 0.005*20
  p.Polm[1:3] ~ ddirch(alpha1.Polm[1:3])
  #
  alpha2.Polm[1] <- 0.985*10
  alpha2.Polm[2] <- 0.01*10
  alpha2.Polm[3] <- 0.005*10
  p.Polm[4:6] ~ ddirch(alpha2.Polm[1:3])
  #
  alpha3.Polm[1] <- 0.2*10
  alpha3.Polm[2] <- 0.64*10
  alpha3.Polm[3] <- 0.16*10
  p.Polm[7:9] ~ ddirch(alpha3.Polm[1:3])
  #
  alpha4.Polm[1] <- 0.12*10
  alpha4.Polm[2] <- 0.79*10
  alpha4.Polm[3] <- 0.19*10
  p.Polm[10:12] ~ ddirch(alpha4.Polm[1:3])
  #
  alpha5.Polm[1] <- 0.6*10
  alpha5.Polm[2] <- 0.32*10
  alpha5.Polm[3] <- 0.08*10
  p.Polm[13:15] ~ ddirch(alpha5.Polm[1:3])
  #
  alpha6.Polm[1] <- 0.85*10
  alpha6.Polm[2] <- 0.05*10
  alpha6.Polm[3] <- 0.1*10
  p.Polm[16:18] ~ ddirch(alpha6.Polm[1:3])  
  #
  alpha7.Polm[1] <- 0.4*10
  alpha7.Polm[2] <- 0.48*10
  alpha7.Polm[3] <- 0.12*10
  p.Polm[19:21] ~ ddirch(alpha7.Polm[1:3])
  #
  alpha8.Polm[1] <- 0.8*10
  alpha8.Polm[2] <- 0.06*10
  alpha8.Polm[3] <- 0.14*10
  p.Polm[22:24] ~ ddirch(alpha8.Polm[1:3])
  #
  alpha9.Polm[1] <- 0.6*10
  alpha9.Polm[2] <- 0.2*10
  alpha9.Polm[3] <- 0.2*10
  p.Polm[25:27] ~ ddirch(alpha9.Polm[1:3]) 
  #
  eta.Polm[1] ~ dnorm(0,1)
  eta.Polm[2] ~ dnorm(0,1)
  eta.Polm[3] ~ dnorm(0,1)
  eta.Polm[4] ~ dnorm(0,1)
        
  }
