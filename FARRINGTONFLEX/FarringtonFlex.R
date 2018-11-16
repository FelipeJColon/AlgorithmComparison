rm(list=ls(all=TRUE))

# FUNCTIONS THAT PRODUCE THE DATA
# DEFINING FUNCTION h

#==============
# 5-day systems
#==============

h1=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2){
     t=1:N
     if(k==0 & k2==0){h1=alpha+beta*t}
     else{
          if(k==0)
          {
               l=1:k2
               h1=rep(0,N)
               for(i in 1:N){
                    h1[i]=alpha+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
               }
          }
          else{
               j=1:k
               l=1:k2
               h1=rep(0,N)
               for(i in 1:N){
                    h1[i]=alpha+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*5))+gama2*sin((2*pi*j*(t[i]+shift2))/(52*5)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/5)+gama4*sin((2*pi*l*(t[i]+shift))/5))
               }
          }
     }
     h1
}

negbinNoise1=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,phi,shift,shift2){
     mu <- exp(h1(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2))
     if(phi==1){yi <- rpois(N,mu)}
     else{
          prob <- 1/phi 
          size <- mu/(phi-1) 
          yi <- rnbinom(N,size=size,prob=prob)
     }
     yi
}


outbreak5=function(currentday,weeklength,wtime,yi,interval,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2,phi,numoutbk,peakoutbk,meanlog,sdlog){
     # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1
     
     N=length(yi)
     t=1:N
     
     mu <- exp(h1(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,shift2))
     s=sqrt(mu*phi)
     
     #wtime = (currentday-49*5+1):currentday         # current outbreaks
     
     # GENERATING OUTBREAKS
     
     # STARTING TIMES OF OUTBREAKS
     
     startoutbk <- sample(wtime, numoutbk, replace = FALSE)
     
     # OUTBREAK SIZE OF CASES
     
     sizeoutbk=rep(0,numoutbk)
     for(i in 1:numoutbk){
          set.seed(i)
          soutbk=1
          sou=1
          while(soutbk<2){
               set.seed(sou)
               soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
               sou=sou+1
          }
          sizeoutbk[i]=soutbk
     }
     
     # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL
     
     outbreak=rep(0,2*N)
     for( j in 1:numoutbk){
          set.seed(j)
          outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
          #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog) 
          #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
          h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
          cases <- h$counts
          weight=rep(0,length(cases))
          duration<-startoutbk:(startoutbk+length(cases)-1)
          dayofweek<-duration%%5 # 0 is friday; 1 is monday; 2 is tuesday etc.
          for(i in 1:length(cases)){
               if(dayofweek[i]==0){weight[i]=1.1}
               if(dayofweek[i]==1){weight[i]=1.5}
               if(dayofweek[i]==2){weight[i]=1.1}
               if(dayofweek[i]==3){weight[i]=1}
               if(dayofweek[i]==4){weight[i]=1}
          }
          cases2 <- cases*weight
          for (l in 1:(length(cases2))){
               outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
          }# l loop 
     }# j loop
     
     #for(v in 1:(currentday-49*5)){if(outbreak[v]>0){outbreak[v]=0}}
     for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
     outbreak=outbreak[1:N]
     
     # ADD NOISE AND OUTBREAKS 
     
     yitot=yi+outbreak
     result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
     #return(result)
}

#==============
# 7-day systems
#==============

h2=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift){
     t=1:N
     if(k==0 & k2==0){h2=alpha+beta*t}
     else{
          if(k==0)
          {
               l=1:k2
               h2=rep(0,N)
               for(i in 1:N){
                    h2[i]=alpha+beta*(t[i]+shift)+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
               }
          }
          else{
               j=1:k
               l=1:k2
               h2=rep(0,N)
               for(i in 1:N){
                    h2[i]=alpha+beta*(t[i]+shift)+sum(gama1*cos((2*pi*j*(t[i]+shift))/(52*7))+gama2*sin((2*pi*j*(t[i]+shift))/(52*7)))+sum(gama3*cos((2*pi*l*(t[i]+shift))/7)+gama4*sin((2*pi*l*(t[i]+shift))/7))
               }
          }
     }
     h2
}


negbinNoise2=function(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,phi,shift){
     mu <- exp(h2(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift))
     if(phi==1){yi <- rpois(N,mu)}
     else{
          prob <- 1/phi 
          size <- mu/(phi-1) 
          yi <- rnbinom(N,size=size,prob=prob)
     }
     yi
}


outbreak7=function(currentday,weeklength,wtime,yi,interval,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift,phi,numoutbk,peakoutbk,meanlog,sdlog){
     # theta, beta, gama1 and gama2 are the parameters of the equation for mu in Section 3.1
     
     N=length(yi)
     t=1:N
     
     mu <- exp(h2(N,k,k2,alpha,beta,gama1,gama2,gama3,gama4,shift))
     s=sqrt(mu*phi)
     
     
     #wtime = (currentday-49*7+1):currentday         # current outbreaks
     # wtime = 350*1:7         # current outbreaks
     
     # GENERATING OUTBREAKS
     
     # STARTING TIMES OF OUTBREAKS
     
     startoutbk <- sample(wtime, numoutbk, replace = FALSE)
     
     # OUTBREAK SIZE OF CASES
     
     sizeoutbk=rep(0,numoutbk)
     for(i in 1:numoutbk){
          set.seed(i)
          soutbk=1
          sou=1
          while(soutbk<2){
               set.seed(sou)
               soutbk=rpois(1,s[startoutbk[i]]*peakoutbk)
               sou=sou+1
          }
          sizeoutbk[i]=soutbk
     }
     
     # DISTRIBUTE THESE CASES OVER TIME USING LOGNORMAL
     
     outbreak=rep(0,2*N)
     for( j in 1:numoutbk){
          set.seed(j)
          outbk <-rlnorm(sizeoutbk[j], meanlog = meanlog, sdlog = sdlog) 
          #outbk <-rnorm(sizeoutbk[j], mean = meanlog2, sd = sdlog) 
          #h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),1),plot=FALSE)
          h<- hist(outbk,breaks=seq(0,ceiling(max(outbk)),interval),plot=FALSE)
          cases <- h$counts
          weight=rep(0,length(cases))
          duration<-startoutbk:(startoutbk+length(cases)-1)
          dayofweek<-duration%%7 # 0 is sunday; 1 is monday; 2 is tuesday etc.
          for(i in 1:length(cases)){
               if(dayofweek[i]==0){weight[i]=2}
               if(dayofweek[i]==1){weight[i]=1}
               if(dayofweek[i]==2){weight[i]=1}
               if(dayofweek[i]==3){weight[i]=1}
               if(dayofweek[i]==4){weight[i]=1}
               if(dayofweek[i]==5){weight[i]=1}
               if(dayofweek[i]==6){weight[i]=2}
          }
          cases2 <- cases*weight
          for (l in 1:(length(cases2))){
               outbreak[startoutbk[j]+(l-1)]= cases2[l]+outbreak[startoutbk[j]+(l-1)]
          }# l loop 
     }# j loop
     
     #for(v in (currentday-49*7):currentday){if(outbreak[v]>0){outbreak[v]=0}}
     for(v in currentday:(currentday+100)){if(outbreak[v]>0){outbreak[v]=0}}
     outbreak=outbreak[1:N]
     
     # ADD NOISE AND OUTBREAKS 
     
     yitot=yi+outbreak
     result=list(yitot=yitot,outbreak=outbreak,startoutbk=startoutbk,sizeoutbk=sizeoutbk,sd=s,mean=mu)
     #return(result)
}

#==========================
# Specify the bank holidays
#==========================

years=7
bankholidays=read.csv("Bankholidays.csv",sep=',')

bankholidays=read.csv("Bankholidays.csv",sep=',')
#fix(bankholidays)
bankhols7=bankholidays$bankhol
bankhols7=as.numeric(bankhols7)
length(bankhols7)
#fix(bankhols7)

bankhols5=bankhols7[-seq(6,length(bankhols7),7)] 
bankhols5=bankhols5[-seq(6,length(bankhols5),6)]
bankhols5=as.numeric(bankhols5)
length(bankhols5)
#fix(bankhols5)

#=======================
# Define the data frames
#=======================

nsim=100

simulateddata1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulateddata17=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))

simulatedtotals1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedtotals17=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))

simulatedoutbreak1=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak2=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak3=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak4=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak5=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak8=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak9=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak10=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak11=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak12=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak13=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak14=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak15=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedoutbreak17=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))


simulatedzseasoutbreak6=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedzseasoutbreak7=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))
simulatedzseasoutbreak16=data.frame(array(rep(0,nsim*52*7*years),dim=c(52*7*years,nsim)))

#################################
#SIMULATE SYNDROMES AND OUTBREAKS
#################################

########################################
# SIMULATIONS WITH VERY SMALL(w=2), SMALL(w=3), MEDIUM(w=5) AND LARGE(w=10) OUTBREAKS
########################################

# VERY SMALL w=2

w=2

#=====================
# 5-day week syndromes
#=====================


days5=5
N=52*days5*years

#sigid6
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50)/10
     #mu=exp(h1(N=N,k=1,k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,shift=-50,shift2=-50))
     
     out1=rep(0,N)
     for(j in 1:years){
          set.seed(j+years*i)
          out=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=((1+(j-1)*days5*52):(20+(j-1)*days5*52)),yi=yt,interval=0.02,k=1,
                        k2=1,alpha=6,beta=0,gama1=0.3,gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days5*80,meanlog=0,sdlog=0.5)
          out1=out1+out$outbreak
     }
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=6,beta=0,gama1=0.3,
                    gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     zseasoutbreak=out2$outbreak +out1
     
     zt=yt +out1 
     zitot=yt + out2$outbreak +out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          yt=append(yt,zeros,after=2*(s-1)+weekend[s])
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
          zseasoutbreak=append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata6[,i]=round(zt)
     simulatedtotals6[,i]=round(zitot)
     simulatedoutbreak6[,i]=round(zoutbreak)
     simulatedzseasoutbreak6[,i]=round(zseasoutbreak)
}

#----------------------------------------------------
# Plot the datasets and outbreaks using the following
#----------------------------------------------------

#plot(1:N,yt,typ='l')
#plot(1:(52*years*7),zt,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(365,728))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(729,1092))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1093,1456))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1457,1820))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1821,2184))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(2185,2548))
#lines(1:(52*years*7),zoutbreak,col='green')

plot(1:(52*years*7),simulatedtotals6[,4],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak6[,4],col='green')
lines(1:(52*years*7),simulatedoutbreak6[,4],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak6[,4],col='green',typ='l')

#sigid7
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=1,beta=0,gama1=0.1,gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50)
     #mu=exp(h1(N=N,k=1,k2=1,alpha=1.5,beta=0,gama1=0.1,gama2=2,gama3=0.1,gama4=0.1,shift=-50,shift2=-50))
     
     out1=rep(0,N)
     for(j in 1:years){
          set.seed(j+years*i)
          out=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=((1+(j-1)*days5*52):(20+(j-1)*days5*52)),yi=yt,interval=0.02,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                        gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days5*50,meanlog=0,sdlog=0.5)
          out1=out1+out$outbreak
     }
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                    gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     zseasoutbreak=out2$outbreak+out1
     
     zt=yt +out1 
     zitot=yt + out2$outbreak  +out1 
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          yt=append(yt,zeros,after=2*(s-1)+weekend[s])
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
          zseasoutbreak=append(zseasoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata7[,i]=round(zt)
     simulatedtotals7[,i]=round(zitot)
     simulatedoutbreak7[,i]=round(zoutbreak)
     simulatedzseasoutbreak7[,i]=round(zseasoutbreak)
}

plot(1:(52*years*7),simulatedtotals7[,7],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak7[,7],col='green')
lines(1:(52*years*7),simulatedoutbreak7[,7],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak7[,4],col='green',typ='l')


#sigid8
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=0,k2=1,alpha=6,beta=0.0001,gama1=0,gama2=0,gama3=0.6,gama4=0.9,phi=1.5,shift=0,shift2=0)/10
     #mu=exp(h1(N=N,k=0,k2=1,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.6,gama4=0.9,shift=0,shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=0,k2=1,alpha=6,beta=0,gama1=0,
                    gama2=0,gama3=0.6,gama4=0.9,phi=1.5,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata8[,i]=round(zt)
     simulatedtotals8[,i]=round(zitot)
     simulatedoutbreak8[,i]=round(zoutbreak)
}


plot(1:(52*years*7),simulateddata8[,1],typ='l',xlim=c(2185,2548),col='green')
lines(1:(52*years*7),simulateddata8[,1]+simulatedoutbreak8[,1])



#sigid9
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=1.5,gama2=0.1,gama3=0.2,gama4=0.3,phi=1,shift=-150,shift2=-150)
     mu=exp(h1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=1.5,gama2=0.1,gama3=0.6,gama4=0.8,shift=-150,shift2=-150))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),interval=0.25,yi=yt,k=1,k2=1,alpha=3,beta=0,gama1=1.5,
                    gama2=0.1,gama3=0.2,gama4=0.3,phi=1,shift=-150,shift2=-150,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata9[,i]=round(zt)
     simulatedtotals9[,i]=round(zitot)
     simulatedoutbreak9[,i]=round(zoutbreak)
}


plot(1:(52*years*7),simulateddata9[,1],typ='l',xlim=c(2185,2548),col='green')
lines(1:(52*years*7),simulateddata9[,1]+simulatedoutbreak9[,1])


#sigid10
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.15,phi=1,shift=-200,shift2=-200)
     #mu=exp(h1(N=N,k=1,k2=1,alpha=3,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.15,shift=-200,shift2=-200))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=3,beta=0,gama1=0.2,
                    gama2=0.1,gama3=0.05,gama4=0.15,phi=1,shift=-200,shift2=-200,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata10[,i]=round(zt)
     simulatedtotals10[,i]=round(zitot)
     simulatedoutbreak10[,i]=round(zoutbreak)
}



#sigid11
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=5,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.1,phi=1,shift=0,shift2=0)
     mu=exp(h1(N=N,k=1,k2=1,alpha=5,beta=0,gama1=0.2,gama2=0.1,gama3=0.05,gama4=0.1,shift=0,shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),interval=0.25,yi=yt,k=1,k2=1,alpha=5,beta=0,gama1=0.2,
                    gama2=0.1,gama3=0.05,gama4=0.1,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata11[,i]=round(zt)
     simulatedtotals11[,i]=round(zitot)
     simulatedoutbreak11[,i]=round(zoutbreak)
}


#sigid12
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,gama2=0,gama3=0.05,gama4=0.15,phi=1,shift=0,shift2=0)
     #mu=exp(h1(N=N,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,gama2=0,gama3=0.05,gama4=0.15,shift=0,shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=2,k2=1,alpha=0.5,beta=0,gama1=0.4,
                    gama2=0,gama3=0.05,gama4=0.15,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata12[,i]=round(zt)
     simulatedtotals12[,i]=round(zitot)
     simulatedoutbreak12[,i]=round(zoutbreak)
}


#sigid13
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise1(N=N,k=1,k2=1,alpha=9,beta=0,gama1=0.5,gama2=0.2,gama3=0.2,gama4=0.5,phi=1,shift=0,shift2=0)/100
     #mu=exp(h1(N=N,k=1,k2=1,alpha=9,beta=0,gama1=0.5,gama2=0.2,gama3=0.2,gama4=0.5,shift=0,shift2=0))
     
     set.seed(i)
     out2=outbreak5(currentday=days5*52*years,weeklength=52*days5*years,wtime=(length(yt)-49*days5+1):length(yt),yi=yt,interval=0.25,k=1,k2=1,alpha=9,beta=0,gama1=0.5,
                    gama2=0.2,gama3=0.2,gama4=0.5,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=w*days5,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     #zitot[(bankhols5==1)]=0
     #zitot[(bankhols5==1)+1]=1.5*zitot[i+1]
     
     for(b in 1:length(zitot)){
          if(bankhols5[b]==1){
               zitot[b]=0
               zitot[b+1]=1.5*zitot[b+1]
          } 
     }
     
     zeros=rep(0,2)
     weekend=seq(days5,days5*years*52,days5)
     #weekend=seq(0,days5*years*52-1,days5)
     for(s in 1:length(weekend)){
          zt=append(zt,zeros,after=2*(s-1)+weekend[s])
          zitot=append(zitot,zeros,after=2*(s-1)+weekend[s])
          zoutbreak=append(zoutbreak,zeros,after=2*(s-1)+weekend[s])
     }
     
     simulateddata13[,i]=round(zt)
     simulatedtotals13[,i]=round(zitot)
     simulatedoutbreak13[,i]=round(zoutbreak)
}

plot(1:length(simulatedtotals13[,1]),simulatedtotals13[,1],typ='l')

plot(1:N,simulatedtotals13[,1],typ='l',xlim=c(2206,2548),col='green')
lines(1:N,simulateddata13[,1],typ='l')


#=====================
# 7-day week syndromes
#=====================

years=7
days7=7
N=52*days7*years


#sigid1
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,gama3=0.5,gama4=0.4,phi=2,shift=29)
     #mu=exp(h2(N=N,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,gama3=0.5,gama4=0.4,shift=29))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,
                    gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata1[,i]=round(zt)
     simulatedtotals1[,i]=round(zitot)
     simulatedoutbreak1[,i]=round(zoutbreak)
}

#sigid3
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,gama2=1.4,gama3=0.5,gama4=0.4,phi=1,shift=-167)
     #mu=exp(h2(N=N,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,gama2=1.4,gama3=0.5,gama4=0.4,shift=-167))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=0.5,beta=0,gama1=1.5,
                    gama2=1.4,gama3=0.5,gama4=0.4,phi=1,shift=-167,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata3[,i]=round(zt)
     simulatedtotals3[,i]=round(zitot)
     simulatedoutbreak3[,i]=round(zoutbreak)
}


#sigid4
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=0,k2=2,alpha=5.5,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1)
     #mu=exp(h2(N=N,k=0,k2=2,alpha=5.5,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*12,wtime=(length(yt)-49*7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=5.5,beta=0,gama1=0,
                    gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata4[,i]=round(zt)
     simulatedtotals4[,i]=round(zitot)
     simulatedoutbreak4[,i]=round(zoutbreak)
}


#sigid5
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=0,k2=2,alpha=2,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1)
     #mu=exp(h2(N=N,k=0,k2=2,alpha=2,beta=0,gama1=0,gama2=0,gama3=0.3,gama4=0.25,shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=2,beta=0,gama1=0,
                    gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata5[,i]=round(zt)
     simulatedtotals5[,i]=round(zitot)
     simulatedoutbreak5[,i]=round(zoutbreak)
}


#sigid14
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=1,k2=2,alpha=2,beta=0.0005,gama1=0.8,gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57)
     #mu=exp(h2(N=N,k=1,k2=2,alpha=2,beta=0,gama1=0.8,gama2=0.8,gama3=0.8,gama4=0.4,shift=57))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=6,beta=0,gama1=0.2,gama2=0.2,
                    gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata14[,i]=round(zt)
     simulatedtotals14[,i]=round(zitot)
     simulatedoutbreak14[,i]=round(zoutbreak)
}


#sigid15
for(i in 1:nsim){
     set.seed(i)
     #yt=0.1*(negbinNoise2(N=N,k=4,k2=1,alpha=1.5,beta=0,gama1=0.1,gama2=0.1,gama3=1.8,gama4=0.1,phi=1,shift=-85)+2)
     
     yt=1*(negbinNoise2(N=N,k=4,k2=1,alpha=0.05,beta=0,gama1=0.01,gama2=0.01,gama3=1.8,gama4=0.1,phi=1,shift=-85)+0)
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=2,beta=0,gama1=0.8,
                    gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata15[,i]=round(zt)
     simulatedtotals15[,i]=round(zitot)
     simulatedoutbreak15[,i]=round(zoutbreak)
}


#plot(1:N,yt,typ='l')
#plot(1:(52*years*7),zt,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1,364))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(365,728))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(729,1092))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1093,1456))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1457,1820))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(1821,2184))
#plot(1:(52*years*7),zitot,typ='l',xlim=c(2185,2548))
#lines(1:(52*years*7),zoutbreak,col='green')

plot(1:(52*years*7),simulatedtotals6[,4],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak6[,4],col='green')
lines(1:(52*years*7),simulatedoutbreak6[,4],col='red')

plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedzseasoutbreak6[,4],col='green',typ='l')


#sigid16
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=1,k2=2,alpha=3,beta=0,gama1=0.8,gama2=0.6,gama3=0.8,gama4=0.4,phi=4,shift=29)
     #mu=exp(h2(N=N,k=1,k2=2,alpha=3,beta=0,gama1=0.8,gama2=0.6,gama3=0.8,gama4=0.4,shift=29))
     
     out1=rep(0,N)
     for(j in 1:years){
          set.seed(j+years*i)
          out=outbreak5(currentday=days7*52*years,weeklength=52*days7*years,wtime=((210+(j-1)*days7*52):(230+(j-1)*days7*52)),yi=yt,interval=0.02,k=1,k2=1,alpha=1,beta=0,gama1=0.1,
                        gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=w*days7*150,meanlog=0,sdlog=0.5)
          out1=out1+out$outbreak
     }
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*days7*years,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=1,k2=2,alpha=3,beta=0,gama1=0.8,
                    gama2=0.6,gama3=0.8,gama4=0.4,phi=4,shift=29,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     zseasoutbreak=out2$outbreak+out1
     
     zt=yt +out1 
     zitot=yt + out2$outbreak +out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata16[,i]=round(zt)
     simulatedtotals16[,i]=round(zitot)
     simulatedoutbreak16[,i]=round(zoutbreak)
     simulatedzseasoutbreak16[,i]=round(zseasoutbreak)
}


plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedtotals16[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak16[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak16[,1],col='red')


plot(1:(52*years*7),simulatedtotals16[,2],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak16[,2],col='green')
lines(1:(52*years*7),simulatedoutbreak16[,2],col='red')


#sigid17
for(i in 1:nsim){
     set.seed(i)
     yt=negbinNoise2(N=N,k=0,k2=2,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.8,gama4=0.4,phi=4,shift=1)
     #mu=exp(h2(N=N,k=0,k2=2,alpha=6,beta=0,gama1=0,gama2=0,gama3=0.8,gama4=0.4,shift=1))
     
     set.seed(i)
     out2=outbreak7(currentday=N,weeklength=52*7*12,wtime=(length(yt)-49*days7+1):length(yt),yi=yt,interval=0.25,k=0,k2=2,alpha=6,beta=0,gama1=0,
                    gama2=0,gama3=0.8,gama4=0.4,phi=4,shift=1,numoutbk=1,peakoutbk=w*days7,meanlog=0,sdlog=0.5)
     
     zoutbreak=out2$outbreak
     
     zt=yt#+out1 
     zitot=yt + out2$outbreak #+out1
     
     for(b in 1:length(zitot)){
          if(bankhols7[b]==1){
               zitot[b]=2*zitot[b]
          } 
     }
     
     simulateddata17[,i]=round(zt)
     simulatedtotals17[,i]=round(zitot)
     simulatedoutbreak17[,i]=round(zoutbreak)
}



###############
#PLOTS###########
################
years=7
i=2
par(mfrow=(c(4,4)))
plot(1:(52*years*7),simulateddata1[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal1')
#plot(1:(52*years*7),simulateddata2[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal2')
plot(1:(52*years*7),simulateddata3[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal2')
plot(1:(52*years*7),simulateddata4[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal3')
plot(1:(52*years*7),simulateddata5[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal4')
plot(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5')
plot(1:(52*years*7),simulateddata7[,i]-simulatedzseasoutbreak7[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6')
plot(1:(52*years*7),simulateddata8[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal7')
plot(1:(52*years*7),simulateddata9[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal8')
plot(1:(52*years*7),simulateddata10[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal9')
plot(1:(52*years*7),simulateddata11[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal10')
plot(1:(52*years*7),simulateddata12[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal11')
plot(1:(52*years*7),simulateddata13[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal12')
plot(1:(52*years*7),simulateddata14[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal13')
plot(1:(52*years*7),simulateddata15[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal14')
plot(1:(52*years*7),simulateddata16[,i]-simulatedzseasoutbreak16[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15')
plot(1:(52*years*7),simulateddata17[,i],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal16')


par(mfrow=(c(2,1)))
plot(1:(52*years*7),simulateddata4[,1],typ='l',xlim=c(1,3*7),xlab='day',ylab='count',main='signal3',xaxt="n")
axis(1,at=1:21,1:21,las=2)
plot(1:(52*years*7),simulateddata8[,1],typ='l',xlim=c(1,3*7),xlab='day',ylab='count',main='signal7',xaxt="n")
axis(1,at=1:21,1:21,las=2)

i=98
par(mfrow=(c(4,4)))
plot(1:(52*years*7),simulateddata1[,i]+simulatedoutbreak1[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal1',col='blue')
lines(1:(52*years*7),simulateddata1[,i],col='black')
#plot(1:(52*years*7),simulateddata2[,i]+simulatedoutbreak2[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal2',col='blue')
#lines(1:(52*years*7),simulateddata2[,i],col='black')
plot(1:(52*years*7),simulateddata3[,i]+simulatedoutbreak3[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal2',col='blue')
lines(1:(52*years*7),simulateddata3[,i],col='black')
plot(1:(52*years*7),simulateddata4[,i]+simulatedoutbreak4[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal3',col='blue')
lines(1:(52*years*7),simulateddata4[,i],col='black')
plot(1:(52*years*7),simulateddata5[,i]+simulatedoutbreak5[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal4',col='blue')
lines(1:(52*years*7),simulateddata5[,i],col='black')
plot(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i]+simulatedoutbreak6[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata6[,i]-simulatedzseasoutbreak6[,i],col='black')
plot(1:(52*years*7),simulateddata7[,i]-simulatedzseasoutbreak7[,i]+simulatedoutbreak7[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata7[,i]-simulatedzseasoutbreak7[,i],col='black')
plot(1:(52*years*7),simulateddata8[,i]+simulatedoutbreak8[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal7',col='blue')
lines(1:(52*years*7),simulateddata8[,i],col='black')
plot(1:(52*years*7),simulateddata9[,i]+simulatedoutbreak9[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal8',col='blue')
lines(1:(52*years*7),simulateddata9[,i],col='black')
plot(1:(52*years*7),simulateddata10[,i]+simulatedoutbreak10[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal9',col='blue')
lines(1:(52*years*7),simulateddata10[,i],col='black')
plot(1:(52*years*7),simulateddata11[,i]+simulatedoutbreak11[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal10',col='blue')
lines(1:(52*years*7),simulateddata11[,i],col='black')
plot(1:(52*years*7),simulateddata12[,i]+simulatedoutbreak12[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal11',col='blue')
lines(1:(52*years*7),simulateddata12[,i],col='black')
plot(1:(52*years*7),simulateddata13[,i]+simulatedoutbreak13[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal12',col='blue')
lines(1:(52*years*7),simulateddata13[,i],col='black')
plot(1:(52*years*7),simulateddata14[,i]+simulatedoutbreak14[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal13',col='blue')
lines(1:(52*years*7),simulateddata14[,i],col='black')
plot(1:(52*years*7),simulateddata15[,i]+simulatedoutbreak15[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal14',col='blue')
lines(1:(52*years*7),simulateddata15[,i],col='black')
plot(1:(52*years*7),simulateddata16[,i]-simulatedzseasoutbreak16[,i]+simulatedoutbreak16[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata16[,i]-simulatedzseasoutbreak16[,i],col='black')
plot(1:(52*years*7),simulateddata17[,i]+simulatedoutbreak17[,i],typ='l',xlim=c(52*years*6,52*years*7),xlab='day',ylab='count',main='signal16',col='blue')
lines(1:(52*years*7),simulateddata17[,i],col='black')



par(mfrow=(c(3,1)))
plot(1:(52*years*7),simulateddata6[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata6[,1]-simulatedzseasoutbreak6[,1],col='black')

plot(1:(52*years*7),simulateddata7[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata7[,1]-simulatedzseasoutbreak7[,1],col='black')

plot(1:(52*years*7),simulateddata16[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata16[,1]-simulatedzseasoutbreak16[,1],col='black')


par(mfrow=(c(3,1)))
plot(1:(52*years*7),simulateddata6[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal5',col='blue')
lines(1:(52*years*7),simulateddata6[,1]-simulatedzseasoutbreak6[,1],col='black')

plot(1:(52*years*7),simulateddata7[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal6',col='blue')
lines(1:(52*years*7),simulateddata7[,1]-simulatedzseasoutbreak7[,1],col='black')

plot(1:(52*years*7),simulateddata16[,1],typ='l',xlim=c(1,52*years*7),xlab='day',ylab='count',main='signal15',col='blue')
lines(1:(52*years*7),simulateddata16[,1]-simulatedzseasoutbreak16[,1],col='black')






plot(1:(52*years*7),zitot,typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),yt,col='blue')
lines(1:(52*years*7),zseasoutbreak,col='green')
lines(1:(52*years*7),zoutbreak,col='red')

plot(1:(52*years*7),simulatedtotals16[,1],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak16[,1],col='green')
lines(1:(52*years*7),simulatedoutbreak16[,1],col='red')


plot(1:(52*years*7),simulatedtotals16[,2],typ='l',xlim=c(1,7*364))
lines(1:(52*years*7),simulatedzseasoutbreak16[,2],col='green')
lines(1:(52*years*7),simulatedoutbreak16[,2],col='red')


#############################
#############################
# FIT FARRINGTON ALGO TO DATA
#############################
#############################

require(surveillance)
require(zoo)

#============================================
# Four functions that implement the algorithm
#============================================

algo.farrington=function (disProgObj, control = list(range = NULL, b = 3, w = 3,
                                                     reweight = TRUE, verbose = FALSE, alpha = 0.01, trend = TRUE,
                                                     limit54 = c(5, 4), powertrans = "2/3", fitFun = "algo.farrington.fitGLM"))
{
     observed <- disProgObj$observed
     
     freq <- disProgObj$freq
     #epochStr <- switch(as.character(freq), `12` = "1 month",
     #    `52` = "1 week", `365` = "1 day")
     if (is.null(control$range)) {
          control$range <- (freq * control$b - control$w):length(observed)
     }
     if (is.null(control$b)) {
          control$b = 5
     }
     if (is.null(control$w)) {
          control$w = 3
     }
     if (is.null(control$reweight)) {
          control$reweight = TRUE
     }
     if (is.null(control$verbose)) {
          control$verbose = FALSE
     }
     if (is.null(control$alpha)) {
          control$alpha = 0.01
     }
     if (is.null(control$trend)) {
          control$trend = TRUE
     }
     if (is.null(control$plot)) {
          control$plot = FALSE
     }
     if (is.null(control$limit54)) {
          control$limit54 = c(5, 4)
     }
     if (is.null(control$powertrans)) {
          control$powertrans = "2/3"
     }
     if (is.null(control$fitFun)) {
          control$fitFun = "algo.farrington.fitGLM"
     }
     #else {
     #    control$fitFun <- match.arg(control$fitFun, c("algo.farrington.fitGLM"))
     #}
     #if (is.null(disProgObj[["epochAsDate", exact = TRUE]])) {
     #    epochAsDate <- FALSE
     #}
     #else {
     #    epochAsDate <- disProgObj[["epochAsDate", exact = TRUE]]
     #}
     if (!((control$limit54[1] >= 0) & (control$limit54[2] > 0))) {
          stop("The limit54 arguments are out of bounds: cases >= 0 and period > 0.")
     }
     # alarmall gives the results for all data sets whereas alarm identifies the
     # ones where there are less than 5 reports in the last 4 weeks
     alarmall <- matrix(data = 0, nrow = length(control$range), ncol = 1)
     alarm<- matrix(data = 0, nrow = length(control$range), ncol = 1)
     upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
     trend <- matrix(data = 0, nrow = length(control$range), ncol = 1)
     pd <- matrix(data = 0, nrow = length(control$range), ncol = 2)
     
     n <- control$b * (2 * control$w + 1)
     
     # identify the leading zeroes in the data and cut them off
     
     #observed1=observed[observed!=NA]
     
     for (k in control$range) {
          #observed<-rollapply(data=observed[1:k], width=7, FUN=sum, ascending = F)
          #observed<-rollapply(data=a, width=7, FUN=sum, ascending = F)
          observed=observed[1:k]
          ob=rep(0,round(length(observed)/7))
          for(w in 0:(round(length(observed)/7)-1)){
               if((k-6*w-6-w)<=0){break}
               ob[w+1]=sum(observed[(k-6*w-w):(k-6*w-6-w)])
          }
          observed=ob[length(ob):1]
          kk=length(observed)
          for(z in 1:length(observed)){
               if(observed[z]==0){observed[z]=NA}
               else{break}
          }
          if (control$verbose) {
               cat("k=", k, "\n")
          }
          # If the remaining data is less than one year's worth,
          # do not fit a model and indicate that by reporting the value 1e+300
          if((kk-z)<(2*round(freq)+control$w+1)){
               upperbound[k - min(control$range) + 1] = 1e+300
               alarmall[k - min(control$range) + 1] = 1e+300
               alarm[k - min(control$range) + 1] = 1e+300
               trend[k - min(control$range) + 1] = 1e+300
          }
          else{
               #if (!epochAsDate) {
               # Assign level factors to the data (I have done that manually given 
               # that the data can be missing...)
               seas1 <- NULL
               seas2 <- NULL
               for (i in control$b:0) {
                    if(i==3){seas1=append(seas1, seq(kk - round(freq * i) - 
                                                          control$w, kk - round(freq * i) +control$w+1, by = 1))}
                    else{seas1=append(seas1, seq(kk - round(freq * i) - 
                                                      control$w, kk - round(freq * i) +control$w, by = 1))}
                    seas2=append(seas2, seq(kk - round(freq * i) - 
                                                 control$w-5, kk - round(freq * i) -control$w-1, by = 1))
               }
               #}
               seas3=seas2-5
               seas4=seas3-5
               seas5=seas4-5
               seas6=seas5-5
               seas7=seas6-5
               seas8=seas7-5
               seas9=seas8-5
               seas10=seas9-5
               seasgroup=rep(0,(length(observed)+control$w))
               seasgroup[seas1]=1
               seasgroup[seas2]=2
               seasgroup[seas3]=3
               seasgroup[seas4]=4
               seasgroup[seas5]=5
               seasgroup[seas6]=6
               seasgroup[seas7]=7
               seasgroup[seas8]=8
               seasgroup[seas9]=9
               seasgroup[seas10]=10
               if((kk-z)<(freq*control$b+control$w+1)){
                    seasgroup=seasgroup[z:(kk-27)]
                    seasgroup=as.factor(seasgroup)
                    wtime= z:(kk-27) 
               }
               else{
                    seasgroup=seasgroup[(kk-control$b*freq-control$w-1):(kk-27)]
                    seasgroup=as.factor(seasgroup)
                    wtime= (kk-control$b*freq-control$w-1):(kk-27) 
               }      
               response <- observed[wtime]   
               if (control$verbose) {
                    print(response)
               }
               
               v=length(response[response>0])
               toosmall=(v<2) # indicator for when (stricktly) less than 2 weeks have non-zero counts
               
               # If stricktly less than 2 weeks have non-zero counts,
               # do not fit a model and indicate that by reporting the value 1e+400
               
               if(toosmall){
                    upperbound[k - min(control$range) + 1] = 1e+400
                    alarmall[k - min(control$range) + 1] = 1e+400
                    alarm[k - min(control$range) + 1] = 1e+400
                    trend[k - min(control$range) + 1] = 1e+400
               }
               else{
                    oneyear=FALSE
                    p=wtime[response>0]
                    oneyear=((p[length(p)]-p[1])<52)
                    # if all non-zero weeks are within one year, do not fit a trend.
                    # This is the only case where we do not fit a trend.
                    if(oneyear){
                         model <- do.call(control$fitFun, args = list(response = response,
                                                                      wtime = wtime, seasgroup=seasgroup,timeTrend = FALSE, reweight = control$reweight))
                    }
                    else{
                         model <- do.call(control$fitFun, args = list(response = response,
                                                                      wtime = wtime, seasgroup=seasgroup,timeTrend = control$trend, reweight = control$reweight))
                    }
                    if(is.null(model)){return(model)}
                    
                    doTrend <- control$trend
                    if (model$T==1) {
                         wp <- summary.glm(model)$coefficients["wtime", 4] # p-value after reweighting
                         # In our new approach, we fit a trend always (unless all non-zero weeks are within one year).
                         # For this reason, we consider the trend to be always significant (wp<=1).
                         # If one wants to fit a trend only if significant, change this to, say, wp<=0.05.
                         significant <- (wp <=1)
                         mu0Hat <- predict.glm(model, data.frame(wtime = c(kk),seasgroup=factor(1)),type = "response")
                         #atLeastThreeYears <- (control$b >= 3)
                         # We remove the noExtrapolation condition
                         #noExtrapolation <- mu0Hat <= max(response)
                         if (!(significant)) {
                              doTrend <- FALSE
                              model <- do.call(control$fitFun, args = list(response = response,
                                                                           wtime = wtime,seasgroup=seasgroup, timeTrend = FALSE, reweight = control$reweight))
                         }
                    }
                    else {
                         doTrend <- FALSE
                    }
                    if(model$phi<1){model$phi=1}
                    pred <- predict.glm(model, data.frame(wtime = c(kk),seasgroup=factor(1)),
                                        dispersion = model$phi, type = "response", se.fit = TRUE)
                    # We use negative binomial quantiles to define the error structure rather 
                    # than the ones based on the transformed Poisson and the Anscombe residuals.         
                    if(model$phi==1){
                         lu<-c(qpois(control$alpha,pred$fit),qpois(1-control$alpha,pred$fit))
                         lu[2]=max(1,lu[2])
                    }
                    else{
                         lu<-c(qnbinom(control$alpha,pred$fit/(model$phi-1),1/model$phi),qnbinom(1-control$alpha,pred$fit/(model$phi-1),1/model$phi))
                         lu[2]=max(1,lu[2])
                    }
                    #if (control$plot) {
                    #        data <- data.frame(wtime = seq(min(wtime), k, length = 1000))
                    #        preds <- predict(model, data, type = "response",
                    #            dispersion = model$phi)
                    #        plot(c(wtime, k), c(response, observed[k]), ylim = range(c(observed[data$wtime],
                    #            lu)), , xlab = "time", ylab = "No. infected",
                    #            main = paste("Prediction at time t=", k, " with b=",
                    #              control$b, ",w=", control$w, sep = ""), pch = c(rep(1,
                    #              length(wtime)), 16))
                    #        lines(data$wtime, preds, col = 1, pch = 2)
                    #        lines(rep(k, 2), lu[1:2], col = 3, lty = 2)
                    #    }
                    
                    enoughCases <- (sum(observed[(kk - control$limit54[2] + 1):kk]) >= control$limit54[1])
                    
                    X <- (observed[kk] - pred$fit)/(lu[2] - pred$fit)
                    upperbound[k - min(control$range) + 1] <- lu[2]
                    alarm[k - min(control$range) + 1] <- (X > 1)
                    # if the last five weeks have less than 4 reports,
                    # indicate that by returning the value 1e-300
                    if(!enoughCases){alarm[k - min(control$range) + 1] = 1e-300}
                    alarmall[k - min(control$range) + 1] <- (X > 1)
                    trend[k - min(control$range) + 1] <- doTrend
               }
          }
          observed <- disProgObj$observed
     }
     control$name <- paste("farrington(", control$w, ",", 0, ",", control$b, ")", sep = "")
     control$data <- paste(deparse(substitute(disProgObj)))
     result <- list(alarm = alarm, alarmall = alarmall, trend = trend,
                    disProgObj = disProgObj, control = control,upperbound=upperbound)
     class(result) <- "survRes"
     return(result)
}

#=========================================

algo.farrington.fitGLM=function (response, wtime,seasgroup, timeTrend = TRUE, reweight = TRUE)
{
     theModel <- as.formula(ifelse(timeTrend, "response~seasgroup+wtime",
                                   "response~seasgroup"))
     model <- glm(theModel, family = quasipoisson(link = "log"))
     # In our new approach, we fit a trend always (unless all non-zero weeks are within one year).
     # For this reason, we consider the trend always significant (p<=1).
     # If one wants to fit a trend only if significant, change this to, say, p<=0.05.
     if (model$converged){
          if (timeTrend) {
               p<- summary.glm(model)$coefficients["wtime", 4] # p-value before reweighting
               if(p<=1){T=1}
               else{
                    T=0
                    model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"))
                    if(!model$converged) {
                         cat("Warning: No convergence without insignificant timeTrend.\n")
                         print(cbind(response, wtime))
                         return(NULL)
                    }
               }
          }
          if(!timeTrend){T=0}
     }
     
     else {
          T=0
          if (timeTrend) {
               model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"))
               cat("Warning: No convergence with timeTrend -- trying without.\n")
          }
          if (!model$converged) {
               cat("Warning: No convergence without timeTrend.\n")
               print(cbind(response, wtime))
               return(NULL)
          }
     }
     
     phi <- max(summary(model)$dispersion, 1)
     if (reweight) {
          s <- anscombe.residuals(model, phi)
          omega <- algo.farrington.assign.weights(s)
          theModel <- as.formula(ifelse(T==1, "response~seasgroup+wtime", "response~seasgroup"))
          model <- glm(theModel, family = quasipoisson(link = "log"),weights = omega)
          if (!model$converged) {
               
               if (T==1) {
                    model <- glm(response ~ seasgroup, family = quasipoisson(link = "log"),weights = omega)
                    cat("Warning: No convergence with weights and timeTrend -- trying without.\n")
               }
               if (!model$converged) {
                    cat("Warning: No convergence with weights but without timeTrend.\n")
                    print(cbind(response, wtime))
                    return(NULL)
               }
               T=0
          }
          
          phi <- max(summary(model)$dispersion, 1)
     }
     model$phi <- phi
     model$T=T
     return(model)
}

#=========================================


algo.farrington.assign.weights=function (s)
{   # we use the cut off point s>2.58 as a compromise between s>2 and s>3 
     gamma <- length(s)/(sum((s^(-2))^(s > 2.58)))
     omega <- numeric(length(s))
     omega[s > 2.58] <- gamma * (s[s > 2.58]^(-2))
     omega[s <= 2.58] <- gamma
     return(omega)
}


#=========================================


anscombe.residuals=function (m, phi) 
{
     y <- m$y
     mu <- fitted.values(m)
     a <- 3/2 * (y^(2/3) * mu^(-1/6) - mu^(1/2))
     a <- a/sqrt(phi * (1 - hatvalues(m)))
     return(a)
}


#==========================================



#=============================
# Define the alarm data frames
#=============================

days=7

nsim=100

alarm1=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm2=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm3=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm4=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm5=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm6=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm7=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm8=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm9=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm10=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm11=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm12=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm13=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm14=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm15=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm16=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarm17=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))



alarmall1=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall2=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall3=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall4=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall5=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall6=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall7=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall8=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall9=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall10=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall11=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall12=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall13=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall14=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall15=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall16=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))
alarmall17=data.frame(array(rep(0,nsim*49*days),dim=c(49*days,nsim)))


#========================================
#Implement the algorithm to data by days
#========================================

#require(surveillance)

b=5
weeks=N/7

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals1[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals1)-49*days+1):nrow(simulatedtotals1),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm1[,i]=a.farrington$alarm[,1]
     alarmall1[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals3[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals3)-49*days+1):nrow(simulatedtotals3),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm3[,i]=a.farrington$alarm[,1]
     alarmall3[,i]=a.farrington$alarmall[,1]
}
s
for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals4[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals4)-49*days+1):nrow(simulatedtotals4),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm4[,i]=a.farrington$alarm[,1]
     alarmall4[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals5[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals5)-49*days+1):nrow(simulatedtotals5),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm5[,i]=a.farrington$alarm[,1]
     alarmall5[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals6[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals6)-49*days+1):nrow(simulatedtotals6),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm6[,i]=a.farrington$alarm[,1]
     alarmall6[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals7[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals7)-49*days+1):nrow(simulatedtotals7),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm7[,i]=a.farrington$alarm[,1]
     alarmall7[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals8[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals8)-49*days+1):nrow(simulatedtotals8),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm8[,i]=a.farrington$alarm[,1]
     alarmall8[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals9[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals9)-49*days+1):nrow(simulatedtotals9),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm9[,i]=a.farrington$alarm[,1]
     alarmall9[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals10[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals10)-49*days+1):nrow(simulatedtotals10),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm10[,i]=a.farrington$alarm[,1]
     alarmall10[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals11[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals11)-49*days+1):nrow(simulatedtotals11),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm11[,i]=a.farrington$alarm[,1]
     alarmall11[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals12[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals12)-49*days+1):nrow(simulatedtotals12),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm12[,i]=a.farrington$alarm[,1]
     alarmall12[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals13[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals13)-49*days+1):nrow(simulatedtotals13),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm13[,i]=a.farrington$alarm[,1]
     alarmall13[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals14[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals14)-49*days+1):nrow(simulatedtotals14),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm14[,i]=a.farrington$alarm[,1]
     alarmall14[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals15[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals15)-49*days+1):nrow(simulatedtotals15),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm15[,i]=a.farrington$alarm[,1]
     alarmall15[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals16[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals16)-49*days+1):nrow(simulatedtotals16),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm16[,i]=a.farrington$alarm[,1]
     alarmall16[,i]=a.farrington$alarmall[,1]
}

for(i in 1:nsim){
     a.disProg=create.disProg(week=weeks,observed=simulatedtotals17[,i],freq=52.18)
     cntrl <- list(range=(nrow(simulatedtotals17)-49*days+1):nrow(simulatedtotals17),w = 3, b = 5, alpha = 0.01,trend=TRUE, fitFun="algo.farrington.fitGLM")
     a.farrington<- algo.farrington(a.disProg, control = cntrl)
     alarm17[,i]=a.farrington$alarm[,1]
     alarmall17[,i]=a.farrington$alarmall[,1]
}


save.image("//penelope/MCSUsers/Staff/an3936/Documents/own research/RAMMIE/R/outx2.Rdata")

#====================================
#====================================
#Summary
#====================================
#====================================

days=7


# FPR false positive rate

fpr=rep(0,17)
fprseas=rep(0,3)

nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall1[nrow(alarmall1)-i+1,j]==1 & simulatedoutbreak1[nrow(simulatedoutbreak1)-i+1,j]==0)+nu
     }
}
a=
     fpr[1]=nu/sum(simulatedoutbreak1[2206:2548,]==0)  



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall2[nrow(alarmall2)-i+1,j]==1 & simulatedoutbreak2[nrow(simulatedoutbreak2)-i+1,j]==0)+nu
     }
}
fpr[2]=nu/sum(simulatedoutbreak2[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall3[nrow(alarmall3)-i+1,j]==1 & simulatedoutbreak3[nrow(simulatedoutbreak3)-i+1,j]==0)+nu
     }
}
fpr[3]=nu/sum(simulatedoutbreak3[2206:2548,]==0)   



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall4[nrow(alarmall4)-i+1,j]==1 & simulatedoutbreak4[nrow(simulatedoutbreak4)-i+1,j]==0)+nu
     }
}
fpr[4]=nu/sum(simulatedoutbreak4[2206:2548,]==0) 




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall5[nrow(alarmall5)-i+1,j]==1 & simulatedoutbreak5[nrow(simulatedoutbreak5)-i+1,j]==0)+nu
     }
}
fpr[5]=nu/sum(simulatedoutbreak5[2206:2548,]==0)   



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedoutbreak6[nrow(simulatedoutbreak6)-i+1,j]==0)+nu
     }
}
fpr[6]=nu/sum(simulatedoutbreak6[2206:2548,]==0)   



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedzseasoutbreak6[nrow(simulatedzseasoutbreak6)-i+1,j]==0)+nu
     }
}
fprseas[1]=nu/sum(simulatedzseasoutbreak6[2206:2548,]==0)   



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedoutbreak7[nrow(simulatedoutbreak7)-i+1,j]==0)+nu
     }
}
fpr[7]=nu/sum(simulatedoutbreak7[2206:2548,]==0)   


nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedzseasoutbreak7[nrow(simulatedzseasoutbreak7)-i+1,j]==0)+nu
     }
}
fprseas[2]=nu/sum(simulatedzseasoutbreak7[2206:2548,]==0) 


nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall8[nrow(alarmall8)-i+1,j]==1 & simulatedoutbreak8[nrow(simulatedoutbreak8)-i+1,j]==0)+nu
     }
}
fpr[8]=nu/sum(simulatedoutbreak8[2206:2548,]==0)   



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall9[nrow(alarmall9)-i+1,j]==1 & simulatedoutbreak9[nrow(simulatedoutbreak9)-i+1,j]==0)+nu
     }
}
fpr[9]=nu/sum(simulatedoutbreak9[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall10[nrow(alarmall10)-i+1,j]==1 & simulatedoutbreak10[nrow(simulatedoutbreak10)-i+1,j]==0)+nu
     }
}
fpr[10]=nu/sum(simulatedoutbreak10[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall11[nrow(alarmall11)-i+1,j]==1 & simulatedoutbreak11[nrow(simulatedoutbreak11)-i+1,j]==0)+nu
     }
}
fpr[11]=nu/sum(simulatedoutbreak11[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall12[nrow(alarmall12)-i+1,j]==1 & simulatedoutbreak12[nrow(simulatedoutbreak12)-i+1,j]==0)+nu
     }
}
fpr[12]=nu/sum(simulatedoutbreak12[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall13[nrow(alarmall13)-i+1,j]==1 & simulatedoutbreak13[nrow(simulatedoutbreak13)-i+1,j]==0)+nu
     }
}
fpr[13]=nu/sum(simulatedoutbreak13[2206:2548,]==0)  



nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall14[nrow(alarmall14)-i+1,j]==1 & simulatedoutbreak14[nrow(simulatedoutbreak14)-i+1,j]==0)+nu
     }
}
fpr[14]=nu/sum(simulatedoutbreak14[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall15[nrow(alarmall15)-i+1,j]==1 & simulatedoutbreak15[nrow(simulatedoutbreak15)-i+1,j]==0)+nu
     }
}
fpr[15]=nu/sum(simulatedoutbreak15[2206:2548,]==0)   




nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedoutbreak16[nrow(simulatedoutbreak16)-i+1,j]==0)+nu
     }
}
fpr[16]=nu/sum(simulatedoutbreak16[2206:2548,]==0)   


nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedzseasoutbreak16[nrow(simulatedzseasoutbreak16)-i+1,j]==0)+nu
     }
}
fprseas[3]=nu/sum(simulatedzseasoutbreak16[2206:2548,]==0)   


nu=0
for(j in 1:nsim){
     for(i in (49*7):1){
          nu=(alarmall17[nrow(alarmall17)-i+1,j]==1 & simulatedoutbreak17[nrow(simulatedoutbreak17)-i+1,j]==0)+nu
     }
}
fpr[17]=nu/sum(simulatedoutbreak17[2206:2548,]==0)   


#--------------------------------------------------------

# POD power of detection

pod=rep(0,17)
podseas=rep(0,3)

mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall1[nrow(alarmall1)-i+1,j]==1 & simulatedoutbreak1[nrow(simulatedoutbreak1)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[1]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall2[nrow(alarmall2)-i+1,j]==1 & simulatedoutbreak2[nrow(simulatedoutbreak2)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[2]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall3[nrow(alarmall3)-i+1,j]==1 & simulatedoutbreak3[nrow(simulatedoutbreak3)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[3]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall4[nrow(alarmall4)-i+1,j]==1 & simulatedoutbreak4[nrow(simulatedoutbreak4)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[4]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall5[nrow(alarmall5)-i+1,j]==1 & simulatedoutbreak5[nrow(simulatedoutbreak5)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[5]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedoutbreak6[nrow(simulatedoutbreak6)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[6]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedzseasoutbreak6[nrow(simulatedzseasoutbreak6)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
podseas[1]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedoutbreak7[nrow(simulatedoutbreak7)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[7]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedzseasoutbreak7[nrow(simulatedzseasoutbreak7)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
podseas[2]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall8[nrow(alarmall8)-i+1,j]==1 & simulatedoutbreak8[nrow(simulatedoutbreak8)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[8]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall9[nrow(alarmall9)-i+1,j]==1 & simulatedoutbreak9[nrow(simulatedoutbreak9)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[9]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall10[nrow(alarmall10)-i+1,j]==1 & simulatedoutbreak10[nrow(simulatedoutbreak10)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[10]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall11[nrow(alarmall11)-i+1,j]==1 & simulatedoutbreak11[nrow(simulatedoutbreak11)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[11]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall12[nrow(alarmall12)-i+1,j]==1 & simulatedoutbreak12[nrow(simulatedoutbreak12)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[12]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall13[nrow(alarmall13)-i+1,j]==1 & simulatedoutbreak13[nrow(simulatedoutbreak13)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[13]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall14[nrow(alarmall14)-i+1,j]==1 & simulatedoutbreak14[nrow(simulatedoutbreak14)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[14]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall15[nrow(alarmall15)-i+1,j]==1 & simulatedoutbreak15[nrow(simulatedoutbreak15)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[15]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedoutbreak16[nrow(simulatedoutbreak16)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[16]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedzseasoutbreak16[nrow(simulatedzseasoutbreak16)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
podseas[3]=mu/nsim


mu=0
for(j in 1:nsim){
     nu=0
     for(i in (49*days):1){
          nu=nu+(alarmall17[nrow(alarmall17)-i+1,j]==1 & simulatedoutbreak17[nrow(simulatedoutbreak17)-i+1,j]>0)
     }
     mu=mu+(nu>0)
}
pod[17]=mu/nsim



#--------------------------------------------------------

# Sensitivity

sensitivity=rep(0,17)
sensitivityseas=rep(0,3)

nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall1[nrow(alarmall1)-i+1,j]==1 & simulatedoutbreak1[nrow(simulatedoutbreak1)-i+1,j]>0)
     }
}
sensitivity[1]=nu/sum(simulatedoutbreak1>0)

nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall2[nrow(alarmall2)-i+1,j]==1 & simulatedoutbreak2[nrow(simulatedoutbreak2)-i+1,j]>0)
     }
}
sensitivity[2]=nu/sum(simulatedoutbreak2>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall3[nrow(alarmall3)-i+1,j]==1 & simulatedoutbreak3[nrow(simulatedoutbreak3)-i+1,j]>0)
     }
}
sensitivity[3]=nu/sum(simulatedoutbreak3>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall4[nrow(alarmall4)-i+1,j]==1 & simulatedoutbreak4[nrow(simulatedoutbreak4)-i+1,j]>0)
     }
}
sensitivity[4]=nu/sum(simulatedoutbreak4>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall5[nrow(alarmall5)-i+1,j]==1 & simulatedoutbreak5[nrow(simulatedoutbreak5)-i+1,j]>0)
     }
}
sensitivity[5]=nu/sum(simulatedoutbreak5>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedoutbreak6[nrow(simulatedoutbreak6)-i+1,j]>0)
     }
}
sensitivity[6]=nu/sum(simulatedoutbreak6>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedzseasoutbreak6[nrow(simulatedzseasoutbreak6)-i+1,j]>0)
     }
}
sensitivityseas[1]=nu/sum(simulatedzseasoutbreak6>0)



nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedoutbreak7[nrow(simulatedoutbreak7)-i+1,j]>0)
     }
}
sensitivity[7]=nu/sum(simulatedoutbreak7>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedzseasoutbreak7[nrow(simulatedzseasoutbreak7)-i+1,j]>0)
     }
}
sensitivityseas[2]=nu/sum(simulatedzseasoutbreak7>0)



nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall8[nrow(alarmall8)-i+1,j]==1 & simulatedoutbreak8[nrow(simulatedoutbreak8)-i+1,j]>0)
     }
}
sensitivity[8]=nu/sum(simulatedoutbreak8>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall9[nrow(alarmall9)-i+1,j]==1 & simulatedoutbreak9[nrow(simulatedoutbreak9)-i+1,j]>0)
     }
}
sensitivity[9]=nu/sum(simulatedoutbreak9>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall10[nrow(alarmall10)-i+1,j]==1 & simulatedoutbreak10[nrow(simulatedoutbreak10)-i+1,j]>0)
     }
}
sensitivity[10]=nu/sum(simulatedoutbreak10>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall11[nrow(alarmall11)-i+1,j]==1 & simulatedoutbreak11[nrow(simulatedoutbreak11)-i+1,j]>0)
     }
}
sensitivity[11]=nu/sum(simulatedoutbreak11>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall12[nrow(alarmall12)-i+1,j]==1 & simulatedoutbreak12[nrow(simulatedoutbreak12)-i+1,j]>0)
     }
}
sensitivity[12]=nu/sum(simulatedoutbreak12>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall13[nrow(alarmall13)-i+1,j]==1 & simulatedoutbreak13[nrow(simulatedoutbreak13)-i+1,j]>0)
     }
}
sensitivity[13]=nu/sum(simulatedoutbreak13>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall14[nrow(alarmall14)-i+1,j]==1 & simulatedoutbreak14[nrow(simulatedoutbreak14)-i+1,j]>0)
     }
}
sensitivity[14]=nu/sum(simulatedoutbreak14>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall15[nrow(alarmall15)-i+1,j]==1 & simulatedoutbreak15[nrow(simulatedoutbreak15)-i+1,j]>0)
     }
}
sensitivity[15]=nu/sum(simulatedoutbreak15>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedoutbreak16[nrow(simulatedoutbreak16)-i+1,j]>0)
     }
}
sensitivity[16]=nu/sum(simulatedoutbreak16>0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedzseasoutbreak16[nrow(simulatedzseasoutbreak16)-i+1,j]>0)
     }
}
sensitivityseas[3]=nu/sum(simulatedzseasoutbreak16>0)



nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall17[nrow(alarmall17)-i+1,j]==1 & simulatedoutbreak17[nrow(simulatedoutbreak17)-i+1,j]>0)
     }
}
sensitivity[17]=nu/sum(simulatedoutbreak17>0)



#--------------------------------------------------------

# Specificity


specificity=rep(0,17)
specificityseas=rep(0,3)


# Specificity
nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall1[nrow(alarmall1)-i+1,j]==0 & simulatedoutbreak1[nrow(simulatedoutbreak1)-i+1,j]==0)
     }
}
specificity[1]=nu/sum(simulatedoutbreak1[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall2[nrow(alarmall2)-i+1,j]==0 & simulatedoutbreak2[nrow(simulatedoutbreak2)-i+1,j]==0)
     }
}
specificity[2]=nu/sum(simulatedoutbreak2[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall3[nrow(alarmall3)-i+1,j]==0 & simulatedoutbreak3[nrow(simulatedoutbreak3)-i+1,j]==0)
     }
}
specificity[3]=nu/sum(simulatedoutbreak3[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall4[nrow(alarmall4)-i+1,j]==0 & simulatedoutbreak4[nrow(simulatedoutbreak4)-i+1,j]==0)
     }
}
specificity[4]=nu/sum(simulatedoutbreak4[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall5[nrow(alarmall5)-i+1,j]==0 & simulatedoutbreak5[nrow(simulatedoutbreak5)-i+1,j]==0)
     }
}
specificity[5]=nu/sum(simulatedoutbreak5[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==0 & simulatedoutbreak6[nrow(simulatedoutbreak6)-i+1,j]==0)
     }
}
specificity[6]=nu/sum(simulatedoutbreak6[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall6[nrow(alarmall6)-i+1,j]==0 & simulatedzseasoutbreak6[nrow(simulatedzseasoutbreak6)-i+1,j]==0)
     }
}
specificityseas[1]=nu/sum(simulatedzseasoutbreak6[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==0 & simulatedoutbreak7[nrow(simulatedoutbreak7)-i+1,j]==0)
     }
}
specificity[7]=nu/sum(simulatedoutbreak7[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall7[nrow(alarmall7)-i+1,j]==0 & simulatedzseasoutbreak7[nrow(simulatedzseasoutbreak7)-i+1,j]==0)
     }
}
specificityseas[2]=nu/sum(simulatedzseasoutbreak7[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall8[nrow(alarmall8)-i+1,j]==0 & simulatedoutbreak8[nrow(simulatedoutbreak8)-i+1,j]==0)
     }
}
specificity[8]=nu/sum(simulatedoutbreak8[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall9[nrow(alarmall9)-i+1,j]==0 & simulatedoutbreak9[nrow(simulatedoutbreak9)-i+1,j]==0)
     }
}
specificity[9]=nu/sum(simulatedoutbreak9[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall10[nrow(alarmall10)-i+1,j]==0 & simulatedoutbreak10[nrow(simulatedoutbreak10)-i+1,j]==0)
     }
}
specificity[10]=nu/sum(simulatedoutbreak10[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall11[nrow(alarmall11)-i+1,j]==0 & simulatedoutbreak11[nrow(simulatedoutbreak11)-i+1,j]==0)
     }
}
specificity[11]=nu/sum(simulatedoutbreak11[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall12[nrow(alarmall12)-i+1,j]==0 & simulatedoutbreak12[nrow(simulatedoutbreak12)-i+1,j]==0)
     }
}
specificity[12]=nu/sum(simulatedoutbreak12[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall13[nrow(alarmall13)-i+1,j]==0 & simulatedoutbreak13[nrow(simulatedoutbreak13)-i+1,j]==0)
     }
}
specificity[13]=nu/sum(simulatedoutbreak13[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall14[nrow(alarmall14)-i+1,j]==0 & simulatedoutbreak14[nrow(simulatedoutbreak14)-i+1,j]==0)
     }
}
specificity[14]=nu/sum(simulatedoutbreak14[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall15[nrow(alarmall15)-i+1,j]==0 & simulatedoutbreak15[nrow(simulatedoutbreak15)-i+1,j]==0)
     }
}
specificity[15]=nu/sum(simulatedoutbreak15[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==0 & simulatedoutbreak16[nrow(simulatedoutbreak16)-i+1,j]==0)
     }
}
specificity[16]=nu/sum(simulatedoutbreak16[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall16[nrow(alarmall16)-i+1,j]==0 & simulatedzseasoutbreak16[nrow(simulatedzseasoutbreak16)-i+1,j]==0)
     }
}
specificityseas[3]=nu/sum(simulatedzseasoutbreak16[2206:2548,]==0)


nu=0
for(j in 1:nsim){
     for(i in (49*days):1){
          nu=nu+(alarmall17[nrow(alarmall17)-i+1,j]==0 & simulatedoutbreak17[nrow(simulatedoutbreak17)-i+1,j]==0)
     }
}
specificity[17]=nu/sum(simulatedoutbreak17[2206:2548,]==0)


#----------------------------------------------

# Timeliness

timeliness=rep(0,17)
timelinessseas=rep(0,3)

n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak1[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak1[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall1[nrow(alarmall1)-i+1,j]==1 & simulatedoutbreak1[nrow(simulatedoutbreak1)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak1)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[1]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak2[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak2[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall2[nrow(alarmall2)-i+1,j]==1 & simulatedoutbreak2[nrow(simulatedoutbreak2)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak2)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[2]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak3[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak3[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall3[nrow(alarmall3)-i+1,j]==1 & simulatedoutbreak3[nrow(simulatedoutbreak3)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak3)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[3]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak4[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak4[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall4[nrow(alarmall4)-i+1,j]==1 & simulatedoutbreak4[nrow(simulatedoutbreak4)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak4)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[4]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak5[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak5[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall5[nrow(alarmall5)-i+1,j]==1 & simulatedoutbreak5[nrow(simulatedoutbreak5)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak5)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[5]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak6[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak6[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedoutbreak6[nrow(simulatedoutbreak6)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak6)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[6]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedzseasoutbreak6[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedzseasoutbreak6[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall6[nrow(alarmall6)-i+1,j]==1 & simulatedzseasoutbreak6[nrow(simulatedzseasoutbreak6)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedzseasoutbreak6)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timelinessseas[1]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak7[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak7[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedoutbreak7[nrow(simulatedoutbreak7)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak7)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[7]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedzseasoutbreak7[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedzseasoutbreak7[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall7[nrow(alarmall7)-i+1,j]==1 & simulatedzseasoutbreak7[nrow(simulatedzseasoutbreak7)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedzseasoutbreak7)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timelinessseas[2]=(ss+n)/nsim



n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak8[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak8[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall8[nrow(alarmall8)-i+1,j]==1 & simulatedoutbreak8[nrow(simulatedoutbreak8)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak8)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[8]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak9[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak9[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall9[nrow(alarmall9)-i+1,j]==1 & simulatedoutbreak9[nrow(simulatedoutbreak9)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak9)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[9]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak10[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak10[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall10[nrow(alarmall10)-i+1,j]==1 & simulatedoutbreak10[nrow(simulatedoutbreak10)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak10)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[10]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak11[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak11[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall11[nrow(alarmall11)-i+1,j]==1 & simulatedoutbreak11[nrow(simulatedoutbreak11)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak11)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[11]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak12[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak12[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall12[nrow(alarmall12)-i+1,j]==1 & simulatedoutbreak12[nrow(simulatedoutbreak12)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak12)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[12]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak13[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak13[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall13[nrow(alarmall13)-i+1,j]==1 & simulatedoutbreak13[nrow(simulatedoutbreak1)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak13)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[13]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak14[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak14[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall14[nrow(alarmall14)-i+1,j]==1 & simulatedoutbreak14[nrow(simulatedoutbreak14)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak14)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[14]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak15[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak15[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall15[nrow(alarmall15)-i+1,j]==1 & simulatedoutbreak15[nrow(simulatedoutbreak15)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak15)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[15]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak16[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak16[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedoutbreak16[nrow(simulatedoutbreak16)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak16)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[16]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedzseasoutbreak16[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedzseasoutbreak16[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall16[nrow(alarmall16)-i+1,j]==1 & simulatedzseasoutbreak16[nrow(simulatedzseasoutbreak16)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedzseasoutbreak16)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timelinessseas[3]=(ss+n)/nsim


n=0
ss=0
for(j in 1:nsim){
     for(i in (52*days*years):(52*days*(years-1)+3*days+1)){
          test=(simulatedoutbreak17[i,j]>0)
          if(test==T){
               r2=i
               break
          }
     }
     for(i in (52*(years-1)*days+3*days+1):(52*years*days)){
          test=(simulatedoutbreak17[i,j]>0)
          if(test==T){
               r1=i
               break
          }
     }
     for(i in (49*days):1){
          test=(alarmall17[nrow(alarmall17)-i+1,j]==1 & simulatedoutbreak17[nrow(simulatedoutbreak17)-i+1,j]>0)
          if(test==T){
               ss=ss+(nrow(simulatedoutbreak17)-i+1-r1)/(r2-r1+1)
               break
          }
     }
     if(i==1 & test!=T){n=n+1}
}
timeliness[17]=(ss+n)/nsim

