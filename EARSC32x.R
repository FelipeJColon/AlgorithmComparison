## ############################################################################
##
## DISCLAIMER: 
## This script has been developed for research purposes only. 
## The script is provided without any warranty of any kind, either express or 
## implied. The entire risk arising out of the use or performance of the sample
## script and documentation remains with you. 
## In no event shall its author, or anyone else involved in the 
## creation, production, or delivery of the script be liable for any damages 
## whatsoever (including, without limitation, damages for loss of business 
## profits, business interruption, loss of business information, or other 
## pecuniary loss) arising out of the use of or inability to use the sample
## scripts or documentation, even if the author has been advised of the
## possibility of such damages. 
##
## ############################################################################
##
## DESCRIPTION
## Simulates outbreaks and analyses them using EARS-C3
##
##
## Written by: Angela Noufaily and Felipe J Colón-González
## For any problems with this code, please contact f.colon@uea.ac.uk
## 
## ############################################################################

rm(list=ls(all=TRUE))

# FUNCTIONS THAT PRODUCE THE DATA
# DEFINING FUNCTION h
require(data.table)
require(dplyr)
require(tidyr)
require(surveillance)
require(lubridate)
require(zoo)

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

myDir <- "/local/zck07apu/Documents/GitLab/rammie_comparison/scripts/C3/2x"
years=7
bankholidays=read.csv(file.path(myDir, "Bankholidays.csv"))

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
                    gama2=2,gama3=0.3,gama4=0.5,phi=1.5,shift=-50,shift2=-50,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0,gama3=0.6,gama4=0.9,phi=1.5,shift=0,shift2=0,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0.1,gama3=0.2,gama4=0.3,phi=1,shift=-150,shift2=-150,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0.1,gama3=0.05,gama4=0.15,phi=1,shift=-200,shift2=-200,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0.1,gama3=0.05,gama4=0.1,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0,gama3=0.05,gama4=0.15,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama2=0.2,gama3=0.2,gama4=0.5,phi=1,shift=0,shift2=0,numoutbk=1,peakoutbk=2*days5,meanlog=0,sdlog=0.5)
     
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
                    gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                    gama2=1.4,gama3=0.5,gama4=0.4,phi=1,shift=-167,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                    gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                    gama2=0,gama3=0.3,gama4=0.25,phi=1,shift=1,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                    gama3=0.5,gama4=0.4,phi=2,shift=29,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                    gama2=0.8,gama3=0.8,gama4=0.4,phi=4,shift=57,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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
                        gama2=2,gama3=0.05,gama4=0.05,phi=1,shift=-50,shift2=-50,numoutbk=1,peakoutbk=3*days7*150,meanlog=0,sdlog=0.5)
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
                    gama2=0,gama3=0.8,gama4=0.4,phi=4,shift=1,numoutbk=1,peakoutbk=2*days7,meanlog=0,sdlog=0.5)
     
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




#=============================
# Define the alarm data frames
#=============================

days=7

nsim=100

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


##############################################
#========================================
#Implement the algorithm to data by days and record the alarms in the above dataframes
#========================================
##############################################
myDates   <- seq(ymd('2010-01-01'), ymd('2016-12-30'), by = '1 day')

dropDays <- as.POSIXct(c('2010-12-31','2011-12-31', '2012-12-31',
                         '2013-12-31', '2014-12-31', '2015-12-31', 
                         '2016-02-29,', '2012-02-29'))

"%ni%" <- Negate("%in%")
myDates <- myDates[myDates %ni% dropDays]

# Convert to 7-day running totals
rolling <- function(x){
     rollapplyr(x, width=7, FUN=sum, na.rm=T, fill=NA)
} 

simdata1 <- apply(simulateddata1, 2, rolling)
# simdata2 <- apply(simulateddata2, 2, rolling)
simdata3 <- apply(simulateddata3, 2, rolling)
simdata4 <- apply(simulateddata4, 2, rolling)
simdata5 <- apply(simulateddata5, 2, rolling)
simdata6 <- apply(simulateddata6, 2, rolling)
simdata7 <- apply(simulateddata7, 2, rolling)
simdata8 <- apply(simulateddata8, 2, rolling)
simdata9 <- apply(simulateddata9, 2, rolling)
simdata10 <- apply(simulateddata10, 2, rolling)
simdata11 <- apply(simulateddata11, 2, rolling)
simdata12 <- apply(simulateddata12, 2, rolling)
simdata13 <- apply(simulateddata13, 2, rolling)
simdata14 <- apply(simulateddata14, 2, rolling)
simdata15 <- apply(simulateddata15, 2, rolling)
simdata16 <- apply(simulateddata16, 2, rolling)
simdata17 <- apply(simulateddata17, 2, rolling)

simtot1 <- apply(simulatedtotals1, 2, rolling)
# simtot2 <- apply(simulatedtotals2, 2, rolling)
simtot3 <- apply(simulatedtotals3, 2, rolling)
simtot4 <- apply(simulatedtotals4, 2, rolling)
simtot5 <- apply(simulatedtotals5, 2, rolling)
simtot6 <- apply(simulatedtotals6, 2, rolling)
simtot7 <- apply(simulatedtotals7, 2, rolling)
simtot8 <- apply(simulatedtotals8, 2, rolling)
simtot9 <- apply(simulatedtotals9, 2, rolling)
simtot10 <- apply(simulatedtotals10, 2, rolling)
simtot11 <- apply(simulatedtotals11, 2, rolling)
simtot12 <- apply(simulatedtotals12, 2, rolling)
simtot13 <- apply(simulatedtotals13, 2, rolling)
simtot14 <- apply(simulatedtotals14, 2, rolling)
simtot15 <- apply(simulatedtotals15, 2, rolling)
simtot16 <- apply(simulatedtotals16, 2, rolling)
simtot17 <- apply(simulatedtotals17, 2, rolling)

# Convert data to sts 
simSts1  <- sts(simdata1, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
# simSts2  <- sts(simdata2, start=c(2010, 1), frequency=364,
# epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts3  <- sts(simdata3, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts4  <- sts(simdata4, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts5  <- sts(simdata5, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts6  <- sts(simdata6, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts7  <- sts(simdata7, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts8  <- sts(simdata8, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts9  <- sts(simdata9, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts10 <- sts(simdata10, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts11 <- sts(simdata11, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts12 <- sts(simdata12, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts13 <- sts(simdata13, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts14 <- sts(simdata14, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts15 <- sts(simdata15, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts16 <- sts(simdata16, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
simSts17 <- sts(simdata17, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)


totSts1  <- sts(simtot1, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
# totSts2 <- sts(simtot2, start=c(2010, 1), frequency=364,
# epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts3  <- sts(simtot3, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts4  <- sts(simtot4, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts5  <- sts(simtot5, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts6  <- sts(simtot6, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts7  <- sts(simtot7, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts8  <- sts(simtot8, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts9  <- sts(simtot9, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts10 <- sts(simtot10, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts11 <- sts(simtot11, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts12 <- sts(simtot12, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts13 <- sts(simtot13, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts14 <- sts(simtot14, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts15 <- sts(simtot15, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts16 <- sts(simtot16, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)
totSts17 <- sts(simtot17, start=c(2010, 1), frequency=364,
                epoch=as.numeric(as.Date(myDates)), epochAsDate=TRUE)

in2016 <- 2206:2548

# Select range of data to monitor, algorithm and prediction interval
control <- list(range=in2016, method="C3", alpha=0.01)


for(sim in seq(nsim)){
     cat("\t", sim)
     
     # Run detection algorithm
     det1  <- earsC(totSts1[,sim], control=control)
     det3  <- earsC(totSts3[,sim], control=control)
     det4  <- earsC(totSts4[,sim], control=control)
     det5  <- earsC(totSts5[,sim], control=control)
     det6  <- earsC(totSts6[,sim], control=control)
     det7  <- earsC(totSts7[,sim], control=control)
     det8  <- earsC(totSts8[,sim], control=control)
     det9  <- earsC(totSts9[,sim], control=control)
     det10 <- earsC(totSts10[,sim], control=control)
     det11 <- earsC(totSts11[,sim], control=control)
     det12 <- earsC(totSts12[,sim], control=control)
     det13 <- earsC(totSts13[,sim], control=control)
     det14 <- earsC(totSts14[,sim], control=control)
     det15 <- earsC(totSts15[,sim], control=control)
     det16 <- earsC(totSts16[,sim], control=control)
     det17 <- earsC(totSts17[,sim], control=control)
     
     # Plot detection results
     dir.create(file.path(myDir, "plots", "totals"), 
                recursive=TRUE)
     
     png(file.path(myDir, "plots", "totals", paste0("Sim_", sim, ".png")),
         width=16,height=14,units="in",res=300)     
     par(mfrow=c(4, 4), oma=c(0, 0, 2, 0))
     plot(det1, main="Dataset 1", legend=NULL)
     plot(det3, main="Dataset 3", legend=NULL)
     plot(det4, main="Dataset 4", legend=NULL)
     plot(det5, main="Dataset 5", legend=NULL)
     plot(det6, main="Dataset 6", legend=NULL)
     plot(det7, main="Dataset 7", legend=NULL)
     plot(det8, main="Dataset 8", legend=NULL)
     plot(det9, main="Dataset 9", legend=NULL)
     plot(det10, main="Dataset 10", legend=NULL)
     plot(det11, main="Dataset 11", legend=NULL)
     plot(det12, main="Dataset 12", legend=NULL)
     plot(det13, main="Dataset 13", legend=NULL)
     plot(det14, main="Dataset 14", legend=NULL)
     plot(det15, main="Dataset 15", legend=NULL)
     plot(det16, main="Dataset 16", legend=NULL)
     plot(det17, main="Dataset 17", legend=NULL)
     title(main=list(paste("Simulation", sim, "Alpha", control$alpha ), 
                     cex=2), outer=TRUE)
     dev.off()
     
     # Retrieve information about alarms
     alarmall1[,sim] <- as.numeric(as.vector(unlist(det1@alarm)))
     alarmall3[,sim] <- as.numeric(as.vector(unlist(det3@alarm)))
     alarmall4[,sim] <- as.numeric(as.vector(unlist(det4@alarm)))
     alarmall5[,sim] <- as.numeric(as.vector(unlist(det5@alarm)))
     alarmall6[,sim] <- as.numeric(as.vector(unlist(det6@alarm)))
     alarmall7[,sim] <- as.numeric(as.vector(unlist(det7@alarm)))
     alarmall8[,sim] <- as.numeric(as.vector(unlist(det8@alarm)))
     alarmall9[,sim] <- as.numeric(as.vector(unlist(det9@alarm)))
     alarmall10[,sim] <- as.numeric(as.vector(unlist(det10@alarm)))
     alarmall11[,sim] <- as.numeric(as.vector(unlist(det11@alarm)))
     alarmall12[,sim] <- as.numeric(as.vector(unlist(det12@alarm)))
     alarmall13[,sim] <- as.numeric(as.vector(unlist(det13@alarm)))
     alarmall14[,sim] <- as.numeric(as.vector(unlist(det14@alarm)))
     alarmall15[,sim] <- as.numeric(as.vector(unlist(det15@alarm)))
     alarmall16[,sim] <- as.numeric(as.vector(unlist(det16@alarm)))
     alarmall17[,sim] <- as.numeric(as.vector(unlist(det17@alarm)))
     
}

# Replace missing values with zero (?)
alarmall1[is.na(alarmall1)] <- 0
alarmall3[is.na(alarmall3)] <- 0
alarmall4[is.na(alarmall4)] <- 0
alarmall5[is.na(alarmall5)] <- 0
alarmall6[is.na(alarmall6)] <- 0
alarmall7[is.na(alarmall7)] <- 0
alarmall8[is.na(alarmall8)] <- 0
alarmall9[is.na(alarmall9)] <- 0
alarmall10[is.na(alarmall10)] <- 0
alarmall11[is.na(alarmall11)] <- 0
alarmall12[is.na(alarmall12)] <- 0
alarmall13[is.na(alarmall13)] <- 0
alarmall14[is.na(alarmall14)] <- 0
alarmall15[is.na(alarmall15)] <- 0
alarmall16[is.na(alarmall16)] <- 0
alarmall17[is.na(alarmall17)] <- 0


# Compare vs data without oubreaks
for(sim in seq(nsim)){
     cat("\t", sim)
     
     det1  <- earsC(simSts1[,sim], control=control)
     det3  <- earsC(simSts3[,sim], control=control)
     det4  <- earsC(simSts4[,sim], control=control)
     det5  <- earsC(simSts5[,sim], control=control)
     det6  <- earsC(simSts6[,sim], control=control)
     det7  <- earsC(simSts7[,sim], control=control)
     det8  <- earsC(simSts8[,sim], control=control)
     det9  <- earsC(simSts9[,sim], control=control)
     det10 <- earsC(simSts10[,sim], control=control)
     det11 <- earsC(simSts11[,sim], control=control)
     det12 <- earsC(simSts12[,sim], control=control)
     det13 <- earsC(simSts13[,sim], control=control)
     det14 <- earsC(simSts14[,sim], control=control)
     det15 <- earsC(simSts15[,sim], control=control)
     det16 <- earsC(simSts16[,sim], control=control)
     det17 <- earsC(simSts17[,sim], control=control)
     
     dir.create(file.path(myDir, "plots", "control"), 
                recursive=TRUE)
     
     png(file.path(myDir, "plots", "control", 
                   paste0("Sim_", sim, ".png")),
         width=16,height=14,units="in",res=300)     
     par(mfrow=c(4, 4), oma=c(0, 0, 2, 0))
     plot(det1, main="Dataset 1", legend=NULL)
     plot(det3, main="Dataset 3", legend=NULL)
     plot(det4, main="Dataset 4", legend=NULL)
     plot(det5, main="Dataset 5", legend=NULL)
     plot(det6, main="Dataset 6", legend=NULL)
     plot(det7, main="Dataset 7", legend=NULL)
     plot(det8, main="Dataset 8", legend=NULL)
     plot(det9, main="Dataset 9", legend=NULL)
     plot(det10, main="Dataset 10", legend=NULL)
     plot(det11, main="Dataset 11", legend=NULL)
     plot(det12, main="Dataset 12", legend=NULL)
     plot(det13, main="Dataset 13", legend=NULL)
     plot(det14, main="Dataset 14", legend=NULL)
     plot(det15, main="Dataset 15", legend=NULL)
     plot(det16, main="Dataset 16", legend=NULL)
     plot(det17, main="Dataset 17", legend=NULL)
     title(main=list(paste("Simulation", sim, "Alpha", control$alpha ), 
                     cex=2), outer=TRUE)
     dev.off()
}



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


#==================================

# Summary=data.frame(fpr,pod,sensitivity,specificity,timeliness)
# row.names(Summary)=c("sigid1","sigid2","sigid3","sigid4","sigid5","sigid6","sigid7","sigid8","sigid9","sigid10","sigid11","sigid12","sigid13","sigid14","sigid15","sigid16","sigid17")
# 
# Summaryseas=data.frame(fprseas,podseas,sensitivityseas,specificityseas,timelinessseas)
# row.names(Summaryseas)=c("sigid6","sigid7","sigid16")
# 
# 
# fix(Summary)
# fix(Summaryseas)
# 
# 
summary1=data.frame(fpr, pod, sensitivity, specificity, timeliness)
row.names(summary1)=c("sigid1", "sigid2", "sigid3", "sigid4", "sigid5",
                      "sigid6", "sigid7", "sigid8", "sigid9", "sigid10",
                      "sigid11", "sigid12", "sigid13", "sigid14", "sigid15",
                      "sigid16","sigid17")


summary2=data.frame(fprseas, podseas, sensitivityseas,
                    specificityseas, timelinessseas)
row.names(summary2)=c("sigid6", "sigid7", "sigid16")

if(!dir.exists(file.path(myDir, "output"))){
     dir.create(file.path(myDir, "output"))
}

write.csv(summary1, file.path(myDir, "output", "summaryC3-18.csv"), 
          row.names=FALSE)


write.csv(summary2, file.path(myDir, "output", "summarySeasC3-18.csv"), 
          row.names=FALSE)
