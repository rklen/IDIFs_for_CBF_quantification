linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
dispersionCorrection<-function(input,tau){
  #cA<-0
  #for(i in 2:length(input)){
  #  cA[i]<-input[i]+tau*(input[i]-input[i-1])
  #}
  #return(cA)
  cA<-0
  for(i in 2:length(input)){
    cA<-c(cA,(1/tau)*input[i-1]+(1-(1/tau))*cA[i-1])
  }
  return(cA)
}
fit_model<-function(param,input,j){
  k1<-abs(param[1])/60
  k2<-abs(param[2])/60
  Va<-1/(1+exp(-param[3]))
  if(length(param)==3){
    cT<-0
    for(i in 2:length(input)){
      if(i-j-1<1){c_i<-0}else{
        if(i-j-1>length(input)){c_i<-input[length(input)]}else{c_i<-input[i-j-1]}
      }
      cT<-c(cT,k1*c_i+(1-k2)*cT[i-1])
    }
    v<-(1-Va)*cT+Va*input
  }
  if(length(param)==4){
    tau<-5/(1+exp(-param[4]))
    input1<-dispersionCorrection(input,tau)
    cT<-0
    for(i in 2:length(input1)){
      if(i-j-1<1){c_i<-0}else{
        if(i-j-1>length(input1)){c_i<-input1[length(input1)]}else{c_i<-input1[i-j-1]}
      }
      cT<-c(cT,k1*c_i+(1-k2)*cT[i-1])
    }
    v<-(1-Va)*cT+Va*input1
  }
  return(v)
}
errorfunction<-function(param,input,organCurve,j){
  return(sum((organCurve-fit_model(param,input,j))^2))
}
meanRelError<-function(modelCurve,organCurve){
  r<-0
  for(i in 2:length(organCurve)){
    if(modelCurve[i]!=organCurve[i]){
      r<-r+abs(modelCurve[i]-organCurve[i])/organCurve[i]
    }
  }
  r<-r/(length(organCurve)-1)
  return(r)
}
optimizationAlg<-function(initialValues,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,j=j,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
fit_input<-function(param,input){
  tau<-5/(1+exp(-param))
  input1<-0
  for(i in 2:length(input)){
    input1<-c(input1,1/tau*input[i-1]+(1-1/tau)*input1[i-1])
  }
  return(input1)
}
errorfunction_input<-function(param,input,input1){
  return(sum((input1-fit_input(param,input))^2))
}
optimizationAlg_input<-function(initialValues,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction_input,initialValues,input=input,
               input1=input1,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}

files<-list.files('C:/Users/oonar/Documents/idif')
studyNumbers<-c()
for(i in 1:length(files)){
  if(substr(files[i],1,7)=='koveri0'){
    studyNumbers<-c(studyNumbers,substr(files[i],1,10))
  }
}
t<-c(0,2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,
     62.5,67.5,75,85,95,110,130,150,175,205,235,265)
times<-c(0:265)

#IDIF, no dispersion
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=8)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[2,])
    organCurve[l]<-linearInter(times[l],t,df[8,])
  }
  plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1]<-abs(param[1])
  df1[i,2]<-abs(param[2])
  df1[i,3]<-1/(1+exp(-param[3]))
  df1[i,5]<-j
  df1[i,6]<-err/1000/1000/length(times)
  df1[i,7]<-length(times)*log(df1[i,6])+2*(length(param)+1)
  df1[i,8]<-meanRelError(fit_model(param,input,j),organCurve)
}
write.csv(df1,file=paste('df_input2_global',sep=''),row.names=F)

#IDIF, dispersion
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=8)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[2,])
    organCurve[l]<-linearInter(times[l],t,df[8,])
  }
  plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9),-log(4))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1]<-abs(param[1])
  df1[i,2]<-abs(param[2])
  df1[i,3]<-1/(1+exp(-param[3]))
  df1[i,4]<-5/(1+exp(-param[4]))
  df1[i,5]<-j
  df1[i,6]<-err/1000/1000/length(times-1)
  df1[i,7]<-length(times-1)*log(df1[i,6])+2*(length(param)+1)
  df1[i,8]<-meanRelError(fit_model(param,input,j),organCurve)
}
write.csv(df1,file=paste('df_input4_global',sep=''),row.names=F)

param_for_tau<-c()
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input<-c()
  input1<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    input1[l]<-linearInter(times[l],t,df[2,])
  }
  plot(times,input1,type='l',ylab=i)
  initialValues<-c(-log(4))
  param<-optimizationAlg_input(initialValues,input,input1)
  points(times,fit_input(param,input),type='l')
  param_for_tau[i]<-param
}
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=8)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input_0<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input_0[l]<-linearInter(times[l],t,df[1,])
    organCurve[l]<-linearInter(times[l],t,df[8,])
  }
  input<-fit_input(param_for_tau[i],input_0)
  plot(times,input,type='l')
  plot(times,organCurve,type='l',ylab=i)
  initialValues<-c(1,1,-log(9))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1]<-abs(param[1])
  df1[i,2]<-abs(param[2])
  df1[i,3]<-1/(1+exp(-param[3]))
  df1[i,4]<-param_for_tau[i]
  df1[i,5]<-j
  df1[i,6]<-err/1000/1000/length(times-1)
  df1[i,7]<-length(times-1)*log(df1[i,6])+2*(length(param)+2)
  df1[i,8]<-meanRelError(fit_model(param,input,j),organCurve)
}
write.csv(df1,file=paste('df_input5_global',sep=''),row.names=F)



print(c(mean(5/(1+exp(-df1[,4]))),sd(5/(1+exp(-df1[,4])))))



print(c(mean(df1[,1]),sd(df1[,1])))#0.45517058 vs 0.6439956 vs 0.45146088
print(c(mean(df1[,2]),sd(df1[,2])))#0.51844758 vs 0.5694702 vs 0.51358781
print(c(mean(df1[,3]),sd(df1[,3])))#0.023339275 vs 0.02622780 vs 0.02053052
print(c(mean(df1[,4]),sd(df1[,4])))# - vs - vs 1.842825
print(c(mean(df1[,5]),sd(df1[,5])))#0.03883495
print(c(mean(df1[,7]),sd(df1[,7])))#4457.5227 vs 4434.135



df_input1<-matrix(data=NA,nrow=length(studyNumbers),ncol=266)
df_input2<-matrix(data=NA,nrow=length(studyNumbers),ncol=266)
df_input3<-matrix(data=NA,nrow=length(studyNumbers),ncol=266)
df_input4<-matrix(data=NA,nrow=length(studyNumbers),ncol=266)
df_input5<-matrix(data=NA,nrow=length(studyNumbers),ncol=266)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input1<-c()
  input2<-c()
  for(l in 1:length(times)){
    input1[l]<-linearInter(times[l],t,df[1,])
    input2[l]<-linearInter(times[l],t,df[2,])
  }
  df_input1[i,]<-input1
  df_input2[i,]<-input2
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input3_global')
  tau<-df1[i,4]
  df_input3[i,]<-dispersionCorrection(input1,tau)
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input4_global')
  tau<-df1[i,4]
  df_input4[i,]<-dispersionCorrection(input2,tau)
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')
  tau<-5/(1+exp(-df1[i,4]))
  df_input5[i,]<-dispersionCorrection(input1,tau)
}
df_input<-df_input5
plot(1,type='n',xlim=c(0,265),ylim=c(0,max(df_input)/1000),
     ylab='Activity concentration (kBq/mL)',
     xlab='Time (s)',
     cex.lab=1.3,cex.axis=1.5,
     main='Hybrid IDIFs')
for(i in 1:length(studyNumbers)){
  points(times,df_input[i,]/1000,type='l',col='steelblue1',lwd=2)
}
points(times,colMeans(df_input)/1000,type='l',col='blue',lwd=2)
colSds<-c()
for(i in 1:266){
  colSds[i]<-sd(df_input[,i])
}
points(times,(colMeans(df_input)-colSds)/1000,type='l',col='blue',lwd=2,lty=2)
points(times,(colMeans(df_input)+colSds)/1000,type='l',col='blue',lwd=2,lty=2)

plot(1,type='n',xlim=c(0,max(df_input1)/1000),
     ylim=c(0,max(df_input2)/1000),
     ylab='Max. activity value in CCA IDIFs (kBq/mL)',
     xlab='Max. activity value in aorta IDIFs (kBq/mL)',
     cex.lab=1.3,cex.axis=1.5)
y1<-c()
y2<-c()
for(i in 1:length(studyNumbers)){
  points(max(df_input1[i,]/1000),max(df_input2[i,]/1000),
         col='steelblue1',lwd=2)
  y1[i]<-max(df_input1[i,]/1000)
  y2[i]<-max(df_input2[i,]/1000)
}
fit<-lm(y2~y1)
abline(fit,col='blue',lty=2,lwd=2)
cor.test(y1,y2,method='spearman')

plot(1,type='n',xlim=c(0,max(df_input2)/1000),
     ylim=c(0,max(df_input5)/1000),
     ylab='Max. activity value in hybrid IDIFs (kBq/mL)',
     xlab='Max. activity value in CCA IDIFs (kBq/mL)',
     cex.lab=1.3,cex.axis=1.5)
y2<-c()
y5<-c()
for(i in 1:length(studyNumbers)){
  points(max(df_input2[i,]/1000),max(df_input5[i,]/1000),
         col='steelblue1',lwd=2)
  y2[i]<-max(df_input2[i,]/1000)
  y5[i]<-max(df_input5[i,]/1000)
}
fit<-lm(y5~y2)
abline(fit,col='blue',lty=2,lwd=2)
cor.test(y2,y5,method='spearman')

plot(1,type='n',xlim=c(0,max(df_input1)/1000),
     ylim=c(0,max(df_input5)/1000),
     ylab='Max. activity value in hybrid IDIFs (kBq/mL)',
     xlab='Max. activity value in aorta IDIFs (kBq/mL)',
     cex.lab=1.3,cex.axis=1.5)
y1<-c()
y5<-c()
for(i in 1:length(studyNumbers)){
  points(max(df_input1[i,]/1000),max(df_input5[i,]/1000),
         col='steelblue1',lwd=2)
  y1[i]<-max(df_input1[i,]/1000)
  y5[i]<-max(df_input5[i,]/1000)
}
fit<-lm(y5~y1)
abline(fit,col='blue',lty=2,lwd=2)
cor.test(y1,y5,method='spearman')

plot(1,type='n',xlim=c(0,max(df_input1[,266])/1000),
     ylim=c(0,max(df_input2[,266])/1000),
     ylab='Last activity value in CCA IDIFs (kBq/mL)',
     xlab='Last activity value in aorta IDIFs (kBq/mL)',
     cex.lab=1.3,cex.axis=1.5)
y1<-c()
y2<-c()
for(i in 1:length(studyNumbers)){
  points(df_input1[i,266]/1000,df_input2[i,266]/1000,
         col='steelblue1',lwd=2)
  y1[i]<-df_input1[i,266]/1000
  y2[i]<-df_input2[i,266]/1000
}
fit<-lm(y2~y1)
abline(fit,col='blue',lty=2,lwd=2)
cor.test(y1,y2,method='spearman')






for(i in 1:length(studyNumbers)){
  if(max(df_input[i,])>500000){
    print(studyNumbers[i])
  }
}





df1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_global')
string<-paste(round(mean(df1[,1]),3),'pmin',round(sd(df1[,1]),3),' & ',
              round(mean(df1[,2]),3),'pmin',round(sd(df1[,2]),3),' & ',
              round(mean(df1[,3]),3),'pmin',round(sd(df1[,3]),3),' & ',
              round(mean(df1[,5]),3),'pmin',round(sd(df1[,5]),3),' & ',
              round(mean(df1[,4]),3),'pmin',round(sd(df1[,4]),3),
              sep='')
print(string)

print(c(mean(5/(1+exp(-df1[,4]))),sd(5/(1+exp(-df1[,4])))))

#df1<-read.csv('C:/Users/oonar/Documents/idif/df_input8_WM')
#df1[,6]<-df1[,6]/1000/1000/length(times)
#df1[,7]<-length(times)*log(df1[,6])+2*5
#write.csv(df1,file=paste('df_input8_WM',sep=''),row.names=F)

df1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')
string<-paste(round(mean(df1[,6]),3),'pmin',round(sd(df1[,6]),3),' & ',
              round(mean(df1[(df1[,8]<1000000),8]),3),
              'pmin',round(sd(df1[(df1[,8]<1000000),8]),3),' & ',
              round(mean(df1[,7]),3),'pmin',round(sd(df1[,7]),3),
              sep='')
print(string)

#4*3*3=36
l<-5
print('Global')
for(i in c(6,8,7)){
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_global')
  y<-df1[,i]
  df1<-read.csv(paste('C:/Users/oonar/Documents/idif/df_input',l,'_global',sep=''))
  y1<-df1[,i]
  pval=wilcox.test(y,y1,paired=TRUE)$p.value
  if(pval*36<=0.05){
    print(c(i,pval*36))
  }
}
print('GM')
for(i in c(6,8,7)){
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_GM')
  y<-df1[,i]
  df1<-read.csv(paste('C:/Users/oonar/Documents/idif/df_input',l,'_GM',sep=''))
  y1<-df1[,i]
  pval=wilcox.test(y,y1,paired=TRUE)$p.value
  if(pval*36<=0.05){
    print(c(i,pval*36))
  }
}
print('WM')
for(i in c(6,8,7)){
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_WM')
  y<-df1[,i]
  df1<-read.csv(paste('C:/Users/oonar/Documents/idif/df_input',l,'_WM',sep=''))
  y1<-df1[,i]
  pval=wilcox.test(y,y1,paired=TRUE)$p.value
  if(pval*36<=0.05){
    print(c(i,pval*36))
  }
}

errorLine1<-rep(0,266)
errorLine2<-rep(0,266)
errorLine3<-rep(0,266)
errorLine4<-rep(0,266)
errorLine5<-rep(0,266)
relErrorLine1<-rep(0,266)
relErrorLine2<-rep(0,266)
relErrorLine3<-rep(0,266)
relErrorLine4<-rep(0,266)
relErrorLine5<-rep(0,266)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/oonar/Documents/idif/',studyNumbers[i],'_idif.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,8)
  input1<-c()
  input2<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input1[l]<-linearInter(times[l],t,df[1,])
    input2[l]<-linearInter(times[l],t,df[2,])
    organCurve[l]<-linearInter(times[l],t,df[8,])
  }
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_global')
  param<-c(df1[i,1],df1[i,2],-log(1/df1[i,3]-1))
  j<-df1[i,5]
  fit1<-fit_model(param,input1,j)
  #plot(times,organCurve,type='l')
  #points(times,fit1,type='l')
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input2_global')
  param<-c(df1[i,1],df1[i,2],-log(1/df1[i,3]-1))
  j<-df1[i,5]
  fit2<-fit_model(param,input2,j)
  #points(times,fit2,type='l')
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input3_global')
  param<-c(df1[i,1],df1[i,2],-log(1/df1[i,3]-1),-log(5/df1[i,4]-1))
  j<-df1[i,5]
  fit3<-fit_model(param,input1,j)
  #points(times,fit3,type='l')
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input4_global')
  param<-c(df1[i,1],df1[i,2],-log(1/df1[i,3]-1),-log(5/df1[i,4]-1))
  j<-df1[i,5]
  fit4<-fit_model(param,input2,j)
  #points(times,fit4,type='l')
  df1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')
  input_new<-fit_input(param_for_tau[i],input1)
  param<-c(df1[i,1],df1[i,2],-log(1/df1[i,3]-1))
  j<-df1[i,5]
  fit5<-fit_model(param,input_new,j)
  #points(times,fit5,type='l')
  errorLine1<-errorLine1+(organCurve/1000-fit1/1000)^2/100
  errorLine2<-errorLine2+(organCurve/1000-fit2/1000)^2/100
  errorLine3<-errorLine3+(organCurve/1000-fit3/1000)^2/100
  errorLine4<-errorLine4+(organCurve/1000-fit4/1000)^2/100
  errorLine5<-errorLine5+(organCurve/1000-fit5/1000)^2/100
  relErrorLine1<-relErrorLine1+abs(organCurve-fit1)/organCurve
  relErrorLine2<-relErrorLine2+abs(organCurve-fit2)/organCurve
  relErrorLine3<-relErrorLine3+abs(organCurve-fit3)/organCurve
  relErrorLine4<-relErrorLine4+abs(organCurve-fit4)/organCurve
  relErrorLine5<-relErrorLine5+abs(organCurve-fit5)/organCurve
}

plot(times,organCurve,type='l',
     ylab='Activity concentration (kBq/mL)',
     xlab='Time (s)',
     cex.lab=1.3,cex.axis=1.5,
     main='Fitted model curves',lwd=2,lty=1,col='steelblue1')
points(times,fit1,type='l',lwd=2,lty=1,col='blue')
points(times,fit2,type='l',lwd=2,lty=1,col='gray')
points(times,fit3,type='l',lwd=2,lty=2,col='blue')
points(times,fit4,type='l',lwd=2,lty=2,col='gray')
points(times,fit5,type='l',lwd=2)
legend('bottomright',
       legend=c('Real brain TAC',
                'Fitted from aorta IDIF',
                'Fitted from CCA IDIF',
                'Fitted from dispersion-corrected 
                aorta IDIF',
                'Fitted from dispersion-corrected 
                CCA IDIF',
                'Fitted from hybrid IDIF'),
       lty=c(1,1,1,2,2,1),
       col=c('steelblue1','blue','gray','blue','gray','black'),
       lwd=c(2,2,2,2,2,2),
       cex=1.3)

plot(times,errorLine2,type='l',
     ylab='MSE',
     xlab='Time (s)',
     cex.lab=1.3,cex.axis=1.5,
     main='MSE with respect to time',lwd=2,lty=1,col='gray')
points(times,errorLine1,type='l',lwd=2,lty=1,col='blue')
points(times,errorLine3,type='l',lwd=2,lty=2,col='blue')
points(times,errorLine4,type='l',lwd=2,lty=2,col='gray')
points(times,errorLine5,type='l',lwd=2)
legend('topright',
       legend=c('Aorta IDIF',
                'CCA IDIF',
                'Dispersion-corrected aorta IDIF',
                'Dispersion-corrected CCA IDIF',
                'Hybrid IDIF'),
       lty=c(1,1,2,2,1),
       col=c('blue','gray','blue','gray','black'),
       lwd=c(2,2,2,2,2),
       cex=1.3)

plot(times[2:266],relErrorLine2[2:266],type='l',
     ylab='MRE (%)',
     xlab='Time (s)',
     cex.lab=1.3,cex.axis=1.5,
     main='MRE with respect to time',lwd=2,lty=1,col='gray')
points(times[2:266],relErrorLine1[2:266],type='l',lwd=2,lty=1,col='blue')
points(times[2:266],relErrorLine3[2:266],type='l',lwd=2,lty=2,col='blue')
points(times[2:266],relErrorLine4[2:266],type='l',lwd=2,lty=2,col='gray')
points(times[2:266],relErrorLine5[2:266],type='l',lwd=2)
legend('topright',
       legend=c('Aorta IDIF',
                'CCA IDIF',
                'Dispersion-corrected aorta IDIF',
                'Dispersion-corrected CCA IDIF',
                'Hybrid IDIF'),
       lty=c(1,1,2,2,1),
       col=c('blue','gray','blue','gray','black'),
       lwd=c(2,2,2,2,2),
       cex=1.3)

df1<-read.csv('C:/Users/oonar/Documents/idif/df_input8_global')
y<-df1[,1]
shapiro.test(y)

#4+3+2+1=10
df1<-read.csv('C:/Users/oonar/Documents/idif/df_input4_global')
y<-df1[,1]
df1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')
y1<-df1[,1]
print(round(cor(y,y1,method='spearman'),3))
print(cor.test(y,y1,method='spearman')$p.value*10)

#3*8=24
y<-read.csv('C:/Users/oonar/Documents/idif/df_input4_global')[,1]
y1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')[,1]
pval<-wilcox.test(y,y1,paired=TRUE)$p.value
if(pval*24<0.05){pval*24}
y<-read.csv('C:/Users/oonar/Documents/idif/df_input4_GM')[,1]
y1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_GM')[,1]
pval<-wilcox.test(y,y1,paired=TRUE)$p.value
if(pval*24<0.05){pval*24}
y<-read.csv('C:/Users/oonar/Documents/idif/df_input4_WM')[,1]
y1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_WM')[,1]
pval<-wilcox.test(y,y1,paired=TRUE)$p.value
if(pval*24<0.05){pval*24}

y1<-read.csv('C:/Users/oonar/Documents/idif/df_input1_global')[,1]
y2<-read.csv('C:/Users/oonar/Documents/idif/df_input2_global')[,1]
y3<-read.csv('C:/Users/oonar/Documents/idif/df_input3_global')[,1]
y4<-read.csv('C:/Users/oonar/Documents/idif/df_input4_global')[,1]
y5<-read.csv('C:/Users/oonar/Documents/idif/df_input5_global')[,1]
boxplot(y1,y2,y3,y4,y5,
        col=c('white','lightblue','white','lightblue','white'),
        names=c('(A)','(B)','(C)','(D)','(E)'),
        lwd=2,cex.lab=1.3,cex.axis=1.5,
        ylab='CBF (mL/min/mL)')
y<-y5
sd(y)/mean(y)

df<-matrix(data=NA,nrow=length(studyNumbers),ncol=13)
df[,1]<-studyNumbers
indexes<-c()
for(i in 1:length(studyNumbers)){
  indexes<-c(indexes,as.numeric(substr(studyNumbers[i],7,10)))
}
library(readxl)
df1<-read_excel('D:/koveri/uusi_metadata_aivot.xlsx')
ids<-as.vector(df1$patient_id)[indexes]
for(i in 1:length(ids)){
  if(as.numeric(substr(ids[i],10,10))%%2==0){
    df[i,3]='Female'
  }else{df[i,3]='Male'}
}
for(i in 1:length(ids)){
  if(is.na(df[i,2])){
    df[i,2]=124-as.numeric(substr(ids[i],5,6))
  }
}
df[,4]<-as.numeric(as.vector(df1$weight)[indexes])
df[,5]<-as.numeric(as.vector(df1$height)[indexes])
df[,6]<-as.numeric(df[,4])/(as.numeric(df[,5])/100)^2
df[,7]<-as.numeric(as.vector(df1$'dose/lepo')[indexes])

i<-7
print(mean(na.omit(as.numeric(df[,i]))))
print(sd(na.omit(as.numeric(df[,i]))))
print(min(na.omit(as.numeric(df[,i]))))
print(max(na.omit(as.numeric(df[,i]))))

print(sum(df[,3]=='Female'))
print(sum(df[,3]=='Male'))

df1<-read.csv('C:/Users/oonar/Documents/idif/df_input5_WM')

print(mean(df1[,6])*265/266)
print(sd(df1[,6])*265/266)

print(mean(265*log(df1[,6]*265/266)+2*5))
print(sd(265*log(df1[,6]*265/266)+2*5))
