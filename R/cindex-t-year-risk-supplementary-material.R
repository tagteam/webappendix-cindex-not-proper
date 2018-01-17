##
## Packages and data 
##
require(data.table)
require(randomForestSRC)
require(survival)
require(rms)
require(riskRegression)
require(pec)
if (file.exists("./recc.RData")){
    recc <- get(load("./recc.RData"))
    setDT(recc)
    recc[,status:=1*(DF.Status!=0)]
    recc <- na.omit(recc)
    recc[,time:=DF.Months/12]
}else{
    stop("You need to download the recc data and save the file recc.RData in the current working directory.")
}
data(pbc,package="survival")
setDT(pbc)
pbc <- na.omit(pbc[,.(time,status,edema,age,bili,protime,albumin)])
pbc[,time:=time/365.25]
pbc[,log.bili:=log(bili)]
pbc[,log.protime:=log(protime)]
pbc[,log.albumin:=log(albumin)]
data(follic) 
setDT(follic)
data(cost,package="pec")
setDT(cost)
cost[,time:=time/365.25]
data(GBSG2,package="pec")
setDT(GBSG2)
GBSG2[,time:=time/365.25]
##
## Table 1
## 
ccc <- function(formula,data,horizon,tau){
    #
    # 0. fit Cox model and evaluate risk prediction at prediction time horizon
    #
    fit <- coxph(formula,data=data,x=TRUE,y=TRUE)
    S <- fit$y
    prediction <- predictRisk(fit,newdata=data,times=horizon)
    #
    # 1. IPCW estimate of AUC(t)
    # 
    AUC.IPCW <- Score(list(cox=prediction),
                      formula=formula,
                      data=data,
                      times=horizon,
                      metrics="AUC",
                      split.method="none",
                      conf.int=FALSE,
                      cens.model="km")$AUC$score[model=="cox",AUC]
    # 
    # 2. IPCW estimate of C(t,tau,infty)
    C.IPCW <- pec::cindex(list(cbind(1-prediction,1-prediction)),
                          formula,
                          data=data,
                          eval.times=c(horizon,tau))$AppC[[1]]
    # 
    # 3. Compute Harrell's C-index with outcome stopped at time horizon 
    #
    harrellC.t <- rcorr.cens(x=1-prediction,S=stopTime(S,horizon))["C Index"]
    # 
    # 4. Compute C-index from coxph
    #
    C.coxph <- summary(fit)$concordance["C"]
    # 
    # 5. Collect results output
    #
    output <- data.table(matrix(round(100*c(AUC.IPCW,C.IPCW,harrellC.t,C.coxph),1),nrow=1))
    output <- cbind("Data set"=as.character(substitute(data)),output,round(tau,1))
    setnames(output,c("Data set",
                      "\\widehat{AUC}_{\\mathrm{IPCW}}(t)",
                      "\\widehat{C}_{\\mathrm{IPCW}}(t)",
                      "\\widehat{C}_{\\mathrm{IPCW}}(\\tau)",
                      "\\widehat{C}(t)",
                      "\\widehat{C}(\\tau)",
                      "\\tau (yrs)"))
    output
}
table1 <- rbindlist(
    list(ccc(formula=Surv(time,status!=0)~edema+age+log.bili+log.protime+log.albumin,
             data=pbc,
             horizon=5,
             tau=max(pbc$time)),
         ccc(Surv(time,status!=0)~age+hgb+clinstg+ch,
             data=follic,
             horizon=5,
             tau=max(follic$time)),
         ccc(Surv(time, cens!=0) ~ horTh+age+menostat+tsize+tgrade+pnodes+progrec+estrec,
             data=GBSG2,
             horizon=5,
             tau=max(GBSG2$time)),
         ccc(Surv(time, status!=0) ~ age+sex+hypTen+ihd+prevStroke+cholest+atrialFib+strokeScore,
             data=cost,
             horizon=5,
             tau=max(cost$time)),
         ccc(Surv(time,status)~Risk.5yr,
             data=recc,
             horizon=5,
             tau=max(recc$time))))
table1
##
## Figure 1
##
predtimes <- seq(from=2,to=10,by=0.5)
## calculate AUC_IPCW(t)
cox.protime <- coxph(Surv(time,status!=0)~protime,data=pbc,x=1L,y=1L)
protime.risk5 <- predict(cox.protime,newdata=pbc,times=5)
AUC.t <- Score(list(protime.risk5),formula=Hist(time,status!=0)~1,
               cause=1,cens.model="marginal",
               contrasts=FALSE,nullModel=FALSE,
               data=pbc,
               times=predtimes,metric="auc")$AUC$score
## calculate C_IPCW(tau)
ipcw.C.tau <- pec::cindex(list(matrix(-protime.risk5,
                                      ncol=1,nrow=nrow(pbc),byrow=FALSE)),
                          formula=Surv(time,status!=0)~protime,
                          data=pbc,eval.times=max(pbc$time))$AppC
## calculate C_IPCW(t)
ipcw.C.t <- pec::cindex(list(matrix(-protime.risk5,ncol=length(predtimes),
                                    nrow=nrow(pbc),byrow=FALSE)),
                        formula=Surv(time,status!=0)~protime,
                        data=pbc,eval.times=predtimes)$AppC
par(mfrow=c(1,1))
par(mai=c(1,1,0.5,0.5))
AUC.t[,{
    plot(times,AUC,type="b",ylim=c(0.4,1),xlim=range(predtimes),
         lwd=2,ylab="Concordance (%)",xlab="time t (years)",axes=FALSE)
    axis(1, at=c(2,4,6,8,10))
    myxat <- c(0.4,0.5,0.627,0.7,0.9,1)
    axis(2,las=2, at=myxat,labels=myxat*100)
    abline(h=0.5,lty=2)}]
lines(predtimes,ipcw.C.t$matrix,lty=1,pch=3,type="b",lwd=2,col="gray40")
lines(predtimes,rep(ipcw.C.tau$matrix,length(predtimes)),lty=1,pch=2,type="b",lwd=2,col="gray60")
legend("topright",c(expression(AUC[t]),expression(C[t]),expression(C[tau])),
       lwd=2,lty=c(1,1),pch=c(1,3,2),col=c("black","gray40","gray60"),ncol=1,bty='n')

##
## Figure 2
##
scale <- 1
tau <- scale+0.2 
sdZ <- 3
meanZ <- 0
theq <- c(0.25,0.4,0.5,0.6,0.75)
genData <- function(n,tau){    
    Z <- rnorm(meanZ,sdZ,n=n)
    expZ <- exp(Z)
    time <- rweibull(n=n,shape=1/expZ,scale=scale)
    status <- rep(1,n)
    status[T>tau] <- 0
    time <- pmin(tau,time)
    data.frame(time=time,status=status,expZ=expZ)
}
set.seed(39)
NMC <- 1000
cindex <- sapply(1:NMC,function(b){
    d <- genData(1000,tau=tau)
    pec::cindex(list(expZ=matrix(-d$expZ)),data=d,formula=Surv(time,status)~1,eval.times=tau)$AppCindex$expZ
})
cindex.t <- sapply(1:NMC,function(b){
    d <- genData(1000,tau=tau)
    pec::cindex(list(expZ=matrix(-d$expZ)),data=d,formula=Surv(time,status)~1,eval.times=1)$AppCindex$expZ
})
shape.values <-  exp(qnorm(p=theq,mean=meanZ,sd=sdZ))
par(mfrow=c(1,2))
# Panel A 
plot(0,0,xlim=c(0,tau),ylim=c(0,1),type="n",xlab="Time (s)",
     ylab=expression(paste(F[s],"(Z)")),axes=FALSE)
for (i in 1:length(shape.values)){
    z <- shape.values[i]
    curve(1-pweibull(x,shape=1/z,scale=scale,lower.tail=FALSE), 
          from=0,to=tau,lty=i,add=TRUE,lwd=2)
}
prodlim::PercentAxis(2,at=seq(0,1,0.25),las=2)
axis(1,at=c(0,1,tau),labels=c("0","t",expression(tau)),xlim=c(0,2.5))
legend(x=1.6,y=1.3,bty="n",
       legend=paste("Z=",round(log(1/shape.values),2),sep=""),       
       lwd=2,lty=1:length(shape.values),xpd=NA,ncol=3)
abline(v=1,lty=3,lwd=3,col="gray50")
text(x=1,y=1.2,xpd=NA,paste0("c-index of Cox-model based on Z= ",round(100*mean(cindex),1),"%"))
text(x=1,y=1.1,xpd=NA,expression(paste("c-index of ",F[t],"(Z)=50%")))
# Panel B
plot(0,0,xlim=c(0,tau),ylim=c(0,6),type="n",xlab="Time (s)",ylab=expression(paste(lambda,"(s|Z)")),axes=FALSE)
for (i in 1:length(shape.values)){
    z <- shape.values[i]
    curve(dweibull(x,shape=1/z,scale=scale)/pweibull(x,shape=1/z,scale=scale,lower.tail=FALSE), 
          from=0,to=tau,lty=i,add=TRUE,lwd=2)
}
axis(2)
axis(1,at=c(0,scale,tau),labels=c("0","t",expression(tau)))
abline(v=1,lty=3,lwd=3,col="gray50")
