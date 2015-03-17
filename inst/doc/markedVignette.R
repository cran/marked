### R code from vignette source 'markedVignette.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: markedVignette.Rnw:21-22
###################################################
  if(exists(".orig.enc")) options(encoding = .orig.enc)


###################################################
### code chunk number 2: markedVignette.Rnw:49-52
###################################################
# See note under Validation and Timing regarding rebuilding this vignette




###################################################
### code chunk number 3: markedVignette.Rnw:151-154
###################################################
if(length(grep("RMark",.packages()))!=0)detach("package:RMark") 
options(width=70)



###################################################
### code chunk number 4: markedVignette.Rnw:157-160
###################################################
library(marked) 
data(dipper)
model=crm(dipper)


###################################################
### code chunk number 5: markedVignette.Rnw:184-185
###################################################
model


###################################################
### code chunk number 6: markedVignette.Rnw:204-206
###################################################
model=cjs.hessian(model)
model


###################################################
### code chunk number 7: markedVignette.Rnw:220-225
###################################################
dipper.proc=process.data(dipper)
dipper.ddl=make.design.data(dipper.proc)
Phi.sex=list(formula=~sex)
model=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=Phi.sex),
          accumulate=FALSE)


###################################################
### code chunk number 8: markedVignette.Rnw:245-259
###################################################
dipper.proc=process.data(dipper)
dipper.ddl=make.design.data(dipper.proc)
fit.models=function()
{
  Phi.sex=list(formula=~sex)
  Phi.time=list(formula=~time)
  p.sex=list(formula=~sex)
  p.dot=list(formula=~1)
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=dipper.proc, ddl=dipper.ddl,
                      external=FALSE,accumulate=FALSE)
  return(results)
}
dipper.models=fit.models()


###################################################
### code chunk number 9: markedVignette.Rnw:264-265
###################################################
dipper.models


###################################################
### code chunk number 10: markedVignette.Rnw:274-276
###################################################
dipper.models[[1]]
dipper.models[["Phi.sex.p.dot"]]


###################################################
### code chunk number 11: markedVignette.Rnw:282-296
###################################################
dipper.proc=process.data(dipper)
dipper.ddl=make.design.data(dipper.proc)
fit.models=function()
{
  Phi.sex=list(formula=~sex)
  Phi.time=list(formula=~time)
  p.sex=list(formula=~sex)
  p.dot=list(formula=~1)
  cml=create.model.list(c("Phi","p"))
  results=crm.wrapper(cml,data=dipper.proc, ddl=dipper.ddl,
                      external=TRUE,accumulate=FALSE)
  return(results)
}
dipper.models=fit.models()


###################################################
### code chunk number 12: markedVignette.Rnw:299-301
###################################################
model=load.model(dipper.models[[1]])
model


###################################################
### code chunk number 13: markedVignette.Rnw:321-341
###################################################
data(dipper)
# Add a dummy weight field which are random values from 1 to 10
set.seed(123)
dipper$weight=round(runif(nrow(dipper),0,9),0)+1
# Add Flood covariate
Flood=matrix(rep(c(0,1,1,0,0,0),each=nrow(dipper)),ncol=6)
colnames(Flood)=paste("Flood",1:6,sep="")
dipper=cbind(dipper,Flood)
# Add td covariate, but exclude first release as a capture
# splitCH and process.ch are functions in the marked package
td=splitCH(dipper$ch)
td=td[,1:6]
releaseocc=process.ch(dipper$ch)$first
releaseocc=cbind(1:length(releaseocc),releaseocc)
releaseocc=releaseocc[releaseocc[,2]<nchar(dipper$ch[1]),]
td[releaseocc]=0
colnames(td)=paste("td",2:7,sep="")
dipper=cbind(dipper,td)
# show names
names(dipper)


###################################################
### code chunk number 14: markedVignette.Rnw:362-372
###################################################
# Process data 
dipper.proc=process.data(dipper)
# Create design data with static and time varying covariates
design.Phi=list(static=c("weight"),time.varying=c("Flood"))
design.p=list(static=c("sex"),time.varying=c("td"),
                        age.bins=c(0,1,20))
design.parameters=list(Phi=design.Phi,p=design.p)
ddl=make.design.data(dipper.proc,parameters=design.parameters)
names(ddl$Phi) 
names(ddl$p) 


###################################################
### code chunk number 15: markedVignette.Rnw:378-382
###################################################
Phi.sfw=list(formula=~Flood+weight)
p.ast=list(formula=~age+sex+td)
model=crm(dipper.proc,ddl,hessian=TRUE,
           model.parameters=list(Phi=Phi.sfw,p=p.ast))


###################################################
### code chunk number 16: markedVignette.Rnw:390-399
###################################################
newdipper=expand.grid(sex=c("Male","Female"),weight=1:10,Flood1=0,
      Flood2=1,Flood3=1,Flood4=0,Flood5=0,Flood6=0,td2=0,td3=c(0,1),
      td4=c(0,1),td5=c(0,1),td6=c(0,1),td7=c(0,1)) 
reals=predict(model,newdata=newdipper,se=TRUE) 
library(ggplot2)
reals$Phi$Flood=factor(reals$Phi$Flood,labels=c("Non-flood","Flood"))
ggplot(reals$Phi,aes(weight,estimate,ymin=lcl,ymax=ucl))+
       geom_errorbar(width=0.2)+geom_point()+geom_line()+
       xlab("\nWeight")+ylab("Survival\n")+facet_grid(Flood~.) 


###################################################
### code chunk number 17: markedVignette.Rnw:410-422
###################################################
data(dipper)
# Add Flood covariate
Flood=matrix(rep(c(0,1,1,0,0,0),each=nrow(dipper)),ncol=6)
colnames(Flood)=paste("Flood",1:6,sep="")
dipper=cbind(dipper,Flood)
design.parameters=list(Phi=list(time.varying="Flood"))
model.parameters=list(Phi=list(formula=~Flood),
                     p=list(formula=~time+sex))
MCMCfit=crm(dipper,model="probitCJS",
            model.parameters=model.parameters,
            design.parameters=design.parameters,
            burnin=1000,iter=5000) 


###################################################
### code chunk number 18: markedVignette.Rnw:429-431
###################################################
# beta estimates
MCMCfit


###################################################
### code chunk number 19: markedVignette.Rnw:434-436
###################################################
# real estimates
MCMCfit$results$reals


###################################################
### code chunk number 20: markedVignette.Rnw:771-777 (eval = FALSE)
###################################################
##  Phiprime=1-Fplus + Phi*Fplus
##  Phistar=t(apply(Phiprime,1,cumprod))
##  pprime=(1-Fplus)+Fplus*(C*p+(1-C)*(1-p))
##  pstar=t(apply(pprime,1,cumprod))
##  pomega=rowSums(L*M*Phistar*pstar)
##  lnl=sum(log(pomega))


