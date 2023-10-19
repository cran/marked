## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----comment="", prompt=TRUE, message=FALSE-----------------------------------
library(marked) 
library(ggplot2)
data(dipper)
model=crm(dipper)

## ----comment="",prompt=TRUE---------------------------------------------------
model

## ----comment="",prompt=TRUE,tidy=FALSE----------------------------------------
model=cjs.hessian(model)
model

## ----comment="",prompt=TRUE,tidy=FALSE,results='hide'-------------------------
dipper.proc=process.data(dipper)
dipper.ddl=make.design.data(dipper.proc)
Phi.sex=list(formula=~sex)
model=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=Phi.sex),
          accumulate=FALSE)

## ----comment="",prompt=TRUE,tidy=FALSE,results='hide',message=FALSE-----------
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

## ----comment="",prompt=TRUE,tidy=FALSE----------------------------------------
dipper.models

## ----comment="",prompt=TRUE,tidy=FALSE----------------------------------------
dipper.models[[1]]
dipper.models[["Phi.sex.p.dot"]]

## ----comment="",prompt=TRUE,tidy=FALSE,results='hide',echo=FALSE, message=FALSE----
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
                      external=TRUE,accumulate=FALSE, replace=TRUE)
  return(results)
}
dipper.models=fit.models()

## ----comment="",prompt=TRUE,tidy=FALSE----------------------------------------
model=load.model(dipper.models[[1]])
model

## ----comment="",prompt=TRUE,tidy=FALSE----------------------------------------
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


## ----comment="",prompt=TRUE,tidy=FALSE,results='hide'-------------------------
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

## ----comment="",prompt=TRUE,tidy=FALSE,results='hide'-------------------------
Phi.sfw=list(formula=~Flood+weight)
p.ast=list(formula=~age+sex+td)
model=crm(dipper.proc,ddl,hessian=TRUE,
           model.parameters=list(Phi=Phi.sfw,p=p.ast))

## ----fig=TRUE,comment="",prompt=TRUE,tidy=FALSE-------------------------------
newdipper=expand.grid(sex=c("Female","Male"),weight=1:10,Flood1=0,
      Flood2=1,Flood3=1,Flood4=0,Flood5=0,Flood6=0,td2=0,td3=c(0,1),
      td4=c(0,1),td5=c(0,1),td6=c(0,1),td7=c(0,1)) 
reals=predict(model,newdata=newdipper,se=TRUE)
reals$Phi$Flood=factor(reals$Phi$Flood,labels=c("Non-flood","Flood"))
ggplot(reals$Phi,aes(weight,estimate,ymin=lcl,ymax=ucl))+
      geom_errorbar(width=0.2)+geom_point()+geom_line()+
      xlab("\nWeight")+ylab("Survival\n")+facet_grid(Flood~.)   

## ----comment="",prompt=TRUE,tidy=FALSE, results='hide'------------------------
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

## ----comment="",prompt=TRUE---------------------------------------------------
# beta estimates
MCMCfit

## ----comment="",prompt=TRUE---------------------------------------------------
# real estimates
MCMCfit$results$reals

## ----comment="",prompt=TRUE, echo=FALSE---------------------------------------
library(kableExtra)

tab <- rbind(
  c("1", " 1"," 2"," 3"," 7"," 8"," 9"),
  c("2", NA, " 4"," 5", NA, "10","11"),
  c("3", NA, NA, " 6", NA, NA, " 12")
)
colnames(tab) <- c("Cohort", "2", "3", "4", "2", "3", "4")

options(knitr.kable.NA = '')
knitr::kable(tab, "html", row.names = FALSE, align = c("l",rep("r",6)), booktabs=TRUE) %>% 
  kableExtra::kable_styling(full_width = FALSE) %>% 
  kableExtra::add_header_above(c(" "=1, "Time" = 3, "Time" = 3), line = FALSE) %>% 
  kableExtra::add_header_above(c(" "=1, "$\\phi$" = 3, "$p$" = 3)) %>% 
  kableExtra::column_spec(c(1,4), border_right = TRUE) %>% 
  kableExtra::column_spec(1:7, width_min="3em", monospace = TRUE)


## ----eval=FALSE---------------------------------------------------------------
#   Phiprime = 1-Fplus + Phi*Fplus
#   Phistar = t(apply(Phiprime,1,cumprod))
#   pprime = (1-Fplus)+Fplus*(C*p+(1-C)*(1-p))
#   pstar = t(apply(pprime,1,cumprod))
#   pomega = rowSums(L*M*Phistar*pstar)
#   lnl = sum(log(pomega))

