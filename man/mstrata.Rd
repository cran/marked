% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marked-package.R
\docType{data}
\name{mstrata}
\alias{mstrata}
\title{Multistrata example data}
\format{
A data frame with 255 observations on the following 2 variables.
\describe{ \item{ch}{a character vector containing the encounter history of
each bird with strata} \item{freq}{the number of birds with that capture
history} }
}
\description{
An example data set which appears to be simulated data that accompanies MARK
as an example analysis using the Multistrata model.
}
\details{
This is a data set that accompanies program MARK as an example for the
Multistrata model and is also in the RMark pacakge. Here I use it to show the 
3 ways models can be fitted to multistrata data. The model MSCJS is not run because it
requires ADMB or the exe constructed from ADMB which is not available if downloaded from CRAN.
}
\examples{
\donttest{
data(mstrata)
ms1=process.data(mstrata,model="MSCJS",strata.labels=c("A","B","C"))
ms1.ddl=make.design.data(ms1)
# following requires ADMB or the exe constructed from ADMB and links set for ADMB
# remove comments if you have admb
#mod1=try(crm(ms1,ms1.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#                                               p=list(formula=~time)),hessian=TRUE))
#mod1
# file.remove("multistate.std")
ms2=process.data(mstrata,model="hmmMSCJS",strata.labels=c("A","B","C"))
ms2.ddl=make.design.data(ms2)
# uses R/Fortran code with hmmMSCJS
mod2=crm(ms2,ms2.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
                                           p=list(formula=~time)),hessian=TRUE)
mod2
# strata.labels for MVMS models must be specified as a list because more than one variable
# can be used
ms3=process.data(mstrata,model="MVMSCJS",strata.labels=list(state=c("A","B","C")))
ms3.ddl=make.design.data(ms3)
ms3.ddl$delta$fix=1
# uses R/Fortran code with MVMSCJS
mod3=crm(ms3,ms3.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
                                           p=list(formula=~time)),hessian=TRUE)
mod3
# requires admb; remove comments if you have admb
#mod4=crm(ms3,ms3.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#                                           p=list(formula=~time)),hessian=TRUE,use.admb=TRUE)
#mod4
#file.remove("mvms.std")
# uses TMB with mvmscjs -remove comment to use
#mod5=crm(ms3,ms3.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum),
#                                          p=list(formula=~time)),hessian=TRUE,use.tmb=TRUE)
#mod5
unlink("*.cpp")
unlink("*.o")
unlink("symbols.rds")
}
}
\keyword{datasets}
