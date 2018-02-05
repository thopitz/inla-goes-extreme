#R code for fitting the models in the following manuscript:
#Opitz, Huser, Bakka, Rue (2018) - 
#"INLA goes extreme: Bayesian tail regression for the estimation of high spatio-temporal quantiles"
#Certain hyperparameters will here be fixed to the optimal values found through cross-validation in the above paper.
#Currently, the generalized Pareto distribution is implemented in the testing version of R-INLA (see www.r-inla.org).
library(INLA)

#set the working directory containing the data files
setwd(".")

#other settings
nthreads=2 #number of threads (=CPU cores) that INLA may use in parallel

#1) Import data ####
data=read.csv("precip.csv",sep=",",header=TRUE,colClasses=c("numeric","numeric","character","numeric"))
coord=read.csv("coord.csv",sep=",",header=TRUE,colClasses=c("numeric","numeric","numeric"))
stat.vec=data[,2] #vector of stations
prcp.vec=data[,4] #vector of precipitation data
#get latitude, longitude:
#(notice: longitude has been shifted by some degrees to "hide" the true station locations in the EVA challenge)
lat=coord[,2]
long=coord[,3]
lonlat=cbind(long,lat)
metric=cbind(long*67.7,lat*111.2) #transform (approximately) to Euclidean coordinates
#plot(metric,asp=1,xlab="x (km)",ylab="y (km)",pch=19,cex=.5)

date.vec=strptime(data[,3],"%m/%d/%Y") #vector of dates (in POSIXlt format)
date.vec2=as.Date(date.vec) #vector of dates (in Date class)
all.stats=unique(stat.vec) 
all.dates=strptime(seq(min(date.vec2),max(date.vec2),by=1),"%Y-%m-%d")
#put data in matrix format (columns correspond to stations)
prcp.mat <- matrix(nrow=length(all.dates),ncol=length(all.stats))
for(j in 1:length(all.stats)){
  prcp.j=prcp.vec[stat.vec==all.stats[j]]
  date.j=date.vec[stat.vec==all.stats[j]]
  ind=match(date.j,all.dates)
  prcp.mat[ind,j]=prcp.j
}
rownames(prcp.mat)=as.character(all.dates)
colnames(prcp.mat)=all.stats
#create vector of station ids
station.id=rep(unique(stat.vec)[1:ncol(prcp.mat)],each=nrow(prcp.mat))
#some data cleaning (see Opitz et al. 2018 for explanation)
prcp.mat[which(as.numeric(all.dates) > as.numeric(all.dates[all.dates$year==77 & all.dates$mon==8 & all.dates$mday==20]) & as.numeric(all.dates) < as.numeric(all.dates[all.dates$year==78 & all.dates$mon==8 & all.dates$mday==20])),c(4,5,7,15,19,22,24,31,34)]=NA
prcp.mat[which(as.numeric(all.dates) < as.numeric(all.dates[all.dates$year==95 & all.dates$mon==4 & all.dates$mday==1])),28]=NA
prcp.mat[,1]=NA

#remove some of the objects that we will not need any more
rm(date.vec,date.vec2)

#2) Construct model components and the regression formula ####

#for an (approximately) "weekly" effect, we divide the year into 52 blocks:
days=all.dates$yday+1
days=rep(days,ncol(prcp.mat))
#select "fractional" days to subdivide the year into weeks (for weekly effect) 
weeks=inla.group(days,method="cut",n=floor(365/7))
#derive the corresponding integer index of the weeks
idx.weeks=match(weeks,sort(unique(weeks)))


#define Matérn model with prefixed range ####
matrange = 50 #fix the Matérn effective range parameter
nu=1 #fix the smoothness parameter
kappa=sqrt(2*nu)/matrange #SPDE parametrization
#fix the Matérn correlation matrix for the 40 stations:
dist=as.matrix(dist(metric))
cormat=as.matrix(2^(1-nu)*(kappa*dist)^nu*besselK(dist*kappa,nu)/gamma(nu))
diag(cormat)=1
prec=solve(cormat) #calculate the corresponding precision matrix
#use PC prior
hyper.pc = list(prec=list(prior="pc.prec",param=c(2,0.01)))

#write inla model formula
form=y~-1+intercept+f(station, model = "generic0", Cmatrix=prec, hyper=hyper.pc,constr=TRUE)+f(week,model="rw2",cyclic=TRUE,hyper=list(prec=list(initial=log(1/.01^2),fixed=TRUE)),constr=TRUE)

#NA value indicators in vector form
is.na=as.vector(is.na(prcp.mat))

#3) Model fitting with INLA####

#3.1) Gamma model for positive precipitation ####

#positive precipitation vector y.prcp
y.prcp=as.numeric(prcp.mat)[prcp.mat>0]
y.prcp=na.omit(y.prcp) #remove NAs
#indicator vector for observed positive precipitation
is.prcp=as.vector(prcp.mat>0)
is.prcp[is.na(is.prcp)]=FALSE

#prepare data frame for observed positive precipiation:
data.inla=data.frame(intercept=1,y=y.prcp,week=weeks[is.prcp&!is.na],idx.week=idx.weeks[is.prcp&!is.na],station=station.id[is.prcp&!is.na])
#add NA values to response for values to predict (52 weeks times 40 stations) 
y.prcp.pred=c(y.prcp,rep(NA,nrow(metric)*52))
week.pred=c(weeks[is.prcp&!is.na],rep(sort(unique(weeks)),each=nrow(metric)))
idx.week.pred=c(idx.weeks[is.prcp&!is.na],rep(1:52,each=nrow(metric)))
station.pred=c(station.id[is.prcp&!is.na],rep(1:nrow(metric),52))
#prepare data frame for prediction of positive precipitation over one year
data.inla.pred=data.frame(intercept=1,y=y.prcp.pred,week=week.pred,idx.week=idx.week.pred,station=station.pred)

#fit gamma model to positive precipitation
fit=inla(form,
         data=data.inla.pred,
         family="gamma",
         control.family=list(hyper=list(theta=list(prior="loggamma",param=c(2,2),initial=log(1)))),
         control.predictor=list(link=1),
         control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         verbose=T,
         num.threads=nthreads
)
save(fit,file="fitprcp.Rdata")


#3.2) Calculate week-station dependent threshold and exceedances based on gamma model ####
pthresh=0.92
pprcp=matrix(pthresh,ncol=nrow(metric),nrow=52)

#threshold from gamma model for positive part
load(file="fitprcp.Rdata")
tmp=fit$summary.fitted.values$mean
prcpfitted=matrix(tmp[(length(tmp)-nrow(metric)*52+1):length(tmp)],ncol=nrow(metric),byrow=TRUE)
shape.gamma=fit$summary.hyperpar$mean[1]
scale.gamma=prcpfitted/shape.gamma
#use fitted gamma parameters for threshold calculation
thr.gamma=qgamma(pprcp,scale=scale.gamma,shape=shape.gamma)

#calculate exceedance above the threshold (for generalized Pareto fit)
y.exc=as.numeric(prcp.mat)
for(i in 1:length(y.exc)){
  thr.i=thr.gamma[idx.weeks[i],station.id[i]]
  y.exc[i]=ifelse(y.exc[i]>thr.i,(y.exc[i]-thr.i)/prcpfitted[idx.weeks[i],station.id[i]],NA)
}
is.exc=!is.na(y.exc)
y.exc=na.omit(y.exc)


#3.3) Bernoulli model for exceedance probability ####

#(proceed as before to prepare response and data for INLA)
#response with observed exceedance indicators
y.pexc=is.exc
y.pexc[is.na]=NA
y.pexc=na.omit(y.pexc)
#add NA values to response for values to predict (52 weeks times 40 stations) 
y.pexc.pred=c(y.pexc,rep(NA,nrow(metric)*52))
week.pred=c(week=weeks[!is.na],rep(sort(unique(weeks)),each=nrow(metric)))
idx.week.pred=c(idx.weeks[!is.na],rep(1:52,each=nrow(metric)))
station.pred=c(station.id[!is.na],rep(1:nrow(metric),52))
#prepare data frame for prediction of positive precipitation over one year
data.inla.pred=data.frame(intercept=1,y=as.integer(y.pexc.pred),week=week.pred,idx.week=idx.week.pred,station=station.pred)

#fit Bernoulli model for exceedance probability
fit=inla(form,
         data=data.inla.pred,
         family="binomial",
         Ntrials=rep(1,length(y.pexc.pred)),
         control.predictor=list(link=1),
         control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         verbose=T,
         num.threads=nthreads
)
save(fit,file="fitprobexcess.Rdata")


#3.4) Generalized Pareto model for exceedance value ####

#add NA values to response for values to predict (52 weeks times 40 stations) 
y.exc.pred=c(y.exc,rep(NA,nrow(metric)*52))
week.pred=c(weeks[is.exc],rep(sort(unique(weeks)),each=nrow(metric)))
idx.week.pred=c(idx.weeks[is.exc],rep(1:52,each=nrow(metric)))
station.pred=c(station.id[is.exc],rep(1:nrow(metric),52))
#prepare data frame for prediction of positive precipitation over one year
data.inla.pred=data.frame(intercept=1,y=y.exc.pred,week=week.pred,idx.week=idx.week.pred,station=station.pred)

#fit generalized Pareto model for exceedance values
#for the generalized Pareto distribution, the linear predictor is linked to the quantile (of probability qprob)
qprob=0.5 #here, we take the median
fit=inla(form,
         data=data.inla.pred,
         family = "gp",
         control.family = list(control.link = list(quantile = alpha)),
         control.predictor=list(link=1),
         control.inla=list(strategy="simplified.laplace",int.strategy="ccd"),
         verbose=T,
         num.threads=nthreads
)
save(fit,file="fitexcess.Rdata")

#4) Some plots of the fitted models ####
load(file="fitexcess.Rdata") #for the GP excess model,
#load(file="fitprcp.Rdata") #for the gamma model for positive precipitation
#load(file="fitprobexcess.Rdata") #for the Bernoulli model for exceedance indicators
summary(fit)

#plot the weekly effect with credible intervals
tmp=fit$summary.random$week[,c("mean","0.025quant","0.975quant")]
ylim=range(tmp)
plot(1:52, fit$summary.random$week$mean,type="l",lwd=3,xlab="week",ylab="linear predictor",ylim=ylim)
lines(1:52, fit$summary.random$week$`0.025quant`,lwd=3,col="blue")
lines(1:52, fit$summary.random$week$`0.975quant`,lwd=3,col="blue")
lines(c(0,360),c(0,0),lty=2,lwd=3,col="gray50")

#plot the spatial effect for observed and predicted stations
library(fields)
cols=tim.colors(n=512)
tmp=seq(-max(abs(fit$summary.random$station$mean))-0.001,max(abs(fit$summary.random$station$mean))+0.001,length=512)
col=c()
for(i in 1:length(fit$summary.random$station$mean)){
  col[i]=cols[which.min(abs(fit$summary.random$station$mean[i]-tmp))]
}
par(mgp=c(2,1,0),mar=c(3.1,3.1,0.2,4.5))
plot(long,lat,pch=20,xlab="Longitude (Shifted)",ylab="Latitude (Shifted)",xlim=c(min(long)-0.1,max(long)+0.1),ylim=c(min(lat)-0.05,max(lat)+0.05),asp=1.3,cex.lab=1.3,cex.axis=1.3,col=col,cex=1.5)
image.plot(x=sort(long),y=sort(lat),z=matrix(c(min(tmp),max(tmp)),nrow=length(long),ncol=length(lat)),legend.only=TRUE,legend.mar=4.5)

#5) Calculate high quantile predictions (for EVA 2017 challenge) ####
ptarget=0.998 #target probability for predicted quantiles
load(file="fitexcess.Rdata")
summary(fit)
tmp=match("Shape", substr(rownames(fit$summary.hyperpar),1,5))
xi.est=fit$summary.hyperpar[tmp,"mean"]
tmp=fit$summary.fitted.values$mean
fitted.exc=matrix(tmp[(length(tmp)-52*40+1):length(tmp)],nrow=52,byrow=TRUE)
#calculate GP scale from fitted GP quantile:
sigma.est=prcpfitted*xi.est*fitted.exc/((1-qprob)^{-xi.est}-1)

load(file="fitprobexcess.Rdata")
tmp=fit$summary.fitted.values$mean
fitted.pexc=matrix(tmp[(length(tmp)-52*40+1):length(tmp)],nrow=52,byrow=TRUE)
pgpd.est=1-(1-ptarget)/fitted.pexc
#calculated estimated quantiles (per week and station):
qtarget.est=thr.gamma+sigma.est/xi.est*((1-pgpd.est)^{-xi.est}-1)

#transform weekly estimates to (approximate) monthly estimates (as required for the EVA2017 challenge):
day2month=function(days){
  monthlen=c(31,28.25,31,30,31,30,31,31,30,31,30,31.75)
  intv=c(0,cumsum(monthlen))
  findInterval(days,intv,all.inside=TRUE)
}

months_from_weeks=day2month(sort(unique(weeks)))
fun=function(vals){
  aggregate(vals,by=list(w2m=months_from_weeks),FUN=mean)$x
}
#final estimates of stationwise monthly quantiles:
qtargetm.est=apply(qtarget.est,2,fun)
