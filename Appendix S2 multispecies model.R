#parameters for generating random species
#these are representative values for growth, survival and canopy allometry from the 
#FORRU database on seedling planting trials: see http://www.forru.org/en/ for more details

survival=c(1.214, 0.781) #survival parameters
growth=c(-9.40, 0.354) #growth parameters
canopy=c(-1.108373,-1.511521 ) #allometry
#set random seed
set.seed(2)

#Function to return a randomly generated name for each species
getRandString<-
  function(n,len=10) {
    spvec<-rep(NA,times=n)
    
    for(i in 1:n) {
      spvec[i]<-(paste(sample(c(rep(0:9,each=5),LETTERS,letters),len,replace=TRUE)
                       ,collapse=''))}
    return(spvec)
  }

#Function to return a dataframe of randomly generated species
#nSP is number of species desitred
#seeds is number of seeds per each species
spGEN<-function(nSP,seeds) {
  Species<-as.character(getRandString(nSP))
  #draws survival from a beta distribution
  mu.c=-log(rbeta(nSP,shape1=exp(survival[1]),shape2=exp(survival[2])))
  #draws growth from a lognormal distribution
  g.c=rlnorm(nSP,meanlog=exp(growth[1]),sdlog=exp(growth[2]))
  #extreme growth values can lead to problems...this truncates growth to a reasonable range
  g.c[which(g.c>2)]<-runif(length(which(g.c>2)),min=0.001,max=2)
  #canopy allometry drawn from a log normal distribution
  c1=rlnorm(nSP,meanlog=exp(canopy[1]),sdlog=exp(canopy[2]))
  #assume understory mortality and growth are less than canopy values by some random amount
  mu.u<-mu.c+runif(nSP,min=-0,max=0.2)*mu.c
  mu.g<-mu.u+runif(nSP,min=0,max=0.5)*mu.u
  g.u<-g.c*runif(nSP,min=0.5,max=.95)
  g.g<-g.u*runif(nSP,min=0.3,max=.95)
  #other parameters (not used in this version of the model)
  c2=rep(1,times=nSP)
  D0=rep(0,times=nSP)
  f1=rep(600,times=nSP)
  f2=rep(0.1,times=nSP)
  stem.weight=rep(1,times=nSP)
  h1=rep(1,times=nSP)
  h2=rep(1,times=nSP)
  #global seed rain
  global<-rpois(lambda=seeds,n=nSP)
  all_yall<-data.frame(Species,mu.u,mu.g,mu.c,g.c,g.g,g.u,c1,c2,h1,h2,f1,f2,global,D0,stem.weight)
  return(all_yall)  
}
  

#use function to create a species frame
  spframe<-spGEN(5,5)
  
  #global parameters
  tsteps<-0.01 #discretized time
  new.dbh<-0 #starting dbh
  plot.area<-1000 #plot area
  grass.ht<-1 #grass height
 
#function to predict canopy cover as a function of species dataframe (datz) and time for a single
  #species
  CA_I<-function(datz,TIME,lagged=0) {
    
    D=exp((log(grass.ht/datz$h1[1]))/datz$h2[1])
    
    
    c.<-grass.ht/datz$g.g+lagged
    
    S<-(datz$global*exp((-datz$mu.g)*D/datz$g.g)*datz$c1)/plot.area
    
    m<-datz$mu.c
    
    g<-datz$g.c
    
    CA_t<-(S*(g+D*m-exp(-((TIME-c.)*m))*(g+D*m+g*(TIME-c.)*m)))/m^2
    
    return(CA_t)
  }
  
#function to optimize time until canopy closure for multiple species
  CA_all=function(spmat,TIME,lagged=0) {
  yqn=sum(CA_I(spmat,TIME,lagged=0))
  yasq=ifelse(yqn<0,0,yqn)
  yasqueen=ifelse(yasq>1,1,yasq)
    return(yasqueen)
  }
  
  
#simulate a list of species frames with different values of seed rain
#number of values of seed rain  
nsams=100
SPs=vector(mode="list",length=nsams)
#seed rain values:  
seseq=exp(seq(from=-5,to=1,length=nsams))
meanseseq=rep(NA,times=length(SPs))

  for(i in 1:length(SPs)){
    SPs[[i]]=spframe
    SPs[[i]]$global=rep(seseq[i],times=nrow(SPs[[i]]))+runif(n=nrow(SPs[[i]]),min=-0.1,max=0.1)
    SPs[[i]]$global=ifelse(SPs[[i]]$global<0,0,SPs[[i]]$global)
    
    meanseseq[i]=mean(SPs[[i]]$global)
    }
num.times=10000
#create empty vector to fill with time until canopy closure  
time.vecALL=rep(NA,times=length(seseq))
CAmat=matrix(NA,nrow=length(seseq),ncol=num.times)
#use optim to optimize time until canopy closure for each of the species frames in the list
  for(i in 1:length(seseq)) {
    for(j in 1:num.times){
CAmat[i,j]=CA_all(SPs[[i]],j)
if(CAmat[i,j]>=1){break}
          }
  }

find.first=function(x){return(which(x>=1)[1])}

canopyclosure=apply(CAmat,1,find.first)*tsteps

seedunits<-seseq/plot.area #seeds in units of meters squared

plot(canopyclosure~seedunits,pch=19,xlab=expression(paste("Seeds ","m"^-2," y"^-1)),
     ylab="Years to canopy closure",cex.axis=1.2,cex.main=1.3,cex.lab=1.3)
lines(canopyclosure[-which(is.na(canopyclosure)==T)]~seedunits[-which(is.na(canopyclosure)==T)],col="darkgreen")
grid()

