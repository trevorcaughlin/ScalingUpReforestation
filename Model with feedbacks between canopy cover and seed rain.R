source("https://raw.githubusercontent.com/trevorcaughlin/ScalingUpReforestation/master/reforestation ANAL FUNCTIONS.R")
source("https://raw.githubusercontent.com/trevorcaughlin/ScalingUpReforestation/master/reforestation%20SIMUL%20functions.R")

#randomly generated species representative of FORRU dataset
framio=structure(list(Species = "9Q3U1d118z", mu.u = 0.589538389965879, 
                      mu.g = 0.837140271057834, mu.c = 0.505895061260364, g.c = 1.91170481390366, 
                      g.g = 0.796126020106151, g.u = 0.995070897643558, c1 = 1.38001718595764, 
                      c2 = 1, h1 = 1, h2 = 1, f1 = 1000, f2 = 0.1, global = 5, 
                      D0 = 0, stem.weight = 1), .Names = c("Species", "mu.u", "mu.g", 
                                                           "mu.c", "g.c", "g.g", "g.u", "c1", "c2", "h1", "h2", "f1", "f2", 
                                                          "global", "D0", "stem.weight"), row.names = c(NA, -1L), class = "data.frame")

#initial conditions for cohort
cohort=structure(list(Species = "early", DBH = 1, CA = 1, Z = 1, N = 0, 
                      age = 1L, layer = "c", cohort.id = 1L), .Names = c("Species", 
                                                                         "DBH", "CA", "Z", "N", "age", "layer", "cohort.id"), row.names = 1L, class = "data.frame")


cohort$Species<-as.character(cohort$Species)
cohort$layer<-as.character(cohort$layer)
  
ori.dim<-dim(cohort)[1]
  
spframe<-framio
spframe$Species<-as.character(spframe$Species)
pars<-names(spframe)[-1]
  
  
  #global pars
  tsteps<-0.05
  new.dbh<-0
  plot.area<-1000
  grass.ht<-1
  
  timez<-500
  
  real.time<-timez*tsteps
  
tvec<-seq(from=tsteps,to=real.time,by=tsteps)  


spframe$global<-5
spframe$f1<-0

  crit<-grass.ht/spframe$g.g
  
BASE<-reTIME(spframe,type="NOframework")

spframe$f1<-1000
spframe$f2<-0.1

par(mfrow=c(1,3))


fram05<-reTIME(spframe,type="framework")
internal05<-reTIME(spframe,type="NOframework")
 
par(mfrow=c(2,2))
par(mar=c(3,3,1,3))
#0.6
  spframe$f2<-0.6

  
  fram05<-reTIME(spframe,type="framework")
  internal05<-reTIME(spframe,type="NOframework")
  
  plot(fram05[[2]]*100~tvec,type="l",xlab="Time (years)",lty=2,
       ylab="Canopy area",main=expression(paste("f"[2],"=0.6")),lwd=4,col="gray60",cex.main=1.4,cex.axis=1.3)
  lines(internal05[[2]]*100~tvec,lwd=4,lty=3,col="gray30")
  
  curve(CA_I(spframe,x)*100,from=crit,add=T,lwd=4)
  
  spframe$f1=1000
  grid()
  #1
  spframe$f2<-1
  
  
  fram05<-reTIME(spframe,type="framework")
  internal05<-reTIME(spframe,type="NOframework")
  
  plot(fram05[[2]]*100~tvec,type="l",xlab="Time (years)",lty=2,
       ylab="Canopy area",main=expression(paste("f"[2],"=0.6")),lwd=4,col="gray60",cex.main=1.4,cex.axis=1.3)
  lines(internal05[[2]]*100~tvec,lwd=4,lty=3,col="gray30")
  
  curve(CA_I(spframe,x)*100,from=crit,add=T,lwd=4)
  
  spframe$f1=1000
  grid()
  
  #1.6
  spframe$f2<-1.66666666666666

  
  fram05<-reTIME(spframe,type="framework")
  internal05<-reTIME(spframe,type="NOframework")
  
  plot(fram05[[2]]*100~tvec,type="l",xlab="Time (years)",lty=2,
       ylab="",main=expression(paste("f"[2],"=1.67")),lwd=4,col="gray60",cex.main=1.4,cex.axis=1.3)
  curve(CA_I(spframe,x)*100,add=T,lwd=4,from=crit)
  lines(internal05[[2]]*100~tvec,lwd=4,lty=3,col="gray30")
  grid()
 
  
  #10
  spframe$f2<-10
  
  
  fram05<-reTIME(spframe,type="framework")
  internal05<-reTIME(spframe,type="NOframework")
  
  plot(fram05[[2]]*100~tvec,type="l",xlab="Time (years)",lty=2,
       ylab="",main=expression(paste("f"[2],"=10")),lwd=4,col="gray60",cex.main=1.4,cex.axis=1.3)
  
  curve(CA_I(spframe,x)*100,add=T,from=crit,lwd=4)
  lines(internal05[[2]]*100~tvec,lwd=4,lty=3,col="gray30")
  grid()
  