
MicMen<-function(val,sat,halfmax) {
  f1<-sat
  f2<-halfmax
  curve((x*f1)/(f2+x))
  return((val*f1)/(f2+val))
}

Sigmoid<-function(val,f1,f2) {
  #f1<-sat
  #f2<-halfmax
  #curve(1/(f1+exp(f2*x)))
  return(1/(f1+exp(f2*val)))
}

plaw<-function(val,a,b) {
  curve(a*x^b)
  return(a*val^b)
  
}

  f.N_t = function(mu,n0,tsteps) n0*exp(-tsteps*mu)          #mortality
  f.DBH_t = function(g,dbh,tsteps) {
  gro<-dbh+tsteps*g
  return(gro) }


  f.CA = function(dbh,c1,c2) c1*dbh^c2            #canopy area
  f.Z = function(dbh,h1,h2) h1*dbh^h2              #height
  f.recruitment<-function(f,dbh) return(rep(f,times=length(dbh)))

  f.Sp=function(global,f1,time) {return(global+f1*time^2)}

  ###MICMEN
  ###MICMEN
  ###MICMEN
  ###MICMEN
  ###MICMEN
  ###MICMEN
  ###MICMEN
  ###MICMEN
  
f.S=function(global,N,dbh,f1,f2,c1,c2) {
  CA<-N*(f.CA(dbh,c1,c2))/plot.area
  #print(CA)
   return(global+sum(f1*CA^f2))}
  
  
  f.Sall=function(global,CA.all,f1,f2) {return(global+f1*CA.all^f2)}


  f.layers = function(cohort) {
      z = cohort[,c('N','Z','CA')]
      z$t_ca = z$N/plot.area*z$CA  #10000#

      ord.coh = order(z$Z,decreasing=T)

      ca = cumsum(z$t_ca[ord.coh])

     cohort$layer[which(cohort$Z<grass.ht)]<-"g"

     cohort$layer[which(cohort$Z>grass.ht)]<-"c"
     if(any(cohort$Z==grass.ht)) {
     tp<-which(cohort$Z==grass.ht)

     tp.g<-cohort[tp,]
     tp.c<-cohort[tp,]

     tp.g$N<-tp.g$N/2
     tp.c$N<-tp.c$N/2

     tp.g$layer<-"g"
     tp.c$layer<-"c"

     cohort2<-cohort[-tp,]

     cohort3<-rbind(cohort2,tp.g,tp.c)
     return(cohort3)
     }
     else{
     return(cohort)
     }
     }



    vitals2<-function(cohort,mu.g,mu.c,g.g,g.c,c1,c2,h1,h2,global,f1,f2,time,CA.all,type) {
    cohort$age<-cohort$age+1
    #reproduction before anything else

    new.cohort<-data.frame(matrix(NA,nrow=1,ncol=length(cohort[1,])))
    colnames(new.cohort)<-colnames(cohort)

    new.cohort$age<-0

    new.cohort$DBH<-new.dbh

    new.cohort$CA<-f.CA(new.dbh,c1,c2)

    new.cohort$Z<-f.Z(new.dbh,h1,h2)

    newz<-new.cohort$Z

     new.cohort$layer<-NA

    #new.cohort$N<-global

    new.cohort$cohort.id<-time+ori.dim

  ###############NO canopy
    if(length(which(cohort$layer=="c"))==0) {
        canopy<-subset(cohort,cohort$layer=="g")

        canopy$DBH<-f.DBH_t(g.g,canopy$DBH,tsteps)

        canopy$N<-f.N_t(mu.g,canopy$N,tsteps)

        updated.cohort<-canopy

        updated.cohort$CA<-rep(0,times=length(updated.cohort$CA))

        updated.cohort$Z<-f.Z(updated.cohort$DBH,h1,h2)

        updated.cohort<-updated.cohort[order(updated.cohort$Z,decreasing=T),]
        
        
        if(type=="framework") {
        new.cohort$N<-f.Sall(global,CA.all,f1,f2)*tsteps
    } 
    else{
        new.cohort$N<-global*tsteps
    }
        #print("UNDERSTORY ONLY")
        return(rbind(updated.cohort,new.cohort))

        }

        ############################NO understory
          if(length(which(cohort$layer=="g"))==0) {
        canopy<-subset(cohort,cohort$layer=="c")

        canopy$DBH<-f.DBH_t(g.c,canopy$DBH,tsteps)

        canopy$N<-f.N_t(mu.c,canopy$N,tsteps)

        updated.cohort<-canopy

        updated.cohort$CA<-f.CA(updated.cohort$DBH,c1,c2)

        updated.cohort$Z<-f.Z(updated.cohort$DBH,h1,h2)

        updated.cohort<-updated.cohort[order(updated.cohort$Z,decreasing=T),]

        
        if(type=="framework") {
          new.cohort$N<-f.Sall(global,CA.all,f1,f2)*tsteps
        } 
        else{
        
        new.cohort$N<-f.S(global,canopy$N,canopy$DBH,f1,f2,c1,c2)*tsteps
        }     
        #print("CANOPY ONLY")
        return(rbind(updated.cohort,new.cohort))
        }



    under.grass<-subset(cohort,cohort$layer=="g")
    canopy<-subset(cohort,cohort$layer=="c")

      
  if(type=="framework") {
    new.cohort$N<-f.Sall(global,CA.all,f1,f2)*tsteps
  } 
  else{
      new.cohort$N<-f.S(global,canopy$N,canopy$DBH,f1,f2,c1,c2)*tsteps
}
    under.grass$DBH<-f.DBH_t(g.g,under.grass$DBH,tsteps)
    canopy$DBH<-f.DBH_t(g.c,canopy$DBH,tsteps)

    under.grass$N<-f.N_t(mu.g,under.grass$N,tsteps)
    canopy$N<-f.N_t(mu.c,canopy$N,tsteps)


    under.grass$CA<-rep(0,times=length(under.grass$CA))
    canopy$CA<-f.CA(canopy$DBH,c1,c2)


    updated.cohort<-rbind(under.grass,canopy)



    updated.cohort$Z<-f.Z(updated.cohort$DBH,h1,h2)

    updated.cohort<-updated.cohort[order(updated.cohort$Z,decreasing=T),]


    fully.up<-rbind(updated.cohort,new.cohort)
    rownames(fully.up)<-NULL


return(fully.up)


    }



  
  species.vitals2<-function(cohort,spframe,time,CA.all,type) {
  new.trees<-vector("list",length=nrow(spframe))
  
  dbh.vec<-rep(NA,length=nrow(spframe))
  
  
  lenzo<-length(pars)+1

  for(j in 1:dim(spframe)[1]) {
  spec.cohort<-subset(cohort,cohort$Species==spframe$Species[j])

  parz<-spframe[j,2:lenzo]

  new.trees[[j]]<-vitals2(spec.cohort,mu.g=parz$mu.g,mu.c=parz$mu.c,g.g=parz$g.g,g.c=parz$g.c,
  c1=parz$c1,c2=parz$c2,h1=parz$h1,h2=parz$h2,global=parz$global,f1=parz$f1,f2=parz$f2,
  time,CA.all,type)
  new.trees[[j]]$Species<-spframe$Species[j]
  
  new<-new.trees[[j]]
  
  cans<-subset(new,new$layer=="c")
  
  dbh.vec[j]<-sum(cans$DBH*cans$N*parz$stem.weight,na.rm=T)

  
  }

  new.mat<-do.call(rbind,new.trees)
  return(list(new.mat,dbh.vec))
  }


    getit<-function(cohort,what) {
    return(sum(cohort[,what],na.rm=T))
    }

      getC<-function(cohort) {
    tots<-nrow(cohort)
    cs<-length(which(cohort$layer=="c"))
    return(cs/tots)
    }


    pout<-function(pp2,...) {


    Ns<-unlist(lapply(pp2,getit,5))
    #print(Ns)
   plot(Ns~c(1:length(Ns)),type="l",xlab="Time",ylab="N",lwd=2,...)
   box()
    #return(Ns)



    par(new=T)

   cprop<-unlist(lapply(pp2,getC))
   plot(cprop~c(1:length(cprop)),type="l",axes=F,col="green",ylab="",xlab="",lwd=2,...)
    axis(4)

   outs<-list(cprop,Ns)
    names(outs)<-c("cprop","Ns")
   return(outs)}


  co.plot<-function(pp2) {
  cohort.mat<-matrix(0,nrow=timez+11,ncol=timey)

  for(i in 1:timey) {
    cohort.now<-pp2[[i]]
    cohort.mat[cohort.now$cohort.id,i]<-cohort.now$N

  }
  
  use<-sample(size=5,x=c(100:timez)) #sampling after initial transient dynamics

   matplot(t(cohort.mat[use,]),type="l")

  }

    getZ<-function(cohort) {
    return(max(cohort$Z,na.rm=T))
    }

 

    reforest<-function(full,type) {

    pp2<-vector("list",length=timez)

    z_vec<-rep(NA,times=timez)

    pp2[[1]]<-cohort
    dbh.all<-rep(NA,times=timez)
    dbh.all[1]<-0
    CA.all<-rep(NA,times=timez)
    CA.all[1]<-0
    
    closed<-rep(0,times=1)

      for(i in 2:timez) {

      layers<-f.layers(pp2[[i-1]])
      the.forest<-layers

      z = layers[,c('N','Z','CA')]
      t_ca = sum(z$N/plot.area*z$CA)  #10000#
      
      CA.all[i]<-t_ca
      
      ord.coh = order(z$Z,decreasing=T)

      ca = cumsum(z$t_ca[ord.coh])

      z_vec[i]<-t_ca
      timeo<-i*tsteps
     
      dyno<-species.vitals2(the.forest,spframe,time=i,CA.all=CA.all[i-1], type=type)
      dynamics<-dyno[[1]]
      dbh.all[i]<-sum(dyno[[2]],na.rm=T)
      #print(dbh.all[i])
      dynamics2<-subset(dynamics,dynamics$N!=0)
      pp2[[i]]<-dynamics2
      }
      if(full=="full") {return(pp2)}
      else{
      return(z_vec)}
      } #pp2

###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
###########TIME TO REFORESTATION ONLY
  
  get.time<-function(listo) {
    return(c(1:length(listo[[2]])*tsteps))
  }
  
  
  reTIME<-function(spframe,type) {
  
    pp2<-vector("list",length=timez)
    
    z_vec<-rep(NA,times=timez)
    
    pp2[[1]]<-cohort
    dbh.all<-rep(NA,times=timez)
    dbh.all[1]<-0
    
    CA.all<-rep(NA,times=timez)
    CA.all[1]<-0
    
    closed<-rep(0,times=1)
    
  
        for(i in 2:timez) {
  #print(i)
        layers<-f.layers(pp2[[i-1]])
        the.forest<-layers
  
        #z = layers[,c('N','Z','CA')]
        t_ca = sum(layers$N/plot.area*layers$CA)  #10000#
  
        ord.coh = order(layers$Z,decreasing=T)
  
        ca = cumsum(t_ca[ord.coh])

          
        
        dyno<-species.vitals2(the.forest,spframe,time=i,CA.all=CA.all[i-1], type=type)
        dynamics<-dyno[[1]]
        dbh.all[i]<-sum(dyno[[2]],na.rm=T)
      
        dynamics2<-subset(dynamics,dynamics$N!=0)
        pp2[[i]]<-dynamics2
        
        z_vec[i]<-t_ca
        CA.all[i]<-t_ca
        #print(CA.all[i])
        timeo<-i*tsteps
        
        
        if(CA.all[i]>1) {
          CA.all[i:timez]=1          
          break}
        else{
        next}
        
        
        }
  
        TIME<-i*tsteps
        
        forest<-pp2
        return(list(TIME=TIME,CA.all=CA.all,forest=forest))
        } #pp2




  #cohort.id in this list is what you need
  get.cohort<-function(datf,what) {
  cots<-subset(datf,datf$cohort.id==what)
  coto<-sum(cots$CA*cots$N/plot.area) #*(cots$N/plot.area))
  return(coto)
  }

  get.CA<-function(dat) {return(sum(dat$N/plot.area*dat$CA))}


#GETTING RESULTS OF SIMULATIONS

louts<-function(co.num,...) {

  ppo<-unlist(lapply(pp2,get.cohort,co.num))

  star.t<-co.num-2
  en.d<-length(ppo)+star.t-1

  tempus<-c(star.t:en.d)*tsteps

  plot(ppo~tempus,type="l",col=co.num,xlab="Time",ylab="Canopy area of cohort",cex.lab=1.5,...) #xlim=c(0,timez),
  return(cumtrapz(c(1:length(ppo)),ppo))} #solve integral using trapezoidal method


  #simulated time to reforestation
refoTint<-function(tvec,datz) {
canopy<-sapply(tvec,CA_int,datz=spframe)
if(max(canopy)<1) {return(NA)}
else{
gops<-tail(which(canopy<1),n=1)
return(tvec[gops])
}
}


