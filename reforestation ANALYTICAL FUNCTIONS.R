
  #
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION
  ###############ANALYTICAL SOLUTION

  #canopy area of a cohort with lag time
  CA_t<-function(datz,TIME,lagged=0) {
    
    D.g=exp((log(grass.ht/datz$h1[1]))/datz$h2[1])
    
    crit.time<-grass.ht/(datz$g.g)+lagged
    
    if(TIME<crit.time) {
      CA_t=0
      return(CA_t)
    }
    else{
      
      CA_t<-datz$global*exp((-datz$mu.g)*D.g/datz$g.g)*exp(-datz$mu.c*(TIME))*datz$c1*(D.g+datz$g.c*(TIME))
      return(CA_t/plot.area)
    }
  }

CA_I<-function(datz,TIME,lagged=0) {
  
  D=exp((log(grass.ht/datz$h1[1]))/datz$h2[1])
  
  
  c.<-grass.ht/datz$g.g+lagged
  
  S<-(datz$global*exp((-datz$mu.g)*D/datz$g.g)*datz$c1)/plot.area
  
  m<-datz$mu.c
  
  g<-datz$g.c
  
  
  CA_t<-(S*(g+D*m-exp(-((TIME-c.)*m))*(g+D*m+g*(TIME-c.)*m)))/m^2
  
  
  #(S*exp(-m*(TIME-c.))*(c.*g*m-D*m-g*m*TIME-g))/(m^2)+
  #(S*(D*m+g))/(m^2)
  
  
  return(CA_t)
}

anal.time<-function(data) {
  
  S<-data$global
  g.g<-data$g.g
  c1<-data$c1
  mu.g<-data$mu.g
  mu.c<-data$mu.c
  h1<-data$h1
  h2<-data$h2
  g.c<-data$g.c
  
  
  
  
  equilibrium<-f.or.no(grass.ht,plot.area,S,g.g,g.c,mu.g,
                       mu.c,c1,h1,h2)
  print(paste("Equilibrium=",equilibrium))
  
  
  
  if(equilibrium<1) {return(Inf)} else{
    
    D=exp((log(grass.ht/h1))/h2)
    
    c.<-grass.ht/g.g #critical time
    
    STUFF<-(S*exp((-mu.g*D)/g.g)*c1)/plot.area
    
    m<-mu.c
    
    g<-g.c
    
  
    t.=c.-(g*W_1((exp(-(g+D*m)/g)*(m*(m-D*STUFF)-g*STUFF))/(g*STUFF))+D*m+g)/(g*m)
    
    
    return(c(t.))
    
  }
}


f.or.no<-function(grass.ht,plot.area,S,g.g,g.c,mu.g,mu.c,c1,h1,h2) {
  
  D.g=exp((log(grass.ht/h1))/h2)
  
  crit.time<-grass.ht/g.g
  
  CA=((S*exp((-mu.g*D.g)/g.g)*c1)*(D.g/mu.c+g.c/mu.c^2))/plot.area
  
  return(CA)
}
  
  
  
 equi.time<-function(data) {
    
    S<-data$global
    g.g<-data$g.g
    c1<-data$c1
    mu.g<-data$mu.g
    mu.c<-data$mu.c
    h1<-data$h1
    h2<-data$h2
    g.c<-data$g.c
    
    equilibrium<-f.or.no(grass.ht,plot.area,S,g.g,g.c,mu.g,
                         mu.c,c1,h1,h2)
    return(equilibrium)
    
    
  }
  
#
