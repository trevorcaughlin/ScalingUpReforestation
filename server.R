shinyServer(function(input, output) {
  
  
  output$text1 <- renderText({ 
    "Dynamics of canopy cover in a reforesting patch. At time 0, the patch is dominated by a grass layer. Over time, 
    tree seeds arrive at a constant rate and expand their crowns above the grass layer. The dashed horizontal line represents a canopy cover of 100%, and the dashed vertical line represents time to canopy closure.
    Infinite time to canopy closure indicates arrested succession."
  })
  
  
  CA_I5<-function(grass.ht=1,plot.area=1000,S,g.g=0.8,g.c=1.9,mu.g=0.837,
                  mu.c=0.50,c1=1.38,h1=1,h2=1,TIME,lagged=0) {
    
    #S=S*1000
    
    D=exp((log(grass.ht/h1[1]))/h2[1])
    
    
    c.<-grass.ht/g.g+lagged
    
    Stuff<-(S*exp((-mu.g)*D/g.g)*c1)/plot.area
    
    m<-mu.c
    
    g<-g.c
    
    
    
    CA_t<-ifelse(TIME>c.,(Stuff*(g+D*m-exp(-((TIME-c.)*m))*(g+D*m+g*(TIME-c.)*m)))/m^2,0)
    
    CA_t<-ifelse(CA_t>1,1,CA_t)
    
    #(S*exp(-m*(TIME-c.))*(c.*g*m-D*m-g*m*TIME-g))/(m^2)+
    #(S*(D*m+g))/(m^2)
    
    
    
    return(CA_t)
  }
  
  
  f.or.no=function(grass.ht=1,plot.area=1000,S,g.g=0.8,g.c=1.9,mu.g=0.837,
                   mu.c=0.50,c1=1.38,h1=1,h2=1) {
    #S=S*1000
    D.g=exp((log(grass.ht/h1))/h2)
    
    crit.time<-grass.ht/g.g
    
    CA=((S*exp((-mu.g*D.g)/g.g)*c1)*(D.g/mu.c+g.c/mu.c^2))/plot.area
    
    return(CA)
  }
  
  SEED.perturb<-function(grass.ht=1,plot.area=1000,S,g.g=0.8,g.c=1.9,mu.g=0.837,
                         mu.c=0.50,c1=1.38,h1=1,h2=1) {
    #S=S*1000
    equilibrium<-f.or.no(grass.ht,plot.area,S,g.g,g.c,mu.g,
                         mu.c,c1,h1,h2)
    D=exp((log(grass.ht/h1))/h2)
    
    c.<-grass.ht/g.g #critical time
    
    STUFF<-(S*exp((-mu.g*D)/g.g)*c1)/plot.area
    
    m<-mu.c
    
    g<-g.c
    
    t.=c.-(g*W_1((exp(-(g+D*m)/g)*(m*(m-D*STUFF)-g*STUFF))/(g*STUFF))+D*m+g)/(g*m)
    
    refo.val=ifelse(equilibrium<1,Inf,t.)
    
    return(c(refo.val))
    
  }
  
  
  
  SEED.perturbTWO<-function(grass.ht=1,plot.area=1000,S,g.g=0.8,g.c=1.9,mu.g=0.837,
                            mu.c=0.50,c1=1.38,h1=1,h2=1) {
    #S=S*1000
    equilibrium<-f.or.no(grass.ht,plot.area,S,g.g,g.c,mu.g,
                         mu.c,c1,h1,h2)
    D=exp((log(grass.ht/h1))/h2)
    
    c.<-grass.ht/g.g #critical time
    
    STUFF<-(S*exp((-mu.g*D)/g.g)*c1)/plot.area
    
    m<-mu.c
    
    g<-g.c
    
    t.=c.-(g*W_1((exp(-(g+D*m)/g)*(m*(m-D*STUFF)-g*STUFF))/(g*STUFF))+D*m+g)/(g*m)
    
    refo.val=ifelse(equilibrium<1,Inf,t.)
    
    return(c(refo.val))
    
  }
  
  
  SEED.perturb3<-function(grass.ht=1,plot.area=1000,S,g.g=0.8,g.c=1.9,mu.g=0.837,
                            mu.c=0.50,c1=1.38,h1=1,h2=1) {
    #S=S*1000
    equilibrium<-f.or.no(grass.ht,plot.area,S,g.g,g.c,mu.g,
                         mu.c,c1,h1,h2)
    D=exp((log(grass.ht/h1))/h2)
    
    c.<-grass.ht/g.g #critical time
    
    STUFF<-(S*exp((-mu.g*D)/g.g)*c1)/plot.area
    
    m<-mu.c
    
    g<-g.c
    
    t.=c.-(g*W_1((exp(-(g+D*m)/g)*(m*(m-D*STUFF)-g*STUFF))/(g*STUFF))+D*m+g)/(g*m)
    
    refo.val=ifelse(equilibrium<1,Inf,t.)
    
    if(is.infinite(refo.val)|refo.val<0) {return("Time to canopy closure= Infinite")}
    else{return(paste("Time to canopy closure=",round(refo.val,2),"yrs",sep=" "))}
    
    
  }
  
  
  
  W_1=function (z) 
  {
    return(W(z, branch = -1))
  }
  
  W=function (z, branch = 0) 
  {
    stopifnot(branch == 0 || branch == -1)
    stopifnot(!any(is.na(z)), !any(is.nan(z)))
    W.z <- rep(NA, length(z))
    if (branch == 0) {
      W.z[z != Inf] <- gsl::lambert_W0(z[z != Inf])
      W.z[z == Inf] <- Inf
    }
    else if (branch == -1) {
      if (any(z == Inf)) {
        stop("Inf is not a valid argument of the non-principal branch W", 
             " (branch = -1).")
      }
      W.z <- gsl::lambert_Wm1(z)
    }
    W.z[is.nan(W.z)] <- NA
    dim(W.z) <- dim(z)
    return(W.z)
  }
  
  output$TIMEplot <- renderPlot({
    par(mar=c(5,5,5,5))
   vecs<-seq(from=0.1,to=input$MAXY,length=100)
    giveittome<-SEED.perturb(S=vecs,g.c=input$GROWTHY,mu.c=-log(input$SURVY))
    if(length(which(is.infinite(giveittome)))==length(giveittome)) 
    {plot(1,axes=F,xlab="",ylab="",pch=""); text(1,1,"Arrested succession: canopy closure will never occur",cex=1.5)
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
     text(1,1,"Arrested succession: canopy closure will never occur",cex=1.9)}
    else{plot(giveittome~vecs,type="l",lwd=7,col="darkgreen",xlab=expression(paste("Seeds ","m"^-2," y"^-1)),
              ylab="Years to canopy closure",cex.lab=1.3,cex.axis=1.3)
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
      lines(giveittome~vecs,lwd=5,col="darkgreen")
      }   
  })
  
  output$CANOPYplot <- renderPlot({
    par(mar=c(5,5,5,5))
    hg=1/input$GROWTHYg
  curve(100*CA_I5(TIME=x,S=input$MAXY,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),g.c=input$GROWTHY,mu.c=-log(input$SURVY)),from=0,to=25,xlab="Time (years)",
        ylab="Tree canopy cover (%)",
  cex.lab=1.5,cex.axis=1.3,lwd=8,col="darkgreen",main=SEED.perturb3(S=input$MAXY,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),g.c=input$GROWTHY,mu.c=-log(input$SURVY)),
  cex.main=1.5,ylim=c(0,100*CA_I5(TIME=25,S=input$MAXY,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),g.c=input$GROWTHY,mu.c=-log(input$SURVY)))) #
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "gray90")
    grid(nx=10,ny=10)
    curve(100*CA_I5(TIME=x,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),S=input$MAXY,g.c=input$GROWTHY,
  mu.c=-log(input$SURVY)),from=0.8,to=25,xlab="Time (years)",ylab="Canopy cover",cex.lab=1.3,cex.axis=1.3,lwd=10,col="green",add=T)
  abline(h=100,lwd=3,col="gray30",lty=2)
  abline(v=SEED.perturb(S=input$MAXY,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),g.c=input$GROWTHY,mu.c=-log(input$SURVY)),lwd=3,col="gray30",lty=2)
  points(SEED.perturb(S=input$MAXY,g.g=input$GROWTHYg,mu.g=-log(input$SURVYg),g.c=input$GROWTHY,mu.c=-log(input$SURVY)),100,pch=19)
 ## text(SEED.perturb(S=input$MAXY,g.c=input$GROWTHY,mu.c=-log(input$SURVY))+5,2,"Canopy closure")
  mtext(side=1,text="Black dot represents time to canopy closure",outer=T)  
  
  
  })
  
})
