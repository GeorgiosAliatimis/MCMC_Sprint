#####
# input a value of the state, and outputs log pi(x)
####
logpi=function(x){ 
  return(-0.5*((0.1*x[1])^2+x[2]^2)) ### 2-d Gaussian -- independent but different scales
}

###function to generate from the proposal distribution -- RWM with Gaussian noise variance h
rprop=function(x,h=1){
  return(rnorm(length(x))*sqrt(h)+x)
}

###simple code for Metropolis Hastings with symmetric proposal

N=20000 ##number of iterations
x.init=c(0,0) #initial values
x.st=matrix(NA,nrow=N,ncol=length(x.init)) #store output from chain
h=1 ##proposal scale

x=x.init ##current state

for(i in 1:N){ ##LOOP over iterations
  
  y=rprop(x,h=h) ##proposal
  log.acc=logpi(y)-logpi(x) ##calculate log of term in acceptance probability
  ##no log q(x|y) terms due to cancellation for symmetric random walk
  
  if(runif(1)<exp(log.acc)){ ##accept-reject
    x=y ##if accept then update state
  }
  
  x.st[i,]=x ##store output
}


#####PLOTS OF THE RESULTS
##TRACE PLOTS
par(mfrow=c(2,1))
plot.ts(x.st[,1])
plot.ts(x.st[,2])

###ACF PLOTS AFTER BURN-IN FOR g(x)=x^{(1)} and g(x)=x^{(2)} and g(x)=x^{(1)}^2x^{(2)}
B=1e3
par(mfrow=c(3,1))
acf(x.st[(B+1):N,1])
acf(x.st[(B+1):N,2])
acf(x.st[(B+1):N,1]^2*x.st[(B+1):N,2])

###ESTIMATE AUTO_CORRELATION TIME!
### we can estimate the auto-correlation time for x^{(2)}

acf.x2=acf(x.st[(B+1):N,2])
##acf seems to be negligible aftera lag fo 20
1+2*sum(acf.x2$acf[2:21]) ##acf.x2$acf starts with lag 0 acf hence we sum vector from entry 2 to 21


