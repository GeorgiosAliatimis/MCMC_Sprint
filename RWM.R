#####
# input a value of the state, and outputs log pi(x)
####
logpi=function(x){ 
	return(sum(dnorm(x,log=T))) ###IID Gaussian
	#return(-0.5*(x[1]^2+x[2]^2-x[1]*x[2])) ### 2-d correlated Gaussian 
}

###function to generate from the proposal distribution -- RWM with Gaussian noise variance h
rprop=function(x,h=1){
	return(rnorm(length(x))*sqrt(h)+x)
}

###simple code for Metropolis Hastings with symmetric proposal

N=100000 ##number of iterations
x.init=c(0,0) #initial values
h=1 ##proposal scale

x=x.init ##current state
x.st=matrix(NA,nrow=N,ncol=length(x.init)) #store output from chain

for(i in 1:N){ ##LOOP over iterations
	
	y=rprop(x,h=h) ##proposal
	log.acc=logpi(y)-logpi(x) ##calculate log of term in acceptance probability
	##no log q(x|y) terms due to cancellation for symmetric random walk
	
	if(runif(1)<exp(log.acc)){ ##accept-reject 
		x=y ##if accept then update state
	}
	
	x.st[i,]=x ##store output
}