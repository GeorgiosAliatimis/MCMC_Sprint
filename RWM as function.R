rm(list=ls())

library(mvtnorm)


# code for the distribution #############################################


mutheta = c(0,0)
sigtheta =  matrix(c(1,0.5,0.5,1),nrow=2,ncol=2)
sigtheta.inv = solve(sigtheta)


logpi = function(th){
  #log target density
  return(as.numeric( -t(th-mutheta)%*%sigtheta.inv%*%(th-mutheta)/2 ))
}


# RWM ###################################################################

rwm = function(x_0,N,h,logpi){
  d = length(x_0)
  
  path = matrix(0,ncol=N,nrow=d)
  path[,1] = x_0
  
  accept = 0
  
  
  for ( i in 2:(N) ){
    x_t = path[,i-1]
    
    x_prop = x_t  + h^0.5*rnorm(d)
    
    loga = logpi(x_prop) - logpi(x_t) 
    
    if ( log(runif(1)) <= loga ){
      accept = accept+1
      x_next = x_prop
    }
    else{
      x_next = x_t
    }
    
    path[,i] = x_next
  }
  
  return (list( path=path, a=accept/N ))
}


##########################################################################

theta_init=c(0,0)

out.rwm = rwm(x_0=theta_init,N=1e4,h = 4.5,logpi = logpi)
out.rwm$a #optimal is 0.234
chain = out.rwm$path


# code for plotting ######################################################

#scatterplot
plot(chain[1,],chain[2,])


#trace plots for each component
plot(1:length(chain[1,]),chain[1,],type='l', main = "traceplot component 1", col='red')
abline(h=0,lty=3)
plot(1:length(chain[2,]),chain[2,],type='l', main = "traceplot component 2", col='red')
abline(h=0,lty=3)


#density plots for each component

true_samples = rmvnorm(1000000,mutheta,sigtheta)

burnin=0

plot(density(chain[1,burnin:length(chain[1,])]), col = "red")
lines(density(true_samples[,1]),col='blue')

plot(density(chain[2,burnin:length(chain[1,])]), col = "red")
lines(density(true_samples[,2]),col='blue')
