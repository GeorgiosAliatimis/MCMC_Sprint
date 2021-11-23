leapfrog <- function(x_init,p_init,K,h,grad_log_pi){
	p <- p_init
	x <- x_init 
	for(i in 1:K){
		p <- p + h/2 * grad_log_pi(x)
		x <- x + h * p
		p <- p + h/2 * grad_log_pi(x)
	}
	c(x,p)
}

HMC <- function(N,x_init,log_pi,K,h,grad_log_pi){
    d <- length(x_init)
	x.st <- matrix(NA,nrow=N,ncol=d)
	x <- x_init
	for(i in 1:N){
		p <- rnorm(d)
		updates <- leapfrog(x,p,K,h,grad_log_pi)
		x.prop <- updates[1:d] 
		p.prop <- updates[(d+1):(2*d)]
		log.acc <- log_pi(x.prop) - log_pi(x) + 
				norm(p,type="2")^2 - norm(p.prop,type="2")^2
		if(runif(1)<exp(log.acc)) x = x.prop 
		x.st[i,] <- x 
	}
	x.st 
}

#Note logpi is defined up to an additive constant
set.seed(1)
log_pi <- function(x) -sum(x^2)/2
grad_log_pi <- function(x) -sum(x) 
N<-1e+5
x_init <- c(0,0)
h <- 1/20
K <- 20
x <- HMC(N,x_init,log_pi,K,h,grad_log_pi)
plot(x[,1],x[,2])