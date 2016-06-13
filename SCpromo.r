# 13/06/2016 - Ramires, T.G
# ---------------------------------------------------------------------------------------
SCp= function (mu.link = "identity", sigma.link = "log", nu.link = "log", tau.link="log") 
{
	mstats <- checklink(   "mu.link", "SCpromotion", substitute(mu.link),    c("inverse", "log", "identity", "own"))
	dstats <- checklink("sigma.link", "SCpromotion", substitute(sigma.link), c("inverse", "log", "identity", "own"))
	vstats <- checklink(   "nu.link", "SCpromotion", substitute(nu.link),    c("inverse", "log", "identity", "own"))
	tstats <- checklink(  "tau.link", "SCpromotion", substitute(tau.link),   c("inverse", "log", "identity", "own")) 
	
	
	structure(list(family = c("SCp", "SCpromotion"),
																parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE, tau=TRUE),
																nopar = 4,
																type = "Continuous",
																mu.link = as.character(substitute(mu.link)),  
																sigma.link = as.character(substitute(sigma.link)), 
																nu.link = as.character(substitute(nu.link)), 
																tau.link = as.character(substitute(tau.link)), 
																
																mu.linkfun = mstats$linkfun, 
																sigma.linkfun = dstats$linkfun, 
																nu.linkfun = vstats$linkfun,
																tau.linkfun = tstats$linkfun,  
																
																mu.linkinv = mstats$linkinv, 
																sigma.linkinv = dstats$linkinv,
																nu.linkinv = vstats$linkinv,
																tau.linkinv = tstats$linkinv, 
																
																mu.dr = mstats$mu.eta, 
																sigma.dr = dstats$mu.eta, 
																nu.dr = vstats$mu.eta,
																tau.dr = tstats$mu.eta, 
																
																###################################################
																########## first and second derivate ##############
																###################################################
																#change log<-log , pi<-pi
																
																
																dldm = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	dldm =-sinh(w)/sigma/cosh(w)+2*nu^2*sinh(w)*cosh(w)/sigma/K+tau/pi*nu*cosh(w)/sigma/K
																	dldm 
																}, 
																d2ldm2 = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau,log = TRUE), 
																																													"mu", delta = 1e-04)
																	dldm = as.vector(attr(nd, "gradient"))
																	d2ldm2 = -dldm * dldm
																	d2ldm2 },
																
																dldd = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	dldd<- -1/sigma-sinh(w)*(log(y)-mu)/sigma^2/cosh(w)+2*nu^2*sinh(w)*cosh(w)*(log(y)-mu)/sigma^2/K+tau/pi*nu*cosh(w)*(log(y)-mu)/sigma^2/K
																	dldd
																}, 
																d2ldd2 = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau,log = TRUE), 
																																													"sigma", delta = 1e-04)
																	dldd = as.vector(attr(nd, "gradient"))
																	d2ldd2 = -dldd * dldd
																	d2ldd2
																}, 
																dldv = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	dldv<- 1/nu-2*nu*sinh(w)^2/K-tau/pi*sinh(w)/K
																	dldv
																}, 
																d2ldv2 = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau,log = TRUE), 
																																													"nu", delta = 1e-04)
																	dldv = as.vector(attr(nd, "gradient"))
																	d2ldv2 = -dldv * dldv
																	d2ldv2
																},
																
																dldt = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	dldt<-1/tau-.5-1/pi*arctan(nu*sinh(w))
																	dldt
																},
																
																d2ldt2 = function(y,mu,sigma,nu,tau){ 
																	d2ldt2=-1/tau^2
																	d2ldt2
																},
																
																
																############################
																# crossed derivatess #######
																############################
																
																d2ldmdd = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu, tau,
																																																		log = TRUE), "mu", delta = 1e-04)
																	dldm = as.vector(attr(nd, "gradient"))
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu, tau,
																																																		log = TRUE), "sigma", delta = 1e-04)
																	dldd = as.vector(attr(nd, "gradient"))           
																	d2ldmdd = -dldm * dldd
																	d2ldmdd },
																
																d2ldmdv = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau, 
																																																		log = TRUE), "mu", delta = 1e-04)
																	dldm = as.vector(attr(nd, "gradient"))
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau, 
																																																		log = TRUE), "nu", delta = 1e-04)
																	
																	dldv = as.vector(attr(nd, "gradient"))
																	d2ldmdv = -dldm * dldv
																	d2ldmdv },
																
																d2ldmdt = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	d2ldmdt <- 1/pi*nu*cosh(w)/sigma/K
																	d2ldmdt
																},
																
																d2ldddv = function(y,mu,sigma,nu,tau){
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau, 
																																																		log = TRUE), "sigma", delta = 1e-04)
																	dldm = as.vector(attr(nd, "gradient"))
																	nd = gamlss:::numeric.deriv(dSCp(y, mu, sigma, nu,tau, 
																																																		log = TRUE), "nu", delta = 1e-04)
																	
																	dldv = as.vector(attr(nd, "gradient"))
																	d2ldmdv = -dldm * dldv
																	d2ldmdv },
																d2ldddt = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	d2ldddt <- 1/pi*nu*cosh(w)*(log(y)-mu)/sigma^2/K
																	d2ldddt 
																},
																d2ldvdt = function(y,mu,sigma,nu,tau){
																	w=(log(y)-mu)/sigma
																	K=nu^2*(sinh(w))^2+1
																	
																	d2ldvdt <- -1/pi*sinh(w)/K
																	d2ldvdt 
																},
																
																G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*dSCp(y, mu = mu, sigma = sigma, nu=nu, tau=tau,log = TRUE), 
																rqres = expression(
																	rqres(pfun="pSCp", type="Continuous",ymin=0, y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)
																), 
																
																mu.initial = expression(     mu <- rep(4,length(y))),   
																sigma.initial = expression(sigma<- rep(0.07, length(y))),
																nu.initial = expression(nu <- rep(0.01, length(y))), 
																tau.initial = expression(tau <-rep(0.5, length(y))), 
																mu.valid = function(mu) TRUE, 
																sigma.valid = function(sigma)  all(sigma > 0),
																nu.valid = function(nu) all(nu > 0.01) , 
																tau.valid = function(tau) all(tau >= 0.01) && all(tau < .99),
																y.valid = function(y)   all(y > 0)
	),
	class = c("gamlss.family","family"))
}

#------------------------------------------------------------------------------------------- # 
dSCp<-function(x, mu = 0.5, sigma = 0.5, nu = 0.1,tau=0, log = FALSE){ 
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
	if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
	w<-(log(x)-mu)/sigma
	fy1 <- tau*((nu)/(x*sigma*pi))*(cosh(w)/(nu^2*(sinh(w))^2+1))*exp(-tau*pSCp(x,mu,sigma,nu,tau))
	if(log==FALSE) fy<-fy1 else fy<-log(fy1)
	fy
}
#------------------------------------------------------------------------------------------ #ok
pSCp <- function(q, mu = 0.5, sigma = 0.5, nu = 0.1, tau=0, lower.tail = TRUE, log.p = FALSE){     
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
	if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
	
	w<-(log(q)-mu)/sigma
	cdf1 <- 1-exp(-tau*(0.5+(1/pi)*atan(nu*sinh(w))))
	if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
	if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
	cdf
}
#------------------------------------------------------------------------------------------ #ok
qESC <- function(p, mu = 0.5, sigma = 0.5, nu = 0.1, tau=0, lower.tail = TRUE, log.p = FALSE){      
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
#	if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
	if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
	if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
	if (log.p==TRUE) p <- exp(p) else p <- p
	if (lower.tail==TRUE) p <- p else p <- 1-p
	q <- exp(mu+sigma*asinh((1/nu)*tan(pi*(-log(1-u)/tau -0.5))))  
	q
}
#------------------------------------------------------------------------------------------ # ok
rESC <- function(n, mu =0.5, sigma =0.5, nu =0.1, tau=0){ 
	if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
	if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
	if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
	if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
	
	uni<- runif(n = n,0,tau)
	r <- qESC(uni,mu =mu, sigma =sigma, nu =  nu, tau=tau)
	r
}
#------------------------------------------------------------------------------------------
