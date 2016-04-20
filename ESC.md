# ELSC
# 01/06/2015 - Ramires, T.G
# ---------------------------------------------------------------------------------------
ESC= function (mu.link = "identity", sigma.link = "log", nu.link = "log", tau.link="log") 
{
    mstats <- checklink(   "mu.link", "ESC", substitute(mu.link),    c("inverse", "log", "identity", "own"))
    dstats <- checklink("sigma.link", "ESC", substitute(sigma.link), c("inverse", "log", "identity", "own"))
    vstats <- checklink(   "nu.link", "ESC", substitute(nu.link),    c("inverse", "log", "identity", "own"))
    tstats <- checklink(  "tau.link", "ESC", substitute(tau.link),   c("inverse", "log", "identity", "own")) 

  
  structure(list(family = c("ESC", "exp-SC"),
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

                 dldm = function(y,mu,sigma,nu,tau){ # ok
                   	dldm = ((-2*nu^2*sinh((y-mu)/sigma)^3-2*sinh((y-mu)/sigma)+4*nu^2*cosh((y-mu)/sigma)^2*sinh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^3*pi+(-pi+2*nu^2*cosh((y-mu)/sigma)^2*pi)*sinh((y-mu)/sigma)-2*nu*cosh((y-mu)/sigma)^2*(tau-1))/sigma/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
                   	dldm 
                 }, 
                 d2ldm2 = function(y,mu,sigma,nu,tau){#ok
                   	dldm = ((-2*nu^2*sinh((y-mu)/sigma)^3-2*sinh((y-mu)/sigma)+4*nu^2*cosh((y-mu)/sigma)^2*sinh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^3*pi+(-pi+2*nu^2*cosh((y-mu)/sigma)^2*pi)*sinh((y-mu)/sigma)-2*nu*cosh((y-mu)/sigma)^2*(tau-1))/sigma/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
            		d2ldm2 = -dldm * dldm
            		d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   	d2ldm2 
                 },
                 dldd = function(y,mu,sigma,nu,tau){ # ok
			dldd<- ((-2*nu^2*(y-mu)*sinh((y-mu)/sigma)^3-2*sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2+4*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*(y-mu)*sinh((y-mu)/sigma)-2*sigma*cosh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*pi*(y-mu)*sinh((y-mu)/sigma)^3-sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2*pi+2*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*pi*(y-mu)*sinh((y-mu)/sigma)-2*(nu*(tau-1)*(y-mu)*cosh((y-mu)/sigma)+1/2*sigma*pi)*cosh((y-mu)/sigma))/sigma^2/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldd
                 }, 
                 d2ldd2 = function(y,mu,sigma,nu,tau){#ok
			dldd<- ((-2*nu^2*(y-mu)*sinh((y-mu)/sigma)^3-2*sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2+4*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*(y-mu)*sinh((y-mu)/sigma)-2*sigma*cosh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*pi*(y-mu)*sinh((y-mu)/sigma)^3-sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2*pi+2*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*pi*(y-mu)*sinh((y-mu)/sigma)-2*(nu*(tau-1)*(y-mu)*cosh((y-mu)/sigma)+1/2*sigma*pi)*cosh((y-mu)/sigma))/sigma^2/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
            		d2ldd2 = -dldd * dldd
            		 d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
			 d2ldd2

                 }, 
                 dldv = function(y,mu,sigma,nu,tau) {#ok
			dldv<- ((-2*nu^2*sinh((y-mu)/sigma)^2+2)*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^2*pi+2*nu*(tau-1)*sinh((y-mu)/sigma)+pi)/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/nu
			dldv
		 }, 
                 d2ldv2 = function(y,mu,sigma,nu,tau) {#ok
			dldv<- ((-2*nu^2*sinh((y-mu)/sigma)^2+2)*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^2*pi+2*nu*(tau-1)*sinh((y-mu)/sigma)+pi)/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/nu
			d2ldv2<- -dldv * dldv
			d2ldv2 <- ifelse(d2ldv2 < -1e-15, d2ldv2,-1e-15)                                    
            		d2ldv2
		 },

	      dldt = function(y,mu,sigma,nu,tau) {#ok
			dldt<- (1+tau*log(1/2+1/pi*atan(nu*sinh((y-mu)/sigma))))/tau
			dldt
                  },

	     d2ldt2 = function(y,mu,sigma,nu,tau) {#ok
			dldt<- (1+tau*log(1/2+1/pi*atan(nu*sinh((y-mu)/sigma))))/tau
			d2ldt2<- -dldt * dldt
   			d2ldt2 <- ifelse(d2ldt2 < -1e-15, d2ldt2,-1e-15)                                    
			d2ldt2
                   },




 
                d2ldmdd = function(y,mu,sigma,nu,tau){

                   	dldm = ((-2*nu^2*sinh((y-mu)/sigma)^3-2*sinh((y-mu)/sigma)+4*nu^2*cosh((y-mu)/sigma)^2*sinh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^3*pi+(-pi+2*nu^2*cosh((y-mu)/sigma)^2*pi)*sinh((y-mu)/sigma)-2*nu*cosh((y-mu)/sigma)^2*(tau-1))/sigma/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldd<- ((-2*nu^2*(y-mu)*sinh((y-mu)/sigma)^3-2*sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2+4*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*(y-mu)*sinh((y-mu)/sigma)-2*sigma*cosh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*pi*(y-mu)*sinh((y-mu)/sigma)^3-sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2*pi+2*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*pi*(y-mu)*sinh((y-mu)/sigma)-2*(nu*(tau-1)*(y-mu)*cosh((y-mu)/sigma)+1/2*sigma*pi)*cosh((y-mu)/sigma))/sigma^2/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
           		 d2ldmdd = -(dldm * dldd)
			 d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
            		 d2ldmdd                 
		   }, 
                 d2ldmdv = function(y,mu,sigma,nu,tau){
                   	dldm = ((-2*nu^2*sinh((y-mu)/sigma)^3-2*sinh((y-mu)/sigma)+4*nu^2*cosh((y-mu)/sigma)^2*sinh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^3*pi+(-pi+2*nu^2*cosh((y-mu)/sigma)^2*pi)*sinh((y-mu)/sigma)-2*nu*cosh((y-mu)/sigma)^2*(tau-1))/sigma/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldv<- ((-2*nu^2*sinh((y-mu)/sigma)^2+2)*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^2*pi+2*nu*(tau-1)*sinh((y-mu)/sigma)+pi)/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/nu

            		 d2ldmdv = -(dldm * dldv)
            		 d2ldmdv				
                   },
 	        d2ldmdt = function(y,mu,sigma,nu,tau){ #---------------- here ok, check errors!
                   	dldm = ((-2*nu^2*sinh((y-mu)/sigma)^3-2*sinh((y-mu)/sigma)+4*nu^2*cosh((y-mu)/sigma)^2*sinh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^3*pi+(-pi+2*nu^2*cosh((y-mu)/sigma)^2*pi)*sinh((y-mu)/sigma)-2*nu*cosh((y-mu)/sigma)^2*(tau-1))/sigma/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldt<- (1+tau*log(1/2+1/pi*atan(nu*sinh((y-mu)/sigma))))/tau
			d2ldmdt <- -(dldm*dldt)
			d2ldmdt
               },

                 d2ldddv = function(y,mu,sigma,nu,tau){
			dldd<- ((-2*nu^2*(y-mu)*sinh((y-mu)/sigma)^3-2*sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2+4*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*(y-mu)*sinh((y-mu)/sigma)-2*sigma*cosh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*pi*(y-mu)*sinh((y-mu)/sigma)^3-sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2*pi+2*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*pi*(y-mu)*sinh((y-mu)/sigma)-2*(nu*(tau-1)*(y-mu)*cosh((y-mu)/sigma)+1/2*sigma*pi)*cosh((y-mu)/sigma))/sigma^2/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldv<- ((-2*nu^2*sinh((y-mu)/sigma)^2+2)*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^2*pi+2*nu*(tau-1)*sinh((y-mu)/sigma)+pi)/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/nu
			d2ldddv = -(dldd * dldv)
            		d2ldddv	
                 },        
      		d2ldddt = function(y,mu,sigma,nu,tau){
			dldd<- ((-2*nu^2*(y-mu)*sinh((y-mu)/sigma)^3-2*sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2+4*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*(y-mu)*sinh((y-mu)/sigma)-2*sigma*cosh((y-mu)/sigma))*atan(nu*sinh((y-mu)/sigma))-nu^2*pi*(y-mu)*sinh((y-mu)/sigma)^3-sigma*cosh((y-mu)/sigma)*nu^2*sinh((y-mu)/sigma)^2*pi+2*(-1/2+nu^2*cosh((y-mu)/sigma)^2)*pi*(y-mu)*sinh((y-mu)/sigma)-2*(nu*(tau-1)*(y-mu)*cosh((y-mu)/sigma)+1/2*sigma*pi)*cosh((y-mu)/sigma))/sigma^2/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/cosh((y-mu)/sigma)
			dldt<- (1+tau*log(1/2+1/pi*atan(nu*sinh((y-mu)/sigma))))/tau
			d2ldddt <- -(dldd*dldt) 
			d2ldddt 
               },
	       d2ldvdt = function(y,mu,sigma,nu,tau){ 
			dldv<- ((-2*nu^2*sinh((y-mu)/sigma)^2+2)*atan(nu*sinh((y-mu)/sigma))-nu^2*sinh((y-mu)/sigma)^2*pi+2*nu*(tau-1)*sinh((y-mu)/sigma)+pi)/(nu^2*sinh((y-mu)/sigma)^2+1)/(pi+2*atan(nu*sinh((y-mu)/sigma)))/nu
			dldt<- (1+tau*log(1/2+1/pi*atan(nu*sinh((y-mu)/sigma))))/tau
			d2ldvdt <- -(dldv*dldt) 
			d2ldvdt 
               },

                 G.dev.incr  = function(y,mu,sigma,nu,tau,...) -2*log((dESC(y, mu = mu, sigma = sigma, nu=nu, tau=tau,log = FALSE))), 
                 rqres = expression(
                   rqres(pfun="pESC", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)
               ), 

			mu.initial = expression(     mu <- rep(4,length(y))),   
			sigma.initial = expression(sigma<- rep(1.5, length(y))),
			nu.initial = expression(nu <- rep(2.5, length(y))), 
			tau.initial = expression(tau <-rep(1, length(y))), 
			mu.valid = function(mu) TRUE, 
			sigma.valid = function(sigma)  all(sigma > 0),
			nu.valid = function(nu) TRUE , 
			tau.valid = function(tau) all(tau > 0), 
			y.valid = function(y)  TRUE 
 		),
  class = c("gamlss.family","family"))
}

#------------------------------------------------------------------------------------------- # ok
dESC<-function(x, mu = 0, sigma = 0.5, nu = 0.1, tau=1, log = FALSE){ 
          if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
          if (any(tau < 0))  stop(paste("tau must be positive", "\n", ""))  
          if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))
w<-(x-mu)/sigma
fy1 <- ((tau*nu)/(sigma*pi))*(cosh(w)/(nu^2*(sinh(w))^2+1))*(0.5+(1/pi)*atan(nu*sinh(w)))^(tau-1)
if(log==FALSE) fy<-fy1 else fy<-log(fy1)
fy
}
#------------------------------------------------------------------------------------------ #ok
pESC <- function(q, mu = 0, sigma = 0.5, nu = 0.1, tau=1, lower.tail = TRUE, log.p = FALSE){     
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          

w<-(q-mu)/sigma
cdf1 <- (0.5+(1/pi)*atan(nu*sinh(w)))^tau
if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
cdf
}
#------------------------------------------------------------------------------------------ #ok
qESC <- function(p, mu = 0, sigma = 0.5, nu = 0.1, tau=1, lower.tail = TRUE, log.p = FALSE){      
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
         if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  	 if (log.p==TRUE) p <- exp(p) else p <- p
  	 if (lower.tail==TRUE) p <- p else p <- 1-p
  q <- mu+sigma*asinh((1/nu)*tan(pi*(p^(1/tau)-0.5)))
  q
}
#------------------------------------------------------------------------------------------ # ok
rESC <- function(n, mu =0, sigma =0.5, nu =0.1, tau=1){ 
         if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
         if (any(tau < 0))  stop(paste("tau must be positive", "\n", "")) 
         if (any(nu < 0))  stop(paste("nu must be positive", "\n", ""))          
         if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  	
	uni<- runif(n = n,0,1)
  	r <- qESC(uni,mu =mu, sigma =sigma, nu =  nu, tau=tau)
        r
}
#------------------------------------------------------------------------------------------

