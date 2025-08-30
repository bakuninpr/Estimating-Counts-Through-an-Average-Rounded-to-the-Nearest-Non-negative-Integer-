##############################################################################
# The open source software R is used to measure Relative Risk of Rounding  
# when Y is a discrete random variable with a binomial, Poisson, or
# negative binomial distribution and U = n[Y/n], where [] means 
# rounding has ocurred.
#   See paper Estimating Counts Through an Average Rounded to the Nearest Non-negative Integer and its
# Theoretical \& Practical Effectsfor details.

# Some functions created include commentary on required arguments 
#
# R o R Studio for Windows can be downloaded from
# https://cran.r-project.org/bin/windows/base/ or https://www.rstudio.com/products/rstudio/download3/. 
# This code runs with R 4.3.0 since August, 2025. In case of error,
# verify your version of R and the data; changes may have occurred.
# Created by R Rivera: January 7, 2025. 

# NOTE: Each line that starts with "#" R detects as a comment.

# The following function will be used to determine the MLE when 
# Y ~ poisson distribution.
        gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
          y<-x   
          if(any(x < 0, na.rm = TRUE)){
            #return(NaN)
            w<-which(x<=0)
            y<-x[-w]
          }
          if(zero.propagate){
            if(any(x <= 0, na.rm = TRUE)){
              return(0)
            }
            exp(mean(log(x), na.rm = na.rm))
          } else {
            exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[x > 0]))
          }
          nm<-ifelse(any(x <= 0, na.rm = TRUE),length(which(x<=0)),0)
          gm<-exp(sum(log(y), na.rm=na.rm) / length(y))
          mlist<-list("gm" = gm, "nm" = nm)
          return(mlist)
        }
        
        
        ###### MLE function for theta
        ###### arguments
        ###### n = sample size used to obtain rounded average
        ###### u = count obtained after rounding average and multiplying by n
        mle_u_p<-function(u,n){
          if(u%%n!=0) return("u has to be a multiple of n")
          if((n %% 2) == 0) {
            #print(paste(n,"is Even"))
            x<-seq(u-n/2,u+n/2-1,by=1)
          } else {
            #print(paste(n,"is Odd"))
            x<-seq(u-n/2+1/2,u+n/2-1/2,by=1)
          }
          if(n==2) {
            ahat<-0
          } else{
            ahat<-as.numeric(gm_mean(x)$gm)
          }
          my_list <- list("t_e" = ahat)
          return(my_list)
        }
        

        
        
        loglike=function(p, u, n, pdf, RoundUp=FALSE) {
          if(u%%n!=0) return("u has to be a multiple of n!")
          y=0:(u+2*n) # all integers that might round to u
          if(RoundUp) z=n*floor((y/n+0.5+1e-10)) # round them, 1/2 up
          else z=n*floor(y/n+.5) # round them, 1/2 down, corrected Wolfgang's original code
          y=y[z==u] #  keep those y rounded to u
          ll=0*p
          #pd=pdf(y, p[i]) #I think this line and the one that follows were added later when having issues with the sum(log(...)) which is supposed to be log(sum(....)) anyway
          #pd=pd[pd>0]
          for(i in seq_along(p)) 
            #ll[i]=sum(log(pd)) 
            ll[i]=log(sum(pdf(y, p[i])))
          # log of sum over y's of pdfs
          ll
        }
        mle = function(interval, u, n, pdf, RoundUp=FALSE) {
          if(u%%n!=0) return("u has to be a multiple of n!")
          optimize(loglike, interval, u=u, n=n, pdf=pdf, 
                   RoundUp=RoundUp, maximum=TRUE)$maximum
        } 
        ci = function(u, n, Model, N, alpha=0.05, RoundUp=FALSE) {
          if(u%%n!=0) return("u has to be a multiple of n!")
          interval=c(0, 1) # optim interval for Binomials
          if(Model=="Poisson") {
            interval=c(max(0, u-5*n), u+5*n)
            pdf=function(x, p) dpois(x, p)
            ci1=c(qchisq(alpha/2, 2*(u + 1))/2, u, qchisq(1-alpha/2, 2*(u + 1))/2)
          }  
          if(Model=="Binomial") {
            pdf=function(x, p) dbinom(x, N, p)  
            ci1=c(binom.test(u, N,conf.level = 1-alpha)$conf.int)
            ci1=c(ci1[1], u/N, ci1[2])
          }  
          if(Model=="Negative Binomial") {
            pdf=function(x, p) dnbinom(x, N, p)  
            ci1=c(0, N/(N+u),0)
          }  
          loglike=function(p, u, n, pdf, RoundUp) {
            y=0:(u+2*n) # all integers that might round to u
            if(RoundUp) z=n*floor(y/n+.5+1e-10) # round them, 1/2 up
            else z=n*floor(y/n+.5) # round them, 1/2 down
            y=y[z==u] #  keep those y rounded to u
            ll=0*p
            for(i in seq_along(p)) 
              #ll[i]=sum(log(pdf(y, p[i])))
              ll[i]=log(sum(pdf(y, p[i])))
            # sum over y's of log(pdf)
            ll
          }
          h=1e-4 # step size for calculation of second derivative
          ml = optimize(loglike, interval, u=u, n=n, pdf=pdf, 
                        RoundUp=RoundUp, maximum=TRUE)$maximum
          tmp = loglike(ml+c(-1,0,1)*ml*h, u, n, pdf, RoundUp)
          hess = (tmp[1]-2*tmp[2]+tmp[3])/(ml*h)^2
          mlci=ml+c(-1,0,1)*qnorm(1-alpha/2)/sqrt(-hess)
          list(ci_y=ci1, ci_u=mlci)
        }

        
        

        # Model = probability distribution; Poisson, Binomial, and Neg Binomial
        # n sample size, may be a vector
        # p or interval value or range of values of parameter. Should be realistic (broad enough)
        # ppoints length of evaluations of ratio
        # lambda2 is for the Skellam dist of Y-X
        #     here Y, X are indep and Y, X have expected values p, lambda2
        # GraphU = to graph ratio of MSEU/MSEMLE; currently just for Poisson
        # logf= if true plot logarithm of mse ratios
        library(ggplot2)
        rror = function(Model, n, interval=c(0.05, 0.95), 
                             N, p, ppoints=250, lambda2,doGraph=FALSE,GraphU=FALSE,logf=FALSE) {  
          if(Model=="Poisson") {
            pdf=function(x, p) dpois(x, p)
            x = 0:qpois(0.999, interval[2])
            #range of likely observations, from 0
            # to 99.9th percentile with largest lambda      
            mle.y = x # mle of a single Poisson rv.
            I=c(0, max(x)+5) #search interval for optimize
          } 
          if(Model=="Binomial") {
            pdf=function(x, p) dbinom(x, N, p)
            x = 0:N
            mle.y = x/N
            I = c(0, 1)
          }
          if(Model=="Negative Binomial") {
            pdf=function(x, p) dnbinom(x, N, p)
            x = 0:qnbinom(0.999, N, interval[1])
            mle.y = N/(N+x)
            I = c(0, 1)
          }
          if(Model=="Skellam") {#added by Rob
            pdf=function(x, p) dskellam(x, p, lambda2)
            
            x = qskellam(0.001, interval[1], lambda2):qskellam(0.999, interval[2], lambda2) #qnbinom(0.999, N, interval[1])
            mle.y = x #not the mle for skellam distribution 
            I=c(min(x)-5, max(x)+5)
          }
          if(missing(p))  p = seq(interval[1], interval[2], length=ppoints)  
          # points to find ratio  
          mseY = matrix(0, length(p), length(n))
          mseUmle = matrix(0, length(p), length(n))
          mseU = matrix(0, length(p), length(n))
          for(j in seq_along(n)) { # loop over n's
            mle.u = 0*x
            v=0*x
            u = 0 # first x value
            if(Model=="Poisson") ml = mle_u_p(u, n[j])$t_e
            else{
              ml = mle(I, u, n[j],  pdf)
            }
            mle.u[1] = ml
            for(i in 2:length(x)) {
              v[i] = n[j]*floor(x[i]/n[j]+.5)#round(x[i]/n[j])
              if(v[i]==u) mle.u[i] = mle.u[i-1] # if u has not changed, mle is the same
              else {
                u = v[i]
                if(Model=="Poisson") mle.u[i] = mle_u_p(u, n[j])$t_e
                else{
                  
                  mle.u[i]=mle(I, u, n[j], pdf)  # need new mle
                }   
              }
            }
            for(i in 1:length(p)) {
              mseUmle[i, j] = sum( (mle.u-p[i])^2*pdf(x, p[i]) ) # mse of u based MLE
              mseU[i, j] = sum( (v-p[i])^2*pdf(x, p[i]) ) # mse of u 
              mseY[i, j] = sum( (mle.y-p[i])^2*pdf(x, p[i]) ) # mse of y
            }
            ratMLE<-c(mseUmle/mseY)
            lratMLE=log(ratMLE)
            ratU<-c(mseU/mseUmle)
            lratU=log(ratU)
            
          }
          if(doGraph) {
            df=data.frame(x=rep(p, length(n)),
                          RatioMLE = ratMLE,RatioU = ratU,logRatioMLE=lratMLE,logRatioU=lratU,
                          n=factor(rep(n, each=length(p))))
            if(logf)
              plt=ggplot2::ggplot(data=df,aes(x,logRatioMLE, colour=n))+
                geom_line()+theme(axis.text=element_text(size=12),
                                  axis.title=element_text(size=14))#,face="bold")) #+ theme(legend.position="none")+theme(axis.title.y=element_blank())#comment last two for Poisson for 3 panel plot
            else{
              plt=ggplot2::ggplot(data=df,aes(x,RatioMLE, colour=n))+
                geom_line()+theme(axis.text=element_text(size=12),
                                  axis.title=element_text(size=14))#,face="bold")) #+ theme(legend.position="none")+theme(axis.title.y=element_blank())#comment out legent for poisson and binomial, keep Y label only for Poisson
            }
            if(GraphU){
              if(logf) 
                plt2=ggplot2::ggplot(data=df,aes(x,logRatioU, colour=n))+
                  geom_line()+theme(axis.text=element_text(size=12),
                                    axis.title=element_text(size=14))#,face="bold"))
              else {
                plt2=ggplot2::ggplot(data=df,aes(x,RatioU, colour=n))+
                  geom_line()+theme(axis.text=element_text(size=12),
                                    axis.title=element_text(size=14))#,face="bold"))
              }
            }
            if(Model=="Poisson") {
              if(logf)
                plt=plt+xlab(expression(theta))+
                  ylab(expression(log(psi[M](theta, n))))+    
                  theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14))#,face="bold"))
              else{
                plt=plt+xlab(expression(theta))+ylab(expression(psi[M](theta, n)))+
                  theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14))#,face="bold"))
              }
              if(GraphU){
                if(logf) 
                  plt2=plt2+xlab(expression(theta))+
                    ylab(expression(log(psi[U](theta, n))))+
                    theme(axis.text=element_text(size=12),
                          axis.title=element_text(size=14))#,face="bold"))
                else{
                  plt2=plt2+xlab(expression(theta))+
                    ylab(expression(psi[U](theta, n)))+theme(axis.text=element_text(size=12),
                                                             axis.title=element_text(size=14))#,face="bold"))
                }
              }
            }
            #for excess deaths application, use: plt=plt+xlab(expression(nu))+ylab(expression(psi(nu, n)))+xlim(20,50)+ylim(0,15)
            else if (Model=="Skellam")  plt=plt+xlab(expression(theta+beta))+ylab(expression(psi[M](theta,beta, n)))
            else 
              if(logf)
                plt=plt+xlab(expression(phi))+
                  ylab(expression(log(psi[M](theta, n))))+theme(axis.text=element_text(size=12),
                                                                axis.title=element_text(size=14))#,face="bold"))
              else{
                plt=plt+xlab(expression(phi))+
                  ylab(expression(psi[M](theta, n)))+theme(axis.text=element_text(size=12),
                                                           axis.title=element_text(size=14))#,face="bold"))
              }
              print(plt)
              if(GraphU){
                print(plt2)
              }
              return(NULL)
          }
          list(p=p, mseY=mseY, mseU=mseU, mseUmle=mseUmle, RatioMLE=round(mseUmle/mseY, 2),RatioU=round(mseU/mseUmle, 2))
        }

        
#### Here are some runs of the function rror
        n=c(1, 5, 10, 25)
        rror("Poisson", n, c(1, 30), doGraph = TRUE,logf=TRUE)
        rror("Poisson", n, c(1, 30), doGraph = TRUE,GraphU=TRUE, logf=TRUE) 
        rror("Binomial", n,  N=100, doGraph = TRUE,logf=TRUE)
        rror("Negative Binomial", n, interval=c(0.2, 0.9), N=100, doGraph = TRUE,logf=TRUE)
        rror("Skellam", n=42,interval=c(0,44),lambda2=24,doGraph=TRUE) #no log needed here
        rror("Poisson", n=42, c(15, 50), doGraph = TRUE)
   
        a=rror("Binomial", 25,  N=100)
        plot(a$p, a$mseUmle, col="blue", type="line", ylim=c(0, max(a$mseUmle, a$mseY)))
        lines(a$p, a$mseY, col="red")

        
 