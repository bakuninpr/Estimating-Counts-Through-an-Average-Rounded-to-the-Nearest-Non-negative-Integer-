##*************
##*
##* Sanity check simulation to check Results.
##* speciffically, this code assesses expected value 
##*
##***************


pu.poisson <- function(n,u,lambda){
  hu <- ceiling(u-n/2)##ifelse(n%%2==0,u-n/2+1,u-n/2)
  gu <- ifelse(u==0,floor(n/2),0)
  q <- 0:(n - (gu+1))
  sum(dpois(x = hu+q+gu,lambda = lambda))
  
}

set.seed(787)
m <- 100
y <- rpois(n = m,lambda=2) ## generate poisson
n <- 3
u1 <- n*(floor(y/n+0.5))

n <- 10
u2 <- n*(floor(y/n+0.5))
## simulation
hist(y,probability = TRUE)
hist(u1,probability = TRUE)
hist(u2,probability = TRUE)


## Plot PMF
y <- seq(0,15,1); lambda <- 2
prob <- dpois(y,lambda) #RRS aqui cambie a lambda para que el codigo sea mas generico
# Plot PMF
plot(y,prob,type="h",xlab="Y",
     ylab="P(Y=y)",main=bquote(lambda == .(lambda)))#main=(sprintf("Y~Poisson(theta = %1.f)",lambda))) #RRS simplifique titulo y corregi de theta a lambda



y <- seq(0,15,1);
n <- 3
u <- cbind(n * round(y/n,digits=0))
prob <- apply(u,1,pu.poisson, n = n,lambda=2)

# Plot PM
plot(u,prob,type="h",xlab="U",ylab="P(U=u)",main=("U~f_u"))

sum(u*prob)
###**************
###* RRS the code below numerically assesses the theoretical expected 
###* value and variance of U when Y ~ Poi(theta)
###* code below replaces lambda with theta to better reflect 
###* paper's notation.
n <- 2
theta <- 5
mean.and.var.u <- function(n,theta){
  r <- ifelse(n%%2==0,1,0.5)
  j <- if(n==1) 0 else 1:(n-1)# RRS replaced 1:(n-1)
  omegaj <- (exp(0+(pi/n)*2i))^(j)#for computational efficiency
  
  if(n%%2==0){
    aj <- (-1)^j
  } else{
    aj <- (-1)^j*(exp(0+(pi/n)*2i))^(j/2)
  }
  sum.complex <- sum(aj* exp(theta/omegaj)/(1-omegaj))
  sum.complex.term2 <- sum(aj*exp(theta/omegaj)/(1-omegaj)*(2*exp(-theta)/(1-omegaj) - exp(-theta)))
  ## expected value of U for Poisson
  mean.u <- theta+0.5*(2*r-1)+ exp(-theta)*sum.complex
  ## variance of u for Poisson
  
  var.u <- theta+(n^2-1)/12- exp(-2*theta)*sum.complex^2-sum.complex.term2
  
  
  c("E(U)" = Re(mean.u),"Var(U)" = Re(var.u),"eComplex.Term.1"=Re(exp(-theta)*sum.complex),"Complex.Term.2"=Re(sum.complex.term2))
}

#### RRS computation below matches paper when n=2
mean.and.var.u(n=2,theta = 0.1)
theta+1/2-exp(-2*theta)/2 #E(U) according to paper when n=2, theta=.1
0.1+0.25-exp(-.1*4)/4#Var(U) according to paper when n=2, theta=.1


## for large n (RRS relative to theta) it seems that the variance is Re(Var) = 2theta 
## it appears that the Re(E(U)) = 0. This happens as long as theta/n=lambda is small enough.... 0.3?


mean.and.var.u(n = 100,theta=5)
##### For n=100, theta=5 the last three terms of Var(U) sum to 5.
(100^2-1)/12-5.5^2-798



mean.and.var.u(n = 1,theta=10) #this confirms that Corollary does not work for n=1 (in that case U has Poisson)

mean.and.var.u(n = 20,theta=25.5) #RRS; Var is rather bigger than theta

##### What about when lambda (=theta/n) >> n
mean.and.var.u(n = 5,theta=1000) #for lambda >> n  Lemma 2 is validated 
