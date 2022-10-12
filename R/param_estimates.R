param_estimates <- function(formula, data,p_shape){

data <- model.frame(formula, data)

  y <- data[,1]
  group <- data[,2]
  n <-tapply(y,group,length)
  a <- length(n)

  p <- c(rep(p_shape,a))
  k <-2*p-3
  nu <- 2*p-1


  ysort <- t1j <-a1 <- b1 <- c1 <- NULL
  beta1j <- alfa1j <-f <- m <- NULL
  muhat <- B11 <- C11 <- D11 <- sigmahat <- muhat <- NULL

for (i in 1:a){
  ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=sort(y[seq(sum(n[1:i-1])+1,sum(n[1:i]))])

  t1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=tinv((1:n[i])/(n[i]+1),nu[i])*sqrt(k[i]/nu[i])

  a1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=1-((1/k[i])*(t1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])^2);

  b1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=1+((1/k[i])*(t1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])^2)


  c1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=(b1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])^2

  beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=a1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]/c1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]
  #
  alfa1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=(2/k[i])*(((t1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])^3)/c1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  #
  f[i]=sum(beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]*ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  m[i]=sum(beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  muhat[i]=f[i]/m[i]
  B11[i]=(2*p[i]/k[i])*sum((ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]-kronecker(muhat[i],ones(1,n[i])))*alfa1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  C11[i]=(2*p[i]/k[i])*sum(((ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]-kronecker(muhat[i],ones(1,n[i])))^2)*beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]);
  #
  D11[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=kronecker(muhat[i],ones(1,n[i]))
 if(C11[i]<0){
   alfa1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=(1/k[i])*(((t1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])^3)/c1[2,seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]);
    beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]=1/c1[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)];

  }
  f[i]=sum(beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]*ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  m[i]=sum(beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  muhat[i]=f[i]/m[i]
  B11[i]=(2*p[i]/k[i])*sum((ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]-kronecker(muhat[i],ones(1,n[i])))*alfa1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)])
  C11[i]=(2*p[i]/k[i])*sum(((ysort[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]-kronecker(muhat[i],ones(1,n[i])))^2)*beta1j[seq(sum(n[1:i-1])+1,sum(n[1:i]),1)]);
  sigmahat[i]=(B11[i]+sqrt(B11[i]^2+4*n[i]*C11[i]))/(2*sqrt(n[i]*(n[i]-1)));

}
list(muhat=muhat, sigmahat=sigmahat, n=n, a=a, m=m)

}







