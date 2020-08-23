

##################### 基于正态总体的比较

############### data generate

data.generate.norm = function(n,mu,sigma){
  return(rnorm(n, mu, sigma))
}

############### 方差已知

###### classical

classical.est1 = function(data,sigma){
  n = length(data)
  est = mean(data)
  se = sigma/sqrt(n)
  return(list(est = est, se = se))
}

classical.main1 = function(n,mu,sigma,N){
  estimate = rep(0, N)
  se.est = rep(0, N)
  cp.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = classical.est1(data1,sigma)
    estimate[i] = result$est
    se.est[i] = result$se
  }
  cp.est[mu>=estimate-qnorm(0.975, 0, 1)*se.est&
           mu<=estimate+qnorm(0.975, 0, 1)*se.est] = 1
  Est = mean(estimate)
  Bias = Est - mu
  SSE = sd(estimate)
  ESE = sqrt(mean(se.est^2))
  MSE = sum((estimate-mu)^2)/N
  CP = mean(cp.est)
  return(list(Est = Est, Real = mu, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE, CP = CP))
}

classical.test1 = function(n,mu,sigma,N,mu0){
  test = rep(NA, length(mu0))
  for(i in 1:length(mu0))
  {
    estimate = rep(0, N)
    se.est = rep(0, N)
    utest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = classical.est1(data1, sigma)
      estimate[j] = result$est
      se.est[j] = result$se
    }
    utest[(estimate-mu0[i])/se.est>=-qnorm(0.95, 0, 1)] = 1
    test[i] = mean(utest)
  }
  return(test)
}

###### bayes

bayes.est1 = function(data,sigma,mu1,sigma1){
  n = length(data)
  sigma0 = sigma/sqrt(n)
  a = 1/(sigma0^2)
  b = 1/(sigma1^2)
  est = (a*mean(data)+b*mu1)/(a+b)
  se = sqrt(1/(a+b))
  return(list(est = est, se = se))
}

bayes.main1 = function(n,mu,sigma,mu1,sigma1,N){
  estimate = rep(0, N)
  se.est = rep(0, N)
  cp.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = bayes.est1(data1,sigma,mu1,sigma1)
    estimate[i] = result$est
    se.est[i] = result$se
  }
  cp.est[mu>=estimate-qnorm(0.975, 0, 1)*se.est&
           mu<=estimate+qnorm(0.975, 0, 1)*se.est] = 1
  Est = mean(estimate)
  Bias = Est - mu
  SSE = sd(estimate)
  ESE = sqrt(mean(se.est^2))
  MSE = sum((estimate-mu)^2)/N
  CP = mean(cp.est)
  return(list(Est = Est, Real = mu, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE, CP = CP))
}

bayes.test1 = function(n,mu,sigma,mu1,sigma1,N,mu0){
  test = rep(NA, length(mu0))
  for(i in 1:length(mu0))
  {
    estimate = rep(0, N)
    se.est = rep(0, N)
    btest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = bayes.est1(data1, sigma, mu1, sigma1)
      estimate[j] = result$est
      se.est[j] = result$se
    }
    btest[(1-pnorm(mu0[i],estimate,se.est))/pnorm(mu0[i],estimate,se.est)>=1]=1
    test[i] = mean(btest)
  }
  return(test)
}

###### result

n = 100
mu0 = seq(2,3,0.01)

set.seed(666)
result1.classical = classical.main1(n,2.5,1.5,1000)
set.seed(666)
test1.classical = classical.test1(n,2.5,1.5,1000,mu0)

mu.history = c(2.15,2.45,2.75,2.25,2.65)
mu1 = mean(mu.history)
sigma1 = sd(mu.history)
set.seed(666)
result1.bayes = bayes.main1(n,2.5,1.5,mu1,sigma1,1000)
set.seed(666)
test1.bayes = bayes.test1(n,2.5,1.5,mu1,sigma1,1000,mu0)

###### compare

result1.compare = cbind(result1.classical,
                       result1.bayes)
write.csv(result1.compare,"norm compare1.csv")
result1.compare = as.data.frame(result1.compare)
cat("campare1")
cat("\n")
print(result1.compare)
#print(xtable::xtable(result1.compare))

plot(mu0, test1.classical, type = "l", lty = 1, lwd = 2,  
     col = "chocolate", xlim = c(mu0[1], mu0[length(mu0)]), ylim = c(0,1),
     xlab = "mu0", ylab = "检验通过率")
lines(mu0, test1.bayes, lty = 2, lwd = 2, col = "darkgreen")
abline(v = 2.5, lty = 2)
abline(h = 0.05, lty = 2)
legend(2.75, 0.9,
       legend = c("Classical","Bayes"),
       lty = 1:5, title = "类型", 
       col = c("chocolate", "darkgreen"),
       cex = 0.6, lwd = 2)

############### 均值已知

####### classical

classical.est2 = function(data,mu){
  n = length(data)
  est = sum((data-mu)^2)/n
  se = sqrt(2*est^4/n)
  return(list(est = est, se = se))
}

classical.main2 = function(n,mu,sigma,N){
  estimate = rep(0, N)
  se.est = rep(0, N)
  cp.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = classical.est2(data1,mu)
    estimate[i] = result$est
    se.est[i] = result$se
  }
  cp.est[sigma^2>=n*estimate/qchisq(0.975, n)&
           sigma^2<=n*estimate/qchisq(0.025, n)] = 1
  Est = mean(estimate)
  Bias = Est - sigma^2
  SSE = sd(estimate)
  ESE = sqrt(mean(se.est^2))
  MSE = sum((estimate-sigma^2)^2)/N
  CP = mean(cp.est)
  return(list(Est = Est, Real = sigma^2, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE, CP = CP))
}

classical.test2 = function(n,mu,sigma,N,sigma0){
  test = rep(NA, length(sigma0))
  for(i in 1:length(sigma0))
  {
    estimate = rep(0, N)
    ctest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = classical.est2(data1, mu)
      estimate[j] = result$est
    }
    ctest[n*estimate/sigma0[i]^2>=qchisq(0.05, n)] = 1
    test[i] = mean(ctest)
  }
  return(test)
}

###### bayes

bayes.est2 = function(data,mu,alpha,lamda){
  n = length(data)
  postalpha = alpha+n/2
  postlamda = lamda+(sum((data-mu)^2)/2)
  est = postlamda/(postalpha-1)
  se = est/sqrt(postalpha-2)
  return(list(est = est, se = se, postalpha = postalpha, postlamda = postlamda))
}

bayes.main2 = function(n,mu,sigma,alpha,lamda,N){
  estimate = rep(0, N)
  se.est = rep(0, N)
  palpha = rep(0, N)
  plamda = rep(0, N)
  cp.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = bayes.est2(data1,mu,alpha,lamda)
    estimate[i] = result$est
    se.est[i] = result$se
    palpha[i] = result$postalpha
    plamda[i] = result$postlamda
  }
  cp.est[sigma^2>=1/qgamma(0.975, shape = palpha, scale = 1/plamda)&
           sigma^2<=1/qgamma(0.025, shape = palpha, scale = 1/plamda)] = 1
  Est = mean(estimate)
  Bias = Est - sigma^2
  SSE = sd(estimate)
  ESE = sqrt(mean(se.est^2))
  MSE = sum((estimate-sigma^2)^2)/N
  CP = mean(cp.est)
  return(list(Est = Est, Real = sigma^2, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE, CP = CP))
}

bayes.test2 = function(n,mu,sigma,alpha,lamda,N,sigma0){
  test = rep(NA, length(sigma0))
  for(i in 1:length(sigma0))
  {
    btest = rep(0, N)
    palpha = rep(0, N)
    plamda = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = bayes.est2(data1, mu, alpha,lamda)
      palpha[j] = result$postalpha
      plamda[j] = result$postlamda
    }
    btest[pgamma(1/sigma0[i]^2,shape=palpha,scale=1/plamda)/
            (1-pgamma(1/sigma0[i]^2,shape=palpha,scale=1/plamda))>=1]=1
    test[i] = mean(btest)
  }
  return(test)
}

###### result

n = 100
sigma0 = seq(1,2,0.01)

set.seed(666)
result2.classical = classical.main2(n,2.5,1.5,1000)
set.seed(666)
test2.classical = classical.test2(n,2.5,1.5,1000,sigma0)

sigma.history = c(1.3,1.6,1.7,1.9,1.1)
alpha = ((mean(sigma.history^2))^2/var(sigma.history^2))+2
lamda = mean(sigma.history^2)*(alpha-1)
set.seed(666)
result2.bayes = bayes.main2(n,2.5,1.5,alpha,lamda,1000)
set.seed(666)
test2.bayes = bayes.test2(n,2.5,1.5,alpha,lamda,1000,sigma0)

###### compare

result2.compare = cbind(result2.classical,
                        result2.bayes)
write.csv(result2.compare,"norm compare2.csv")
result2.compare = as.data.frame(result2.compare)
cat("campare2")
cat("\n")
print(result2.compare)
#print(xtable::xtable(result2.compare))

plot(sigma0, test2.classical, type = "l", lty = 1, lwd = 2,  
     col = "chocolate", xlim = c(sigma0[1], sigma0[length(sigma0)]), ylim = c(0,1),
     xlab = "sigma0", ylab = "检验通过率")
lines(sigma0, test2.bayes, lty = 2, lwd = 2, col = "darkgreen")
abline(v = 1.5, lty = 2)
abline(h = 0.05, lty = 2)
legend(1.75, 0.9,
       legend = c("Classical","Bayes"),
       lty = 1:5, title = "类型", 
       col = c("chocolate", "darkgreen"),
       cex = 0.6, lwd = 2)

############ 均值与方差均未知

###### classical

classical.est3 = function(data){
  n = length(data)
  mu.est = mean(data)
  sqsigma.est = var(data)
  mu.se = sqsigma.est/sqrt(n)
  sqsigma.se = sqrt(2*sqsigma.est^2/(n-1))
  return(list(mu.est = mu.est, sqsigma.est = sqsigma.est,
              mu.se = mu.se, sqsigma.se = sqsigma.se))
}

classical.main3 = function(n,mu,sigma,N){
  mu.estimate = rep(0, N)
  sqsigma.estimate = rep(0, N)
  mu.serror = rep(0, N)
  sqsigma.serror = rep(0, N)
  mu.cp = rep(0, N)
  sqsigma.cp = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = classical.est3(data1)
    mu.estimate[i] = result$mu.est
    sqsigma.estimate[i] = result$sqsigma.est
    mu.serror[i] = result$mu.se
    sqsigma.serror[i] = result$sqsigma.se
  }
  mu.cp[mu>=mu.estimate-qt(0.975, n-1)*mu.serror&
          mu<=mu.estimate+qt(0.975, n-1)*mu.serror] = 1
  sqsigma.cp[sigma^2>=(n-1)*sqsigma.estimate/qchisq(0.975,n-1)&
             sigma^2<=(n-1)*sqsigma.estimate/qchisq(0.025,n-1)] = 1
  mu.Est = mean(mu.estimate)
  mu.Bias = mu.Est - mu
  mu.SSE = sd(mu.estimate)
  mu.ESE = sqrt(mean(mu.serror^2))
  mu.MSE = sum((mu.estimate-mu)^2)/N
  mu.CP = mean(mu.cp)
  sqsigma.Est = mean(sqsigma.estimate)
  sqsigma.Bias = sqsigma.Est - sigma^2
  sqsigma.SSE = sd(sqsigma.estimate)
  sqsigma.ESE = sqrt(mean(sqsigma.serror^2))
  sqsigma.MSE = sum((sqsigma.estimate-sigma^2)^2)/N
  sqsigma.CP = mean(sqsigma.cp)
  return(list(mu.Est = mu.Est, mu.Bias = mu.Bias, mu.SSE = mu.SSE, mu.ESE = mu.ESE,
              mu.MSE = mu.MSE, mu.CP = mu.CP, sqsigma.Est = sqsigma.Est, 
              sqsigma.Bias = sqsigma.Bias, sqsigma.SSE = sqsigma.SSE, 
              sqsigma.ESE = sqsigma.ESE, sqsigma.MSE = sqsigma.MSE,
              sqsigma.CP = sqsigma.CP))
}

classical.test3 = function(n,mu,sigma,N,mu0){
  mu.test = rep(NA, length(mu0))
  for(i in 1:length(mu0))
  {
    mu.estimate = rep(0, N)
    sqsigma.estimate = rep(0, N)
    ctest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = classical.est3(data1)
      mu.estimate[j] = result$mu.est
      sqsigma.estimate[j] = result$sqsigma.est
    }
    ctest[sqrt(n)*(mu.estimate-mu0[i])/sqrt(sqsigma.estimate)>=qt(0.05, n-1)] = 1
    mu.test[i] = mean(ctest)
  }
  return(mu.test)
}

classical.test4 = function(n,mu,sigma,N,sigma0){
  sigma.test = rep(NA, length(sigma0))
  for(i in 1:length(sigma0))
  {
    sqsigma.estimate = rep(0, N)
    ctest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = classical.est3(data1)
      sqsigma.estimate[j] = result$sqsigma.est
    }
    ctest[(n-1)*sqsigma.estimate/(sigma0[i]^2)>=qchisq(0.05,n-1)] = 1
    sigma.test[i] = mean(ctest)
  }
  return(sigma.test)
}

###### bayes

bayes.est3 = function(data,mu1,k1,alpha,lamda){
  n = length(data)
  postmu1 = (k1*mu1+n*mean(data))/(k1+n)
  postk1 = k1+n
  postalpha = alpha+n/2
  postlamda = lamda+((n-1)*var(data)+k1*n*(mu1-mean(data))^2/(k1+n))/2
  mu.est = postmu1
  sqsigma.est = postlamda/(postalpha-1)
  mu.se = sqrt(postlamda/((postalpha-1)*postk1))
  sqsigma.se = sqsigma.est/sqrt(postalpha-2)
  return(list(postmu1 = postmu1, postk1 = postk1, postalpha = postalpha, 
              postlamda = postlamda, mu.est = mu.est, sqsigma.est = sqsigma.est,
              mu.se = mu.se, sqsigma.se = sqsigma.se))
}

bayes.main3 = function(n,mu,sigma,mu1,k1,alpha,lamda,N){
  mu.estimate = rep(0, N)
  sqsigma.estimate = rep(0, N)
  mu.serror = rep(0, N)
  sqsigma.serror = rep(0, N)
  mu.cp = rep(0, N)
  sqsigma.cp = rep(0, N)
  pmu1 = rep(0, N)
  pk1 = rep(0, N)
  palpha = rep(0, N)
  plamda = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.norm(n, mu, sigma)
    result = bayes.est3(data1,mu1,k1,alpha,lamda)
    mu.estimate[i] = result$mu.est
    sqsigma.estimate[i] = result$sqsigma.est
    mu.serror[i] = result$mu.se
    sqsigma.serror[i] = result$sqsigma.se
    pmu1[i] = result$postmu1
    pk1[i] = result$postk1
    palpha[i] = result$postalpha
    plamda[i] = result$postlamda
  }
  mu.cp[mu>=pmu1-qt(0.975,2*palpha)*sqrt(plamda/(palpha*pk1))] = 1
  sqsigma.cp[sigma^2>=1/qgamma(0.975, shape = palpha, scale = 1/plamda)&
               sigma^2<=1/qgamma(0.025, shape = palpha, scale = 1/plamda)] = 1
  mu.Est = mean(mu.estimate)
  mu.Bias = mu.Est - mu
  mu.SSE = sd(mu.estimate)
  mu.ESE = sqrt(mean(mu.serror^2))
  mu.MSE = sum((mu.estimate-mu)^2)/N
  mu.CP = mean(mu.cp)
  sqsigma.Est = mean(sqsigma.estimate)
  sqsigma.Bias = sqsigma.Est - sigma^2
  sqsigma.SSE = sd(sqsigma.estimate)
  sqsigma.ESE = sqrt(mean(sqsigma.serror^2))
  sqsigma.MSE = sum((sqsigma.estimate-sigma^2)^2)/N
  sqsigma.CP = mean(sqsigma.cp)
  return(list(mu.Est = mu.Est, mu.Bias = mu.Bias, mu.SSE = mu.SSE, mu.ESE = mu.ESE,
              mu.MSE = mu.MSE, mu.CP = mu.CP, sqsigma.Est = sqsigma.Est, 
              sqsigma.Bias = sqsigma.Bias, sqsigma.SSE = sqsigma.SSE, 
              sqsigma.ESE = sqsigma.ESE, sqsigma.MSE = sqsigma.MSE,
              sqsigma.CP = sqsigma.CP))
}

bayes.test3 = function(n,mu,sigma,mu1,k1,alpha,lamda,N,mu0){
  mu.test = rep(NA, length(mu0))
  for(i in 1:length(mu0))
  {
    pmu1 = rep(0, N)
    pk1 = rep(0, N)
    palpha = rep(0, N)
    plamda = rep(0, N)
    btest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = bayes.est3(data1,mu1,k1,alpha,lamda)
      pmu1[j] = result$postmu1
      pk1[j] = result$postk1
      palpha[j] = result$postalpha
      plamda[j] = result$postlamda
    }
    btest[(1-pt((mu0[i]-pmu1)/sqrt(plamda/(palpha*pk1)),2*palpha))/
            pt((mu0[i]-pmu1)/sqrt(plamda/(palpha*pk1)),2*palpha)>=1] = 1
    mu.test[i] = mean(btest)
  }
  return(mu.test)
}

bayes.test4 = function(n,mu,sigma,mu1,k1,alpha,lamda,N,sigma0){
  sigma.test = rep(NA, length(sigma0))
  for(i in 1:length(sigma0))
  {
    palpha = rep(0, N)
    plamda = rep(0, N)
    btest = rep(0, N)
    for(j in 1:N){
      data1 = data.generate.norm(n, mu, sigma)
      result = bayes.est3(data1,mu1,k1,alpha,lamda)
      palpha[j] = result$postalpha
      plamda[j] = result$postlamda
    }
    btest[pgamma(1/sigma0[i]^2,shape=palpha,scale=1/plamda)/
            (1-pgamma(1/sigma0[i]^2,shape=palpha,scale=1/plamda))>=1] = 1
    sigma.test[i] = mean(btest)
  }
  return(sigma.test)
}

###### result

n = 100
mu0 = seq(2,3,0.01)
sigma0 = seq(1,2,0.01)

set.seed(666)
result3.classical = classical.main3(n,2.5,1.5,1000)
set.seed(666)
test3.classical = classical.test3(n,2.5,1.5,1000,mu0)
set.seed(666)
test4.classical = classical.test4(n,2.5,1.5,1000,sigma0)

mu.history = c(2.15,2.45,2.75,2.25,2.65)
sigma.history = c(1.3,1.6,1.7,1.9,1.1)
mu1 = mean(mu.history)
alpha = ((mean(sigma.history^2))^2/var(sigma.history^2))+2
lamda = mean(sigma.history^2)*(alpha-1)
k1 = lamda/((alpha-1)*var(mu.history))
set.seed(666)
result3.bayes = bayes.main3(n,2.5,1.5,mu1,k1,alpha,lamda,1000)
set.seed(666)
test3.bayes = bayes.test3(n,2.5,1.5,mu1,k1,alpha,lamda,1000,mu0)
set.seed(666)
test4.bayes = bayes.test4(n,2.5,1.5,mu1,k1,alpha,lamda,1000,sigma0)

###### compare

result3.compare = cbind(result3.classical,
                        result3.bayes)
write.csv(result3.compare,"norm compare3.csv")
result3.compare = as.data.frame(result3.compare)
cat("campare3")
cat("\n")
print(result3.compare)
#print(xtable::xtable(result3.compare))

plot(mu0, test3.classical, type = "l", lty = 1, lwd = 2,  
     col = "chocolate", xlim = c(mu0[1], mu0[length(mu0)]), ylim = c(0,1),
     xlab = "mu0", ylab = "检验通过率")
lines(mu0, test3.bayes, lty = 2, lwd = 2, col = "darkgreen")
abline(v = 2.5, lty = 2)
abline(h = 0.05, lty = 2)
legend(2.75, 0.9,
       legend = c("Classical","Bayes"),
       lty = 1:5, title = "类型", 
       col = c("chocolate", "darkgreen"),
       cex = 0.6, lwd = 2)

plot(sigma0, test4.classical, type = "l", lty = 1, lwd = 2,  
     col = "chocolate", xlim = c(sigma0[1], sigma0[length(sigma0)]), ylim = c(0,1),
     xlab = "sigma0", ylab = "检验通过率")
lines(sigma0, test4.bayes, lty = 2, lwd = 2, col = "darkgreen")
abline(v = 1.5, lty = 2)
abline(h = 0.05, lty = 2)
legend(1.75, 0.9,
       legend = c("Classical","Bayes"),
       lty = 1:5, title = "类型", 
       col = c("chocolate", "darkgreen"),
       cex = 0.6, lwd = 2)