

################# 基于随机参数的比较

########### data generate

data.generate.exp = function(n){
  data = rep(NA,n)
  u = runif(n,0,1)
  data[u<=1/3] = rexp(sum(u<=1/3),1/2)
  data[u>1/3] = rexp(sum(u>1/3),1/3)
  return(data)
}

data.generate.rate2 = function(n){
  data = rep(NA,n)
  u = runif(n,0,1)
  data[u<=1/3] = rbinom(sum(u<=1/3),1,1/3)
  data[u>1/3] = rbinom(sum(u>1/3),1,2/3)
  return(data)
}

############ exp

###### classical

classical.exp.est = function(data){
  n = length(data)
  est = mean(data)
  mse = est^2/n
  return(list(est = est, mse = mse))
}

classical.exp.main = function(n,N){
  estimate = rep(0,N)
  mse.est = rep(0,N)
  cp1.est = rep(0,N)
  cp2.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp(n)
    result = classical.exp.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  cp1.est[n*estimate/2>=qgamma(0.025,n,1)&
            n*estimate/2<=qgamma(0.975,n,1)] = 1
  cp2.est[n*estimate/3>=qgamma(0.025,n,1)&
            n*estimate/3<=qgamma(0.975,n,1)] = 1
  Est = mean(estimate)
  SSE = sd(estimate)
  ESE = mean(sqrt(mse.est))
  CP1 = mean(cp1.est)
  CP2 = mean(cp2.est)
  return(list(Est = Est, SSE = SSE, ESE = ESE, CP1 = CP1, CP2 = CP2))
}

classical.exp.main2 = function(n,N){
  estimate = rep(0,N)
  mse.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp2(n)
    result = classical.exp.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  Est = mean(estimate)
  SSE = sd(estimate)
  ESE = mean(sqrt(mse.est))
  return(list(Est = Est, SSE = SSE, ESE = ESE))
}

###### bayes

bayes.exp.est = function(data){
  n = length(data)
  postalpha = n
  postlamda = sum(data)
  est = sum(data)/(n-1)
  mse = est^2/(n-2)
  return(list(postalpha = postalpha, postlamda = postlamda,
              est = est, mse = mse))
}

bayes.exp.main = function(n,N){
  estimate = rep(0,N)
  mse.est = rep(0,N)
  palpha = rep(0,N)
  plamda = rep(0,N)
  cp1.est = rep(0,N)
  cp2.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp(n)
    result = bayes.exp.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
    palpha[i] = result$postalpha
    plamda[i] = result$postlamda
  }
  cp1.est[2>=1/qgamma(0.975, palpha, plamda)&
            2<=1/qgamma(0.025, palpha, plamda)] = 1
  cp2.est[3>=1/qgamma(0.975, palpha, plamda)&
            3<=1/qgamma(0.025, palpha, plamda)] = 1
  Est = mean(estimate)
  SSE = sd(estimate)
  ESE = mean(sqrt(mse.est))
  palpha.est = mean(palpha)
  plamda.est = mean(plamda)
  palpha.se = sd(palpha)
  plamda.se = sd(plamda)
  CP1 = mean(cp1.est)
  CP2 = mean(cp2.est)
  return(list(Est = Est, SSE = SSE, ESE = ESE, palpha.est = palpha.est, 
              plamda.est = plamda.est, palpha.se = palpha.se,
              plamda.se = plamda.se, CP1 = CP1, CP2 = CP2))
}

bayes.exp.main2 = function(n,N){
  pai1.est = rep(0,N)
  pai2.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp(n)
    result = bayes.exp.est(data1)
    pai1.est[i] = result$pai1
    pai2.est[i] = result$pai2
  }
  p1.est = mean(pai1.est)
  p2.est = mean(pai2.est)
  p1.se = sd(pai1.est)
  p2.se = sd(pai2.est)
  return(list(p1.est = p1.est, p2.est = p2.est,
              p1.se = p1.se, p2.se = p2.se))
}

###### result

n = 100

set.seed(666)
result.exp.classical = classical.exp.main(n,1000)
result.exp.bayes = bayes.exp.main(n,1000)
result2.exp.bayes = bayes.exp.main2(40,1000)

###### compare

result.exp.compare = cbind(result.exp.classical,
                        list(result.exp.bayes$Est,result.exp.bayes$SSE,
                             result.exp.bayes$ESE,result.exp.bayes$CP1,
                             result.exp.bayes$CP2))
write.csv(result.exp.compare,"expcompare.csv")
result.exp.compare = as.data.frame(result.exp.compare)
print(result.exp.compare)
#print(xtable::xtable(result.exp.compare))

x = 1/rgamma(1000, result.exp.bayes$palpha.est, result.exp.bayes$plamda.est)
hist(x, probability = T)
lines(density(x))

############## rate

#source("ratemain.R")

###### classical

classical.rate.main = function(n,N){
  estimate = rep(0, N)
  mse.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.rate2(n)
    result = classical.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  Est = mean(estimate)
  SSE = sd(estimate)
  ESE = mean(sqrt(mse.est))
  return(list(Est = Est, SSE = SSE, ESE = ESE))
}

###### bayes


bayes.rate.main = function(n,N){
  estimate = rep(0, N)
  mse.est = rep(0, N)
  pa = rep(0,N)
  pb = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.rate2(n)
    result = bayes.est(data1,1,1)
    estimate[i] = result$est
    mse.est[i] = result$mse
    pa[i] = result$posta
    pb[i] = result$postb
  }
  Est = mean(estimate)
  SSE = sd(estimate)
  ESE = mean(sqrt(mse.est))
  pa.est = mean(pa)
  pb.est = mean(pb)
  pa.se = sd(pa)
  pb.se = sd(pb)
  return(list(Est = Est, SSE = SSE, ESE = ESE, pa.est = pa.est, 
              pb.est = pb.est, pa.se = pa.se, pb.se = pb.se))
}

###### result

n = 100

set.seed(666)
result.rate.classical = classical.rate.main(n,1000)
result.rate.bayes = bayes.rate.main(n,1000)

###### compare

result.rate.compare = cbind(result.rate.classical,
                           list(result.rate.bayes$Est,result.rate.bayes$SSE,
                                result.rate.bayes$ESE))
write.csv(result.rate.compare,"ratecompare.csv")
result.rate.compare = as.data.frame(result.rate.compare)
print(result.rate.compare)
#print(xtable::xtable(result.rate.compare))

x = rbeta(1000, result.rate.bayes$pa.est, result.rate.bayes$pb.est)
hist(x, probability = T)
lines(density(x))