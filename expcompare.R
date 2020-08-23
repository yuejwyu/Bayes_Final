
################# 基于指数分布的比较

########### data generate

data.generate.exp = function(n){
  data = rexp(n,1/3)
  return(data)
}

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
  cp.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp(n)
    result = classical.exp.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  cp.est[n*estimate/3>=qgamma(0.025,n,1)&
            n*estimate/3<=qgamma(0.975,n,1)] = 1
  Est = mean(estimate)
  Bias = Est - 3 
  SSE = sd(estimate)
  ESE = sqrt(mean(mse.est))
  MSE = sum((estimate-3)^2)/N
  CP = mean(cp.est)
  return(list(Est = Est, Bias = Bias, SSE = SSE, ESE = ESE, 
              MSE = MSE, CP = CP))
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
  cp.est = rep(0,N)
  for(i in 1:N){
    data1 = data.generate.exp(n)
    result = bayes.exp.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
    palpha[i] = result$postalpha
    plamda[i] = result$postlamda
  }
  cp.est[3>=1/qgamma(0.975, palpha, plamda)&
            3<=1/qgamma(0.025, palpha, plamda)] = 1
  Est = mean(estimate)
  Bias = Est - 3 
  SSE = sd(estimate)
  ESE = sqrt(mean(mse.est))
  MSE = sum((estimate-3)^2)/N
  palpha.est = mean(palpha)
  plamda.est = mean(plamda)
  palpha.se = sd(palpha)
  plamda.se = sd(plamda)
  CP = mean(cp.est)
  return(list(Est = Est, Bias = Bias, SSE = SSE, ESE = ESE, 
              MSE = MSE, palpha.est = palpha.est, 
              plamda.est = plamda.est, palpha.se = palpha.se,
              plamda.se = plamda.se, CP = CP))
}

###### result

n = 100

set.seed(666)
result.exp.classical = classical.exp.main(n,1000)
result.exp.bayes = bayes.exp.main(n,1000)

###### compare

result.exp.compare = cbind(result.exp.classical,
                           list(result.exp.bayes$Est, result.exp.bayes$Bias,
                                result.exp.bayes$SSE, result.exp.bayes$ESE,
                                result.exp.bayes$MSE, result.exp.bayes$CP))
write.csv(result.exp.compare,"expcompare.csv")
result.exp.compare = as.data.frame(result.exp.compare)
print(result.exp.compare)
#print(xtable::xtable(result.exp.compare))

x = 1/rgamma(1000, result.exp.bayes$palpha.est, result.exp.bayes$plamda.est)
hist(x, probability = T)
lines(density(x))