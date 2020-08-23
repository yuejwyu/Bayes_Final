

############### 基于成功概率的比较

###### data generate

data.generate.rate = function(n,theta){
  return(rbinom(n,1,theta))
}

###### classical

classical.est = function(data){
  n = length(data)
  est = mean(data)
  mse = est*(1-est)/n
  return(list(est = est, mse = mse))
}

classical.main = function(n,theta,N){
  estimate = rep(0, N)
  mse.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.rate(n, theta)
    result = classical.est(data1)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  Est = mean(estimate)
  Bias = Est - theta
  SSE = sd(estimate)
  ESE = sqrt(mean(mse.est))
  MSE = sum((estimate-theta)^2)/N
  return(list(Est = Est, Real = theta, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE))
}

###### bayes

bayes.est = function(data,a,b){
  n = length(data)
  x = sum(data)
  est = (a+x)/(a+b+n)
  mse = (a+x)*(b+n-x)/(((a+b+n)^2)*(a+b+n+1))
  posta = a+x
  postb = b+n-x
  return(list(est = est, mse = mse, posta = posta, postb = postb))
}

bayes.main = function(n,theta,a,b,N){
  estimate = rep(0, N)
  mse.est = rep(0, N)
  for(i in 1:N){
    data1 = data.generate.rate(n, theta)
    result = bayes.est(data1,a,b)
    estimate[i] = result$est
    mse.est[i] = result$mse
  }
  Est = mean(estimate)
  Bias = Est - theta
  SSE = sd(estimate)
  ESE = sqrt(mean(mse.est))
  MSE = sum((estimate-theta)^2)/N
  return(list(Est = Est, Real = theta, Bias = Bias, SSE = SSE, 
              ESE = ESE, MSE = MSE))
}

###### result

n = 300 #10,30,300

set.seed(666)
result.classical = classical.main(n,0.8,1000)

set.seed(666)
result.bayes1 = bayes.main(n,0.8,1,1,1000)

theta.history = c(0.72,0.77,0.78,0.83,0.87)
theta.mean = mean(theta.history)
theta.var = var(theta.history)
a2 = theta.mean*(((1-theta.mean)*theta.mean/theta.var)-1)
b2 = (1-theta.mean)*(((1-theta.mean)*theta.mean/theta.var)-1)
set.seed(666)
result.bayes2 = bayes.main(n,0.8,a2,b2,1000)

a3 = 2
b3 = 3
set.seed(666)
result.bayes3 = bayes.main(n,0.8,a3,b3,1000)

###### compare

result.compare = cbind(result.classical,
                       result.bayes1,
                       result.bayes2,
                       result.bayes3)
write.csv(result.compare,"rate compare.csv")
result.compare = as.data.frame(result.compare)
print(result.compare)
#print(xtable::xtable(result.compare))