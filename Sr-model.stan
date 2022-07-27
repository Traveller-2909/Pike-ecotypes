data {
  int <lower=0> n; //number of data points
  int <lower=0> k; //number of x-variables
  real <lower=0> Salcapt;
  real <lower=0> sigma;
  matrix[n,k] X;
}
parameters {
  real a;
  vector [k] B;
}
transformed parameters {
  vector <lower=0> [n] mu;
  mu = exp(a+X*B);
}
model {
  Salcapt~normal(mu,sigma);
  a~normal(0,10);
  B~gamma(0.01, 0.01);
}
generated quantities {
  real r_predict [n];
  for(i in 1:n){
    r_predict[i] = poisson_rng(mu[i]);
  }
}