data {
  int<lower=0> N1;
  int<lower=0> N2;
  real time_int[N1 + N2];
  int<lower=0> n[N1];
  int<lower=0> y[N1];
  
  real<lower=0> ig_alpha;
  real<lower=0> ig_beta;
  vector[N1] x1;
  vector[N1] x2;
  vector[N1] x3;
}

transformed data {
  real delta = 1e-9;
  int<lower=1> N = N1 + N2;
}

parameters {
  real<lower=0> elleq;
  real<lower=0> sigma;
  vector[N] eta;
  real beta1;
  real beta2;
  real beta3;
}

transformed parameters {
  vector[N] z;
  {
    matrix[N, N] L_K;
    matrix[N, N] K = cov_exp_quad(time_int, sigma, elleq);

    for (i in 1:N)
      K[i, i] = K[i, i] + delta;

    L_K = cholesky_decompose(K);
    z = L_K * eta;
  }
}

model {
  elleq ~ inv_gamma(ig_alpha, ig_beta);
  sigma ~ normal(0, 0.5);
  eta ~ std_normal();
  beta1 ~ normal(-2, 1);
  beta2 ~ normal(-2, 1);
  beta3 ~ normal(-2, 1);
  for (i in 1:N1){
      y[i] ~ binomial(n[i], Phi(z[i] + beta1*x1[i] + beta2 * x2[i] + beta3 * x3[i])); 
  }
}

generated quantities {
  vector[N1] log_lik;
  for (i in 1:N1)
    log_lik[i] = binomial_lpmf(y[i]| n[i], Phi(z[i] + beta1*x1[i] + beta2 * x2[i] + beta3 * x3[i]));
}
