data {
  int<lower=0> N; // nro de observaciones
  int<lower=0> S; // semillas sembradas (constante)
  array[N] int y;

  vector[N] D;
  vector[N] E;
  vector[N] H;
}

parameters {
  real alpha;
  real beta_D;
  real beta_E;
  real beta_H;
  real log_phi;
}

transformed parameters {
  real phi = exp(log_phi);
  vector[N] mu = inv_logit(
    alpha + beta_D * D + beta_E * E + beta_H * H
  );

  // Parámetros de la Beta-Binomial (parte beta)
  vector[N] a = mu * phi;
  vector[N] b = (1 - mu) * phi;
}

model {
  // Previas
  alpha ~ normal(0, 10);
  beta_D ~ normal(0, 5);
  beta_E ~ normal(0, 5);
  beta_H ~ normal(0, 5);
  log_phi ~ normal(log(7), 1.5);

  // Verosimilitud
  y ~ beta_binomial(S, a, b);
}

generated quantities {
  // Point-wise log verosimilitud (para comparar modelos y otras cosas)
  vector[N] log_lik;
  for (n in 1:N) log_lik[n] = beta_binomial_lpmf(y[n] | S, a, b);
}