data {
  int<lower=0> N; // nro de observaciones
  int<lower=0> S; // semillas sembradas (constante)
  array[N] int y; // nro germinadas

  vector[N] D; // daño por fuego
}

parameters {
  real alpha;
  real beta;
}

transformed parameters {
  vector[N] mu = inv_logit(alpha + beta * D); // plogis en R
}

model {
  // Previas
  alpha ~ normal(0, 10);
  beta ~ normal(0, 5);

  // Verosimilitud
  y ~ binomial(S, mu);
}
