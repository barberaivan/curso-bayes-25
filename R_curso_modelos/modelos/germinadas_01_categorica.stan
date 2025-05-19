data {
  int<lower=0> N; // nro de observaciones
  int<lower=0> S; // semillas sembradas (constante)
  array[N] int y; // nro germinadas

  vector[N] DM; // daño por fuego medio
  vector[N] DA; // daño por fuego alto

}

parameters {
  real alpha;
  real betaM;
  real betaA;
}

transformed parameters {
  vector[N] mu = inv_logit(alpha + betaM * DM + betaA * DA); // plogis en R
}

model {
  // Previas
  alpha ~ normal(0, 10);
  betaM ~ normal(0, 5);
  betaA ~ normal(0, 5);

  // Verosimilitud
  y ~ binomial(S, mu);
}
