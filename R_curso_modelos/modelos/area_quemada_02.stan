data {
  int<lower=0> N;  // observaciones
  int<lower=0> K;  // tipos de vegetación
  int<lower=0> B;  // filas con ocurrencia de fuego

  vector[N] prop;
  array[N] int fire;
  vector[N] fwiz;
  array[N] int veg;

  array[B] int rows_fire; // id de las filas con fuego
}

/*
  Los parámetros serán a = intercepto, b = pendiente;
  1 indica modelo de ocurrencia (bernoulli, fire = 1);
  2 indica proporción quemada
*/
parameters {
  // Parámetros de cada tipo de vegetación
  vector[K] a1;
  vector[K] b1;
  vector[K] a2;
  vector[K] b2;
  vector[K] log_phi; // precisión de la beta

  // Hiperparámetros
  real a1_mu;
  real b1_mu;
  real<lower=0> a1_sd;
  real<lower=0> b1_sd;

  real a2_mu;
  real b2_mu;
  real<lower=0> a2_sd;
  real<lower=0> b2_sd;

  real log_phi_mu;
  real<lower=0> log_phi_sd;
}

transformed parameters {
  vector[K] phi = exp(log_phi); // precisión de la beta

  vector[N] p1;  // probabilidad de presencia
  vector[N] p2;  // proporción quemada esperada
  vector[N] phi_vec; // dispersión de la beta

  for (n in 1:N) {
    p1[n] = inv_logit(a1[veg[n]] + b1[veg[n]] * fwiz[n]);
    p2[n] = inv_logit(a2[veg[n]] + b2[veg[n]] * fwiz[n]);
    phi_vec[n] = phi[veg[n]];
  }
}

model {
  // Previas
  // Los _sd tienen previa truncada en cero porque se definen con <lower=0>
  a1 ~ normal(a1_mu, a1_sd);
  a1_mu ~ normal(0, 10);
  a1_sd ~ normal(0, 3);

  a2 ~ normal(a2_mu, a2_sd);
  a2_mu ~ normal(0, 10);
  a2_sd ~ normal(0, 3);

  b1 ~ normal(b1_mu, b1_sd);
  b1_mu ~ normal(0, 5);
  b1_sd ~ normal(0, 3);

  b2 ~ normal(b2_mu, b2_sd);
  b2_mu ~ normal(0, 5);
  b2_sd ~ normal(0, 3);

  log_phi ~ normal(log_phi_mu, log_phi_sd);
  log_phi_mu ~ normal(log(10), 1.5);
  log_phi_sd ~ normal(0, 2);

  // Verosimilitud

  // Para la presencia de fuego:
  fire ~ bernoulli(p1);

  // Para la proporción quemada, sólo cuando algo se quema
  prop[rows_fire] ~ beta_proportion(p2[rows_fire], phi_vec[rows_fire]);
}