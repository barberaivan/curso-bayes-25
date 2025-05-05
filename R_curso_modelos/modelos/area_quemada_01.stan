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
  vector[K] a1;
  vector[K] b1;
  vector[K] a2;
  vector[K] b2;
  vector<lower=0>[K] phi; // precisión de la beta
}

transformed parameters {
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
  // Previas (muy poco informativas)
  a1 ~ normal(0, 10);
  a2 ~ normal(0, 10);

  b1 ~ normal(0, 5);
  b2 ~ normal(0, 5);

  phi ~ lognormal(log(10), 1.5);

  // Verosimilitud

  // Para la presencia de fuego:
  fire ~ bernoulli(p1);

  // Para la proporción quemada, sólo cuando algo se quema
  prop[rows_fire] ~ beta_proportion(p2[rows_fire], phi_vec[rows_fire]);
}