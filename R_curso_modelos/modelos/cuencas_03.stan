data {
  int nobs;
  int nc;
  vector[nobs] caudal;
  array[nobs] int cuenca;
  vector[nobs] semana;    // mejor que sea real
  vector[nc] rocaz;
}

transformed data {
  vector[nc] rocaz2 = rocaz ^ 2;
}

parameters {
  // Parámetros a nivel de cuenca estandarizados (raw)
  vector[nc] log_a_raw;
  vector[nc] log_l_raw;
  vector[nc] log_k_raw;
  vector[nc] log_phi_raw;     // dispersión

  // Hiperparámetros: media (intercept), pendientes (b1 y b2) y desvío de cada uno
  real log_a_mu;
  real log_l_mu;
  real log_k_mu;
  real log_phi_mu;

  real log_a_b1;
  real log_l_b1;
  real log_k_b1;
  real log_phi_b1;

  real log_a_b2;
  real log_l_b2;
  real log_k_b2;
  real log_phi_b2;

  real<lower=0> log_a_sd;
  real<lower=0> log_l_sd;
  real<lower=0> log_k_sd;
  real<lower=0> log_phi_sd;
}

transformed parameters {
  // media de los parámetros a nivel de cuenca (regresión con rocaz)
  vector[nc] log_a_mean = log_a_mu + log_a_b1 * rocaz + log_a_b2 * rocaz2;
  vector[nc] log_l_mean = log_l_mu + log_l_b1 * rocaz + log_l_b2 * rocaz2;
  vector[nc] log_k_mean = log_k_mu + log_l_b1 * rocaz + log_k_b2 * rocaz2;
  vector[nc] log_phi_mean = log_phi_mu + log_phi_b1 * rocaz + log_phi_b2 * rocaz2;

  // Verdaderos parámetros a nivel de cuenca, no estandarizados, en log
  vector[nc] log_a = log_a_raw * log_a_sd + log_a_mean;
  vector[nc] log_l = log_l_raw * log_l_sd + log_l_mean;
  vector[nc] log_k = log_k_raw * log_k_sd + log_k_mean;
  vector[nc] log_phi = log_phi_raw * log_phi_sd + log_phi_mean;

  // y en la escala de interés
  vector[nc] a = exp(log_a);
  vector[nc] l = exp(log_l);
  vector[nc] k = exp(log_k) + 1; // límite inferior es 1
  vector[nc] phi = exp(log_phi);

  // media a nivel de observación
  vector<lower=0>[nobs] mu;
  for (i in 1:nobs) {
    mu[i] = a[cuenca[i]] * exp(-(semana[i] / l[cuenca[i]]) ^ k[cuenca[i]]);
  }
}

model {
  // Previas para parámetros estandarizados.
  log_a_raw ~ std_normal();
  log_l_raw ~ std_normal();
  log_k_raw ~ std_normal();
  log_phi_raw ~ std_normal();
  /*
    Parametrización no centrada:
    definir
      log_a_raw ~ std_normal(); = normal(0, 1)
    con
      log_a = log_a_raw * sigma + mu;
    implica
      log_a ~ normal(mu, sigma);
    pero mejora la forma de al posterior, facilitando el muestreo.
  */

  // Previas sobre los hiperparámetros
  log_a_mu ~ normal(log(15), 5);
  log_l_mu ~ normal(log(2.5), 1.5);
  log_k_mu ~ normal(log(2), 0.7);
  log_phi_mu ~ normal(log(0.3), 1.5);

  // (pendientes)
  log_a_b1 ~ normal(0, 1);
  log_l_b1 ~ normal(0, 1);
  log_k_b1 ~ normal(0, 1);
  log_phi_b1 ~ normal(0, 1);

  log_a_b2 ~ normal(0, 1);
  log_l_b2 ~ normal(0, 1);
  log_k_b2 ~ normal(0, 1);
  log_phi_b2 ~ normal(0, 1);

  // como log_*_sd está definido con <lower=0>, estas son normales truncadas
  // en cero
  log_a_sd ~ normal(0, 5);
  log_l_sd ~ normal(0, 1.5);
  log_k_sd ~ normal(0, 0.7);
  log_phi_sd ~ normal(0, 1.5);

  // Verosimilitud
  {
    real alpha;
    real beta;

    for (i in 1:nobs) {
      alpha = 1 / phi[cuenca[i]];
      beta = 1 / (mu[i] * phi[cuenca[i]]);
      caudal[i] ~ gamma(alpha, beta);
    }
  }
}
