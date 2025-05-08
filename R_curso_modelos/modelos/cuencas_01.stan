data {
  int nobs;
  int nc;
  vector[nobs] caudal;
  array[nobs] int cuenca;
  vector[nobs] semana;    // mejor que sea real
}

parameters {
  vector[nc] log_a;
  vector[nc] log_l;
  vector[nc] log_k;
  vector[nc] log_phi;     // dispersión
}

transformed parameters {
  vector[nc] a = exp(log_a);
  vector[nc] l = exp(log_l);
  vector[nc] k = exp(log_k) + 1; // límite inferior es 1
  vector[nc] phi = exp(log_phi);

  vector<lower=0>[nobs] mu;
  for (i in 1:nobs) {
    mu[i] = a[cuenca[i]] * exp(-(semana[i] / l[cuenca[i]]) ^ k[cuenca[i]]);
  }
}

model {
  // Previas
  log_a ~ normal(log(15), 5);
  log_l ~ normal(log(2.5), 1.5);
  log_k ~ normal(log(2), 0.7);
  log_phi ~ normal(log(0.3), 1.5);

  // Verosimilitud
  {
    //vector[nobs] alpha = rep_vector(phi ^ (-1), nobs);
    //vector[nobs] beta = (mu * phi) ^ (-1);
    real alpha;
    real beta;

    for (i in 1:nobs) {
      alpha = 1 / phi[cuenca[i]];
      beta = 1 / (mu[i] * phi[cuenca[i]]);
      caudal[i] ~ gamma(alpha, beta);
    }
    //caudal ~ gamma(alpha, beta);
  }
}
