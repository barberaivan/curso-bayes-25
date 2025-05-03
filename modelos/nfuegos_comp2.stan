data {
  int N;
  array[N] int y; // Número de incendios
  vector[N] fwi;
  vector[N] pp;
  vector[N] temp;
}

// Parámetros a muestrear
parameters {
  real alpha;
  real beta1;
  real beta2;
  real beta3;
  real<lower=0> phi; // Precisión de la binomial negativa (infinito = poisson)
}

// Calculamos las cantidades derivadas
transformed parameters {
  vector[N] lambda = exp(
    alpha + beta1 * fwi + beta2 * pp + beta3 * temp
  );
}

// Acá definimos la log densidad posterior (o la que sea)
model {
  // Densidad previa
  alpha ~ normal(log(10), 2); // Centrada cerca de la media de los datos
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  beta3 ~ normal(0, 1);
  phi ~ gamma(2, 0.1); // Recommended by ChatGPT

  // Verosimilitud
  y ~ neg_binomial_2(lambda, phi);
}

generated quantities {
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = neg_binomial_2_lpmf(y[i] | lambda[i], phi);
  }
}