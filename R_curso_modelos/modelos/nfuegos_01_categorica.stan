data {
  int N;
  array[N] int y; // Número de incendios
  vector[N] M;    // FWI medio
  vector[N] A;    // FWI alto
}

// Parámetros a muestrear
parameters {
  real alpha;
  real betaM;
  real betaA;
}

// Calculamos las cantidades derivadas
transformed parameters {
  vector[N] lambda = exp(alpha + betaM * M + betaA * A);
}

// Acá definimos la log densidad posterior (o la que sea)
model {
  // Densidad previa
  alpha ~ normal(0, 10);
  betaM ~ normal(0, 10);
  betaA ~ normal(0, 10);

  // Verosimilitud
  y ~ poisson(lambda);
}
