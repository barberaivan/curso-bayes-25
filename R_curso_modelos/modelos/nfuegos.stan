data {
  int N;
  array[N] int y; // Número de incendios
  vector[N] x;    // FWI
}

// Parámetros a muestrear
parameters {
  real alpha;
  real beta;
}

// Calculamos las cantidades derivadas
transformed parameters {
  vector[N] lambda = exp(alpha + beta * x);
}

// Acá definimos la log densidad posterior (o la que sea)
model {
  // Densidad previa
  alpha ~ normal(0, 1);
  beta ~ normal(0, 0.1);

  // Verosimilitud
  y ~ poisson(lambda);
}