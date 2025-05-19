data {
  int N;          // número de observaciones (plazas), 50
  vector[N] y;    // número de hill q1
  vector[N] x;    // cobertura vegetal
  // array[N] real x; un vector es un array de reales
}

transformed data {
  real log_mean_y = log(mean(y));
}

// Parámetros a muestrear
parameters {
  real alpha;
  real beta;
  real<lower=0> phi;
}

// Calculamos las cantidades derivadas
transformed parameters {
  vector[N] mu = exp(alpha + beta * x);

  // Reparametrizamos la gamma con a y b
  real a = 1 / phi;
  vector[N] b = (mu * phi) ^ (-1);
}

// Acá definimos la log densidad posterior (o la que sea)
model {
  // Acá va todo lo que tiene virgulilla

  // Densidad previa
  alpha ~ normal(log_mean_y, 10);
  beta ~ normal(0, 2);
  phi ~ gamma(2, 0.1);

  // Verosimilitud
  y ~ gamma(a, b);
}
