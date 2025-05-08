data {
  int N;
  array[N] int y; // Número de incendios
  vector[N] fwi;
  vector[N] pp;
  vector[N] temp;

  // Datos para validación cruzada
  int O; // número de observaciones dejadas afuera
  array[N-O] int train_ids; // id de las observaciones de entrenamiento
  array[O] int test_ids; // id de las observacione de prueba
}

// Parámetros a muestrear
parameters {
  real alpha;
  real beta1;
  real<lower=0> phi; // Precisión de la binomial negativa (infinito = poisson)
}

// Calculamos las cantidades derivadas
transformed parameters {
  vector[N] lambda = exp(
    alpha + beta1 * fwi
  );
}

// Acá definimos la log densidad posterior (o la que sea)
model {
  // Densidad previa
  alpha ~ normal(log(10), 2); // Centrada cerca de la media de los datos
  beta1 ~ normal(0, 1);
  phi ~ gamma(2, 0.1); // Recommended by ChatGPT

  // Verosimilitud de las observaciones de entrenamiento
  for (i in train_ids) {
    y[i] ~ neg_binomial_2(lambda[i], phi);
  }

}

generated quantities {
  vector[O] log_lik;
  for (i in 1:O) {
    int j = test_ids[i];
    log_lik[i] = neg_binomial_2_lpmf(y[j] | lambda[j], phi);
  }
}
