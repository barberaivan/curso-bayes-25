data {
  int<lower=0> N;
  int<lower=0> K;  // n efectos fijos
  int<lower=0> np; // n plots

  vector[N] g;

  vector[N] r;
  vector[N] b;
  vector[N] rb;

  vector[N] h;
  vector[N] a;
  vector[N] f;

  array[N] int p; // plot id
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[K] coef;
  vector[np] delta;
  real<lower=0> sigma_ranef; // desvío del efecto aleatorio de la parcela
  real<lower=0> sigma_g;     // desvío de la variable respuesta
  real<lower=0> nu;          // grados de libertad de la t de student
}

transformed parameters {
  vector[N] fixef =
    coef[1] +                                  // intercept
    coef[2] * r + coef[3] * b + coef[4] * rb + // pred categóricas
    coef[5] * h + coef[6] * a + coef[7] * f;   // pred continuas

  vector[N] ranef;
  // efectos aleatorios asignados a cada observación, según su parcela

  vector[N] mu;

  for (n in 1:N) {
    ranef[n] = delta[p[n]];
  }

  mu = fixef + ranef;
}

model {
  // Previas

  // Efectos fijos
  coef ~ normal(0, 1000);

  // Efectos aleatorios
  delta ~ normal(0, sigma_ranef);
  sigma_ranef ~ normal(0, 1000);

  // Desvío residual
  sigma_g ~ normal(0, 1000);

  // Grados de libertad
  nu ~ gamma(2, 0.1);
  // Previa sugerida en
  // https://statmodeling.stat.columbia.edu/2015/05/17/do-we-have-any-recommendations-for-priors-for-student_ts-degrees-of-freedom-parameter/

  // Verosimilitud
  g ~ student_t(nu, mu, sigma_g);
}