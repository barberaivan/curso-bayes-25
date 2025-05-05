rm(list = ls()) # Limpiamos el entorno

# Paquetes ----------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(posterior)
library(viridis)
library(truncnorm)

# Funciones ---------------------------------------------------------------

# ggplot custom theme
nice_theme <- function() {
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major =  element_line(),

    axis.line = element_line(linewidth = 0.35),

    axis.text = element_text(size = 9),
    axis.title = element_text(size = 11),

    strip.text = element_text(size = 11),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

theme_set(theme_classic())

# Funciones para resumir posteriores
ci <- function(x, width = 0.95, name = "mu") {
  # definimos probabilidad acumulada según la amplitud del intervalo (ci)
  low <- (1 - width) / 2
  high <- width + low

  # calculamos cuantiles asociados
  qq <- quantile(x, probs = c(low, high), method = 8)

  # damos nombres y ordenamos
  nn <- paste(name, c("_lwr", "_upr"), sep = "")
  names(qq) <- nn

  return(qq)
}

mean_ci <- function(x, width = 0.95, name = "mu") {
  # media
  mm <- mean(x)
  names(mm) <- paste(name, "_mean", sep = "")

  # cuantiles
  qq <- ci(x, width, name)

  return(c(mm, qq))
}


# Cargamos y preparamos datos ---------------------------------------------

# Registros de caudal por cuenca y fecha
dcaudal <- read.csv(file.path("datos", "cingolani_data_caudales_fechas.csv"))

# Variables predictoras por cuenca
dvar <- read.csv(file.path("datos", "cingolani_data_caudales_variables.csv"))

# Extraemos las columnas con los registros por fecha (fechas en columnas)
date_names_0 <- names(dcaudal)[grep("F", names(dcaudal))]
date_names <- date_names_0[grep("mm", date_names_0, invert = T)]
dsub <- dcaudal[, c("cca", date_names)]

# número de semanas (t), cuencas (c) y observaciones (o)
nt <- length(date_names)
nc <- nrow(dsub)
nobs <- nt * nc

# Alargamos
dlong <- pivot_longer(
  dsub, all_of(date_names), names_to = "fecha_char", values_to = "caudal"
)

# Extraemos el día como numérico, arrancando en cero
dlong$dia <- sapply(1:nrow(dlong), function(i) {
  strsplit(dlong$fecha_char[i], split = "F")[[1]][2] |> as.numeric()
})

# Calculamos semana, para facilitar interpretación
dlong$semana <- dlong$dia / 7

# Cuenca como factor, servirá para graficar
dlong$cca_fac <- factor(as.character(dlong$cca),
                        levels = as.character(unique(dsub$cca)))

# Extraemos rocosidad total
dlong$roca <- rep(dvar$RT, each = nt)

# Gráfico exploratorio ----------------------------------------------------

ggplot(dlong, aes(semana, caudal, color = roca, group = cca_fac)) +
  # geom_line() +
  geom_smooth(se = F, linewidth = 0.3) +
  geom_point() +
  scale_color_viridis(option = "A") +
  nice_theme()

# Exploramos función de decaimiento de Weibul -----------------------------

a <- 1   # máximo
l <- 3   # cuánto tarda en bajar
k <- 1   # forma (de exponencial a sigmoidea). limitamos a >1
curve(a * exp(-(x / l) ^ k), to = 10, ylim = c(0, 1))

# Simulamos k
cc <- rgb(0, 0, 0, 0.1)

a <- 1      # máximo
l <- 2      # cuánto tarda en bajar
k <- 1   # forma (de exponencial a sigmoidea)
curve(a * exp(-(x / l) ^ k), to = 10, col = cc, ylim = c(0, a))
for (i in 1:100) {
  k <- rtruncnorm(1, a = 1, mean = 1.5, sd = 3)
  curve(a * exp(-(x / l) ^ k), add = T, col = cc)
}

# Simulamos k y l
a <- 1      # máximo
l <- 2      # cuánto tarda en bajar
k <- 1   # forma (de exponencial a sigmoidea)
curve(a * exp(-(x / l) ^ k), to = 10, col = cc, ylim = c(0, a))
for (i in 1:100) {
  k <- rtruncnorm(1, a = 1, mean = 1.5, sd = 2)
  l <- rtruncnorm(1, a = 0, mean = 2.5, sd = 2)
  curve(a * exp(-(x / l) ^ k), add = T, col = cc)
}

# a es el valor en la primera fecha.
mean(dsub$F0)
range(dsub$F0)

# a puede estar centrada en
mu_a <- mean(log(dsub$F0))
sd_a <- 2
curve(dlnorm(x, mu_a, sd_a), to = 80)

# k >= 1, con 10 siendo un valor muy extremo
mu_k <- log(2)
sd_k <- 1.5
curve(dlnorm(1+x, mu_k, sd_k), to = 10)

# l debería estar medio lejos de cero, con 6 siendo un valor alto
mu_l <- log(2.5)
sd_l <- 1.5
curve(dlnorm(x, mu_l, sd_l), to = 10)

# A las medias les pongo estas previas,
# y los log-sds, les pongo normal truncada con sd_param * 2


# Previa para phi
mu <- 15
log_phi <- rnorm(1, log(0.3), 1)
phi <- exp(log_phi)
curve(dgamma(x, shape = 1 / phi, rate = 1 / (mu * phi)), to = 80, col = cc,
      ylim = c(0, 0.15))
for (i in 1:100) {
  log_phi <- rnorm(1, log(0.3), 1.5)
  phi <- exp(log_phi)

  phi <- rtruncnorm(1, 0, mean = 0.5, sd = 2)
  curve(dgamma(x, shape = 1 / phi, rate = 1 / (mu * phi)), col = cc, add = T)
}


# Datos para Stan ---------------------------------------------------------

# Reemplazamos ceros por algo pequeño
dlong$caudal2 <- dlong$caudal
dlong$caudal2[dlong$caudal == 0] <- 0.0001

stan_data <- list(
  nc = nc, nobs = nobs,
  semana = dlong$semana,
  caudal = dlong$caudal2,
  cuenca = dlong$cca
)

# Modelo 01: no pooling = efectos fijos ---------------------------------

# Compilamos
model1 <- cmdstan_model(file.path("modelos", "cuencas_01.stan"))

opt1 <- model1$optimize(data = stan_data)

# Uniform initialization function.
# Tenía problemas para arrancar con el default (uniforme(-2, 2))
lim <- 0.1
init_fun1 <- function() {
  list(
    log_a   = runif(nc, -lim, lim),
    log_l   = runif(nc, -lim, lim),
    log_k   = runif(nc, -lim, lim),
    log_phi = runif(nc, -lim, lim)
  )
}

# iteraciones y cadenas
nchain <- 4
nwarmup <- 1000
nsampling <- 1000

# Muestreamos
fit1 <- model1$sample(
  data = stan_data,
  iter_warmup = nwarmup,
  iter_sampling = nsampling,
  chains = nchain,
  parallel_chains = 4,
  init = init_fun1,
  adapt_delta = 0.95
)

fit1$cmdstan_diagnose()
fit1$diagnostic_summary()
# pocas transiciones divergentes, estamos ok.
# podemos eliminarlas con más warmup o con mayor adapt_delta.

## Graficamos todas las curvas
semana_seq <- seq(0, 6, by = 0.1)
nrep <- length(semana_seq)

# data.frame con variables para calcular media y simular datos.
pred <- data.frame(
  semana = rep(semana_seq, nc),
  cuenca = rep(1:nc, each = nrep)
)

# Matrices con media predicha y datos simulados, para obtener resumen de la media
# y resumen de la distribución predictiva posterior.
np <- nrow(pred)
S <- nsampling * nchain
mumat <- ymat2 <- matrix(NA, np, S)

# extraemos muestras, en una matriz separada por parámetro. los d van por "draws"
# estas matrices de muestras tienen las iteraciones en filas, y las cuencas en
# columnas
da <- fit1$draws("a", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dl <- fit1$draws("l", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dk <- fit1$draws("k", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dphi <- fit1$draws("phi", format = "draws_array") |> as_draws_matrix() |> as.matrix()
# Y phi es un vector
# dphi <- fit1$draws("phi", format = "draws_array") |> as.numeric()

# Llenamos las matrices
for (i in 1:S) {
  print(i) # para ansiosos

  # calculamos mu y simulamos loopeando sobre cuencas
  for (cc in 1:nc) {
    # filtramos filas de la cuenca "cc"
    rr <- pred$cuenca == cc
    mumat[rr, i] <- da[i, cc] * exp(-(pred$semana[rr] / dl[i, cc]) ^ dk[i, cc])
    ymat2[rr, i] <- rgamma(
      nrep, 1 / dphi[i, cc], 1 / (mumat[rr, i] * dphi[i, cc])
    )
  }
}

# Resumimos. Utilizaremos las funciones mean_ci y ci definidas a arriba, para
# aplicarlas fila por fila en las matrices mumat y ymat2, con apply.
mu_summ <- apply(mumat, 1, mean_ci) |> t() |> as.data.frame()
y_summ <- apply(ymat2, 1, ci, name = "y") |> t() |> as.data.frame()

# Unimos resúmenes con los datos usados (semana y cuenca)
pred2 <- cbind(pred, mu_summ, y_summ)

# nombramos "cuenca" a la columna en los datos
dlong$cuenca <- dlong$cca

# Graficamos
ggplot(pred2, aes(semana, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # intervalo de predicción
  geom_ribbon(
    aes(semana, ymin = y_lwr, ymax = y_upr),
    fill = "blue", color = NA, alpha = 0.3
  ) +
  # intervalo de credibilidad para la media
  geom_ribbon(
    fill = "orange", color = NA, alpha = 0.6
  ) +
  # media de la media
  geom_line() +
  # datos
  geom_point(
    data = dlong, mapping = aes(x = semana, y = caudal), inherit.aes = F,
    alpha = 0.7
  ) +
  # cosmética
  facet_wrap(vars(cuenca), axes = "all", scales = "free_y") +
  nice_theme() +
  ylab("Caudal") +
  xlab("Semana")

# Sobreestima la variabilidad. Veamos previa vs. posterior de phi
plot(density(dphi), xlab = expression(phi), main = NA)
curve(dlnorm(x, log(0.3), 1.5), col = 4, add = T)


# Modelo 02: jerárquico = partial pooling --------------------------------

# Compilamos
model2 <- cmdstan_model(file.path("modelos", "cuencas_02.stan"))

# Uniform initialization function.
# Tenía problemas para arrancar con el default (uniforme(-2, 2))
lim <- 0.1
init_fun2 <- function() {
  list(
    # parámetros estandarizados a nivel de cuenca
    log_a_raw   = runif(nc, -lim, lim),
    log_l_raw   = runif(nc, -lim, lim),
    log_k_raw   = runif(nc, -lim, lim),
    log_phi_raw = runif(nc, -lim, lim),

    # hiperparámetros
    log_a_mu = runif(1, -lim, lim),
    log_l_mu = runif(1, -lim, lim),
    log_k_mu = runif(1, -lim, lim),
    log_phi_mu = runif(1, -lim, lim),

    log_a_sd = runif(1, 0.1, 0.2),
    log_l_sd = runif(1, 0.1, 0.2),
    log_k_sd = runif(1, 0.1, 0.2),
    log_phi_sd = runif(1, 0.1, 0.2)
  )
}

# iteraciones y cadenas
nchain <- 4
nwarmup <- 1000
nsampling <- 1000

# Muestreamos (ojo, toma un rato)
fit2 <- model2$sample(
  data = stan_data,
  iter_warmup = nwarmup,
  iter_sampling = nsampling,
  chains = nchain,
  parallel_chains = 4,
  init = init_fun2,
  adapt_delta = 0.95,
  max_treedepth = 15
)

fit2$cmdstan_diagnose()
fit2$diagnostic_summary()
# pocas transiciones divergentes, estamos ok.
# podemos eliminarlas con más warmup o con mayor adapt_delta.

## Graficamos todas las curvas
semana_seq <- seq(0, 6, by = 0.1)
nrep <- length(semana_seq) # largo de secuencia temporal para graficar (de 0 a 7 semanas)

# data.frame con variables para calcular media y simular datos.
pred <- data.frame(
  semana = rep(semana_seq, nc),
  cuenca = rep(1:nc, each = nrep)
)

# Matrices con media predicha y datos simulados, para obtener resumen de la media
# y resumen de la distribución predictiva posterior.
np <- nrow(pred)
S <- nsampling * nchain
mumat2 <- ymat2 <- matrix(NA, np, S)

# extraemos muestras, en una matriz separada por parámetro. los d van por "draws"
# estas matrices de muestras tienen las iteraciones en filas, y las cuencas en
# columnas
da <- fit2$draws("a", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dl <- fit2$draws("l", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dk <- fit2$draws("k", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dphi <- fit2$draws("phi", format = "draws_array") |> as_draws_matrix() |> as.matrix()
# Y phi es un vector
# dphi <- fit1$draws("phi", format = "draws_array") |> as.numeric()

# Llenamos las matrices
for (i in 1:S) {
  print(i) # para ansiosos

  # calculamos mu y simulamos loopeando sobre cuencas
  for (cc in 1:nc) {
    # filtramos filas de la cuenca "cc"
    rr <- pred$cuenca == cc
    mumat2[rr, i] <- da[i, cc] * exp(-(pred$semana[rr] / dl[i, cc]) ^ dk[i, cc])
    ymat2[rr, i] <- rgamma(
      nrep, 1 / dphi[i, cc], 1 / (mumat2[rr, i] * dphi[i, cc])
    )
  }
}

# Resumimos. Utilizaremos las funciones mean_ci y ci definidas a arriba, para
# aplicarlas fila por fila en las matrices mumat2 y ymat2, con apply.
mu_summ <- apply(mumat2, 1, mean_ci) |> t() |> as.data.frame()
y_summ <- apply(ymat2, 1, ci, name = "y") |> t() |> as.data.frame()

# Unimos resúmenes con los datos usados (semana y cuenca)
pred2 <- cbind(pred, mu_summ, y_summ)

# nombramos "cuenca" a la columna en los datos
dlong$cuenca <- dlong$cca

# Graficamos
ggplot(pred2, aes(semana, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # intervalo de predicción
  geom_ribbon(
    aes(semana, ymin = y_lwr, ymax = y_upr),
    fill = "blue", color = NA, alpha = 0.3
  ) +
  # intervalo de credibilidad para la media
  geom_ribbon(
    fill = "orange", color = NA, alpha = 0.6
  ) +
  # media de la media
  geom_line() +
  # datos
  geom_point(
    data = dlong, mapping = aes(x = semana, y = caudal), inherit.aes = F,
    alpha = 0.7
  ) +
  # cosmética
  facet_wrap(vars(cuenca), axes = "all", scales = "free_y") +
  nice_theme() +
  ylab("Caudal") +
  xlab("Semana")

# Modelo 03: jerárquico con predictora a nivel de cuenca ------------------
# la media de los parámetros es una función lineal de una predictora.

# Datos incluyendo predictora a nivel de cuenca (rocosidad total)
roca_mean <- mean(dvar$RT)
roca_sd <- sd(dvar$RT)
rocaz <- (dvar$RT - roca_mean) / roca_sd

stan_data3 <- list(
  nc = nc, nobs = nobs,
  semana = dlong$semana,
  caudal = dlong$caudal2,
  cuenca = dlong$cca,
  rocaz = rocaz
)

# Compilamos
model3 <- cmdstan_model(file.path("modelos", "cuencas_03.stan"))

# Uniform initialization function.
# Tenía problemas para arrancar con el default (uniforme(-2, 2))
lim <- 0.1
init_fun3 <- function() {
  list(
    # parámetros estandarizados a nivel de cuenca
    log_a_raw   = runif(nc, -lim, lim),
    log_l_raw   = runif(nc, -lim, lim),
    log_k_raw   = runif(nc, -lim, lim),
    log_phi_raw = runif(nc, -lim, lim),

    # hiperparámetros
    # interceptos
    log_a_mu = runif(1, -lim, lim),
    log_l_mu = runif(1, -lim, lim),
    log_k_mu = runif(1, -lim, lim),
    log_phi_mu = runif(1, -lim, lim),

    # pendientes (linear y cuadrática)
    log_a_b1 = runif(1, -lim, lim),
    log_l_b1 = runif(1, -lim, lim),
    log_k_b1 = runif(1, -lim, lim),
    log_phi_b1 = runif(1, -lim, lim),

    log_a_b2 = runif(1, -lim, lim),
    log_l_b2 = runif(1, -lim, lim),
    log_k_b2 = runif(1, -lim, lim),
    log_phi_b2 = runif(1, -lim, lim),

    # desvíos
    log_a_sd = runif(1, 0.1, 0.2),
    log_l_sd = runif(1, 0.1, 0.2),
    log_k_sd = runif(1, 0.1, 0.2),
    log_phi_sd = runif(1, 0.1, 0.2)
  )
}

# iteraciones y cadenas
nchain <- 4
nwarmup <- 1000
nsampling <- 1000

# Muestreamos (ojo, toma un rato)
fit3 <- model3$sample(
  data = stan_data3,
  iter_warmup = nwarmup,
  iter_sampling = nsampling,
  chains = nchain,
  parallel_chains = 4,
  init = init_fun3,
  adapt_delta = 0.95,
  max_treedepth = 15
)

fit3$cmdstan_diagnose()
fit3$diagnostic_summary()

## Graficamos todas las curvas
semana_seq <- seq(0, 6, by = 0.1)
nrep <- length(semana_seq) # largo de secuencia temporal para graficar (de 0 a 7 semanas)

# data.frame con variables para calcular media y simular datos.
pred <- data.frame(
  semana = rep(semana_seq, nc),
  cuenca = rep(1:nc, each = nrep)
)

# Matrices con media predicha y datos simulados, para obtener resumen de la media
# y resumen de la distribución predictiva posterior.
np <- nrow(pred)
S <- nsampling * nchain
mumat3 <- ymat3 <- matrix(NA, np, S)

# extraemos muestras, en una matriz separada por parámetro. los d van por "draws"
# estas matrices de muestras tienen las iteraciones en filas, y las cuencas en
# columnas
da <- fit3$draws("a", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dl <- fit3$draws("l", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dk <- fit3$draws("k", format = "draws_array") |> as_draws_matrix() |> as.matrix()
dphi <- fit3$draws("phi", format = "draws_array") |> as_draws_matrix() |> as.matrix()
# Y phi es un vector
# dphi <- fit1$draws("phi", format = "draws_array") |> as.numeric()

# Llenamos las matrices
for (i in 1:S) {
  print(i) # para ansiosos

  # calculamos mu y simulamos loopeando sobre cuencas
  for (cc in 1:nc) {
    # filtramos filas de la cuenca "cc"
    rr <- pred$cuenca == cc
    mumat3[rr, i] <- da[i, cc] * exp(-(pred$semana[rr] / dl[i, cc]) ^ dk[i, cc])
    ymat3[rr, i] <- rgamma(
      nrep, 1 / dphi[i, cc], 1 / (mumat3[rr, i] * dphi[i, cc])
    )
  }
}

# Resumimos. Utilizaremos las funciones mean_ci y ci definidas a arriba, para
# aplicarlas fila por fila en las matrices mumat3 y ymat3, con apply.
mu_summ <- apply(mumat3, 1, mean_ci) |> t() |> as.data.frame()
y_summ <- apply(ymat3, 1, ci, name = "y") |> t() |> as.data.frame()

# Unimos resúmenes con los datos usados (semana y cuenca)
pred2 <- cbind(pred, mu_summ, y_summ)

# nombramos "cuenca" a la columna en los datos
dlong$cuenca <- dlong$cca

# Graficamos
ggplot(pred2, aes(semana, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # intervalo de predicción
  geom_ribbon(
    aes(semana, ymin = y_lwr, ymax = y_upr),
    fill = "blue", color = NA, alpha = 0.3
  ) +
  # intervalo de credibilidad para la media
  geom_ribbon(
    fill = "orange", color = NA, alpha = 0.6
  ) +
  # media de la media
  geom_line() +
  # datos
  geom_point(
    data = dlong, mapping = aes(x = semana, y = caudal), inherit.aes = F,
    alpha = 0.7
  ) +
  # cosmética
  facet_wrap(vars(cuenca), axes = "all", scales = "free_y") +
  nice_theme() +
  ylab("Caudal") +
  xlab("Semana")


# Modelo 03: predicciones heterodoxas -------------------------------------

# Derivada promedio en las 7 fechas en función de la rocosidad.
# Promedio de las 7 fechas en función de la rocosidad.
# Idealmente, calcular también para lo observado.

semanas <- 0:6
dx <- 1e-3
semanas_fore <- semanas + dx
# para calcular derivadas. Sólo hacia adelante, porque la función no soporta
# tiempo negativo.

# Generamos nrr valores de rocosidad, de mínimo a máximo, sobre las cuales
# calcularemos la derivada media y la media media.
nrr <- 50
rocaz_seq <- seq(min(rocaz), max(rocaz), length.out = nrr)
rocaz_seq2 <- rocaz_seq ^ 2

# matrices de nrr * S para guardar la media y la derivada
meanmat <- derivmat <- matrix(NA, nrr, S)

# nro de muestras de efectos aleatorios para tomar por cada iteración de la
# posterior
nran <- 100

# Extraemos muestras de la posterior

# Hiperparámetros, para simular parámetros
fixmu <- fit3$draws(
  c("log_a_mu", "log_l_mu", "log_k_mu", "log_phi_mu"),
  format = "draws_df"
)

fixb1 <- fit3$draws(
  c("log_a_b1", "log_l_b1", "log_k_b1", "log_phi_b1"),
  format = "draws_df"
)

fixb2 <- fit3$draws(
  c("log_a_b2", "log_l_b2", "log_k_b2", "log_phi_b2"),
  format = "draws_df"
)

fixsd <- fit3$draws(
  c("log_a_sd", "log_l_sd", "log_k_sd", "log_phi_sd"),
  format = "draws_df"
)

# Matrices temporal para guardar las curvas (en 0:6 y en su diff hacia adelante)
mtemp_focal <- mtemp_fore <- matrix(NA, 7, nran)

for (i in 1:S) {
  print(i)
  # Loopoeamos sobre los valores de la predictora rocaz
  for (j in 1:nrr) {
    # Simulamos parámetros para nran cuencas nuevas
    log_a_mean <-
      fixmu$log_a_mu[i] +
      fixb1$log_a_b1[i] * rocaz_seq[j] +
      fixb2$log_a_b2[i] * rocaz_seq2[j]
    a <- exp(rnorm(nran, log_a_mean, fixsd$log_a_sd[i]))

    log_l_mean <-
      fixmu$log_l_mu[i] +
      fixb1$log_l_b1[i] * rocaz_seq[j] +
      fixb2$log_l_b2[i] * rocaz_seq2[j]
    l <- exp(rnorm(nran, log_l_mean, fixsd$log_l_sd[i]))

    log_k_mean <-
      fixmu$log_k_mu[i] +
      fixb1$log_k_b1[i] * rocaz_seq[j] +
      fixb2$log_k_b2[i] * rocaz_seq2[j]
    k <- exp(rnorm(nran, log_k_mean, fixsd$log_k_sd[i])) + 1

    # evaluamos nran curvas
    for (s in 1:nran) {
      mtemp_focal[, s] <- a[s] * exp(-(semanas / l[s]) ^ k[s])
      mtemp_fore[, s] <- a[s] * exp(-(semanas_fore / l[s]) ^ k[s])
    }

    # las resumimos en media y derivada promedio
    # media:
    meanmat[j, i] <- mean(mtemp_focal)

    # derivada:
    dy_dx <- (mtemp_fore - mtemp_focal) / dx
    derivmat[j, i] <- mean(dy_dx)
  }
}

# derivada relativizada a la media
derivmat_rel <- derivmat / meanmat

# Resumimos
pred_summ0 <- data.frame(
  rocosidad = rocaz_seq * roca_sd + roca_mean # desestandarizamos para graficar
)

msumm <- apply(meanmat, 1, mean_ci, name = "mu") |> t() |> as.data.frame()
dsumm <- apply(derivmat, 1, mean_ci, name = "d") |> t() |> as.data.frame()
dsumm_rel <- apply(derivmat_rel, 1, mean_ci, name = "dr") |> t() |> as.data.frame()

pred_summ <- cbind(pred_summ0, msumm, dsumm, dsumm_rel)

## Antes de graficar, calculamos las mismas métricas estimadas para las cuencas
## observadas.
meanmat_obs <- derivmat_obs <- matrix(NA, nc, S)

# Llenamos las matrices
for (i in 1:S) {
  print(i) # para ansiosos
  # calculamos mu y simulamos loopeando sobre cuencas
  for (cc in 1:nc) {
    mu_focal <- da[i, cc] * exp(-(semanas / dl[i, cc]) ^ dk[i, cc])
    mu_fore <- da[i, cc] * exp(-(semanas_fore / dl[i, cc]) ^ dk[i, cc])

    meanmat_obs[cc, i] <- mean(mu_focal)
    derivmat_obs[cc, i] <- mean((mu_fore - mu_focal) / dx)
  }
}

# derivada relativizada a la media
derivmat_rel_obs <- derivmat_obs / meanmat_obs

# Resumimos
pred_obs0 <- data.frame(
  rocosidad = rocaz * roca_sd + roca_mean # desestandarizamos para graficar
)

msumm_obs <- apply(meanmat_obs, 1, mean_ci, name = "mu") |> t() |> as.data.frame()
dsumm_obs <- apply(derivmat_obs, 1, mean_ci, name = "d") |> t() |> as.data.frame()
dsumm_rel_obs <- apply(derivmat_rel_obs, 1, mean_ci, name = "dr") |> t() |> as.data.frame()

pred_obs <- cbind(pred_obs0, msumm_obs, dsumm_obs, dsumm_rel_obs)



ggplot(pred_summ, aes(rocosidad, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  geom_ribbon(color = NA, fill = "orange", alpha = 0.6) +
  geom_line() +
  geom_point(data = pred_obs, size = 2, alpha = 0.8) +
  geom_linerange(data = pred_obs, alpha = 0.7) +
  labs(x = "Rocosidad (%)", y = "Caudal medio")


ggplot(pred_summ, aes(rocosidad, d_mean, ymin = d_lwr, ymax = d_upr)) +
  geom_ribbon(color = NA, fill = "orange", alpha = 0.6) +
  geom_line() +
  geom_point(data = pred_obs, size = 2, alpha = 0.8) +
  geom_linerange(data = pred_obs, alpha = 0.7) +
  labs(x = "Rocosidad (%)", y = "d Caudal / d t")

ggplot(pred_summ, aes(rocosidad, dr_mean, ymin = dr_lwr, ymax = dr_upr)) +
  geom_ribbon(color = NA, fill = "orange", alpha = 0.6) +
  geom_line() +
  geom_point(data = pred_obs, size = 2, alpha = 0.8) +
  geom_linerange(data = pred_obs, alpha = 0.7) +
  labs(x = "Rocosidad (%)", y = "Derivada / media")

# El promedio predicho, en particular para el caudal medio, sobreestima las
# estimaciones de las cuencas reales. Esto puede deberse a que los parámetros
# de la curva tienen cierta correlación, pero nosotros asumimos que son
# independientes. En general, ignorar las correlaciones infla la incertidumbre.