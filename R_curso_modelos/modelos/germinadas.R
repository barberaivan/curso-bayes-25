rm(list = ls()) # Limpiamos el entorno

# Eliminamos modelos compilados de stan, para evitar incompatibilidades
file.remove(list.files("modelos", pattern = "^[^.]+$", full.names = TRUE))

# Paquetes ----------------------------------------------------------------

library(tidyverse)  # gráficos et al
library(cmdstanr)   # ajustar modelos con Stan
library(posterior)  # resumir posteriores y diagnósticos de MCMC
library(bayesplot)  # visualizar posteriores
library(DHARMa)     # probabilidad acumulada
library(extraDistr) # rbbinom, la beta-binomial

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

# Datos -------------------------------------------------------------------

datos <- read.csv(file.path("datos", "alinari_data_espinillo.csv"))

# Conservamos los que tienen datos de germinación
datos <- subset(datos, !is.na(germinadas) & sembradas > 0)

# Estandarizamos las predictoras para tener efectos comparables y para
# facilitar la interpretación de las previas
D_mean <- mean(datos$danio)
D_sd <- sd(datos$danio)

E_mean <- mean(datos$altitud)
E_sd <- sd(datos$altitud)

H_mean <- mean(datos$altura)
H_sd <- sd(datos$altura)

datos$D <- (datos$danio - D_mean) / D_sd
datos$E <- (datos$altitud - E_mean) / E_sd
datos$H <- (datos$altura - H_mean) / H_sd

# Datos para Stan
stan_data <- list(
  N = nrow(datos),
  y = datos$germinadas,
  S = 120,

  D = datos$D, # En modelo 1 sólo se usa esta predictora
  E = datos$E,
  H = datos$H
)

# Un enfoque más familiar -------------------------------------------------

# Ajustamos el mismo modelo pero con enfoque frecuentista
fit_freq <- glm(
  cbind(germinadas, sembradas - germinadas) ~ D,
  family = "binomial", data = datos
)
summary(fit_freq)

# Modelo 01: binomial, una predictoras ------------------------------------

# Ajuste del modelo -------------------------------------------------------

# Compilamos
model1 <- cmdstan_model(here::here("modelos", "germinadas_01.stan"))

# Muestreamos
fit1 <- model1$sample(data = stan_data, parallel_chains = 4)

# Chequeamos MCMC
fit1$cmdstan_diagnose()

# Resultados del ajuste ---------------------------------------------------

summ <- fit1$summary(c("alpha", "beta"))
print(summ)

# Previa vs. posterior ----------------------------------------------------

# Extraemos muestras de los parámetros que necesitamos
d <- fit1$draws(c("alpha", "beta"), format = "draws_df")
colnames(d) <- c("a", "b")
S <- nrow(d)
N <- nrow(datos)

# alpha (intercept)
plot(density(d$a), main = NA, xlab = expression(alpha), xlim = c(-2, 2))
curve(dnorm(x, 0, 10), add = T, col = 4) # previa definida en Stan

# beta (efecto del daño)
plot(density(d$b), main = NA, xlab = expression(beta), xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# Probabilidad de que el daño posfuego sobre la planta madre disminuya la
# capacidad germinativa de sus semillas
sum(d$b < 0) / S # total certeza

# Verificación predictiva posterior ---------------------------------------

# Verificamos el ajuste simulando de la predictiva posterior. Para simular de
# la beta-binomial en R usaremos el paquete extraDistr, que incluye la función
# rbbinom
ysim <- matrix(NA, N, S)
for (i in 1:S) {
  mu <- plogis(d$a[i] + d$b[i] * datos$D)
  # plogis equivale a inv_logit = logit inverso

  ysim[, i] <- rbinom(N, size = datos$sembradas, prob = mu)
}

ysim_t <- t(ysim)

# Miramos con bayesplot
ppc_intervals(
  datos$germinadas, yrep = ysim_t, x = datos$danio,
  prob = 0.5, prob_outer = 0.9
)

use <- sample(1:S, 100)
ppc_dens_overlay(datos$germinadas, yrep = ysim_t[use, ])

# DHARMa
res <- createDHARMa(ysim, datos$germinadas, integerResponse = T)
plot(res)
plotResiduals(res, form = datos$danio, rank = F)

# R2 ----------------------------------------------------------------------

# Vector para almacenar el R2
R2 <- numeric(S)

# Por comodidad, definimos las sembradas
SS <- datos$sembradas[1]

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # En la beta-binomial, la media es mu * size,
  # siendo size la cantidad de sembradas

  # predichos para cada observación
  mu <- plogis(d$a[i] + d$b[i] * datos$D)
  # plogis equivale a inv_logit = logit inverso

  media <- mu * SS # 120 sembradas

  # varianza de los predichos = varianza explicada
  var_fit <- var(media)

  # varianza residual para cada observación. En la binomial, la
  # varianza se define así:
  var_res <- mu * (1 - mu) * SS # N * p * (1 - p)

  # y tomamos su media entre observaciones
  var_res_mean <- mean(var_res)

  # R2
  R2[i] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")
# muy alto

# Predicciones ------------------------------------------------------------

# Predicción parcial para visualizar el efecto del daño por fuego (D).

# Graficaremos curvas basadas en nrep valores de daño en el rango observado
nrep <- 200
danio_seq <- seq(min(datos$danio), max(datos$danio), length.out = nrep)
# Estandarizamos la secuencia, porque los parámetros viven en escala estandarizada.
D_seq <- (danio_seq - D_mean) / D_sd

# Creamos data.frame de predicciones con las predictoras.
pred <- data.frame(D = D_seq)
pred$danio <- danio_seq # la agregamos en la escala original, para visualizar

# Matriz para llenar con predicciones de mu
mu_mat <- matrix(NA, nrep, S)

# Y matriz para llenar con valores simulados de semillas germinadas
y_mat <- matrix(NA, nrep, S) # se parece a ysim, pero tiene nrep filas, no N.

# Las llenamos
for (i in 1:S) {
  mu <- plogis(d$a[i] + d$b[i] * pred$D)
  # equivale a inv_logit = logit inverso

  mu_mat[, i] <- mu

  # simulamos datos
  y_mat[, i] <- rbinom(nrep, size = datos$sembradas[1], prob = mu)
}

# Resumimos la distribución posterior de mu fila por fila
pred$mu_mean <- rowMeans(mu_mat)
pred$mu_lwr <- apply(mu_mat, 1, quantile, probs = 0.025)
pred$mu_upr <- apply(mu_mat, 1, quantile, probs = 0.975)

# Convertimos la distribución de y en proporciones, para que esté en la misma
# escala que mu
y_mat_p <- y_mat / datos$sembradas[1]

# y resumimos
pred$y_lwr <- apply(y_mat_p, 1, quantile, probs = 0.025)
pred$y_upr <- apply(y_mat_p, 1, quantile, probs = 0.975)

# Calculamos la proporción observada en los datos
datos$prop <- datos$germinadas / datos$sembradas

# Graficamos
ggplot(pred, aes(danio, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # Intervalo de predicción (posterior predictive)
  geom_ribbon(aes(danio, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
              color = NA, alpha = 0.3, fill = "blue") +
  # Intervalo para la media
  geom_ribbon(color = NA, alpha = 0.7, fill = "orange") +
  # Media de la media
  geom_line() +
  # datos
  geom_point(aes(danio, prop), data = datos, inherit.aes = F,
             size = 2, alpha = 0.7) +
  labs(
    y = "Proporción de germinación",
    x = "Daño (%)"
  )

# Modelo 02: binomial, tres predictoras -----------------------------------

# Ajuste del modelo -------------------------------------------------------

# Compilamos
model2 <- cmdstan_model(here::here("modelos", "germinadas_02.stan"))

# Muestreamos
fit2 <- model2$sample(data = stan_data, parallel_chains = 4)

# Chequeamos MCMC
fit2$cmdstan_diagnose()


# Verificación predictiva posterior ---------------------------------------

# Extraemos muestras de los parámetros que necesitamos
d <- fit2$draws(
  c("alpha", "beta_D", "beta_E", "beta_H"), format = "draws_df"
)
colnames(d) <- c("a", "bD", "bE", "bH")
S <- nrow(d)
N <- nrow(datos)

# Verificamos el ajuste simulando de la predictiva posterior. Para simular de
# la beta-binomial en R usaremos el paquete extraDistr, que incluye la función
# rbbinom
ysim <- matrix(NA, N, S)
for (i in 1:S) {
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
    d$bD[i] * datos$D + d$bE[i] * datos$E + d$bH[i] * datos$H
  )
  ysim[, i] <- rbinom(N, size = datos$sembradas, prob = mu)
}

ysim_t <- t(ysim)

# Miramos con bayesplot
ppc_intervals(
  datos$germinadas, yrep = ysim_t, x = datos$danio,
  prob = 0.5, prob_outer = 0.9
)

use <- sample(1:S, 100)
ppc_dens_overlay(datos$germinadas, yrep = ysim_t[use, ])

# DHARMa
res <- createDHARMa(ysim, datos$germinadas, integerResponse = T)
plot(res)
plotResiduals(res, form = datos$danio, rank = F)
plotResiduals(res, form = datos$altitud, rank = F)
plotResiduals(res, form = datos$altura, rank = F)


# Previa vs. posterior ----------------------------------------------------

# beta_D (daño)
plot(density(d$bD), main = NA, xlab = "beta_D", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# beta_E (altitud)
plot(density(d$bE), main = NA, xlab = "beta_E", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# beta_H (altura)
plot(density(d$bH), main = NA, xlab = "beta_H", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan


# R2 ----------------------------------------------------------------------

# Vector para almacenar el R2
R2 <- numeric(S)

# Por comodidad, definimos las sembradas
SS <- datos$sembradas[1]

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # En la beta-binomial, la media es mu * size,
  # siendo size la cantidad de sembradas

  # predichos para cada observación
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
      d$bD[i] * datos$D + d$bE[i] * datos$E + d$bH[i] * datos$H
  )
  media <- mu * SS # 120 sembradas

  # varianza de los predichos = varianza explicada
  var_fit <- var(media)

  # varianza residual para cada observación. En la binomial, la
  # varianza se define así:
  var_res <- mu * (1 - mu) * SS # N * p * (1 - p)

  # y tomamos su media entre observaciones
  var_res_mean <- mean(var_res)

  # R2
  R2[i] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")
# muy alto

# Predicciones ------------------------------------------------------------

# Predicción parcial para visualizar el efecto del daño por fuego (D).

# Graficaremos curvas basadas en nrep valores de daño en el rango observado
nrep <- 200
danio_seq <- seq(min(datos$danio), max(datos$danio), length.out = nrep)
# Estandarizamos la secuencia, porque los parámetros viven en escala estandarizada.
D_seq <- (danio_seq - D_mean) / D_sd

# Creamos data.frame de predicciones con las predictoras. Las predictoras no
# focales se fijan en su media (0 porque están estandarizadas), por eso es una
# predicción parcial.
pred <- data.frame(
  D = D_seq,
  E = 0,
  H = 0
)
pred$danio <- danio_seq # la agregamos en la escala original, para visualizar

# Matriz para llenar con predicciones de mu
mu_mat <- matrix(NA, nrep, S)

# Y matriz para llenar con valores simulados de semillas germinadas
y_mat <- matrix(NA, nrep, S) # se parece a ysim, pero tiene nrep filas, no N.

# Las llenamos
for (i in 1:S) {
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
      d$bD[i] * pred$D + d$bE[i] * pred$E + d$bH[i] * pred$H
    # los dos últimos términos pueden obviarse, pero los dejamos para
    # explicitar.
  )
  mu_mat[, i] <- mu

  # simulamos datos
  y_mat[, i] <- rbinom(nrep, size = datos$sembradas[1], prob = mu)
}

# Resumimos la distribución posterior de mu fila por fila
pred$mu_mean <- rowMeans(mu_mat)
pred$mu_lwr <- apply(mu_mat, 1, quantile, probs = 0.025)
pred$mu_upr <- apply(mu_mat, 1, quantile, probs = 0.975)

# Convertimos la distribución de y en proporciones, para que esté en la misma
# escala que mu
y_mat_p <- y_mat / datos$sembradas[1]

# y resumimos
pred$y_lwr <- apply(y_mat_p, 1, quantile, probs = 0.025)
pred$y_upr <- apply(y_mat_p, 1, quantile, probs = 0.975)

# Calculamos la proporción observada en los datos
datos$prop <- datos$germinadas / datos$sembradas

# Graficamos
ggplot(pred, aes(danio, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # Intervalo de predicción (posterior predictive)
  geom_ribbon(aes(danio, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
              color = NA, alpha = 0.3, fill = "blue") +
  # Intervalo para la media
  geom_ribbon(color = NA, alpha = 0.7, fill = "orange") +
  # Media de la media
  geom_line() +
  # datos
  geom_point(aes(danio, prop), data = datos, inherit.aes = F,
             size = 2, alpha = 0.7) +
  labs(
    y = "Proporción de germinación",
    x = "Daño (%)"
  )


# Modelo 03: beta-binomial, tres predictoras ------------------------------

# Ajuste del modelo -------------------------------------------------------

# Compilamos
model3 <- cmdstan_model(here::here("modelos", "germinadas_03.stan"))

# Muestreamos
fit3 <- model3$sample(data = stan_data, parallel_chains = 4)

# Chequeamos MCMC
fit3$cmdstan_diagnose()


# Verificación predictiva posterior ---------------------------------------

# Extraemos muestras de los parámetros que necesitamos
d <- fit3$draws(
  c("alpha", "beta_D", "beta_E", "beta_H", "phi"), format = "draws_df"
)
colnames(d) <- c("a", "bD", "bE", "bH", "phi")
S <- nrow(d)
N <- nrow(datos)

# Verificamos el ajuste simulando de la predictiva posterior. Para simular de
# la beta-binomial en R usaremos el paquete extraDistr, que incluye la función
# rbbinom
ysim <- matrix(NA, N, S)
for (i in 1:S) {
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
      d$bD[i] * datos$D + d$bE[i] * datos$E + d$bH[i] * datos$H
  )
  a <- mu * d$phi[i]
  b <- (1 - mu) * d$phi[i]
  ysim[, i] <- rbbinom(N, datos$sembradas, a, b)
}

ysim_t <- t(ysim)

# Miramos con bayesplot
ppc_intervals(
  datos$germinadas, yrep = ysim_t, x = datos$danio,
  prob = 0.5, prob_outer = 0.9
)

use <- sample(1:S, 100)
ppc_dens_overlay(datos$germinadas, yrep = ysim_t[use, ])

# DHARMa
res <- createDHARMa(ysim, datos$germinadas, integerResponse = T)
plot(res)
plotResiduals(res, form = datos$danio, rank = F)
plotResiduals(res, form = datos$altitud, rank = F)
plotResiduals(res, form = datos$altura, rank = F)


# Previa vs. posterior ----------------------------------------------------

# beta_D (daño)
plot(density(d$bD), main = NA, xlab = "beta_D", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# beta_E (altitud)
plot(density(d$bE), main = NA, xlab = "beta_E", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# beta_H (altura)
plot(density(d$bH), main = NA, xlab = "beta_H", xlim = c(-1, 1))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# phi
plot(density(d$phi, from = 0), main = NA, xlab = expression(phi),
     xlim = c(0, 10))
curve(dlnorm(x, log(7), 1.5), add = T, col = 4) # previa definida en Stan


# R2 ----------------------------------------------------------------------

# Vector para almacenar el R2
R2 <- numeric(S)

# Por comodidad, definimos las sembradas
SS <- datos$sembradas[1]

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # En la beta-binomial, la media es mu * size,
  # siendo size la cantidad de sembradas

  # predichos para cada observación
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
      d$bD[i] * datos$D + d$bE[i] * datos$E + d$bH[i] * datos$H
  )
  media <- mu * SS # 120 sembradas

  # varianza de los predichos = varianza explicada
  var_fit <- var(media)

  # varianza residual para cada observación. En la beta-binomial, la
  # varianza se define así:
  var_res <- SS * mu * (1 - mu) * ((SS + d$phi[i]) / (1 + d$phi[i]))

  # y tomamos su media entre observaciones
  var_res_mean <- mean(var_res)

  # R2
  R2[i] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")


# Predicciones ------------------------------------------------------------

# Predicción parcial para visualizar el efecto del daño por fuego (D).

# Graficaremos curvas basadas en nrep valores de daño en el rango observado
nrep <- 200
danio_seq <- seq(min(datos$danio), max(datos$danio), length.out = nrep)
# Estandarizamos la secuencia, porque los parámetros viven en escala estandarizada.
D_seq <- (danio_seq - D_mean) / D_sd

# Creamos data.frame de predicciones con las predictoras. Las predictoras no
# focales se fijan en su media (0 porque están estandarizadas), por eso es una
# predicción parcial.
pred <- data.frame(
  D = D_seq,
  E = 0,
  H = 0
)
pred$danio <- danio_seq # la agregamos en la escala original, para visualizar

# Matriz para llenar con predicciones de mu
mu_mat <- matrix(NA, nrep, S)

# Y matriz para llenar con valores simulados de semillas germinadas
y_mat <- matrix(NA, nrep, S) # se parece a ysim, pero tiene nrep filas, no N.

# Las llenamos
for (i in 1:S) {
  mu <- plogis( # equivale a inv_logit = logit inverso
    d$a[i] +
      d$bD[i] * pred$D + d$bE[i] * pred$E + d$bH[i] * pred$H
    # los dos últimos términos pueden obviarse, pero los dejamos para
    # explicitar.
  )
  mu_mat[, i] <- mu

  # calculamos alpha y beta para simular
  a <- mu * d$phi[i]
  b <- (1 - mu) * d$phi[i]
  y_mat[, i] <- rbbinom(nrep, datos$sembradas[1], a, b)
}

# Resumimos la distribución posterior de mu fila por fila
pred$mu_mean <- rowMeans(mu_mat)
pred$mu_lwr <- apply(mu_mat, 1, quantile, probs = 0.025)
pred$mu_upr <- apply(mu_mat, 1, quantile, probs = 0.975)

# Convertimos la distribución de y en proporciones, para que esté en la misma
# escala que mu
y_mat_p <- y_mat / datos$sembradas[1]

# y resumimos
pred$y_lwr <- apply(y_mat_p, 1, quantile, probs = 0.025)
pred$y_upr <- apply(y_mat_p, 1, quantile, probs = 0.975)

# Calculamos la proporción observada en los datos
datos$prop <- datos$germinadas / datos$sembradas

# Graficamos
ggplot(pred, aes(danio, mu_mean, ymin = mu_lwr, ymax = mu_upr)) +
  # Intervalo de predicción (posterior predictive)
  geom_ribbon(aes(danio, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
              color = NA, alpha = 0.3, fill = "blue") +
  # Intervalo para la media
  geom_ribbon(color = NA, alpha = 0.7, fill = "orange") +
  # Media de la media
  geom_line() +
  # datos
  geom_point(aes(danio, prop), data = datos, inherit.aes = F,
             size = 2, alpha = 0.7) +
  labs(
    y = "Proporción de germinación",
    x = "Daño (%)"
  )