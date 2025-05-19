rm(list = ls()) # Limpiamos el entorno

# Eliminamos modelos compilados de stan, para evitar incompatibilidades
file.remove(list.files("modelos", pattern = "^[^.]+$", full.names = TRUE))

# Paquetes ----------------------------------------------------------------

library(tidyverse)  # gráficos et al
library(cmdstanr)   # ajustar modelos con Stan
library(posterior)  # resumir posteriores y diagnósticos de MCMC
library(bayesplot)  # visualizar posteriores
library(DHARMa)     # probabilidad acumulada

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

datos <- read.csv(file.path("datos", "barbera_data_fire_total_climate.csv"))

# Datos ordenados para Stan:
# Creamos una lista nombrada para pasarle los datos. Los nombres de cada
# elemento de la lista tienen que ser exactamente los que definimos en
# la sección data {}
stan_data <- list(
  N = nrow(datos),
  y = datos$fires,
  x = datos$fwi
)


# Un enfoque más familiar -------------------------------------------------

# Ajustamos el mismo modelo pero con enfoque frecuentista
fit_freq <- glm(fires ~ fwi, family = "poisson", data = datos)
summary(fit_freq)
car::Anova(fit_freq)

# graficamos
plot(fires ~ fwi, datos, pch = 19)
curve(exp(0.0842 + 0.16532 * x), add = TRUE, col = 4)

# Modelo 01: poisson -----------------------------------------------------

# Compilamos
model1 <- cmdstan_model(file.path("modelos", "nfuegos_01.stan"))

# Y muestreamos la posterior
fit1 <- model1$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit1$cmdstan_diagnose()

# Resultados del ajuste ---------------------------------------------------

# Resumen de las marginales (lo que le pidamos en "variables")
summ <- fit1$summary(variables = c("alpha", "beta"))
print(summ)

# Extraemos las muestras para calculas más cosas
d <- fit1$draws(variables = c("alpha", "beta"), format = "draws_df")

# calcular percentiles 2.5 y 97.5 %, el intervalo de credibilidad del 95 %
quantile(d$alpha, probs = c(0.025, 0.975))
quantile(d$beta, probs = c(0.025, 0.975))

# Algunos enunciados probabilísticos que pueden interesar (dependiendo del
# contexto)

S <- nrow(d) # cantidad de muestras de la posterior

# Pr(beta > 0)
sum(d$beta > 0) / S

# Pr(exp(beta) > 1.2)
ebeta <- exp(d$beta)
sum(ebeta > 1.2) / S

# Pr(lambda > 20 | FWI = 18)
lambda <- exp(d$alpha + d$beta * 18)
sum(lambda > 20) / S

# Pr(beta < 0.01 & beta > -0.01)
sum(d$beta < 0.01 & d$beta > -0.01) / S


# Comparación de previa contra posterior.
# Usamos density() para obtener la densidad empírica a partir de muestras,
# y curve() para dibujar la previa.

# alpha
plot(density(d$alpha), main = NA, xlab = expression(alpha), xlim = c(-5, 5))
curve(dnorm(x, 0, 1), col = 4, add = T)

# beta
plot(density(d$beta), main = NA, xlab = expression(beta), xlim = c(-0.4, 0.4))
curve(dnorm(x, 0, 0.1), col = 4, add = T)


# Verificación predictiva posterior  --------------------------------------

# Extraemos las muestras de alpha y beta, en formato data.frame
d <- fit1$draws(variables = c("alpha", "beta"), format = "draws_df")
colnames(d) <- c("a", "b")

# Ahora podemos simular desde la distribución predictiva posterior.
# Para cada observación simularemos S = 4000 muestras de la predictiva posterior.
# Guardaremos las simulaciones en una matriz de N * S:

S <- nrow(d)     # número de muestras de la posterior
N <- nrow(datos) # número de observaciones

ysim <- matrix(NA, N, S)

# Simulamos
for (s in 1:S) {
  # calculamos lambda para la muestra "s"
  lambda <- exp(d$a[s] + d$b[s] * datos$fwi)
  # es un vector de largo N, porque datos$fwi también lo es.

  # llenamos una columna entera en ysim, que equivale a un set de datos simulado
  ysim[, s] <- rpois(N, lambda)
}


# visualizamos ejemplos
fila <- sample(1:N, 1)
fila
mm <- max(c(ysim[fila, ], datos$fires[fila]))
hist(ysim[fila, ], xlim = c(0, mm * 1.05))
abline(v = datos$fires[fila], col = "red", lwd = 3, lty = 2)

# Ahora cada fila tiene muestras de la distribución predictiva posterior de y
# para el valor correspondiente de fwi

res <- createDHARMa(
  simulatedResponse = ysim,
  observedResponse = datos$fires,
  integerResponse = TRUE # se da cuenta solo, pero por las dudas
)

plot(res)
hist(res$scaledResiduals)

# Podemos ver si hay alguna tendencia en función de una predictora
plotResiduals(res, form = datos$fwi)

# Agregamos los PIT calculados por DHARMa a la tabla de datos
datos$pit <- res$scaledResiduals

ggplot(datos, aes(pp, pit)) +
  geom_point(size = 1.5) +
  labs(y = "Prob. acumulada",  # o "residuos dharma"
       x = "Precipitación (mm)") +
  geom_smooth(method = "gam")

# No parece haber una tendencia, lo cual sugiere que no ganaríamos mucho
# incluyendo la precipitación.

# Trasponemos la matriz de datos simulados porque a bayesplot le gustan las
# iteraciones de la posterior como filas (las teníamos en las columnas)
ysim_t <- t(ysim)

# posterior predictive intervals en función de algo
ppc_intervals(
  y = datos$fires, # observaciones
  yrep = ysim_t,   # simulaciones de la predictiva posterior
  x = datos$fwi,   # alguna predictora de interés
  prob = 0.5       # probabilidad de los intervalos
) +
  labs(x = "FWI", y = "Número de incendios")

# Y la densidad marginal de la respuesta, según lo observado (y) y los sets de
# datos simulados (yrep). Como graficará una curva por set de datos simulado,
# será mejor elegir al azar unos pocos (menos que 4000)
use <- sample(1:S, 50) # elige 50 al azar
ppc_dens_overlay(y = datos$fires, yrep = ysim_t[use, ])

# Predicciones ------------------------------------------------------------

nrep <- 200 # largo de la secuencia de FWI
fwi_seq <- seq(min(datos$fwi), max(datos$fwi), length.out = nrep)

# Guardaremos la función evaluada sobre la secuencia de fwi en una matriz,
# donde cada columna corresponderá a una muestra de la posterior
lambda_mat <- matrix(NA, nrep, S)

for (s in 1:S) {
  lambda_mat[, s] <- exp(d$a[s] + d$b[s] * fwi_seq)
}

# Podemos graficar algunas curvas ("spaghetti plot").
# Ordenaremos la matriz para graficarla con ggplot
ncurves <- 50
use <- sample(1:S, size = ncurves, replace = F) # elegimos 150 curvas al azar
lambda_mat_sub <- as.data.frame(lambda_mat[, use])
lambda_mat_sub$fwi <- fwi_seq

# Elongamos la matriz
llong <- pivot_longer(
  lambda_mat_sub, cols = all_of(1:ncurves),
  names_to = "sim", values_to = "lambda"
)

ggplot(llong, aes(fwi, lambda, group = sim)) +
  geom_line(alpha = 0.1) +
  labs(x = "FWI", y = expression(lambda))

# Lambda es una cantidad derivada de alpha, beta y el fwi. Todo lo que se calcula
# a partir de los parámetros tiene su propia distribución posterior. Y acá
# tenemos, además, una distribución posterior de lambda para cada valor de
# fwi que evaluemos.

# Otra forma de graficar la función es resumir la distribución de cada lambda,
# es decir, la distribución de lambda marginal a alpha y beta para cada x.

pred <- data.frame(
  fwi = fwi_seq,
  lambda_mean = rowMeans(lambda_mat),
  lambda_lwr = apply(lambda_mat, MARGIN = 1, quantile, probs = 0.025),
  lambda_upr = apply(lambda_mat, 1, quantile, probs = 0.975)
)

# apply() es una forma concisa de loopear. equivale a hacer
for (i in 1:N) {
  pred$lambda_lwr[i] <- quantile(lambda_mat[i, ], probs = 0.025)
}
# Es decir, aplica la función quantile, con argumento probs = 0.025, sobre el
# array llamado "lambda_mat", iterando sobre la dimensión 1 (filas)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(x = fwi, y = lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  geom_ribbon(alpha = 0.4, color = NA, fill = "orange") +
  geom_line() +
  labs(x = "FWI", y = expression(lambda)) +
  # Agregamos los datos
  geom_point(aes(fwi, fires), data = datos, inherit.aes = F,
             alpha = 0.7, size = 2)


# La cinta muestra los percentiles 2.5 y 97.5 % de la distribución posterior de
# lambda. Pero también podemos explorar los percentiles extremos de la
# distribución predictiva posterior de y.
# Para eso, simulamos nuevas observaciones sobre la secuencia de FWI.
# (Antes, para el chequeo posterior, habíamos simulado nuevos valores de y
# usando los valores observados de FWI, no sobre una secuencia.)

ysim_seq <- matrix(NA, nrep, S)
for (s in 1:S) {
  ysim_seq[, s] <- rpois(nrep, lambda = lambda_mat[, s])

  # plot(lambda_mat[, s] ~ fwi_seq, type = "l", ylim = c(0, 50))
  # points(ysim_seq[, s] ~ fwi_seq, pch = 19, col = rgb(1, 0, 0, 0.1))
}

# Obtenemos los cuantiles extremos:
pred$y_lwr <- apply(ysim_seq, 1, quantile, probs = 0.025)
pred$y_upr <- apply(ysim_seq, 1, quantile, probs = 0.975)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(fwi, lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  # Intervalo para y
  geom_ribbon(aes(fwi, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
              alpha = 0.3, color = NA, fill = "blue") +
  # Intervalo para lambda
  geom_ribbon(alpha = 0.6, color = NA, fill = "orange") +
  geom_line() +
  labs(x = "FWI", y = "Número de incendios") +
  # Agregamos los datos
  geom_point(aes(fwi, fires), data = datos, inherit.aes = F,
             alpha = 0.7, size = 2)

# La cinta más ancha es donde esperamos ver los datos con 95 % de probabilidad;
# la cinta más fina es donde esperamos ver la media de muchos datos, con 95 %
# de probabilidad.
# El más ancho se suele llamar intervalo de predicción; el más fino, intervalo de
# credibilidad para la media.


# R2 ----------------------------------------------------------------------

# Vector para almacenar el R2
R2 <- numeric(S)

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # predichos para cada observación
  mu <- exp(d$a[i] + d$b[i] * datos$fwi) # es lambda!

  # varianza de los predichos = varianza explicada
  var_fit <- var(mu)

  # varianza residual para cada observación. Esto depende de la
  # distribución. En la poisson, la varianza es igual al media
  var_res <- mu # sólo para ser explícitos
  # y tomamos su media entre observaciones
  var_res_mean <- mean(var_res)

  # R2
  R2[i] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")

# Modelo 02: binomial negativa -------------------------------------------

# Compilamos
model2 <- cmdstan_model(file.path("modelos", "nfuegos_02.stan"))

# Y muestreamos la posterior
fit2 <- model2$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Verificamos el HMC
fit2$cmdstan_diagnose()

# Resultados del ajuste ---------------------------------------------------

# Resumen de las marginales (lo que le pidamos en "variables")
summ <- fit2$summary(variables = c("alpha", "beta", "phi"))
print(summ)

# Extraemos las muestras para calculas más cosas
d <- fit2$draws(variables = c("alpha", "beta", "phi"), format = "draws_df")

# calcular percentiles 2.5 y 97.5 %, el intervalo de credibilidad del 95 %
quantile(d$alpha, probs = c(0.025, 0.975))
quantile(d$beta, probs = c(0.025, 0.975))
quantile(d$phi, probs = c(0.025, 0.975))

# Algunos enunciados probabilísticos que pueden interesar (dependiendo del
# contexto)

S <- nrow(d) # cantidad de muestras de la posterior

# Pr(beta > 0)
sum(d$beta > 0) / S

# Pr(exp(beta) > 1.2)
ebeta <- exp(d$beta)
sum(ebeta > 1.2) / S

# Pr(lambda > 20 | FWI = 18)
lambda <- exp(d$alpha + d$beta * 18)
sum(lambda > 20) / S

# Comparación de previa contra posterior.
# Usamos density() para obtener la densidad empírica a partir de muestras,
# y curve() para dibujar la previa.

# alpha
plot(density(d$alpha), main = NA, xlab = expression(alpha), xlim = c(-5, 5))
curve(dnorm(x, 0, 1), col = 4, add = T)

# beta
plot(density(d$beta), main = NA, xlab = expression(beta), xlim = c(-0.4, 0.4))
curve(dnorm(x, 0, 0.1), col = 4, add = T)

# phi
plot(density(d$phi), main = NA, xlab = expression(phi), xlim = c(0, 100))
curve(dgamma(x, 2, 0.1), col = 4, add = T)


# Verificación predictiva posterior  --------------------------------------

# Extraemos las muestras de alpha y beta, en formato data.frame
d <- fit2$draws(variables = c("alpha", "beta", "phi"), format = "draws_df")
colnames(d) <- c("a", "b", "phi")

# Ahora podemos simular desde la distribución predictiva posterior.
# Para cada observación simularemos S = 4000 muestras de la predictiva posterior.
# Guardaremos las simulaciones en una matriz de N * S:

S <- nrow(d)     # número de muestras de la posterior
N <- nrow(datos) # número de observaciones

ysim <- matrix(NA, N, S)

# Simulamos
for (s in 1:S) {
  # calculamos lambda para la muestra "s"
  lambda <- exp(d$a[s] + d$b[s] * datos$fwi)
  # es un vector de largo N, porque datos$fwi también lo es.

  # llenamos una columna entera en ysim, que equivale a un set de datos simulado
  ysim[, s] <- rnbinom(N, mu = lambda, size = d$phi[s])
}

# Ahora cada fila tiene muestras de la distribución predictiva posterior de y
# para el valor correspondiente de fwi

res <- createDHARMa(
  simulatedResponse = ysim,
  observedResponse = datos$fires,
  integerResponse = TRUE # se da cuenta solo, pero por las dudas
)

plot(res)
hist(res$scaledResiduals)

# Podemos ver si hay alguna tendencia en función de una predictora
plotResiduals(res, form = datos$fwi)

# Agregamos los PIT calculados por DHARMa a la tabla de datos
datos$pit <- res$scaledResiduals

ggplot(datos, aes(pp, pit)) +
  geom_point(size = 1.5) +
  labs(y = "Prob. acumulada",  # o "residuos dharma"
       x = "Precipitación (mm)") +
  geom_smooth(method = "gam")

# No parece haber una tendencia, lo cual sugiere que no ganaríamos mucho
# incluyendo la precipitación.

# Trasponemos la matriz de datos simulados porque a bayesplot le gustan las
# iteraciones de la posterior como filas (las teníamos en las columnas)
ysim_t <- t(ysim)

# posterior predictive intervals en función de algo
ppc_intervals(
  y = datos$fires, # observaciones
  yrep = ysim_t,   # simulaciones de la predictiva posterior
  x = datos$fwi,   # alguna predictora de interés
  prob = 0.9       # probabilidad de los intervalos
) +
  labs(x = "FWI", y = "Número de incendios")

# Y la densidad marginal de la respuesta, según lo observado (y) y los sets de
# datos simulados (yrep). Como graficará una curva por set de datos simulado,
# será mejor elegir al azar unos pocos (menos que 4000)
use <- sample(1:S, 50) # elige 50 al azar
ppc_dens_overlay(y = datos$fires, yrep = ysim_t[use, ])

# Predicciones ------------------------------------------------------------

nrep <- 200 # largo de la secuencia de FWI
fwi_seq <- seq(min(datos$fwi), max(datos$fwi), length.out = nrep)

# Guardaremos la función evaluada sobre la secuencia de fwi en una matriz,
# donde cada columna corresponderá a una muestra de la posterior
lambda_mat <- matrix(NA, nrep, S)

for (s in 1:S) {
  lambda_mat[, s] <- exp(d$a[s] + d$b[s] * fwi_seq)
}

# Podemos graficar algunas curvas ("spaghetti plot").
# Ordenaremos la matriz para graficarla con ggplot
ncurves <- 50
use <- sample(1:S, ncurves) # elegimos 150 curvas al azar
lambda_mat_sub <- as.data.frame(lambda_mat[, use])
lambda_mat_sub$fwi <- fwi_seq

# Elongamos la matriz
llong <- pivot_longer(
  lambda_mat_sub, cols = all_of(1:ncurves),
  names_to = "sim", values_to = "lambda"
)

ggplot(llong, aes(fwi, lambda, group = sim)) +
  geom_line(alpha = 0.1) +
  labs(x = "FWI", y = expression(lambda))

# Lambda es una cantidad derivada de alpha, beta y el fwi. Todo lo que se calcula
# a partir de los parámetros tiene su propia distribución posterior. Y acá
# tenemos, además, una distribución posterior de lambda para cada valor de
# fwi que evaluemos.

# Otra forma de graficar la función es resumir la distribución de cada lambda,
# es decir, la distribución de lambda marginal a alpha y beta para cada x.

pred <- data.frame(
  fwi = fwi_seq,
  lambda_mean = rowMeans(lambda_mat),
  lambda_lwr = apply(lambda_mat, 1, quantile, probs = 0.025),
  lambda_upr = apply(lambda_mat, 1, quantile, probs = 0.975)
)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(fwi, lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  geom_ribbon(alpha = 0.4, color = NA, fill = "orange") +
  geom_line() +
  labs(x = "FWI", y = expression(lambda)) +
  # Agregamos los datos
  geom_point(aes(fwi, fires), data = datos, inherit.aes = F,
             alpha = 0.7, size = 2)


# La cinta muestra los percentiles 2.5 y 97.5 % de la distribución posterior de
# lambda. Pero también podemos explorar los percentiles extremos de la
# distribución predictiva posterior de y.
# Para eso, simulamos nuevas observaciones sobre la secuencia de FWI.
# (Antes, para el chequeo posterior, habíamos simulado nuevos valores de y
# usando los valores observados de FWI, no sobre una secuencia.)

ysim_seq <- matrix(NA, nrep, S)
for (s in 1:S) {
  ysim_seq[, s] <- rnbinom(nrep, mu = lambda_mat[, s], size = d$phi[s])
}

# Obtenemos los cuantiles extremos:
pred$y_lwr <- apply(ysim_seq, 1, quantile, probs = 0.025)
pred$y_upr <- apply(ysim_seq, 1, quantile, probs = 0.975)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(fwi, lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  # Intervalo para y
  geom_ribbon(aes(fwi, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
              alpha = 0.3, color = NA, fill = "blue") +
  # Intervalo para lambda
  geom_ribbon(alpha = 0.6, color = NA, fill = "orange") +
  geom_line() +
  labs(x = "FWI", y = expression(lambda)) +
  # Agregamos los datos
  geom_point(aes(fwi, fires), data = datos, inherit.aes = F,
             alpha = 0.7, size = 2)

# La cinta más ancha es donde esperamos ver los datos con 95 % de probabilidad;
# la cinta más fina es donde esperamos ver la media de muchos datos, con 95 %
# de probabilidad.
# El más ancho se suele llamar intervalo de predicción; el más fino, intervalo de
# credibilidad para la media.


# R2 ----------------------------------------------------------------------

# Vector para almacenar el R2
R2 <- numeric(S)

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # predichos para cada observación
  mu <- exp(d$a[i] + d$b[i] * datos$fwi) # es lambda!

  # varianza de los predichos = varianza explicada
  var_fit <- var(mu)

  # varianza residual para cada observación. Esto depende de la
  # distribución. En la binomial negativa, esta es la ecuación:
  var_res <- mu + mu ^ 2 / d$phi[i]

  # y tomamos su media entre observaciones
  var_res_mean <- mean(var_res)

  # R2
  R2[i] <- var_fit / (var_fit + var_res_mean)
}

plot(density(R2, from = 0, to = 1), main = NA, xlab = "R2")