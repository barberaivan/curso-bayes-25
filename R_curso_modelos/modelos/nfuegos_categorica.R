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

# Para ejemplificar el uso de una predictora categórica, categorizamos el FWI
# en bajo, medio y alto. En general esto no es una práctica deseable, sólo lo
# hacemos para tener un ejemplo.

# Creamos predictora categórica (normalmente, ya la tendríamos en la tabla;
# este paso no sería necesario).
datos$fwi_cat <- cut(
  datos$fwi, quantile(datos$fwi, probs = c(0, 0.33, 0.66, 1)),
  include.lowest = T
)
# renombramos las categorías
datos$fwi_cat <- factor(datos$fwi_cat, labels = c("bajo", "medio", "alto"))

# Para codificar un modelo (log) lineal, necesitamos tantas dummies como niveles
# menos 1. En este caso, tomamos "bajo" como el nivel de referencia, por lo que
# creamos las dummies M y A (para medio y alto):
datos$M <- as.numeric(datos$fwi_cat == "medio")
datos$A <- as.numeric(datos$fwi_cat == "alto")

# (incluirlas como columnas de "datos" es opcional)

# Datos ordenados para Stan:
# Creamos una lista nombrada para pasarle los datos. Los nombres de cada
# elemento de la lista tienen que ser exactamente los que definimos en
# la sección data {}
stan_data <- list(
  N = nrow(datos),
  y = datos$fires,
  M = datos$M,
  A = datos$A
)


# Ajuste del modelo (MCMC) ------------------------------------------------

# Compilamos
model1 <- cmdstan_model(file.path("modelos", "nfuegos_01_categorica.stan"))

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
summ <- fit1$summary(variables = c("alpha", "betaM", "betaA"))
print(summ)

# Extraemos las muestras para calculas más cosas
d <- fit1$draws(variables = c("alpha", "betaM", "betaA"), format = "draws_df")

# calcular percentiles 2.5 y 97.5 %, el intervalo de credibilidad del 95 %
quantile(d$alpha, probs = c(0.025, 0.975))
quantile(d$betaM, probs = c(0.025, 0.975))
quantile(d$betaA, probs = c(0.025, 0.975))

# Comparación de previa contra posterior.
# Usamos density() para obtener la densidad empírica a partir de muestras,
# y curve() para dibujar la previa.

# alpha
plot(density(d$alpha), main = NA, xlab = expression(alpha), xlim = c(-5, 5))
curve(dnorm(x, 0, 10), col = 4, add = T)

# betaM
plot(density(d$betaM), main = NA, xlab = expression(beta~M), xlim = c(-5, 5))
curve(dnorm(x, 0, 10), col = 4, add = T)

# betaA
plot(density(d$betaA), main = NA, xlab = expression(beta~A), xlim = c(-5, 5))
curve(dnorm(x, 0, 10), col = 4, add = T)


# betaM y betaA son diferencias en escala log con respecto a FWI bajo:

# lambda(FWI = bajo) = exp(alpha)
# lambda(FWI = medio) = exp(alpha + betaM)
# lambda(FWI = alto) = exp(alpha + betaA)

# Calculamos los lambdas de cada categoría de esa manera:
d$lambda_bajo = exp(d$alpha)
d$lambda_medio = exp(d$alpha + d$betaM)
d$lambda_alto = exp(d$alpha + d$betaA)

# Los lambdas son el número esperado de incendios según la categoría de FWI.

# Algunos enunciados probabilísticos que pueden interesar (dependiendo del
# contexto)
S <- nrow(d) # cantidad de muestras de la posterior

# Cuán probable es que con FWI medio tenga mayor número de incendios promedio
# que con FWI bajo?
# Pr(lambda_medio > lambda_bajo)
sum(d$lambda_medio > d$lambda_bajo) / S

# Cuán probable es que con FWI alto tenga mayor número de incendios promedio
# que con FWI medio?
# Pr(lambda_alto > lambda_medio)
sum(d$lambda_alto > d$lambda_medio) / S

# Cuán probable es que con FWI alto tenga un número de incendios promedio que
# al menos triplique el número de incendios promedio de FWI bajo?
# Pr[(lambda_alto / lambda_bajo) >= 3]
q <- d$lambda_alto / d$lambda_bajo
sum(q >= 3) / S

# Este tipo de comparaciones pueden hacerse entre todos los pares de niveles
# de un factor. Serían análogos a las pruebas de comparaciones múltiples
# (e.g., Tukey), pero más flexibles.


# Verificación predictiva posterior  --------------------------------------

# Extraemos las muestras de alpha y beta, en formato data.frame
d <- fit1$draws(variables = c("alpha", "betaM", "betaA"), format = "draws_df")
colnames(d) <- c("a", "bM", "bA")

# Ahora podemos simular desde la distribución predictiva posterior.
# Para cada observación simularemos S = 4000 muestras de la predictiva posterior.
# Guardaremos las simulaciones en una matriz de N * S:

S <- nrow(d)     # número de muestras de la posterior
N <- nrow(datos) # número de observaciones

ysim <- matrix(NA, N, S)

# Simulamos
for (s in 1:S) {
  # calculamos lambda para la muestra "s"
  lambda <- exp(d$a[s] + d$bM[s] * datos$M + d$bA * datos$A)
  # es un vector de largo N, porque datos$M y datos$A también lo son.

  # llenamos una columna entera en ysim, que equivale a un set de datos simulado
  ysim[, s] <- rpois(N, lambda)
}

# Ahora cada fila tiene muestras de la distribución predictiva posterior de y
# para el valor correspondiente de fwi categorizado.
res <- createDHARMa(
  simulatedResponse = ysim,
  observedResponse = datos$fires,
  integerResponse = TRUE # se da cuenta solo, pero por las dudas
)

plot(res)
hist(res$scaledResiduals)

# Podemos ver si hay alguna tendencia en función de una predictora
plotResiduals(res, form = datos$fwi_cat)

# Trasponemos la matriz de datos simulados porque a bayesplot le gustan las
# iteraciones de la posterior como filas (las teníamos en las columnas)
ysim_t <- t(ysim)

# posterior predictive intervals en función de algo.
# Cuando la predictora es categórica, no podemos usar ppc_intervals con la
# predictora en x. Por eso mejor usemos "obs_id", la fila de cada dato.
datos$obs_id <- 1:nrow(datos)
ppc_intervals(
  y = datos$fires, # observaciones
  yrep = ysim_t,   # simulaciones de la predictiva posterior
  x = datos$obs_id,   # alguna predictora de interés
  prob = 0.5       # probabilidad de los intervalos
) +
  labs(x = "Observación (ID)", y = "Número de incendios")

# Y la densidad marginal de la respuesta, según lo observado (y) y los sets de
# datos simulados (yrep). Como graficará una curva por set de datos simulado,
# será mejor elegir al azar unos pocos (menos que 4000)
use <- sample(1:S, 50) # elige 50 al azar
ppc_dens_overlay(y = datos$fires, yrep = ysim_t[use, ])

# Predicciones ------------------------------------------------------------

# Para graficar la predicción con respecto a una variable categórica calcularemos
# la media predicha y el intervalo de predicción para todos sus valores posibles
# (bajo, medio, alto)
fwi_seq <- c("bajo", "medio", "alto")
nrep <- length(fwi_seq)

# Creamos las dummies asociadas
fwi_seq_M <- as.numeric(fwi_seq == "medio")
fwi_seq_A <- as.numeric(fwi_seq == "alto")

# Guardaremos la función evaluada sobre la secuencia de fwi en una matriz,
# donde cada columna corresponderá a una muestra de la posterior
lambda_mat <- matrix(NA, nrep, S)

for (s in 1:S) {
  lambda_mat[, s] <- exp(d$a[s] + d$bM[s] * fwi_seq_M + d$bA[s] * fwi_seq_A)
}

# Resumimos la distribución posterior de cada lambda.
pred <- data.frame(
  fwi = fwi_seq,
  lambda_mean = rowMeans(lambda_mat),
  lambda_lwr = apply(lambda_mat, MARGIN = 1, quantile, probs = 0.025),
  lambda_upr = apply(lambda_mat, MARGIN = 1, quantile, probs = 0.975)
)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(x = fwi, y = lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  geom_linerange(color = "orange", size = 2) +
  geom_point(size = 3, color = "black", fill = "red", shape = 21) +
  labs(x = "FWI", y = expression(lambda)) +
  # Agregamos los datos
  geom_point(
    aes(fwi_cat, fires), data = datos, inherit.aes = F,
    alpha = 0.7, size = 2,
    position = position_jitter(width = 0.07)
  )

# En negro, los datos, con un jitter horizontal para que no se apilen.
# En naranja, el intervalo de credibilidad del 95 % para lambda, y en rojo,
# la media de la distribución posterior de lambda.

# También podemos graficar percentiles de la distribución predictiva posterior
# de y. Para eso, simulamos nuevas observaciones sobre la secuencia de FWI.
# (Antes, para el chequeo posterior, habíamos simulado nuevos valores de y
# usando los valores observados de FWI, no sobre una secuencia.)

ysim_seq <- matrix(NA, nrep, S)
for (s in 1:S) {
  ysim_seq[, s] <- rpois(nrep, lambda = lambda_mat[, s])
}

# Obtenemos los cuantiles extremos:
pred$y_lwr <- apply(ysim_seq, 1, quantile, probs = 0.025)
pred$y_upr <- apply(ysim_seq, 1, quantile, probs = 0.975)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(fwi, lambda_mean, ymin = lambda_lwr, ymax = lambda_upr)) +
  # Intervalo para y
  geom_linerange(
    aes(fwi, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
    alpha = 0.4, color = "blue", size = 3
  ) +
  # Intervalo para lambda
  geom_linerange(color = "orange", size = 2) +
  geom_point(size = 3, color = "black", fill = "red", shape = 21) +
  labs(x = "FWI", y = "Número de incendios") +
  # Agregamos los datos
  geom_point(
    aes(fwi_cat, fires), data = datos, inherit.aes = F,
    alpha = 0.7, size = 2,
    position = position_jitter(width = 0.07)
  )

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
  mu <- exp(d$a[i] + d$bM[i] * datos$M + d$bA[i] * datos$A) # es lambda!

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
