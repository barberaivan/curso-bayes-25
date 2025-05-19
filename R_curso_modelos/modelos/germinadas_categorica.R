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

# Para ejemplificar el uso de una predictora categórica, categorizamos el daño
# por fuego en bajo, medio y alto. En general esto no es una práctica deseable,
# sólo lo hacemos para tener un ejemplo.

# Creamos predictora categórica (normalmente, ya la tendríamos en la tabla;
# este paso no sería necesario).
datos$danio_cat <- cut(
  datos$danio, quantile(datos$danio, probs = c(0, 0.33, 0.66, 1)),
  include.lowest = T
)
# renombramos las categorías
datos$danio_cat <- factor(datos$danio_cat, labels = c("bajo", "medio", "alto"))

# Para codificar un modelo (logit) lineal, necesitamos tantas dummies como niveles
# menos 1. En este caso, tomamos "bajo" como el nivel de referencia, por lo que
# creamos las dummies M y A (para medio y alto):
datos$DM <- as.numeric(datos$danio_cat == "medio")
datos$DA <- as.numeric(datos$danio_cat == "alto")

# Datos para Stan
stan_data <- list(
  N = nrow(datos),
  y = datos$germinadas,
  S = 120,

  DM = datos$DM,
  DA = datos$DA
)

# Ajuste del modelo (MCMC) ------------------------------------------------

# Compilamos
model1 <- cmdstan_model(here::here("modelos", "germinadas_01_categorica.stan"))

# Muestreamos
fit1 <- model1$sample(data = stan_data, parallel_chains = 4)

# Chequeamos MCMC
fit1$cmdstan_diagnose()

# Resultados del ajuste ---------------------------------------------------

summ <- fit1$summary(c("alpha", "betaM", "betaA"))
print(summ)

# Previa vs. posterior ----------------------------------------------------

# Extraemos muestras de los parámetros que necesitamos
d <- fit1$draws(c("alpha", "betaM", "betaA"), format = "draws_df")
colnames(d) <- c("a", "bM", "bA")
S <- nrow(d)
N <- nrow(datos)

# alpha (intercept)
plot(density(d$a), main = NA, xlab = expression(alpha), xlim = c(-2, 2))
curve(dnorm(x, 0, 10), add = T, col = 4) # previa definida en Stan

# betaM
plot(density(d$bM), main = NA, xlab = expression(beta~M), xlim = c(-2, 2))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan

# betaA
plot(density(d$bA), main = NA, xlab = expression(beta~A), xlim = c(-2, 2))
curve(dnorm(x, 0, 5), add = T, col = 4) # previa definida en Stan


# Interpretación de los parámetros ----------------------------------------

# betaM y betaA son diferencias en escala logit con respecto a daño bajo. Como
# son poco informativas, es mejor obtener la probabilidad de germinación a partir
# de ellos y compararla entre niveles de daño.
# En Stan llamabamos mu a la probabilidad de germinación; acá seremos más
# explícitos y le llamaremos p.
# En R, la inversa del logit se llama plogis (en Stan se llama inv_logit)

# p(danio = bajo) = plogis(alpha)
# p(danio = medio) = plogis(alpha + betaM)
# p(danio = alto) = plogis(alpha + betaA)

# Calculamos los lambdas de cada categoría de esa manera:
d$p_bajo = plogis(d$a)
d$p_medio = plogis(d$a + d$bM)
d$p_alto = plogis(d$a + d$bA)

# Los p_ son la probabilidad de germinación según la categoría de daño por fuego.

# Algunos enunciados probabilísticos que pueden interesar (dependiendo del
# contexto)
S <- nrow(d) # cantidad de muestras de la posterior

# Cuán probable es que con daño medio tenga menor probabilidad de germinación
# que con daño bajo?
# Pr(p_medio < p_bajo)
sum(d$p_medio < d$p_bajo) / S

# Este tipo de comparaciones pueden hacerse entre todos los pares de niveles
# de un factor. Serían análogos a las pruebas de comparaciones múltiples
# (e.g., Tukey), pero más flexibles.

# Verificación predictiva posterior ---------------------------------------

# Verificamos el ajuste simulando de la predictiva posterior
ysim <- matrix(NA, N, S)
for (i in 1:S) {
  mu <- plogis(d$a[i] + d$bM[i] * datos$DM + d$bA[i] * datos$DA)
  # plogis equivale a inv_logit = logit inverso

  ysim[, i] <- rbinom(N, size = datos$sembradas, prob = mu)
}

ysim_t <- t(ysim)

# Miramos con bayesplot

# Como no nos deja poner una predictora categórica en el x, utilizamos el ID de
# cada observación
datos$obs_id <- 1:nrow(datos)

ppc_intervals(
  datos$germinadas, yrep = ysim_t, x = datos$obs_id,
  prob = 0.7, prob_outer = 0.9
) +
  labs(
    y = "Número de semillas germinadas",
    x = "Observación (ID)"
  )

# Densidad marginal predicha y observada:
use <- sample(1:S, 100)
ppc_dens_overlay(datos$germinadas, yrep = ysim_t[use, ])

# DHARMa
res <- createDHARMa(ysim, datos$germinadas, integerResponse = T)
plot(res)
plotResiduals(res, form = datos$danio_cat, rank = F)

# Predicciones ------------------------------------------------------------

# Para graficar la predicción con respecto a una variable categórica calcularemos
# la media predicha y el intervalo de predicción para todos sus valores posibles
# (bajo, medio, alto)
danio_seq <- c("bajo", "medio", "alto")
nrep <- length(danio_seq)

# Creamos las dummies asociadas
danio_seq_M <- as.numeric(danio_seq == "medio")
danio_seq_A <- as.numeric(danio_seq == "alto")

# Guardaremos la función evaluada sobre la secuencia de daño en una matriz,
# donde cada columna corresponderá a una muestra de la posterior
p_mat <- matrix(NA, nrep, S)

for (s in 1:S) {
  p_mat[, s] <- plogis(d$a[s] + d$bM[s] * danio_seq_M + d$bA[s] * danio_seq_A)
}

# Resumimos la distribución posterior de cada lambda.
pred <- data.frame(
  danio = danio_seq,
  p_mean = rowMeans(p_mat),
  p_lwr = apply(p_mat, MARGIN = 1, quantile, probs = 0.025),
  p_upr = apply(p_mat, MARGIN = 1, quantile, probs = 0.975)
)

# Ahora graficamos la predicción media y su incertidumbre. Pero para agregar
# los datos, que son conteos, los convertimos a proporción:
datos$germ_prop <- datos$germinadas / datos$sembradas

ggplot(pred, aes(x = danio, y = p_mean, ymin = p_lwr, ymax = p_upr)) +
  geom_linerange(color = "orange", size = 2) +
  geom_point(size = 3, color = "black", fill = "red", shape = 21) +
  labs(x = "Daño por fuego", y = "Probabilidad de germinación") +
  # Agregamos los datos
  geom_point(
    aes(danio_cat, germ_prop), data = datos, inherit.aes = F,
    alpha = 0.7, size = 2,
    position = position_jitter(width = 0.07)
  )

# En negro, los datos, con un jitter horizontal para que no se apilen.
# En naranja, el intervalo de credibilidad del 95 % para mu o p, y en rojo,
# la media de la distribución posterior de mu o p

# También podemos graficar percentiles de la distribución predictiva posterior
# de y. Para eso, simulamos nuevas observaciones sobre la secuencia de daño.
# (Antes, para el chequeo posterior, habíamos simulado nuevos valores de y
# usando los valores observados de daño, no sobre una secuencia.)

ysim_seq <- matrix(NA, nrep, S)
for (s in 1:S) {
  ysim_seq[, s] <- rbinom(nrep, size = datos$sembradas[1], prob = p_mat[, s])
}

# Para graficar junto a la probabilidad, dividimos los datos simulados (conteos
# de germinadas) por lo sembrado
ysim_seq_prop <- ysim_seq / datos$sembradas[1] # todas las filas de
# datos$sembradas tienen el mismo valor, por eso tomamos sólo uno (no importa cuál)

# Obtenemos los cuantiles extremos:
pred$y_lwr <- apply(ysim_seq_prop, 1, quantile, probs = 0.025)
pred$y_upr <- apply(ysim_seq_prop, 1, quantile, probs = 0.975)

# Ahora graficamos la predicción media y su incertidumbre
ggplot(pred, aes(danio, p_mean, ymin = p_lwr, ymax = p_upr)) +
  # Intervalo para y
  geom_linerange(
    aes(danio, ymin = y_lwr, ymax = y_upr), inherit.aes = F,
    alpha = 0.4, color = "blue", size = 3
  ) +
  # Intervalo para p (o mu)
  geom_linerange(color = "orange", size = 2) +
  geom_point(size = 3, color = "black", fill = "red", shape = 21) +
  labs(x = "Daño por fuego", y = "Proporción de semillas germinadas") +
  # Agregamos los datos
  geom_point(
    aes(danio_cat, germ_prop), data = datos, inherit.aes = F,
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

# Por comodidad, definimos las sembradas
SS <- datos$sembradas[1]

# Loop sobre muestras de la posterior
for (i in 1:S) {
  # En la beta-binomial, la media es mu * size,
  # siendo size la cantidad de sembradas

  # predichos para cada observación
  mu <- plogis(d$a[i] + d$bM[i] * datos$DM + d$bA[i] * datos$DA)
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
