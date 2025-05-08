## ## TP 5: Estimación de parámetros e inferencia bayesiana ## ##

# [Codificación: UTF-8]

# Paquetes ----------------------------------------------------------------

library(tidyverse) # ggplot y pivot_longer
library(viridis)   # colores
library(dplyr)
library(ggplot2)
library(scales)    # leyenda de gráfico
library(patchwork) # unir gráficos

# Funciones ---------------------------------------------------------------

# Para normalizar densidades o likelihoods (con fines gráficos).
# area es el tamaño del segmento o celda en la aproximación discreta.
# En 1D, es la distancia entre valores de la secuencia de parámetros.
# En 2D, es el área de cada celda.
# Devuelve un vector de densidad o likelihood, tal que
# sum(normalized_density * area) ~= 1
normalize_dens <- function(dens, area) {
  dens / sum(dens * area)
}

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

# Actividad 1 -------------------------------------------------------------

# 1.1 ---------------------------------------------------------------------
# Elegir un set de datos y definir un modelo sencillo (1 o 2 parámetros).
# Puede ser estadístico o no (la función de costo no necesariamente tiene que
# ser una -log verosimilitud).

# y_i ~ Poisson(lambda)

# (no hay predictoras)

# 1.2 ---------------------------------------------------------------------
# Evaluar qué restricciones debería imponerse sobre los parámetros para que el
# modelo sea sensato.

# lambda >= 0, ¿quizás no superior a 200?

# 1.3 ---------------------------------------------------------------------
# Graficar la función de verosimilitud del modelo (o la función de costo) y
# encontrá el máximo utilizando el método de la grilla. Si el modelo tiene
# D > 2 parámetros, se deberán fijar D - 2 para poder graficar. Graficar la
# predicción en función de alguna variable predictora (si es que la hay)
# utilizando el valor de la grilla con mayor verosimilitud o menor costo.

like_fire <- function(lambda, log = T, filas = 1:nrow(datos)) {
  # "like" por likelihood (sería raro decirle "vero");
  # "fire" porque es la función para estos datos de fuego, no cualquier
  # verosimilitud.

  # usando lambda evaluamos la verosimilitud de cada observación,
  # en escala log:
  like_pointwise <- dpois(datos$fires[filas], lambda, log = T)

  # la log-verosimilitud conjunta es la suma de las log-verosimilitudes por
  # observación:
  like <- sum(like_pointwise)

  if (log) return(like)
  else return(exp(like))
}

# Creamos grilla de valores de lambda.
side <- 500 # cantidad de valores por lado
grilla <- data.frame(
  lambda = seq(0, 30, length.out = side)
)

# Evaluamos la verosimilitud para cada valor
grilla$loglike <- NA
for (i in 1:side) {
  grilla$loglike[i] <- like_fire(grilla$lambda[i])
}
grilla$like <- exp(grilla$loglike)

# Elegimos el punto con mayor verosimilitud
row <- which.max(grilla$loglike)
lambda <- grilla$lambda[row]

# Curva de verosimilitud a partir de la grilla (1D)
p1 <- ggplot(grilla, aes(x = lambda, y = like)) +
  geom_line() +
  # geom_vline(xintercept = lambda, linetype = "dashed") +
  # scale_y_continous(trans="log10")+
  labs(x = expression(lambda), y = "Verosimilitud") +
  nice_theme()

# Plot de puntos + curva
p2 <- ggplot(datos, aes(x = fwi, y = fires)) +
  geom_point(color = rgb(0, 0, 0, 0.8), size = 2) +
  stat_function(fun = function(x) lambda,
                color = "blue", linewidth = 1) +
  ylab("y") + xlab("FWI") +
  nice_theme()

# Unimos
(p1 | p2)

# 1.4 ---------------------------------------------------------------------
# Comparár la estimación basada en la grilla con la obtenida usando optim.

opt <- optim(
  par = 1,                     # Valor inicial de lambda
  fn = like_fire,              # Función a optimizar
  method = "Brent",            # Optmización en 1D
  lower = 0, upper = 1000,     # Limites
  control = list(fnscale = -1) # Para que maximice en vez de minimizar
)

opt$par
grilla[row, c("lambda")]

# 1.5 ---------------------------------------------------------------------
# ¿Existen otros valores de los parámetros que también presentan un
# ajuste razonable a los datos, aunque sea subóptimo? ¿Cómo se podría
# escoger un set de valores razonables?

# 1.6 ---------------------------------------------------------------------
# Repetir los puntos 1.3 y 1.5 pero utilizando el 20 % y el 50 % de los
# datos disponibles, y compará los resultados.

# Ejemplo con N = 20 %

filas20 <- sample(1:nrow(datos), size = round(0.2 * nrow(datos)),
                  replace = F)

# Evaluamos la verosimilitud para cada valor
grilla$loglike20 <- NA
for (i in 1:side) {
  grilla$loglike20[i] <- like_fire(grilla$lambda[i], filas = filas20)
}
grilla$like20 <- exp(grilla$loglike20)

row20 <- which.max(grilla$loglike20)
lambda20 <- grilla$lambda[row20]

# Curva de verosimilitud a partir de la grilla (1D)
p1 <- ggplot(grilla, aes(x = lambda, y = like20)) +
  geom_line() +
  # geom_vline(xintercept = lambda, linetype = "dashed") +
  labs(x = expression(lambda), y = "Verosimilitud") +
  nice_theme()

# Plot de puntos + curva
p2 <- ggplot(datos[filas20, ], aes(x = fwi, y = fires)) +
  geom_point(color = rgb(0, 0, 0, 0.8), size = 2) +
  stat_function(fun = function(x) lambda20,
                color = "blue", linewidth = 1) +
  ylab("y") + xlab("FWI") +
  nice_theme()

# Unimos usando patchwork
(p1 | p2)


# Actividad 2 -------------------------------------------------------------

# 2.1 ---------------------------------------------------------------------
# Utilizando el mismo set de datos y modelo escogido previamente (ejercicio 1),
# definir previas para los parámetros, convirtiéndolo en un modelo bayesiano.

# y_i ~ Poisson(lambda)
# lambda ~ Gama(1 / tau, 1 / (mu * tau))
# tau = 1; mu = 10

# (tau es un parámetro de dispersión, mu es la media)

# 2.2 ---------------------------------------------------------------------
# Graficar la función de verosimilitud, la densidad previa y la densidad
# posterior.

# Creamos grilla de valores de lambda.
side <- 200 # cantidad de valores por lado

lambda_seq <- seq(0, 30, length.out = side)

grilla <- expand.grid(
  lambda = lambda_seq,
  loglike = NA,   # para llenar luego
  like = NA,
  prior = NA,
  post = NA
)

# Evaluamos la verosimilitud para cada valor
for (i in 1:side) {
  grilla$loglike[i] <- like_fire(grilla$lambda[i])
}
grilla$like <- exp(grilla$loglike)

# Definimos y evaluamos la previa (gama)
mu <- 10 # media
tau <- 1 # dispersion
a <- 1 / tau; b <- 1 / (tau * mu)

for (i in 1:side) {
  grilla$prior[i] <- dgamma(grilla$lambda[i], a, b)
}

# Posterior
grilla$post <- grilla$like * grilla$prior

# Normalizamos para poder visualizar
g <- grilla
area <- diff(lambda_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla[, v], area)
}

# Elongamos g
cols <- which(names(g) %in% c("like", "prior", "post"))
glong <- pivot_longer(g, all_of(cols), names_to = "variable",
                      values_to = "density")
glong$variable <- factor(
  glong$variable,
  levels = c("like", "prior", "post"),
  labels = c("Verosimilitud", "Previa", "Posterior")
)

# Graficamos
ggplot(glong, aes(x = lambda, y = density, ymax = density, ymin = 0,
                  color = variable, fill = variable)) +
  geom_line() +
  geom_ribbon(alpha = 0.4, color = NA) +
  scale_color_viridis(discrete = T, option = "C", end = 0.6) +
  scale_fill_viridis(discrete = T, option = "C", end = 0.6) +
  theme(legend.title = element_blank()) +
  xlab(expression(lambda)) +
  ylab("Densidad o Verosimilitud") +
  nice_theme()


# 2.3 ---------------------------------------------------------------------
# Explorar distintas previas, ya sea variando los parámetros de una
# misma distribución o utilizando otra distribución. Intentar crear
# previas que influyan en el resultado y otras que no.

grilla2 <- grilla

# Redefinimos y evaluamos la previa
mu <- 20 # media
tau <- 0.01 # dispersion
a <- 1 / tau; b <- 1 / (tau * mu)

for (i in 1:side) {
  grilla2$prior[i] <- dgamma(grilla$lambda[i], a, b)
}

# Posterior
grilla2$post <- grilla2$like * grilla2$prior

# Normalizamos para poder visualizar
g <- grilla2
area <- diff(lambda_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla2[, v], area)
}

# Elongamos g
cols <- which(names(g) %in% c("like", "prior", "post"))
glong <- pivot_longer(g, all_of(cols), names_to = "variable",
                      values_to = "density")
glong$variable <- factor(
  glong$variable,
  levels = c("like", "prior", "post"),
  labels = c("Verosimilitud", "Previa", "Posterior")
)

# Graficamos
ggplot(glong, aes(x = lambda, y = density, ymax = density, ymin = 0,
                  color = variable, fill = variable)) +
  geom_line() +
  geom_ribbon(alpha = 0.4, color = NA) +
  scale_color_viridis(discrete = T, option = "C", end = 0.6) +
  scale_fill_viridis(discrete = T, option = "C", end = 0.6) +
  theme(legend.title = element_blank()) +
  xlab(expression(lambda)) +
  ylab("Densidad o Verosimilitud") +
  nice_theme()

## Para probar con otra distribución, editar estas líneas:

# Redefinimos y evaluamos la previa
mu <- 20    # media
tau <- 0.01 # dispersion
a <- 1 / tau; b <- 1 / (tau * mu)

for (i in 1:side) {
  grilla2$prior[i] <- dgamma(grilla$lambda[i], a, b)
}


# (El resto del código, igual)

# 2.4 ---------------------------------------------------------------------
# Repetír la comparación entre posteriores pero utilizando únicamente 2 o 3
# observaciones. Las previas que no influían en el resultado usando todos los
# datos, ¿siguen sin influir?

# Creamos grilla de valores de lambda.
side <- 200 # cantidad de valores por lado

lambda_seq <- seq(0, 30, length.out = side)

grilla <- expand.grid(
  lambda = lambda_seq,
  loglike = NA,   # para llenar luego
  like = NA,
  prior = NA,
  post = NA
)

# Evaluamos la verosimilitud para cada valor
# ACÁ DEFINIMOS QUÉ FILAS USAR!
for (i in 1:side) {
  grilla$loglike[i] <- like_fire(grilla$lambda[i], filas = c(2, 10))
}
grilla$like <- exp(grilla$loglike)

# Definimos y evaluamos la previa (gama)
mu <- 10 # media
tau <- 1 # dispersion
a <- 1 / tau; b <- 1 / (tau * mu)

for (i in 1:side) {
  grilla$prior[i] <- dgamma(grilla$lambda[i], a, b)
}

# Posterior
grilla$post <- grilla$like * grilla$prior

# Normalizamos para poder visualizar
g <- grilla
area <- diff(lambda_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla[, v], area)
}

# Elongamos g
cols <- which(names(g) %in% c("like", "prior", "post"))
glong <- pivot_longer(g, all_of(cols), names_to = "variable",
                      values_to = "density")
glong$variable <- factor(
  glong$variable,
  levels = c("like", "prior", "post"),
  labels = c("Verosimilitud", "Previa", "Posterior")
)

# Graficamos
ggplot(glong, aes(x = lambda, y = density, ymax = density, ymin = 0,
                  color = variable, fill = variable)) +
  geom_line() +
  geom_ribbon(alpha = 0.4, color = NA) +
  scale_color_viridis(discrete = T, option = "C", end = 0.6) +
  scale_fill_viridis(discrete = T, option = "C", end = 0.6) +
  theme(legend.title = element_blank()) +
  xlab(expression(lambda)) +
  ylab("Densidad o Verosimilitud") +
  nice_theme()


## y con una previa más informativa:

# Redefinimos y evaluamos la previa (gama)
mu <- 20    # media
tau <- 0.01 # dispersion
a <- 1 / tau; b <- 1 / (tau * mu)

for (i in 1:side) {
  grilla$prior[i] <- dgamma(grilla$lambda[i], a, b)
}

# Posterior
grilla$post <- grilla$like * grilla$prior

# Normalizamos para poder visualizar
g <- grilla
area <- diff(lambda_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla[, v], area)
}

# Elongamos g
cols <- which(names(g) %in% c("like", "prior", "post"))
glong <- pivot_longer(g, all_of(cols), names_to = "variable",
                      values_to = "density")
glong$variable <- factor(
  glong$variable,
  levels = c("like", "prior", "post"),
  labels = c("Verosimilitud", "Previa", "Posterior")
)

# Graficamos
ggplot(glong, aes(x = lambda, y = density, ymax = density, ymin = 0,
                  color = variable, fill = variable)) +
  geom_line() +
  geom_ribbon(alpha = 0.4, color = NA) +
  scale_color_viridis(discrete = T, option = "C", end = 0.6) +
  scale_fill_viridis(discrete = T, option = "C", end = 0.6) +
  theme(legend.title = element_blank()) +
  xlab(expression(lambda)) +
  ylab("Densidad o Verosimilitud") +
  nice_theme()

# 2.5 ---------------------------------------------------------------------
# ¿Cómo se podría caracterizar la distribución posterior de una forma más
# compacta que un gráfico? ¿Cómo se podrían describir resultados
# sobre un parámetro en particular si el modelo tiene más de un parámetro?