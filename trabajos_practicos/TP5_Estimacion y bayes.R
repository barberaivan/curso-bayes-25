## ## TP 5: Estimación de parámetros e inferencia bayesiana ## ##

# [Codificación: UTF-8]

# Paquetes ----------------------------------------------------------------

library(tidyverse) # ggplot y pivot_longer
library(viridis)   # colores
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

# 01 ----------------------------------------------------------------------
# Elegí un set de datos y definí un modelo sencillo (1 o 2 parámetros).
# Puede ser estadístico o no (la función de costo no necesariamente tiene que
# ser una -log verosimilitud).

# y_i ~ Poisson(lambda_i)
# lambda_i = exp(alpha + beta * x)

# 02 ----------------------------------------------------------------------
# Evaluá qué restricciones debería imponerse sobre los parámetros para que el
# modelo sea sensato.

# beta > 0, y no muy grande

# 03 ----------------------------------------------------------------------
# Graficá la función de verosimilitud del modelo (o la función de costo) y
# encontrá el máximo utilizando el método de la grilla. Si el modelo tiene
# D > 2 parámetros, se deberán fijar D - 2 para poder graficar. Graficá la
# predicción en función de alguna variable predictora (si es que la hay)
# utilizando el valor de la grilla con mayor verosimilitud o menor costo.

like_fire <- function(alpha, beta, log = T, filas = 1:nrow(datos)) {
  # "like" por likelihood (sería raro decirle "vero");
  # "fire" porque es la función para estos datos de fuego, no cualquier
  # verosimilitud.

  # Calculamos la media, lambda:
  lambda <- exp(alpha + beta * datos$fwi[filas])

  # usando lambda evaluamos la verosimilitud de cada observación,
  # en escala log:
  like_pointwise <- dpois(datos$fires[filas], lambda, log = T)

  # la log-verosimilitud conjunta es la suma de las log-verosimilitudes por
  # observación:
  like <- sum(like_pointwise)

  if (log) return(like)
  else return(exp(like))
}

# Creamos grilla de valores de alpha y beta.
side <- 100 # cantidad de valores por lado
grilla <- expand.grid(
  alpha = seq(log(0.01), log(30), length.out = side),
  beta = seq(0, 0.3, length.out = side)
)

# Evaluamos la verosimilitud para cada valor
grilla$loglike <- NA
size <- side ^ 2 # cantidad de valores en la grilla
for (i in 1:size) {
  grilla$loglike[i] <- like_fire(grilla$alpha[i], grilla$beta[i])
}
grilla$like <- exp(grilla$loglike)

# Elegimos el punto con mayor verosimilitud
row <- which.max(grilla$loglike)
alpha <- grilla$alpha[row]
beta <- grilla$beta[row]

# Grilla
p1 <- ggplot(grilla, aes(x = alpha, y = beta, fill = like)) +
  geom_tile() +
  geom_point(data = grilla[row, , drop = F]) +
  scale_fill_viridis_c(option = "viridis", name = "Verosimilitud",
                       label = scientific_format(digits = 2)) +
  labs(x = expression(alpha), y = expression(beta)) +
  theme_minimal(base_size = 14)

# Plot de puntos + curva
p2 <- ggplot(datos, aes(x = fwi, y = fires)) +
  geom_point(color = rgb(0, 0, 0, 0.8), size = 2) +
  stat_function(fun = function(x) exp(alpha + beta * x),
                color = "blue", linewidth = 1) +
  ylab("y") + xlab("FWI") +
  nice_theme()

# Unimos
(p1 | p2)

# 04 ----------------------------------------------------------------------
# Compará la estimación basada en la grilla con la obtenida usando optim.

like_fire_opt <- function(x) {
  alpha <- x[1]
  beta <- x[2]
  return(-like_fire(alpha, beta))
}

opt <- optim(
  par = c(0, 0.1),    # vector inicial de parámetros, donde comienza la búsqueda
  fn = like_fire_opt, # función a optimizar
  method = "BFGS"
)

opt$par
grilla[row, c("alpha", "beta")]

# 05 ----------------------------------------------------------------------
# ¿Existen otros valores de los parámetros que también presentan un
# ajuste razonable a los datos, aunque sea subóptimo? ¿Cómo se podría
# escoger un set de valores razonables?

# 06 ----------------------------------------------------------------------
# Repetí los puntos 3 y 5 pero utilizando el 20 % y el 50 % de los
# datos disponibles, y compará los resultados.

# Ejemplo con N = 20 %

filas20 <- sample(1:nrow(datos), size = round(0.2 * nrow(datos)))

# Evaluamos la verosimilitud para cada valor
grilla$loglike20 <- NA
for (i in 1:size) {
  grilla$loglike20[i] <- like_fire(grilla$alpha[i], grilla$beta[i],
                                   filas = filas20)
}
grilla$like20 <- exp(grilla$loglike20)

row20 <- which.max(grilla$loglike20)
alpha20 <- grilla$alpha[row20]
beta20 <- grilla$beta[row20]

# Grilla
p1 <- ggplot(grilla, aes(x = alpha, y = beta, fill = like20)) +
  geom_tile() +
  geom_point(data = grilla[row20, , drop = F]) +
  scale_fill_viridis_c(option = "viridis", name = "Verosimilitud",
                       label = scientific_format(digits = 2)) +
  labs(x = expression(alpha), y = expression(beta)) +
  theme_minimal(base_size = 14)

# Plot de puntos + curva
p2 <- ggplot(datos[filas20, ], aes(x = fwi, y = fires)) +
  geom_point(color = rgb(0, 0, 0, 0.8), size = 2) +
  stat_function(fun = function(x) exp(alpha20 + beta20 * x),
                color = "blue", linewidth = 1) +
  ylab("y") + xlab("FWI") +
  nice_theme()

# Unimos usando patchwork
(p1 | p2)



# 07 ----------------------------------------------------------------------
# Utilizando el mismo set de datos y modelo escogido previamente (ejercicio 1),
# definí previas para los parámetros, convirtiéndolo en un modelo bayesiano.

# y_i ~ Poisson(lambda_i)
# lambda_i = exp(alpha + beta * x)
# alpha ~ Normal(0, 1)
# beta ~ Normal(0, 0.5)

# 08 ----------------------------------------------------------------------
# Graficá la función de verosimilitud, la densidad previa y la densidad
# posterior.

# Creamos grilla de valores de alpha y beta.
side <- 100 # cantidad de valores por lado

alpha_seq <- seq(log(0.01), log(30), length.out = side)
beta_seq <- seq(-0.5, 0.5, length.out = side)

grilla <- expand.grid(
  alpha = alpha_seq,
  beta = beta_seq,
  loglike = NA,   # para llenar luego
  like = NA,
  prior = NA,
  post = NA
)

# Evaluamos la verosimilitud para cada valor
size <- side ^ 2 # cantidad de valores en la grilla
for (i in 1:size) {
  grilla$loglike[i] <- like_fire(grilla$alpha[i], grilla$beta[i])
}
grilla$like <- exp(grilla$loglike)

# Definimos y evaluamos la previa
mu_a <- 0; sigma_a <- 1
mu_b <- 0; sigma_b <- 0.5

for (i in 1:size) {
  grilla$prior[i] <-
    dnorm(grilla$alpha[i], mu_a, sigma_a) *
    dnorm(grilla$beta[i], mu_b, sigma_b)
}

# Posterior
grilla$post <- grilla$like * grilla$prior

# Normalizamos para poder visualizar
g <- grilla
a <- diff(alpha_seq)[1] * diff(beta_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla[, v], a)
}

# Graficamos por separado, para apreciar las distintas escalas
p1 <-
  ggplot(g, aes(x = alpha, y = beta, fill = like)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Verosimilitud") +
  nice_theme()

p2 <-
  ggplot(g, aes(x = alpha, y = beta, fill = prior)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Previa") +
  nice_theme()

p3 <-
  ggplot(g, aes(x = alpha, y = beta, fill = post)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Posterior") +
  nice_theme()

# Unimos usando patchwork
(p1 | p2 | p3) + plot_layout(nrow = 2)

# 09 ----------------------------------------------------------------------
# Explorá distintas previas, ya sea variando los parámetros de una
# misma distribución o utilizando otra distribución. Intentá crear
# previas que influyan en el resultado y otras que no.

grilla2 <- grilla

# Redefinimos y evaluamos la previa
mu_a <- 1; sigma_a <- 0.1
mu_b <- -0.1; sigma_b <- 0.05

for (i in 1:size) {
  grilla2$prior[i] <-
    dnorm(grilla2$alpha[i], mu_a, sigma_a) *
    dnorm(grilla2$beta[i], mu_b, sigma_b)
}

# Posterior
grilla2$post <- grilla2$like * grilla2$prior

# Normalizamos para poder visualizar
g2 <- grilla2
a <- diff(alpha_seq)[1] * diff(beta_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g2[, v] <- normalize_dens(grilla2[, v], a)
}

# Graficamos por separado, para apreciar las distintas escalas
p1 <-
  ggplot(g2, aes(x = alpha, y = beta, fill = like)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Verosimilitud") +
  nice_theme()

p2 <-
  ggplot(g2, aes(x = alpha, y = beta, fill = prior)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Previa") +
  nice_theme()

p3 <-
  ggplot(g2, aes(x = alpha, y = beta, fill = post)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Posterior") +
  nice_theme()

# Unimos usando patchwork
(p1 | p2 | p3) + plot_layout(nrow = 2)

## Para probar con otra distribución, editar estas líneas:

# Redefinimos y evaluamos la previa
mu_a <- 1; sigma_a <- 0.1
mu_b <- 0.3; sigma2_b <- 0.1 # beta tendrá una previa gamma

# obtenemos a y b, que son los parámetros que requiere dgamma
a <- mu_b ^ 2 / sigma2_b # sigma2_b es la varianza de la previa para b
b <- mu_b / sigma2_b

for (i in 1:size) {
  grilla2$prior[i] <-
    dnorm(grilla2$alpha[i], mu_a, sigma_a) * # normal
    dgamma(grilla2$beta[i], a, b)            # gamma
}

# (El resto del código, igual)

# 10 ----------------------------------------------------------------------
# Repetí la comparación entre posteriores pero utilizando únicamente 2 o 3
# observaciones. Las previas que no influían en el resultado usando todos los
# datos, ¿siguen sin influir?

# Creamos grilla de valores de alpha y beta.
side <- 100 # cantidad de valores por lado

alpha_seq <- seq(log(0.01), log(30), length.out = side)
beta_seq <- seq(-0.5, 0.5, length.out = side)

grilla <- expand.grid(
  alpha = alpha_seq,
  beta = beta_seq,
  loglike = NA,   # para llenar luego
  like = NA,
  prior = NA,
  post = NA
)

# Evaluamos la verosimilitud para cada valor
size <- side ^ 2 # cantidad de valores en la grilla
for (i in 1:size) {
  grilla$loglike[i] <- like_fire(grilla$alpha[i], grilla$beta[i],
                                 filas = filas20)
}
grilla$like <- exp(grilla$loglike)

# Definimos y evaluamos la previa
mu_a <- 1; sigma_a <- 0.1
mu_b <- 0.3; sigma2_b <- 0.1 # beta tendrá una previa gamma

# obtenemos a y b, que son los parámetros que requiere dgamma
a <- mu_b ^ 2 / sigma2_b # sigma2_b es la varianza de la previa para b
b <- mu_b / sigma2_b

for (i in 1:size) {
  grilla$prior[i] <-
    dnorm(grilla$alpha[i], mu_a, sigma_a) * # normal
    dgamma(grilla$beta[i], a, b)            # gamma
}

# Posterior
grilla$post <- grilla$like * grilla$prior

# Normalizamos para poder visualizar
g <- grilla
a <- diff(alpha_seq)[1] * diff(beta_seq)[1] # para normalizar
for (v in c("like", "prior", "post")) {
  g[, v] <- normalize_dens(grilla[, v], a)
}

# Graficamos por separado, para apreciar las distintas escalas
p1 <-
  ggplot(g, aes(x = alpha, y = beta, fill = like)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Verosimilitud") +
  nice_theme()

p2 <-
  ggplot(g, aes(x = alpha, y = beta, fill = prior)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Previa") +
  nice_theme()

p3 <-
  ggplot(g, aes(x = alpha, y = beta, fill = post)) +
  geom_tile() +
  scale_fill_viridis() +
  theme(legend.title = element_blank()) +
  xlab(expression(alpha)) +
  ylab(expression(beta)) +
  ggtitle("Posterior") +
  nice_theme()

# Unimos usando patchwork
(p1 | p2 | p3) + plot_layout(nrow = 2)

# 11 ----------------------------------------------------------------------
# ¿Cómo podrías caracterizar la distribución posterior de una forma más
# compacta que un gráfico? ¿Cómo lo harías si quisieras describir resultados
# sobre un parámetro en particular si el modelo tiene más de un parámetro?