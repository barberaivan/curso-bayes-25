## ## TP 2E: MODELO LINEAL CON PARÁMETROS VARIABLES ## ##

# Notas ----
# - Especifica modelo lineal con incertidumbre en variable dependiente
#   y en parámetros
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Definiciones ----
x <- seq(-5, 25, 2)
a <- 2 #ordenada al origen
b <- .7 #pendiente
desvio_media <- 1.2 #desvío estándar de la media del modelo
N <- length(x)
media_a <- 2    #media de 'a'
desvio_a <- .5  #desvío estándar de 'a'
media_b <- .7   #media de 'b'
desvio_b <- .2  #desvío estándar de 'b'

# Simular datos ----
y_media <- a + b * x

y_media2 <- rnorm(1, media_a, desvio_a) + rnorm(1, media_b, desvio_b) * x

y_obs <- rnorm(N, y_media2, desvio_media)

# Organizar en un dataframe los datos simulados ----
d_sim <- cbind.data.frame(
  x = x,
  y_media = y_media,
  y_media2 = y_media2,
  y_obs = y_obs
)

# Visualizar datos simulados ----
ggplot(d_sim) +
  geom_line(aes(x, y_media), color = "grey50") +
  geom_point(aes(x, y_media), color = "grey50") +
  geom_point(aes(x, y_media2),color = "orange") +
  geom_line(aes(x, y_media2),color = "orange") +
  geom_point(aes(x, y_obs), color = "blue") +
  coord_cartesian(ylim = c(-4, 24)) +
  labs(x = "x", y = "y", title = "TP2E: Modelo lineal con parámetros variables")
