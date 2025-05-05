## ## TP 2C: ERROR EN MODELO LINEAL ## ##

# Notas ----
# - Especifica modelo lineal con incertidumbre ("error") en la variable dependiente
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Librerías ----
library(tidyverse)
library(ggridges)

# Definiciones ----
x <- seq(-5, 25, 2)
a <- 2 #ordenada al origen
b <- .7 #pendiente
media_error <- 0 #media del error del modelo
desvio_error <- 1.2 #desvío estándar del error del modelo
N <- length(x) #cantidad de valores de 'x' considerados

# Simular datos ----
y_media <- a + b * x

error <- rnorm(N, media_error, desvio_error)
y_obs <- a + b * x + error

# y_obs <- rnorm(N, y_media, desvio.error) #expresión alternativa

# Organizar en un dataframe los datos simulados ----
d_sim <- cbind.data.frame(
  x = x,
  y_media = y_media,
  y_obs = y_obs
)

# Visualizar datos simulados ----
ggplot(d_sim) +
  geom_line(aes(x, y_media), color = "grey50") +
  geom_point(aes(x, y_media), color = "grey50") +
  geom_point(aes(x, y_obs), color = "orange") +
  labs(x = "x", y = "y", title = "TP2C: Modelo lineal con error")
