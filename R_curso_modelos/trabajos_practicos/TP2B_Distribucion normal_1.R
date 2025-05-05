## ## TP 2B: DISTRIBUCIÓN NORMAL ## ##

# Notas ----
# - Simula de datos según distribución normal
# - Calcula pdf y cdf de distribución normal
# - Compara distribución normal con datos de fuego
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)

# Generar muestras de distribución normal ----
## Definiciones ----
path.datos <- file.path("datos", "barbera_data_fire_total_climate.csv")

x <- seq(-5, 10, .1) #valores de 'x'
media <- 3 #media de la distribución normal
desvio <- 2 #desvío estándar de la distribución normal
N <- 100 #cantidad de muestras a simular de la normal

# Abrir datos ----
d <- read.csv(path.datos)

## Simular datos ----
muestras <- rnorm(N, media, desvio) #obtiene 'N' muestras de una normal con media 'media' y desvío 'desvio'

## Visualizar simulación ----
hist(muestras) #histograma

# Explorar pdf y cdf de normal ----
## Simular datos ----
pdf <- dnorm(x, media, desvio) #función de densidad de probabilidad
cdf <- pnorm(x, media, desvio) #función de distribución acumulada
x_q <- qnorm(cdf, media, desvio) #cuantiles

## Organizar en un dataframe los datos simulados----
d_sim <- cbind.data.frame(
  x = x,
  cdf = cdf,
  pdf = pdf,
  x_q = x_q
)

## Visualizar simulación ----
ggplot(d_sim) +
  geom_line(aes(x, cdf), color = "blue") +
  geom_line(aes(x, pdf), color = "red") +
  geom_vline(aes(xintercept = media), linetype = 2) +
  labs(x = "x", y = "cdf, pdf", title = "Distribución normal")

# Calcular normal para la precipitación en los datos de fuego ----
## Calcular parámetros de la normal
media_pp <- mean(d$pp) #media de los datos de precipitación
desvio_pp <- sd(d$pp) #desvío estándar de los datos de precipitación

## Calcular pdf de la normal
x_pp <- seq(0, 350, 5) #secuencia de valores posibles de precipitación
pdf_pp <- dnorm(x_pp, media_pp, desvio_pp)
d_pp <- cbind.data.frame(
  x = x_pp,
  pdf = pdf_pp
)

## Graficar datos de precipitación y distribución normal
ggplot() +
  geom_density(data = d, aes(pp)) +
  geom_line(data = d_pp, aes(x, pdf), color = "blue")

# Calcular normal para el log10 del área quemada en los datos de fuego ----
## Calcular parámetros de la normal
media_area <- mean(log10(d$area_ha))
desvio_area <- sd(log10(d$area))

## Calcular pdf de la normal
x_area <- seq(0, 5, .1)
pdf_area <- dnorm(x_area, media_area, desvio_area)
d_area <- cbind.data.frame(
  x = x_area,
  pdf = pdf_area
)

## Graficar datos de precipitación y distribución normal
ggplot() +
  geom_density(data = d, aes(log10(area_ha))) +
  geom_line(data = d_area, aes(x, pdf), color = "blue")
