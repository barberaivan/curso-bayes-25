## ## TP 2D: ERROR EN MODELO LINEAL - PARTE 2 ## ##

# Notas ----
# - Especifica modelo lineal con incertidumbre en la variable dependiente
# - Simula datos muchas veces para representar incertidumbre
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Paquetes ----
library(tidyverse)
library(ggridges)

# Definiciones ----
x <- seq(-5, 25, 2)
# x <- c(-5, -4, -3, -1, seq(9, 25, 2)) #para probar después
a <- 2 #ordenada al origen
b <- .7 #pendiente
media_error <- 0 #media del error del modelo
desvio_error <- 1.2 #desvío estándar del error del modelo
N <- length(x) #cantidad de valores de 'x' considerados
veces <- 100 #cantidad de veces que se repite la simulación

# Simular datos ----
y_media <- a + b * x

y_obs <- rnorm(N * veces, y_media, desvio_error)

# Organizar en un dataframe los datos simulados ----
d_sim <- cbind.data.frame(
  x = x,
  y_media = y_media,
  y_obs = y_obs
)

# Visualizar datos simulados ----
p <- ggplot(d_sim) +
  geom_point(aes(x, y_obs), color = "orange2", alpha = .1) +
  geom_vridgeline(aes(x, y_obs, group = x, width = after_stat(density)),
                  stat = "ydensity", alpha = 0.2, scale = 4, trim = FALSE,
                  fill = "orange2", color = "orange2") +
  labs(x = "x", y = "y", title = "TP2D: Modelo lineal con error")
p

## Resumir datos simulados
rsm_d_sim <- d_sim %>%            #resume 'd_sim'
  group_by(x) %>%                 #para cada valor de 'x'
  summarise(                      #calcula...
    mediana = median(y_obs),      #mediana de 'y_obs'
    q_05 = quantile(y_obs, .05),  #cuantil del 5% de 'y_obs'
    q_95 = quantile(y_obs, .95)   #cuantil del 95% de 'y_obs'
  )

p2 <- p +
  geom_ribbon(data = rsm_d_sim, aes(x, ymin = q_05, ymax = q_95),
              fill = "orange", color = NA, alpha = .4) +
  geom_line(data = rsm_d_sim, aes(x, mediana), color = "orange")
p2

## Distribución univariada de 'y_obs'
hist(d_sim$y_obs) #histograma
plot(density(d_sim$y_obs)) #densidad (una forma de suavizar el histograma)
