## ## TP 1B: MODELO LINEAL SIMPLE ## ##

# Notas ----
# - Script para especificar un modelo lineal, simular datos y graficarlos
# - El modelo es: y = a + bx
# - Compara datos simulados con datos de área quemada en función de la temperatura
# - Codificación: UTF-8

# Entorno ----
rm(list = ls()) # vaciar entorno de trabajo

# Paquetes ----
library(tidyverse) # cargar paquete 'tidyverse'

# Definiciones ----
path_datos <- file.path("datos", "barbera_data_fire_total_climate.csv") # ubicación del csv

x <- seq(-5, 25, 2) # definir una secuencia que va desde -5 a 20, cada 2; es un vector
a <- 2              # crear y asignar valor a 'a' (que será un escalar)
b <- .7             # crear y asignar valor a 'b' (que será un escalar)

# Simular datos ----
y <- a + b * x      # calcular línea recta; 'y' es un vector (porque 'x' es un vector)

# Organizar en un dataframe los datos simulados ----
d_sim <- cbind.data.frame( # generar un dataframe al unir columnas (cbind)
  x = x,                   # la primera columna se llama 'x' y contiene los valores de 'x'
  y = y                    # la segunda columna se llama 'x' y contiene los valores de 'x'
)

# Visualizar datos simulados ----
ggplot(d_sim, aes(x, y)) +   # graficar en eje 'x' columna 'x' y en eje 'y' columna 'y' de 'd_sim'
  geom_line() +              # con línea para 'y' en función de 'x'
  geom_point(color = "blue") # con puntos en color azul para 'y' en función de 'x'

# Comparar con datos de fuego ----
## Abrir datos
d <- read.csv(path_datos) ## d <- leer_csv(ubicación_y_nombre_archivo_de_datos)

## Transformar área quemada
d$log10_area_ha <- log10(d$area_ha) # crear columna 'log10_area_ha' en 'd' y asignar log10 de columna 'area_ha' de 'd'

## Graficar modelo y datos de fuego
ggplot(d_sim, aes(x, y)) +    # graficar en eje 'x' columna 'x' y en eje 'y' columna 'y' de 'd_sim'
  geom_line(color = "blue") + # con línea para 'y' en función de 'x'
  geom_point(                 # con puntos...
    data = d,                 # de datos 'd'
    aes(temp, log10_area_ha)) # con columna 'temp' en eje 'x' y columna 'log10_area_ha' en eje 'y'
