## ## EXPLORACIÓN MANEJO DE R ## ##

# Notas ----
# - Script para explorar manejo y comprensión de R para el curso
# - Se simulan datos para un modelo lineal y se comparan gráficamente con datos reales
# - Codificación: UTF-8 (si no se leen bien los acentos, usar 'Reopen with Encoding...')

# Entorno ----
rm(list = ls()) # vaciar entorno de trabajo

# Paquetes ----
library(tidyverse) # cargar paquete 'tidyverse'

# Definiciones ----
path_datos <- file.path("datos", "barbera_data_fire_total_climate.csv") # ubicación del csv

x <- seq(-5, 25, 2) # definir una secuencia que va desde -5 a 20, cada 2; es un vector
a <- 2              # crear y asignar valor a 'a' (que será un escalar)
b <- .7             # crear y asignar valor a 'b' (que será un escalar)
veces <- 10         # veces que se repetirá tomar muestras al azar
N <- length(x)      # cantidad de observaciones en 'x'

# Simular datos ----
mu <- a + b * x      # calcular línea recta; 'y' es un vector (porque 'x' es un vector)

d_sim <- NULL                # dataframe para guardar datos simulados
for(i in 1:veces){           # para cada vez (repetición)
  
  y_tmp <- rnorm(N, mu)      # calcular 'y' según una distribución normal
  
  ## Guardar simulación
  d_sim <- rbind.data.frame( # unir filas
    d_sim,
    cbind.data.frame(        # unir columnas
      vez = i,
      x = x,
      mu = mu,
      y = y_tmp
    )
  )
}

# Visualizar datos simulados ----
ggplot(d_sim) +   # graficar en eje 'x' columna 'x' y en eje 'y' columna 'y' de 'd_sim'
  geom_line(aes(x, mu)) +               # con línea para 'y' en función de 'x'
  geom_point(aes(x, y),                 # con puntos de 'y' en función de 'x'
             color = "blue", alpha= .4) # en color azul y transparencia .4

# Comparar con datos de fuego ----
## Abrir datos
d <- read.csv(path_datos) ## d <- leer_csv(ubicación_y_nombre_archivo_de_datos)

## Transformar área quemada
d$log10_area_ha <- log10(d$area_ha) # crear columna 'log10_area_ha' en 'd' y asignar log10 de columna 'area_ha' de 'd'

## Graficar modelo y datos de fuego
ggplot(d_sim, aes(x, mu)) +   # graficar en eje 'x' columna 'x' y en eje 'y' columna 'mu' de 'd_sim'
  geom_line(color = "blue") + # con línea azul para 'y' en función de 'x'
  geom_point(                 # con puntos...
    data = d,                 # de datos 'd'
    aes(temp, log10_area_ha)) # con columna 'temp' en eje 'x' y columna 'log10_area_ha' en eje 'y'
