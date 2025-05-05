## ## TP 1A: VISUALIZACIÓN DE DATOS ## ##

# Notas ----
# - Abre, reorganiza y visualiza datos
# - Codificación: UTF-8

# Entorno ----
# vaciar entorno de trabajo
rm(list = ls()) ## remover(lista_de_objetos = objetos_existentes) ##

# Paquetes ----
# cargar paquete 'tidyverse'
library(tidyverse) ## cargar_libreria(libreria_a_cargar) ##

# Definiciones ----
# generar cadena de texto (string) con nombre y ubicación del archivo de datos
path_datos <- file.path("datos", "barbera_data_fire_total_climate.csv")

# Abrir datos ----
# cargar datos
d <- read.csv(path_datos) ## d <- leer_csv(ubicación_y_nombre_archivo_de_datos)

# Inspeccionar datos ----
is(d)         # muestra qué tipo de objeto 'd'
str(d)        # muestra estructura de 'd'
summary(d)    # resumen numérico de cada columna de 'd'
View(d)       # abre 'd' como tabla en nueva pestaña
table(d$year) # resume en una tabla la frecuencia de cada valor de la columna 'year' de 'd'

# Reorganizar datos ----
# pasar 'd' a formato largo (long)
d_lng <- d %>%                  # aplica a 'd'
  pivot_longer(                 # pasar a formato long (una variable debajo de la otra, no al lado)
    c(-year, -area_ha, -fires), # excepto a las columnas 'year', 'area_ha', 'fires'
    names_to = "var",           # nombrar 'var' la columna que indicará el nombre de cada variable
    values_to = "valor") %>%    # nombrar 'valor' la columna que contendrá los valores de cada variable
  mutate(lg10_area_ha = log10(area_ha)) # agregar la columna 'lg10_area_ha' que es el log10 de 'area_ha'
str(d_lng)    # muestra estrucutra tiene 'd_lng'

# Visualizar datos ----
## Gráfico simple
ggplot() +                                  # genera un gráfico tipo 'ggplot'
  geom_line(data = d, aes(year, area_ha)) + # con línea de 'area_ha' en función de 'year' de los datos 'd'
  geom_point(data = d, aes(year, area_ha))  # con puntos de 'area_ha' en función de 'year' de los datos 'd'

## Gráfico simple (compactando código y agregando rótulos)
# lo mismo que antes pero más compacto el código y agregando color y rótulos
ggplot(d, aes(year, area_ha)) +     # gráfico 'ggplot' para datos 'd', con 'year' en eje 'x' y 'area_ha' en eje 'y'
  geom_line(color = "blue") +       # con línea color azul
  geom_point() +                    # con puntos (según lo especificado en 'ggplot')
  labs(                             # con rótulos
    x = "Año",                      # rótulo para el eje 'x'
    y = "Área quemada (ha)",        # rótulo para el eje 'y'
    title = "Área quemada por año") # título del gráfico

## Gráfico según variables
# graficar 'area_ha' en función de otras variables
ggplot(d_lng) +                         # gráfico 'ggplot' para datos 'd_lng'
  geom_point(aes(valor, area_ha)) +     # con puntos, con 'valor' en eje 'x' y 'area_ha' en eje 'y'
  scale_y_continuous(trans = "log10") + # transformando valores de 'y' con log base 10
  facet_wrap(                           # mostrando subpaneles (facets)
    vars(var),                          # para cada valor de 'var'
    scales = "free") +                  # con escalas libres entre facets para ambos ejes
  labs(x = "Valor",                     # agregar rótulos de ejes y título de gráfico
       y = "Área quemada (ha)",
       title = "Área quemada según variable")
