## ## TP 2A: DISTRIBUCIÓN BINOMIAL ## ##

# Notas ----
# - Simulación de datos según distribución binomial
# - Codificación: UTF-8

# Entorno ----
rm(list = ls())

# Explorar función 'sample' ----
## Definiciones
elementos1 <- c(0, 1) #elementos que se sortearán

## Simular datos
muestra <- sample(elementos1, 10, replace=T) #sortea 10 veces elementos de 'elementos1' con reposición

## Explorar simulación
print(muestra)
table(muestra)
mean(muestra)

# Calcular distribución binomial ----
## Definiciones
elementos2 <- c(0, 1) #elementos que se sortearán
veces <- 1000 #veces que se sortearán elementos
muestras <- integer(veces) #vector de enteros, con longitud 'veces'

## Simular datos
for (i_vez in 1:veces) { #para cada vez
  muestras[i_vez] <- sum(sample(elementos2, 10, replace=T)) #sortea 10 veces elementos de 'elementos2' y suma
}

## Explorar simulación
table(muestras)
mean(muestras)
sd(muestras)

## Visualizar datos
hist(muestras) #histograma
