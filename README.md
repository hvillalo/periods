periods: Extracción y modelación de periodicidades en series de tiempo
regulares
================

## Instalación

El paquete **periods** es la implementación en R del método presentado
en González-Rodríguez et al. (2015). Se puede instalar desde github con
ayuda de **devtools**, el cual a su vez se instala de la manera habitual
en caso de no estar ya disponible.

``` r
install.packages("devtools") # si no está ya instalado

library(devtools)
install_github("hvillalo/periods")
```

## Ejemplos de uso

### Serie simulada

Se generó una serie de tiempo (n = 220) con cuatro componentes
harmónicos definidos por: periodos = 25, 10, 16 y 73; amplitudes 0 40,
20, 10 y 5; fases = 2, 5, 1 y 0; media = 0; sin tendencia lineal y 10 %
de ruido aleatorio.

``` r
library(periods)

## Serie simulada ----
data(sim)

# Plot
plot(sim, type = "l")
```

![](README_files/figure-commonmark/unnamed-chunk-2-1.png)
