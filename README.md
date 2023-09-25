# linear_algebra

Librería de algebra lineal en C++ basada en los videos del canal de YouTube `QuantitativeBytes` https://www.youtube.com/@QuantitativeBytes.

Gracias a QuantitativeBytes por su excelente trabajo en la construcción de esta libreria, la cual me ha
permitido aprender más sobre métodos numéricos para el álgebra y a prácticamente programar en C++ desde cero.


## Funcionalidades

### EIG.h

Funciones para el cálculo del valores y vectores propios.

## LinSolve.h

Funciones para resolver sistemas lineales representados a través de matrices mediante la eliminación Gaussiana.

## LSQ.h

Dado un sistema de equaciones lineales en forma de matriz se calcula la solución por mínimos cuadrados.

## Matrix.h

Classe principal para la gestión de las matrices. Tiene implementadas funciones para las principales operaciones a realizar en matrices como por ejemplo la traspuesta, la inversa, determinante, rango...

## PCA.h

Funciones para realizar análisis de componentes principales.

### QR.h

Funciones para aplicar una descomposición QR a una matriz dada, donde Q es una matriz ortogonal y R triangular superior. El método utilizado es el de Householder.
