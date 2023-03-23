import numpy as np
import matplotlib.pyplot as plt

def fractional_deriv(x, y, alpha, h):
    """
    Implementa el método de Grünwald-Letnikov para la derivada fraccionaria.

    Parámetros:
    x -- arreglo de valores de la variable independiente.
    y -- arreglo de valores de la variable dependiente.
    alpha -- orden de la derivada fraccionaria.
    h -- tamaño del paso.

    Retorna:
    deriv -- arreglo de valores de la derivada fraccionaria.
    """

    n = len(y)
    m = int(alpha)

    # Calcula las diferencias fraccionarias
    frac_diff = np.zeros(n)

    for i in range(m, n):
        sum = 0.0
        for j in range(m+1):
            sign = (-1) ** j
            coeff = np.math.factorial(m) / (np.math.factorial(j) * np.math.factorial(m-j))
            sum += sign * coeff * y[i-j]
        frac_diff[i] = sum / (h ** alpha)

    # Interpola los puntos faltantes
    deriv = np.zeros(n)
    deriv[:m] = frac_diff[m]
    deriv[n-m:] = frac_diff[n-m-1]

    for i in range(m, n-m):
        sum = 0.0
        for j in range(1, m+1):
            sum += (-1) ** j * np.math.factorial(alpha) / (np.math.factorial(j) * np.math.factorial(alpha-j)) * frac_diff[i+j]
        deriv[i] = sum / h ** alpha

    return deriv

    # Definimos una función de ejemplo
def f(x):
    return np.sin(x)

# Definimos los parámetros para el cálculo de la derivada fraccionaria
alpha = 0.5
h = 0.1

# Creamos los arreglos de valores de la variable independiente y dependiente
x = np.arange(0, 2*np.pi, h)
y = f(x)

# Calculamos la derivada fraccionaria de la función
deriv = fractional_deriv(x, y, alpha, h)

# Graficamos la función y su derivada fraccionaria
plt.plot(x, f(x), label="Función original")
plt.plot(x, deriv, label="Derivada fraccionaria")
plt.legend()
plt.show()
