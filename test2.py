#===============================================================================
#   Información
#   Programa que utiliza la definición del operador diferencial de Grünwald-Letnikov
#   para aproximar la solución numérica de un sistema de ecuaciones de orden
#   fraccionario alpha con valores iniciales.
#   Se utiliza el sistema caótico de Lorenz para calibrar el método numérico.
#===============================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
#===============================================================================
#Funciones
#===============================================================================
def lorenz_frac(x):
    # Definición de parámetros del sistema
    sigma = 10.0
    beta = 8./3.
    rho = 28.0
    xdot = np.zeros_like(x)
    xdot[0] = sigma * (x[1] - x[0])
    xdot[1] = rho * x[0] - x[1] - x[0] * x[2]
    xdot[2] = -beta * x[2] + x[0] * x[1]
    return xdot
# Definición del orden fraccionario
alpha = 0.9929
# Definición del intervalo de tiempo y el paso de integración
t0 = 0.0
tf = 100.0
h = 0.001
h_alpha = h**alpha
# Definición del número de puntos en la suma de Grünwald-Letnikov
k = int((tf-t0)/h)+1

# Inicialización del vector de estado y tiempo
x = np.zeros((int((tf-t0)/h)+1, 3))
t = np.linspace(t0, tf, len(x))
xt = np.zeros_like(x)
xt[:,0] = t
xt[:,1] = t
xt[:,2] = t

sum_x = np.zeros(3)
print(xt)
