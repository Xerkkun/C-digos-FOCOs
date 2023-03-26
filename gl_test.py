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
#Ecuaciones del sistema de Lorenz
def f(x):
    xdot = np.zeros_like(x)
    xdot = x**2
    return xdot

#coeficientes binomiales
def binomial_coef(alpha,k,nu):
    c = np.zeros((k+nu,1))
    c[0] = 1
    for j in range(1,k+nu):
        c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)
        #c[j] = int(c*10**10)/10**10
    return c

# Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario
def grunwald_letnikov(x,h_alpha,k,alpha,xt,nu):
    # Integración numérica del sistema de Lorenz con Grünwald-Letnikov

    # Iteraciones del método
    sum_x = 0

    c = binomial_coef(alpha,k,nu)
    for i in range(1,k+1):
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        for j in range(nu,k+1):
            sum_x += c[j]*x[k-j,0]
            # Las variables x,y,z se evaluan con el vector de tiempo
        x[i,0] = f(x[i-1,0])*h_alpha - sum_x
        print(i,x[i-1,0])
        sum_x = 0

    return x

# Definición del orden fraccionario
alpha = 0.9

# Definición del vector de estado inicial
x0 = 1.0

# Definición del intervalo de tiempo y el paso de integración
t0 = 0.0
tf = 1.0
h = 0.01
h_alpha = h**alpha

# Longitud de memoria
Lm = 10

# Número de coeficientes binomiales
m = Lm/h

# Definición del número de puntos en la suma de Grünwald-Letnikov
k = int((tf-t0)/h)

#Principio de memoria corta
if k<=m:
    nu = 1
else:
    nu = k - m

# Inicialización del vector de estado y tiempo
x = np.zeros((k+1, 1))
t = np.linspace(t0, tf, len(x))
xt = np.zeros_like(x)
xt[:,0] = t

# Condición inicial
x[0,0] = x0

x = grunwald_letnikov(x,h_alpha,k,alpha,xt,nu)

# Gráfica de la solución
# fig = plt.figure(figsize=(10, 8))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x[:,0], x[:,1], x[:,2])
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.show()
