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

#coeficientes binomiales
def binomial_coef(alpha, k):
    c = np.zeros((k+1,1))
    c[0] = 1
    for j in range(1,k+1):
        c[j] = c[j-1] - c[j-1]*(1+alpha)/j
        #c[j] = int(c*10**10)/10**10
        #print(c[j])
    return c

# Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario
def grunwald_letnikov(x,h_alpha,k,alpha,xt):
    # Integración numérica del sistema de Lorenz con Grünwald-Letnikov

    # Iteraciones del método
    sum_x = np.zeros(3)
    c = binomial_coef(alpha,k)

    for i in range(1,k):
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        for j in range(1,i):
            sum_x += c[j]*x[k-j,:]
            # Las variables x,y,z se evaluan con el vector de tiempo
        print(sum_x)
        x[i,:] = lorenz_frac(xt[k-1,:])*h_alpha - sum_x
        #print(x[i,0],x[i,1],x[i,2])
    return x

# Definición del orden fraccionario
alpha = 0.993

# Definición del vector de estado inicial
x0 = np.array([1.0, 1.0, 1.0])

# Definición del intervalo de tiempo y el paso de integración
t0 = 0.0
tf = 10.0
h = 0.001
h_alpha = h**alpha
# Definición del número de puntos en la suma de Grünwald-Letnikov
k = int((tf-t0)/h)

# Inicialización del vector de estado y tiempo
x = np.zeros((int((tf-t0)/h)+1, 3))
t = np.linspace(t0, tf, len(x))
xt = np.zeros_like(x)
xt[:,0] = t
xt[:,1] = t
xt[:,2] = t
#print(lorenz_frac(xt[100,:]))
# Condición inicial
x[0,:] = x0

x = grunwald_letnikov(x,h_alpha,k,alpha,xt)

# Gráfica de la solución
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
