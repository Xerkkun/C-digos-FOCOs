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
from pylab import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
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
    xdot[1] = (rho * x[0]) - x[1] - (x[0] * x[2])
    xdot[2] = (-beta * x[2]) + (x[0] * x[1])
    return xdot

def ho2_system(x):
    # Definición de parámetros del sistema
    b = 0.1
    a = 0.2
    xdot = np.zeros_like(x)
    xdot[0] = x[1] * x[2]
    xdot[1] = x[0] - x[1] - a*x[3]
    xdot[2] = 1 - x[0]*x[0]
    xdot[3] = b * x[1]
    return xdot

#coeficientes binomiales
def binomial_coef(alpha,k,nu):
    c = np.zeros((k+nu,1))
    c[0] = 1
    for j in range(1,k+nu):
        c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)
        #c[j] = int(c*10**10)/10**10
        #print(c[j])
    return c

# Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario
def grunwald_letnikov(x,h,h_alpha,k,alpha,xt,nu,d):
    # Integración numérica del sistema de Lorenz con Grünwald-Letnikov

    # Iteraciones del método
    sum_x = np.zeros(d)
    c = binomial_coef(alpha,k,nu)

    for i in range(1,k+1):

        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        for j in range(nu,i):
            sum_x += c[j]*x[i-j,:]
            # Las variables x,y,z se evaluan con el vector de tiempo
        x[i,:] = ho2_system(x[i-1,:])*h_alpha - sum_x

        if i%1000 == 0 :
            #print(sum_x)
            #print(i,c[k-i+1],x[i-1,:])
            print(xt[i,0],x[i,0],x[i,1],x[i,2],x[i,3])
            #arch.write('%.5f' % t + '\t' + '%.5f' % x + '\t' + '%.5f' % y + '\t' + '%.5f' % z + '\t' + '%.5f' % w + '\n')
        sum_x = np.zeros(d)

    return x

# Definición del orden fraccionario
alpha = 0.995

# Definición del vector de estado inicial
x0 = np.array([0.1, 0.1, 0. , 0.])

# Definición del intervalo de tiempo y el paso de integración
t0 = 0.0
tf = 250.0
h = 0.005
h_alpha = h**alpha

# Longitud de memoria
Lm = 3000

# Número de coeficientes binomiales
m = Lm/h

# Definición del número de puntos en la suma de Grünwald-Letnikov
k = int((tf-t0)/h)

#Principio de memoria corta
if k<=m:
    nu = 1
else:
    nu = k - m

#dimensiones del sistema
d = 4

# Inicialización del vector de estado y tiempo
x = np.zeros((k+1, d))
t = np.linspace(t0, tf, len(x))
xt = np.zeros((k+1,1))
xt[:,0] = t

#print(xt)
# Condición inicial
x[0,:] = x0

x = grunwald_letnikov(x,h,h_alpha,k,alpha,xt,nu,d)
#print(lorenz_frac(x0)*h_alpha)

t, xx, y, z, w = xt[:,0], x[:,0], x[:,1], x[:,2], x[:,3]

subplot(221)
p1,=plot(xx,y,"m",lw=0.3)
xlabel("x")
ylabel("y")
#plt.axis('square')

subplot(222)
p2,=plot(y,z,"m",lw=0.3)
xlabel("y")
ylabel("z")

subplot(223)
p3,=plot(xx,z,"m",lw=0.3)
xlabel("x")
ylabel("z")

subplot(224)
p4,=plot(xx,w,"m",lw=0.3)
xlabel("x")
ylabel("w")

plt.savefig("atractores.pdf",dpi=300,bbox_inches= 'tight')
