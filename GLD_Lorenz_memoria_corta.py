"""   Información
   Programa que utiliza la definición del operador diferencial de Grünwald-Letnikov
   para aproximar la solución numérica de un sistema de ecuaciones de orden
   fraccionario alpha con valores iniciales.
   Se utiliza el sistema caótico de Lorenz para calibrar el método numérico."""
# ===============================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import style
#===============================================================================
#Funciones
#===============================================================================

def lorenz_frac(x):
    """ Ecuaciones del sistema caótico de Lorenz """
    sigma = 10.0
    beta = 8./3.
    rho = 28.0
    xdot = np.zeros_like(x)
    xdot[0] = sigma * (x[1] - x[0])
    xdot[1] = (rho * x[0]) - x[1] - (x[0] * x[2])
    xdot[2] = (-beta * x[2]) + (x[0] * x[1])
    return xdot

def ho2_system(x):
    """ Cálculo de los coeficientes binomiales """
    b = 0.1
    a = 0.2
    xdot = np.zeros_like(x)
    xdot[0] = x[1] * x[2]
    xdot[1] = x[0] - x[1] - a*x[3]
    xdot[2] = 1 - x[0]*x[0]
    xdot[3] = b * x[1]
    return xdot

#coeficientes binomiales
def binomial_coef(alpha,mm,decimal):
    """ Cálculo de los coeficientes binomiales """
    c = np.zeros((mm,1))
    c[0] = 1
    arch = open("lorenz_coeficientes" + ".rnd","w")
    for j in range(1,mm+1):
        c[j] = c[j-1] - (c[j-1]*(1+alpha)/j)
        c[j] = int(c*10**decimal)/10**decimal #para cortar la cantidad de bits de los valores decimales

        arch.write('%.15f' % c[j] + '\n')
    arch.close
    return c

def grunwald_letnikov(x,h,h_alpha,k,alpha,x_t,nu,d,mm,decimal):
    """Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario"""

    # Iteraciones del método
    sum_x = np.zeros(d)
    c = binomial_coef(alpha,mm,decimal)
    arch_1 = open("lorenz_variables" + ".rnd","w") #"wb" para escribir archivos con formato binario
    arch_2 = open("lorenz_sumatoria" + ".rnd","w")

    for i in range(1,k+1):
        # mm = k, mm = m
        # Se calculan las sumas de los coeficientes binomiales en cada iteracion
        if nu == 1:
            for j in range(nu,i):
                sum_x += c[j]*x[k-j,:]
            x[i,:] = lorenz_frac(x[i-1,:])*h_alpha - sum_x
        else:
#============================================================
            for j in range(1,nu):
                sum_x += c[j]*x[k-j,:]
            if i <= m: # Se llenan los valores de x hasta la iteracion m (maximo de memoria)
                x[i,:] = lorenz_frac(x[i-1,:])*h_alpha - sum_x
            # Luego se van desechando los primeros valores de x que se habian calculado con anterioridad para ir guardando los nuevos
            # Si el  numero de iteracion (i) ya excede el valor de la memoria (m)
            # Desplazar todos los valores del vecto x una posicion hacia arriba y dejar vacio el último elemento del vector
            else:
                new_x = lorenz_frac(x[i-1,:])*h_alpha - sum_x
                
#============================================================
        if d==3:
            if i%100 == 0 :
                print(x_t[i,0],x[i,0],x[i,1],x[i,2])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.10f' % x[i,0] + '\t' + '%.10f' % x[i,1] + '\t' + '%.10f' % x[i,2] + '\n')
            arch_2.write('%.15f' % sum_x[0] + '\t' + '%.15f' % sum_x[1] + '\t' + '%.15f' % sum_x[2] + '\n')

        else:
            if i%100 == 0 :
                print(x_t[i,0],x[i,0],x[i,1],x[i,2],x[i,3])
            arch_1.write('%.3f' % x_t[i,0] + '\t' + '%.10f' % x[i,0] + '\t' + '%.10f' % x[i,1] + '\t' + '%.10f' % x[i,2] + '%.10f' % x[i,3] + '\n')
            arch_2.write('%.15f' % sum_x[0] + '\t' + '%.15f' % sum_x[1] + '\t' + '%.15f' % sum_x[2] + '%.15f' % sum_x[3] + '\n')

        sum_x = np.zeros(d)
    arch_1.close
    arch_2.close
    return x

def grafica3d(x,y,z,t):
    """Función para graficar los atractores en el espacio de fases y guardar en pdf"""
    subplot(221)
    p1,=plot(x,y,"m",lw=0.3)
    xlabel("x")
    ylabel("y")
    #plt.axis('square')

    subplot(222)
    p2,=plot(x,z,"m",lw=0.3)
    xlabel("x")
    ylabel("z")

    subplot(223)
    p3,=plot(y,z,"m",lw=0.3)
    xlabel("y")
    ylabel("z")
    plt.savefig("lorenz_atractores.pdf",dpi=300,bbox_inches= 'tight')
    clf()

def grafica4d(x,y,z,w,t):
    """Función para graficar los atractores en el espacio de fases y guardar en pdf"""
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


def main():
    """Integración numérica del sistema de Lorenz con Grünwald-Letnikov"""
    # Definición de la parte decimal
    decimal = 10

    # Definición del orden fraccionario
    alpha = 0.8

    # Definición del vector de estado inicial
    x_0 = np.array([0.1, 0.1, 0.1])

    # Definición del intervalo de tiempo y el paso de integración
    t_0 = 0.0
    t_f = 100.0
    h = 0.01
    h_alpha = h**alpha

    # Longitud de memoria
    Lm = 100

    # Número de coeficientes binomiales
    m = Lm/h #2000

    # Definición del número de puntos en la suma de Grünwald-Letnikov
    k = int((t_f-t_0)/h) #10000

    #Principio de memoria corta
    if k<m:
        nu,mm = 1,k
    else:
        nu,mm = m,m

    #dimensiones del sistema
    d = 3

    # Inicialización del vector de estado y tiempo
    x = np.zeros((k+1, d))
    t = np.linspace(t_0, t_f, len(x))
    x_t = np.zeros((k+1,1))
    x_t[:,0] = t

    # Condición inicial
    x[0,:] = x_0

    x = grunwald_letnikov(x,h,h_alpha,k,alpha,x_t,nu,d,mm,decimal)

    if d==3:
        t, xx, y, z = x_t[:,0], x[:,0], x[:,1], x[:,2]
        grafica3d(xx,y,z,t)
    else:
        t, xx, y, z, w = x_t[:,0], x[:,0], x[:,1], x[:,2], x[:,3]
        grafica4d(xx,y,z,w,t)


if __name__ == '__main__':
    main()
