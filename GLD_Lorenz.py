#===============================================================================
#   Información
#   Programa que utiliza la definición del operador diferencial de Grünwald-Letnikov
#   para aproximar la solución numérica de un sistema de ecuaciones de orden
#   fraccionario q con valores iniciales.
#   Se utiliza el sistema caótico de Lorenz para calibrar el método numérico.
#===============================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
#===============================================================================
#Funciones
#===============================================================================
#Ecuaciones del sistema de Lorenz
def Fx1(x,y,z):
    a = 10.
    Fx = a*(y-x)
    return Fx

def Fx2(x,y,z):
    c = 28.
    Fy = x*(c-z)-y
    return Fy

def Fx3(x,y,z):
    b = 8./3.
    Fz = x*y-b*z
    return Fz
#===============================================================================
#   Métodos numéricos
#===============================================================================
def GLD(xn,yn,zn,h,q,k):
    #   Grünwald-Letnikov
    # h: ancho de paso
    # q: orden fraccionario
    # k: número de coeficientes binomiales
    h_alpha = h**(q)
    #coeficientes binomiales
    wjx,wjy,wjz = 1,1,1
    sum_x,sum_y,sum_z = 0,0,0
    #revisar la suma por que los productos van en orden diferente para la variable y el coeficiente binomial
    for j in range(0,k):
        sum_x,sum_y,sum_z = sum_x + wjx*Fx1(xn,yn,zn),sum_y + wjy*Fx2(xn,yn,zn), sum_z+wjz*Fx3(xn,yn,zn)
        wjx,wjy,wjz = wjx - (wjx*(1+q))/j,wjy - (wjy*(1+q))/j,wjz - (wjz*(1+q))/j

    x,y,z = sum_x + h_alpha*Fx1(xn,yn,zn),sum_y + h_alpha*Fx2(xn,yn,zn),sum_z + h_alpha*Fx3(xn,yn,zn)
