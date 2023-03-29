import numpy as np
import math
import random
#===============================================================================
#Funciones
#===============================================================================
#Ecuaciones del oscilador HO1
def Fx1(x,y,z,w):
    a,b = 4., 6.
    Fx = (a*x) - (b*y*z) - 10.
    return Fx

def Fx2(x,y,z,w):
    c,k = 10.,2.5
    Fy = (-c*y) + (x*z) + (k*w)
    return Fy

def Fx3(x,y,z,w):
    d,e = 5.,2.
    Fz = (-d*z) + (e*x*y)
    return Fz

def Fx4(x,y,z,w):
    f = 0.05
    Fw = f*(x+z)
    return Fw
#===============================================================================
#Métodos numéricos
def forward_euler(xn,yn,zn,wn,h):
    x,y,z,w = xn+h*Fx1(xn,yn,zn,wn),yn+h*Fx2(xn,yn,zn,wn),zn+h*Fx3(xn,yn,zn,wn),wn+h*Fx4(xn,yn,zn,wn)
    return x,y,z,w
#-------------------------------------------------------------------------------
def backward_euler(xn,yn,zn,wn,h):
    xn1,yn1,zn1,wn1 = forward_euler(xn,yn,zn,wn,h)
    x,y,z,w = xn+h*Fx1(xn1,yn1,zn1,wn1),yn+h*Fx2(xn1,yn1,zn1,wn1),zn+h*Fx3(xn1,yn1,zn1,wn1),wn+h*Fx4(xn1,yn1,zn1,wn1)
    return x,y,z,w
#-------------------------------------------------------------------------------
def runge_kutta4(xn,yn,zn,wn,h):
    k1 = h*Fx1(xn,yn,zn,wn)
    k2 = Fx1(xn + h*0.5*k1,yn + h*0.5*k1,zn + h*0.5*k1,wn + h*0.5*k1)
    k3 = Fx1(xn + h*0.5*k2,yn + h*0.5*k2,zn + h*0.5*k2,wn + h*0.5*k2)
    k4 = Fx1(xn + h*k3,yn + h*k3,zn + h*k3,wn + h*k3)
    x = xn + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)

    k1 = h*Fx2(xn,yn,zn,wn)
    k2 = Fx2(xn + h*0.5*k1,yn + h*0.5*k1,zn + h*0.5*k1,wn + h*0.5*k1)
    k3 = Fx2(xn + h*0.5*k2,yn + h*0.5*k2,zn + h*0.5*k2,wn + h*0.5*k2)
    k4 = Fx2(xn + h*k3,yn + h*k3,zn + h*k3,wn + h*k3)
    y = yn + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)

    k1 = h*Fx3(xn,yn,zn,wn)
    k2 = Fx3(xn + h*0.5*k1,yn + h*0.5*k1,zn + h*0.5*k1,wn + h*0.5*k1)
    k3 = Fx3(xn + h*0.5*k2,yn + h*0.5*k2,zn + h*0.5*k2,wn + h*0.5*k2)
    k4 = Fx3(xn + h*k3,yn + h*k3,zn + h*k3,wn + h*k3)
    z = zn + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)

    k1 = h*Fx4(xn,yn,zn,wn)
    k2 = Fx4(xn + h*0.5*k1,yn + h*0.5*k1,zn + h*0.5*k1,wn + h*0.5*k1)
    k3 = Fx4(xn + h*0.5*k2,yn + h*0.5*k2,zn + h*0.5*k2,wn + h*0.5*k2)
    k4 = Fx4(xn + h*k3,yn + h*k3,zn + h*k3,wn + h*k3)
    w = wn + (h/6.)*(k1 + 2*k2 + 2*k3 + k4)

    return x,y,z,w
#-------------------------------------------------------------------------------
def adams_bashforth6(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,
                        xn3,yn3,zn3,wn3,xn4,yn4,zn4,wn4,xn5,yn5,zn5,wn5,h):

    x = (xn + (h/1440.) * (4277*Fx1(xn,yn,zn,wn) - 7923*Fx1(xn1,yn1,zn1,wn1) +
    9982*Fx1(xn2,yn2,zn2,wn2)-7298*Fx1(xn3,yn3,zn3,wn3) + 2877*Fx1(xn4,yn4,zn4,wn4) -
    475*Fx1(xn5,yn5,zn5,wn5)))

    y = (yn + (h/1440.) * (4277*Fx2(xn,yn,zn,wn) - 7923*Fx2(xn1,yn1,zn1,wn1) +
    9982*Fx2(xn2,yn2,zn2,wn2)-7298*Fx2(xn3,yn3,zn3,wn3) + 2877*Fx2(xn4,yn4,zn4,wn4) -
    475*Fx2(xn5,yn5,zn5,wn5)))

    z = (zn + (h/1440.) * (4277*Fx3(xn,yn,zn,wn) - 7923*Fx3(xn1,yn1,zn1,wn1) +
    9982*Fx3(xn2,yn2,zn2,wn2)-7298*Fx3(xn3,yn3,zn3,wn3) + 2877*Fx3(xn4,yn4,zn4,wn4) -
    475*Fx3(xn5,yn5,zn5,wn5)))

    w = (wn + (h/1440.) * (4277*Fx4(xn,yn,zn,wn) - 7923*Fx4(xn1,yn1,zn1,wn1) +
    9982*Fx4(xn2,yn2,zn2,wn2)-7298*Fx4(xn3,yn3,zn3,wn3) + 2877*Fx4(xn4,yn4,zn4,wn4) -
    475*Fx4(xn5,yn5,zn5,wn5)))

    return x,y,z,w
#-------------------------------------------------------------------------------
def adams_moulton4(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,h):

    x1,y1,z1,w1 = forward_euler(xn,yn,zn,wn,h)

    x = xn + (h/24.) * (9*Fx1(x1,y1,z1,w1) + 19*Fx1(xn,yn,zn,wn) -
    5*Fx1(xn1,yn1,zn1,wn1) + Fx1(xn2,yn2,zn2,wn2))

    y = yn + (h/24.) * (9*Fx2(x1,y1,z1,w1) + 19*Fx2(xn,yn,zn,wn) -
    5*Fx2(xn1,yn1,zn1,wn1) + Fx2(xn2,yn2,zn2,wn2))

    z = zn + (h/24.) * (9*Fx3(x1,y1,z1,w1) + 19*Fx3(xn,yn,zn,wn) -
    5*Fx3(xn1,yn1,zn1,wn1) + Fx3(xn2,yn2,zn2,wn2))

    w = wn + (h/24.) * (9*Fx4(x1,y1,z1,w1) + 19*Fx4(xn,yn,zn,wn) -
    5*Fx4(xn1,yn1,zn1,wn1) + Fx4(xn2,yn2,zn2,wn2))

    return x,y,z,w
#===============================================================================
def gear4(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,xn3,yn3,zn3,wn3,h):
    x1,y1,z1,w1 = forward_euler(xn,yn,zn,wn,h)

    x = (48./25.)*xn - (36./25.)*xn1 + (16./25.)*xn2 - (3./25.)*xn3 + (12./25.)*h*Fx1(x1,y1,z1,w1)
    y = (48./25.)*yn - (36./25.)*yn1 + (16./25.)*yn2 - (3./25.)*yn3 + (12./25.)*h*Fx2(x1,y1,z1,w1)
    z = (48./25.)*zn - (36./25.)*zn1 + (16./25.)*zn2 - (3./25.)*zn3 + (12./25.)*h*Fx3(x1,y1,z1,w1)
    w = (48./25.)*wn - (36./25.)*wn1 + (16./25.)*wn2 - (3./25.)*wn3 + (12./25.)*h*Fx4(x1,y1,z1,w1)

    return x,y,z,w
#===============================================================================
#Inicio

#parámetros de entrada
# n: pasos por corridas (1e6)
# met: método numérico de resolución (FE,BE,RK4,AB6,AM4,G4)

nn,met = input('Parámetros de entrada: ').split()
n = int(float(nn))

c = (50,30,30,2)

xo = (random.random()-0.5)*c[0]
yo = (random.random()-0.5)*c[1]
zo = (random.random()-0.5)*c[2]
wo = (random.random()-0.5)*c[3]

print("Número de pasos:", n)
print("Método numérico: ", met)

hh = (0.001,0.01,0.001,0.001,0.005,0.005) #Ancho de paso para cada método
arch = open(met + ".rnd","w") #"wb" para escribir archivos con formato binario
i,t = 0,0

xn1,yn1,zn1,wn1 = 0,0,0,0
xn2,yn2,zn2,wn2 = 0,0,0,0
xn3,yn3,zn3,wn3 = 0,0,0,0
xn4,yn4,zn4,wn4 = 0,0,0,0
xn5,yn5,zn5,wn5 = 0,0,0,0

x,y,z,w = xo,yo,zo,wo #condiciones iniciales

for i in range(0,n):
    xn,yn,zn,wn = x,y,z,w

    if met == 'FE':
        h = hh[0]
        x,y,z,w = forward_euler(xn,yn,zn,wn,h)
    elif met == 'BE':
        h = hh[1]
        x,y,z,w = backward_euler(xn,yn,zn,wn,h)
    elif met == 'RK4':
        h = hh[2]
        x,y,z,w = runge_kutta4(xn,yn,zn,wn,h)
    elif met == 'AB6':
        h = hh[3]
        x,y,z,w = adams_bashforth6(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,
                                xn3,yn3,zn3,wn3,xn4,yn4,zn4,wn4,xn5,yn5,zn5,wn5,h)
    elif met == 'AM4':
        h = hh[4]
        x,y,z,w = adams_moulton4(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,h)
    elif met == 'G4':
        h = hh[5]
        x,y,z,w = gear4(xn,yn,zn,wn,xn1,yn1,zn1,wn1,xn2,yn2,zn2,wn2,xn3,yn3,zn3,wn3,h)

    xn5,yn5,zn5,wn5 = xn4,yn4,zn4,wn4
    xn4,yn4,zn4,wn4 = xn3,yn3,zn3,wn3
    xn3,yn3,zn3,wn3 = xn2,yn2,zn2,wn2
    xn2,yn2,zn2,wn2 = xn1,yn1,zn1,wn1
    xn1,yn1,zn1,wn1 = xn,yn,zn,wn

    arch.write('%.5f' % t + '\t' + '%.5f' % x + '\t' + '%.5f' % y + '\t' + '%.5f' % z + '\t' + '%.5f' % w + '\n')
    t = t + h

arch.close
