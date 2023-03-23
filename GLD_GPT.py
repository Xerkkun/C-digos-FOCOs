import numpy as np
import matplotlib.pyplot as plt

# Definición de parámetros del sistema
sigma = 10.0
beta = 8./3.
rho = 28.0

# Definición del orden fraccionario
q = 0.5

def binomial_coef(n, k):
    if k == 0 or k == n:
        return 1
    else:
        return binomial_coef(n-1, k-1) + binomial_coef(n-1, k)

# Definición de la función del sistema de Lorenz
def lorenz_frac(t, x):
    xdot = np.zeros_like(x)
    xdot[0] = sigma * (x[1] - x[0])
    xdot[1] = rho * x[0] - x[1] - x[0] * x[2]
    xdot[2] = -beta * x[2] + x[0] * x[1]
    return xdot

# Definición del método de Grünwald-Letnikov para la derivada de orden fraccionario
def grunwald_letnikov(f, x0, q, h, N):
    gl = 0.0
    for k in range(N):
        gl += (-1)**k * binomial_coef(q, k) * f(x0 - k*h)
    gl /= h**q
    return gl

# Definición del intervalo de tiempo y el paso de integración
t0 = 0.0
tf = 100.0
dt = 0.01

# Definición del vector de estado inicial
x0 = np.array([1.0, 1.0, 1.0])

# Definición del número de puntos en la suma de Grünwald-Letnikov
N = 10

# Inicialización del vector de estado y tiempo
x = np.zeros((int((tf-t0)/dt)+1, 3))
t = np.linspace(t0, tf, len(x))

# Condición inicial
x[0,:] = x0

# Integración numérica del sistema de Lorenz con Grünwald-Letnikov
for i in range(len(x)-1):
    # Cálculo de la derivada de orden fraccionario con Grünwald-Letnikov
    dxdt = np.zeros(3)
    for j in range(3):
        f = lambda x: lorenz_frac(x, x0)[j]
        dxdt[j] = grunwald_letnikov(f, t[i], q, dt, N)

    # Cálculo del siguiente estado del sistema
    x[i+1,:] = x[i,:] + dt * dxdt

# Gráfica de la solución
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x[:,0], x[:,1], x[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
