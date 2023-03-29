# importar todas las funciones de pylab
from pylab import *

# importar el módulo pyplot
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import axes3d
from matplotlib import style

t, x, y, z, w = np.loadtxt("FE.dat",unpack=True)  #Asigna a las variables los datos de cada columna del archivo de datos

subplot(411)
p1,=plot(t,x,"m",lw=0.7)
ylabel("x")
plt.xticks([])

subplot(412)
p2,=plot(t,y,"m",lw=0.7)
ylabel("y")
plt.xticks([])

subplot(413)
p2,=plot(t,z,"m",lw=0.7)
ylabel("z")
plt.xticks([])

subplot(414)
p3,=plot(t,w,"m",lw=0.7)
ylabel("w")
xlabel("Tiempo")
plt.savefig("time.pdf",dpi=300,bbox_inches= 'tight')

clf()

subplot(221)
p1,=plot(x,y,"m",lw=0.3)
xlabel("x")
ylabel("y")
#plt.axis('square')

subplot(222)
p2,=plot(y,z,"m",lw=0.3)
xlabel("y")
ylabel("z")

subplot(223)
p3,=plot(x,z,"m",lw=0.3)
xlabel("x")
ylabel("z")

subplot(224)
p2,=plot(x,w,"m",lw=0.3)
xlabel("x")
ylabel("w")

plt.savefig("atractores.pdf",dpi=300,bbox_inches= 'tight')
clf()

#fig = plt.figure()
#ax1 = fig.add_subplot(111,projection='3d')
#ax1.plot(x,y,z,"m",lw=0.2)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
#ax1.set_zlabel('z')
#plt.savefig("atractor3d.pdf",dpi=300,bbox_inches= 'tight')
#clf()

p1,=plot(x,y,"m",lw=0.3)
xlabel("x")
ylabel("y")
plt.savefig("x-y.pdf",dpi=300,bbox_inches= 'tight')
clf()

p2,=plot(y,z,"m",lw=0.3)
xlabel("y")
ylabel("z")
plt.savefig("y-z.pdf",dpi=300,bbox_inches= 'tight')
clf()

p3,=plot(x,z,"m",lw=0.3)
xlabel("x")
ylabel("z")
plt.savefig("x-z.pdf",dpi=300,bbox_inches= 'tight')
clf()

p4,=plot(x,w,"m",lw=0.3)
xlabel("x")
ylabel("w")
plt.savefig("x-w.pdf",dpi=300,bbox_inches= 'tight')
clf()
#plt.show()

metodos = ["FE","BE","RK4","AB6","AM4_FE","GEAR4"]
colores = ["m","g","r","c","b","y"]
ancho = ["0.001", "0.01", "0.001", "0.001", "0.005", "0.005"]

j = -1
h = -1
for i in metodos:
	archivo = i + ".dat"
	t, x, y, z, w = np.loadtxt(archivo,unpack=True)
	j = j + 1
	h = h + 1
	etiquetas = i + ", h=" + ancho[h]

	subplot(411)
	p1,=plot(t,x,colores[j],lw=0.3,label=etiquetas)
	ylabel("x")
	plt.xticks([])
	plt.xlim([0,5])

	subplot(412)
	p2,=plot(t,y,colores[j],lw=0.3,label=etiquetas)
	ylabel("y")
	plt.xticks([])
	plt.xlim([0,5])
	plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	subplot(413)
	p3,=plot(t,z,colores[j],lw=0.3,label=etiquetas)
	ylabel("z")
	plt.xticks([])
	plt.xlim([0,5])

	subplot(414)
	p4,=plot(t,w,colores[j],lw=0.3,label=etiquetas)
	ylabel("w")
	xlabel("Tiempo")
	plt.xlim([0,5])


	#plt.legend(loc='upper center', bbox_to_anchor=(0.5, -.25),
         # fancybox=True, shadow=True, ncol=5)

plt.savefig("times.pdf",dpi=300,bbox_inches= 'tight')
