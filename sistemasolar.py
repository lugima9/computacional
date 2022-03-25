#//////////////////////////////////////////////////////////////
# Programa que simula la dinámica del Sistema Solar.
# OUTPUT: Periodo de traslación de cada planeta y comparación el tabulado,
#         energía cinética y momento angular total del sistema (conservado?), 
#         animación de las órbitas en modelo heliocéntrico y geocéntrico
#         (restando la posición de la Tierra a la de cada planeta)         
# Autor: Luis Gil Martín
# Universidad de Granada
#//////////////////////////////////////////////////////////////

import math as m
import matplotlib.pyplot as plt
import numpy as np
import random

#-------------------------------------CLASES--------------------------------------------
class Planetas:
    def __init__(self, nombre, masa, distMedia, vel, *args):
        self.nombre = nombre
        self.masa = masa
        self.distMedia = distMedia
        self.vel = vel

        self.args = args

#------------------------------------FUNCIONES------------------------------------------

# Función que calcula la aceleración (x,y) de un planeta en un instante por acción del Sol y el resto de planetas
# Recibe: índice del planeta, array de posiciones de cada planeta en (x,y), vector de masas de cada planeta
# Devuelve: array 2D con aceleración (x,y) del planeta
def calcAcelPlaneta(estePlaneta, posiciones, masas):

    # Posición (x,y) con respecto al Sol del planeta y distancia cubo
    difPos = posiciones[estePlaneta]

    distCubo = np.linalg.norm(difPos)**3

    # Aceleración (x,y) del planeta debido al Sol
    acel = -difPos/distCubo

    # Aceleración (x,y) del planeta debido al resto de planetas
    for i in range(len(posiciones)):
        if i != estePlaneta:
            difPos = posiciones[estePlaneta] - posiciones[i]      #Posición relativa entre planetas
            distCubo = np.linalg.norm(difPos)**3
            acel -= masas[i]*difPos/distCubo                      #Sumamos la aceleración por el resto de planetas a la del Sol

    return acel


# Función que calcula un paso del algoritmo de Verlet para obtener la posición y velocidad de cada planeta en el tiempo
# Recibe: paso, masas, posiciones, velocidades
# Devuelve: array de posiciones, array de velocidades
def pasoVerlet(h, m, r, v):

    #Creo dos arrays 3D auxiliares de 8 (nº planetas) x 2 (x,y)
    acel, w = np.zeros(8, 2), np.zeros(8, 2)

    # Para cada uno de los 8 planetas...
    for i in range(len(r)):
        #Calcula la aceleración (x,y) de cada planeta       
        acel[i] = calcAcelPlaneta(i, r, m) #OJO: se mete el array r entero
        
        #Cálculo intermedio (x,y) de cada planeta
        w[i] = v[i] + 0.5*h*acel[i]

        #Calcula la nueva posición (x,y) de cada planeta
        r[i] = r[i] + h*w[i]

        #Calcula la nueva aceleración (x,y) de cada planeta
        acel[i] = calcAcelPlaneta(i, r, m)

        #Calcula la nueva velocidad (x,y) de cada planeta
        v[i] = w[i] + 0.5*h*acel[i]

    return r, v


# Función que obtiene los datos de los planetas de un fichero de texto
# Recibe:  paso, nombre del fichero con una fila por planeta, bool de alineación
# Devuelve: lista de planetas con sus atributos actualizados
def condIni(paso, nombreFichero, alinear):
    # Abrimos el fichero en modo lectura
    f = open(nombreFichero, "r")

    # Leemos el archivo saltando el header y separamos en lista por planeta
    next(f)
    datos = f.read().split('\n')
    f.close()

    # Guardamos cada dato de cada planeta en un elemento distinto de cada array
    nombres, masas, distancias, velocidades = [], [], [], []
    for i in range(len(datos)):
        d = datos[i].split()
        nombres[i] = d[0]
        masas[i] = d[1]
        distancias[i] = d[2]
        velocidades[i] = d[3]
    
    # Reescalamos los datos
    paso, masas, distancias, velocidades = reescala(paso, masas, distancias, velocidades)

    # Iniciamos posiciones y velocidades

    if alinear == True:
        r, v = inicioAlineado(distancias, velocidades)
    else:
        r, v = inicioAleatorio(distancias, velocidades)

    return paso, nombres, masas, r, v


# Función que realiza el cambio de unidades de reescalado
# Recibe: paso, array de masas, array de posiciones iniciales y array de velocidades iniciales
# Devuelve: paso, array de masas, array de posiciones iniciales y array de velocidades iniciales (reescalado)
def reescala(paso, masas, distancias, velocidades):

    c = 299792458           # [m s-1]
    masaSol = 1.98847e30    # [kg]
    G = 6.67430e-11         # [m3 kg-1 s-2]

    paso = paso*(G*masaSol/c**3)**(0.5)
    masas = masas/masaSol
    distancias = distancias/c
    velocidades = velocidades*(c/(G*masaSol))**(0.5)

    return paso, masas, distancias, velocidades

# Función que aleatoriza el inicio de la simulación
# Recibe: distancias y velocidades iniciales
# Devuelve: array de posiciones y array de velocidades con formato (x,y) para cada elemento-planeta
def inicioAleatorio(distancias, velocidades):

    # Para cada planeta (tantos como elementos en distancias[])
    r, v = np.zeros(8, 2), np.zeros(8, 2)
    for i in range(len(distancias)):

        # Genero una coordenada polar aleatoria
        angulo = random.randrange(0., 2.*np.pi)
        coseno = m.cos(angulo)
        seno = m.sin(angulo)

        # Calculo una posición inicial en cartesianas
        r[i,1] = distancias[i]*coseno
        r[i,2] = distancias[i]*seno

        # Calculo una velocidad inicial perpendicular al vector r en cartesianas
        v[i,1] = -velocidades[i]*seno
        v[i,2] = velocidades[i]*coseno        

    return r, v


# Función que alinea ápsides planetas al inicio de la simulación
# Recibe: distancias y velocidades iniciales
# Devuelve: array de posiciones y array de velocidades con formato (x,y) para cada elemento-planeta
def inicioAlineado(distancias, velocidades):

    # Para cada planeta (tantos como elementos en distancias[])
    r, v = np.zeros(8, 2), np.zeros(8, 2)
    for i in range(len(distancias)):

        # Fijo posición inicial en (x,0)
        r[i,1] = distancias[i]
        r[i,2] = 0

        # Fijo velocidad inicial en (0, v)
        v[i,1] = 0
        v[i,2] = velocidades[i]       

    return r, v

#------------------------------------FUNCIÓN PRINCIPAL------------------------------------------

# Leemos las condiciones iniciales y las guardamos en los correspondientes arrays (ya reescaladas)
h = float(input("Introduzca el paso: \n"))

alinear = ""
while alinear != "Y" and alinear != "N":
    alinear = bool("¿Desea que los planetas comiencen alineados? (Y/N): \n")

nombreArchivo = "datos_planetas.txt"
h, nombres, masas, r, v = condIni(h, nombreArchivo, alinear)

# Simulamos el Sistema Solar durante un tiempo fijado
t=0
t_max = 100

while t < t_max:

    r, v = pasoVerlet(h, m, r, v)

    # DEBO IR ESCRIBIENDO LOS DATOS EN UN FICHERO EN EL FORMATO QUE INDICA EL SCRIPT DE JUANI

    t += h









