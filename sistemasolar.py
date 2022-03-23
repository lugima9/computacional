#//////////////////////////////////////////////////////////////
# Programa que simula la dinámica del Sistema Solar.
# OUTPUT: Periodo de traslación de cada planeta y comparación el tabulado,
#         energía cinética y momento angular total del sistema (conservado?), 
#         animación de las órbitas en modelo heliocéntrico y geocéntrico
#         (restando la posición de la Tierra a la de cada planeta)         
# Autor: Luis Gil Martín
# Universidad de Granada
#//////////////////////////////////////////////////////////////

#Nota: Para el cálculo del periodo no podemos comprobar que pase dos veces por exactamente el mismo punto.
#Podemos comprobar que el vector posición cambia de cuadrante con desigualdades u otra cosa que se nos ocurra.

#Nota: Instalar ffmpeg extensión para convertir vídeo

import math as m
import matplotlib.pyplot as plt
import numpy as np

#Clase de planetas (PROBAR CON LOS **ARGS)
class Planetas:
    def __init__(self, nombre, masa, distMedia, vel):
        self.nombre = nombre
        self.masa = masa
        self.distMedia = distMedia
        self.vel = vel

# AYUDA PARA VISUALIZACIÓN
# --------------
# Array de posiciones instantáneas:
# r = [[x_0, y_0], ..., [x_7, y_7]]
# Si quiero la posición del primer planeta: posicion = r[0]
# --------------
# Array de velocidades instantáneas:
# (ANÁLOGO)
# --------------

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
# Recibe: paso, instante, posiciones, velocidades
# Devuelve: nuevos instante, array de posiciones, array de velocidades

def pasoVerlet(h, t, r, v):

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
    
    t += h

    return t, r, v

# Función que obtiene los datos de los planetas de un fichero de texto
# Recibe:  nombre del fichero con una fila por planeta
# Devuelve: lista con string de datos para cada planeta 

def leerCondIni(nombreFichero):
    # Abrimos el fichero en modo lectura
    f = open(nombreFichero, "r")

    # Leemos el archivo saltando el header y separamos en lista por planeta
    next(f)
    datos = f.read().split('\n')

    # Asignamos a cada planeta sus atributos
    p = []
    for i in range(len(datos)):
        d = datos[i].split()
        p[i] = Planetas(d[0], d[1], d[2], d[3])

    f.close()

    return p

# Función que realiza el cambio de unidades de reescalado
# Recibe: paso, array de posiciones iniciales, array de masas
# Devuelve: lista con paso, array de posiciones iniciales, array de masas (reescalado)
def reescala(h, r, m):
    c = 299792458           #[m s-1]
    masaSol = 1.98847e30    # [kg]
    G = 6.67430e-11         # [m3 kg-1 s-2]

    h = sqrt(G*masaSol/c**3)*h
    r = r/c
    m = m/masaSol

    return (h, r, m)

# FUNCIÓN PRINCIPAL

#

t=0
h = 0.01
t_max = 100
