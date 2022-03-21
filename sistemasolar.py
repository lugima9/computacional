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

class Planetas:
    def __init__(self, nombre, masa, distMedia, periodo, posicion, velAng):
        self.nombre = nombre
        self.masa = masa
        self.distMedia = distMedia
        self.periodo = periodo

        self.posicion = posicion
        self.velAng = velAng

# Constantes globales

global masaSol = 1.98847e30 # [kg]
global G = 6.67430e-11 # [m3 kg-1 s-2]

# Lectura de datos de planetas

# Asignación de condiciones iniciales

# Algoritmo de Verlet

def Verlet(h, t, r, v):
    a

    t = t+h
    for i in range(0, 7, 1)
    r[i][0] = 


    def VerletPendulo(g, t_0, h, ang_0, vAng_0):
    #Valores iniciales de las listas
    t, ang, vAng, aAng, w = []*101, []*101, []*101, []*101, []*101

    t.insert(0, t_0)
    ang.insert(0, ang_0)
    vAng.insert(0, vAng_0)
    aAng.insert(0, -g*math.sin(ang_0))
    for i in range(0, 100, 1):
        t.insert(i+1, t[i] + h)
        w.insert(i, vAng[i] + 0.5*h*aAng[i])
        ang.insert(i+1, ang[i] + h*w[i])
        aAng.insert(i+1, -g*math.sin(ang[i+1]))
        vAng.insert(i+1, w[i] + 0.5*h*aAng[i+1])

    return (t, ang, vAng)


