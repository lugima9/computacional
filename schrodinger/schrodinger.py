"//////////////////////////////////////////////////////////////"
# Programa que resuelve la ecuación de Schrödinger para un pozo
# de potencial cuadrado.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import math as m
from math import pi
import numpy as np
import random as rd
from numba import njit

"-------------------------------------FUNCIONES------------------------------------------"

# Función que devuelve una función de onda compleja de tipo gaussiano normalizada a la unidad
# Recibe: N, k0
# Devuelve: array espacial con función de onda compleja gaussiana normalizada a 1
@njit
def iniGauss(N, k0):

    # Condiciones de contorno: phi(0) = phi(N-1) = 0]
    waveFunc = np.zeros(N, dtype=np.complex_)

    # phi = 1/norm * e^{i k0 j} e^{(x-x0)^2/2sigma^2}}
    norm = 0.
    for p in range(1, N-1):
        waveFunc[p] = m.exp(-8*((4*p - N)/N)**2) * (m.cos(k0*p) + m.sin(k0*p)*1j)
        aux = np.abs(waveFunc[p])
        norm += aux*aux

    norm = norm**0.5
    waveFunc /= norm

    return waveFunc

# Función que define una barrera de potencial cuadrada
# Recibe: lambda, k0
# Devuelve: array espacial con barrera de potential
@njit
def definePotential(L, K0):
    
    V = L*k0*k0

    pot = np.zeros(N, dtype=np.complex_)
    potWidth = N/5
    potCenter = N/2

    for j in range(N):
        if j > (potCenter - potWidth/2) and j < (potCenter +  potWidth/2):
            pot[j] = V
    
    return pot

# Función que computa la evolución temporal de la función de onda dada inicial
# Recibe: tiempo de simulación (pasos), número de celdas espaciales, función de onda, potencial
# Devuelve: función de onda actualizada en cada instante
@njit
def timeEvolution(tSteps, N, phi, pot):

    # CÁLCULOS PREVIOS
    #-----------------

    # Definimos los array A
    A_minus = np.ones(N, dtype=np.complex_)
    A_plus = A_minus

    A_zero = np.zeros(N, dtype=np.complex_)
    A_zero += -2 + (2/s_tilde)*1j - pot

    # Calculamos el array alpha por recurrencia (alpha(N-1) = 0)
    alpha = np.zeros(N, dtype=np.complex_)
    gamma = np.zeros(N, dtype=np.complex_)

    for j in range(N-2, 0, -1):
        gamma[j] = 1/(A_zero[j] + A_plus[j]*alpha[j])
        alpha[j-1] = -A_minus[j]*gamma[j]

    gamma[0] = 1/(A_zero[0] + A_plus[0]*alpha[0])

    # EVOLUCIÓN TEMPORAL
    #-------------------

    beta = np.zeros(N, dtype=np.complex_)
    b = np.zeros(N, dtype=np.complex_)
    chi = np.zeros(N, dtype=np.complex_)

    for n in range(tSteps-1):

        # Calculamos el array b a partir de la función de onda en este instante
        b = (4/s_tilde)*1j*phi[n]

        # Calculamos el array beta por recurrencia (beta(n, N-1) = 0)
        for j in range(N-2, 0, -1):
            beta[j-1] = gamma[j]*(b[j] - A_plus[j]*beta[j])

        # Calculamos el array chi por recurrencia
        for j in range(0, N-1):
            chi[j+1] = alpha[j]*chi[j] + beta[j]

        # Calculamos el array phi en el nuevo instante
        for j in range(1, N-1):
            phi[n+1, j] = chi[j] - phi[n, j]

    return phi

# Función que escribe los resultados en un fichero en el formato de animación.
# Recibe: función de onda en cada instante, barrera de potencial, elección de representación,
#         tiempo de simulación (pasos)
# Devuelve: -
def waveToFile(phi, pot, choice, tSteps, h):
    
    file = open("results.dat", "w")
    norm = np.zeros(tSteps)
    
    pot = np.real(pot)

    if choice == 1:
        phiReal = np.real(phi)
        phiMod = np.abs(phi)

        for n in range(tSteps):
            norm[n] = np.sum(np.square(phiMod[n]))
            for j in range(N):
                file.write(f"{h*j}, \t {phiReal[n,j]}, \t {norm[n]}, \t {pot[j]} \n")

            file.write("\n")

    elif choice == 2:
        phiImag = np.imag(phi)
        phiMod = np.abs(phi)

        for n in range(tSteps):
            norm[n] = np.sum(np.square(phiMod[n]))
            for j in range(N):
                file.write(f"{h*j}, \t {phiImag[n,j]}, \t {norm[n]}, \t {pot[j]} \n")

            file.write("\n")

    elif choice == 3:
        phiMod = np.abs(phi)

        for n in range(tSteps):
            norm[n] = np.sum(np.square(phiMod[n]))
            for j in range(N):
                file.write(f"{h*j}, \t {phiMod[n,j]}, \t {norm[n]}, \t {pot[j]} \n")

            file.write("\n")

    else:
        phiReal = np.real(phi)
        phiMod = np.abs(phi)

        for n in range(tSteps):
            norm[n] = np.sum(np.square(phiMod[n]))
            for j in range(N):
                file.write(f"{h*j}, \t {phiReal[n,j]}, \t {phiMod[n,j]}, \t {pot[j]} \n")

            file.write("\n")

    file.close()


"-------------------------------------MAIN------------------------------------------"

#----------------------------------------------------
# PARÁMETROS Y CONDICIONES INICIALES DE LA SIMULACIÓN
#----------------------------------------------------

# Discretizamos espacio y tiempo
N = 200                 
tSteps = 1000           
h = 1                   

# Parámetros iniciales
L = 1.3
cycles = N/16
k0 = 2.*pi*cycles/(N*h)
s_tilde = 1/(4.*k0*k0)

# Definimos la función de onda en el instante inicial con condiciones de
# contorno (pozo infinito de ancho N*h) y forma gaussiana
# phi(n, 0) = phi(n, N) = 0     n = 0, 1, ..., tSteps

phi = np.zeros((tSteps, N), dtype=np.complex_)     # Matriz función de onda en cada celda y paso temporal

phi[0] = iniGauss(N, k0)

#-----------------------------------
# CÁLCULO DE LA BARRERA DE POTENCIAL
#-----------------------------------
pot = definePotential(L, k0)

#------------------------------------
# SIMULACIÓN DE LA EVOLUCIÓN TEMPORAL
#------------------------------------
phi = timeEvolution(tSteps, N, phi, pot)

#--------------------------------------
# GUARDADO DE LOS RESULTADOS EN FICHERO
#--------------------------------------

# Elegimos si guardar parte real, imaginaria, módulo o real-módulo
choice = int(input("Real (1), Imag. (2), Modulo (3), Real y Modulo (4): "))

waveToFile(phi, pot, choice, tSteps, h)

print("El archivo se ha creado con exito en la localizacion del script.")