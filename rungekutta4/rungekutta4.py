"//////////////////////////////////////////////////////////////"
# Programa que simula el movimiento de una nave que despega de la
# Tierra bajo el efecto gravitatorio de la Tierra y la Luna mediante RK4.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

from cmath import pi
from pyexpat.errors import XML_ERROR_PARTIAL_CHAR
import numpy as np
import matplotlib.pyplot as plt
from numba import njit


"-------------------------------------FUNCIONES------------------------------------------"
# Función que calcula las condiciones iniciales del problema reescaladas.
# Recibe: -
# Devuelve: r_ini, phi_ini, p_r_ini, p_phi_ini (reescalados)
@njit
def condIni():

    # Parámetros
    G = 6.67e-11        #[N m^2 kg^-2]
    M_T = 5.9736e24     #[kg]
    R_T = 6.378160e6    #[m]
    M_L = 0.07349e24    #[kg]
    R_L = 1.7374e6      #[m]
    d_TL = 3.844e8      #[m]
    omega = 2.6617e-6   #[rad/s]

    # Punto de lanzamiento: (r = R_T, longitud = 90º)
    r_ini = R_T / d_TL
    phi_ini = 90.*pi/180.
    # Lanzamiento vertical con una fracción de la velocidad de escape terrestre
    vEscape = (2.*G*M_T/R_T)**0.5 / d_TL
    p_r_ini = 0.99213*vEscape
    p_phi_ini = 0

    return r_ini, phi_ini, p_r_ini, p_phi_ini


# Función que calcula el coeficiente k_i_j para cada variable j
# Recibe: r, phi, p_r, p_phi, instante, paso, omega, nabla, mu
# Devuelve: step*array j de coeficientes k_i_j
@njit
def evalEquation(r, phi, p_r, p_phi, t, step, omega, nabla, mu):

    k = np.zeros(4)

    phase = phi - omega*t
    r2 = r*r
    r3 = r2*r
    r3_prime = (1 + r2 - 2*r*np.cos(phase))**(3/2)

    k[0] = p_r
    k[1] = p_phi/r2
    k[2] = p_phi*p_phi/r3 - nabla*(1/r2 + mu/r3_prime * (r - np.cos(phase)))
    k[3] = -nabla*mu*r/r3_prime * np.sin(phase)

    return step*k


# Función que computa la simulación RK4 del sistema
# Recibe: paso, duración de la simulación
# Devuelve: r, phi, p_r, p_phi
@njit
def rk4(step, duration):

    # Parámetros
    G = 6.67e-11        #[N m^2 kg^-2]
    M_T = 5.9736e24     #[kg]
    R_T = 6.378160e6    #[m]
    M_L = 0.07349e24    #[kg]
    R_L = 1.7374e6      #[m]
    d_TL = 3.844e8      #[m]
    omega = 2.6617e-6   #[rad/s]

    nabla = G*M_T/(d_TL**3)
    mu = M_L/M_T

    # Variables 
    iters = duration//step

    r, phi = np.zeros(iters+1), np.zeros(iters+1)
    p_r, p_phi = np.zeros(iters+1), np.zeros(iters+1) 

    k1, k2, k3, k4 = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)

    # Introducimos condiciones iniciales
    r[0], phi[0], p_r[0], p_phi[0] = condIni()

    # Simulación en el tiempo
    for i in range(iters):
        
        t = i*step  # Instante actual

        # Evaluamos k1
        k1 = evalEquation(r[i], phi[i], p_r[i], p_phi[i], t, step, omega, nabla, mu)

        # Evaluamos k2
        taux = t + step*0.5
        a = r[i] + k1[0]*0.5
        b = phi[i] + k1[1]*0.5
        c = p_r[i] + k1[2]*0.5
        d = p_phi[i] + k1[3]*0.5

        k2 = evalEquation(a, b, c, d, taux, step, omega, nabla, mu)

        # Evaluamos k3
        a = r[i] + k2[0]*0.5
        b = phi[i] + k2[1]*0.5
        c = p_r[i] + k2[2]*0.5
        d = p_phi[i] + k2[3]*0.5

        k3 = evalEquation(a, b, c, d, taux, step, omega, nabla, mu)

        # Evaluamos k4
        taux = t + step
        a = r[i] + k3[0]
        b = phi[i] + k3[1]
        c = p_r[i] + k3[2]
        d = p_phi[i] + k3[3]

        k4 = evalEquation(a, b, c, d, taux, step, omega, nabla, mu)

        # Actualizamos el valor de las variables solución
        r[i+1] = r[i] + (k1[0] + 2.*k2[0] + 2.*k3[0] + k4[0])/6.               
        phi[i+1] = phi[i] + (k1[1] + 2.*k2[1] + 2.*k3[1] + k4[1])/6. 
        p_r[i+1] = p_r[i] + (k1[2] + 2.*k2[2] + 2.*k3[2] + k4[2])/6. 
        p_phi[i+1] = p_phi[i] + (k1[3] + 2.*k2[3] + 2.*k3[3] + k4[3])/6. 

    return r, phi, p_r, p_phi


# Función que computa la órbita de la Luna asumiendo MCU
# Recibe: paso, duración de la simulación
# Devuelve: r, phi (Luna)
@njit
def moonTraj(step, duration):

    omega = 2.6617e-6   #[rad/s]

    iters = duration//step

    # MCU reescalado con fase inicial 0 (empieza sobre eje X)
    r = np.ones(iters)
    
    phi = np.zeros(iters)
    for i in range(iters):
        t = i*step
        phi[i] = t
    phi *= omega

    return r, phi


# Función que realiza el cambio de coordenadas polares a cartesianas
# Recibe: r, phi
# Devuelve: x, y
@njit
def cilToCart(r, phi):

    x = r*np.cos(phi)
    y = r*np.sin(phi)

    return x, y

"-------------------------------------MAIN------------------------------------------"

# Condiciones de simulación
step = 1
duration = int(1e6)

# Nave
r_ship, phi_ship, p_r, p_phi = rk4(step, duration)
x_ship, y_ship = cilToCart(r_ship, phi_ship)

# Luna
r_moon, phi_moon = moonTraj(step, duration)
x_moon, y_moon = cilToCart(r_moon, phi_moon)

# Escribimos los resultados en formato animacion.py
fileName = "resultados.dat"
with open(fileName, "w") as file:

    iters = duration//step

    mostrar = 1000
    for i in range(iters):

        # Escribimos las posiciones en el fichero (1 de cada 50)
        if mostrar == 1000:
            file.write(f"{0}, {0}\n") # Tierra
            file.write(f"{x_moon[i]},  {y_moon[i]}\n") # Luna           
            file.write(f"{x_ship[i]},  {y_ship[i]}\n") # Nave

            file.write("\n")

            mostrar = 0
        
        mostrar += 1

print("El archivo ha sido creado con exito en la localizacion del script.")
