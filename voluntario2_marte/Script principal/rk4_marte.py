"//////////////////////////////////////////////////////////////"
# Programa que simula el envío de una nave desde la Tierra hasta
# una órbita baja en Marte mediante una órbita de transferencia
# de Hohmann.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

from cmath import pi
import numpy as np
from numba import njit # Compilador de Python que mejora el rendimiento en la ejecución del código

"-------------------------------------GLOBALES-------------------------------------------"

# Parámetros
G = 6.67e-11                #[N m^2 kg^-2]
M_S = 1.98847e30            #[kg]

M_T = 5.9736e24             #[kg]
R_T = 6.378160e6            #[m]

M_M = 6.39e23               #[kg]
R_M = 3.3895e6              #[m]

d_ST = 1.49597870700e11     #[m]
d_SM = 2.27940000000e11     #[m]

omegas = np.array([(G*M_S/d_ST**3)**0.5, (G*M_S/d_SM**3)**0.5])          #[rad/s] (Tierra, Marte)


"-------------------------------------FUNCIONES------------------------------------------"
# Función que calcula las condiciones iniciales del problema reescaladas.
# Recibe: radio de órbita terrestre de partida de la Nave, factor de corrección
# Devuelve: r_N, phi_N, p_r_N, p_phi_N, delta_v1, v_ini, v_p
@njit
def condIni(r_orb_T, correction):

    # Radio de la SOI terrestre
    R_SOI = (M_T/M_S)**(2/5) * d_ST

    # Semieje mayor de la órbita de transferencia
    a_hoh = (d_ST + d_SM)/2 * correction

    # Exceso de velocidad hiperbólico al llegar a la SOI terrestre
    v_inf = (G*M_S/d_ST)**0.5 * ((2 - d_ST/a_hoh)**0.5 - 1)

    # Semieje mayor de la hipérbola de escape
    a_hyp = 1/(2/R_SOI - v_inf**2/(G*M_T))

    # Velocidad inicial para escape desde órbita de aparcamiento terrestre
    v_p = (G*M_T*(2/r_orb_T - 1/a_hyp))**0.5

    # Delta-v de salida
    v_ini = (G*M_T/r_orb_T)**0.5
    delta_v1 = np.abs(v_p - v_ini)

    """"""""""""""""""""""""
    print("Delta v1 (m/s):")
    print(delta_v1)
    """"""""""""""""""""""""

    # Ángulo entre perigeo de la hipérbola y vector velocidad de la Tierra
    beta = np.arccos(1 / (1 + r_orb_T*v_inf**2 / (G*M_T)))

    # Posición / velocidad inicial en órbita de aparcamiento con respecto a la Tierra
    r_NT = r_orb_T
    phi_NT = -(pi/2 - beta)
    x_NT, y_NT = polToCart(r_NT, phi_NT)

    vr_NT = 0
    omega_NT = v_p / r_orb_T
    vx_NT, vy_NT = polToCartVel(r_NT, phi_NT, vr_NT, omega_NT)

    # Posición / velocidad inicial en órbita de aparcamiento con respecto al Sol
    # La Tierra comienza en y = 0 => x = r, vy = omega_T * d_ST
    x_N, y_N = (x_NT + d_ST), y_NT 
    r_N, phi_N = cartToPol(x_N, y_N)

    vx_N, vy_N = vx_NT, (vy_NT + omegas[0]*d_ST)
    vr_N, omega_N = cartToPolVel(x_N, y_N, vx_N, vy_N)

    """"""""""""""""""""""""
    print("r_ini de la nave (m):")
    print(r_N)
    print("phi_ini de la nave (rad)")
    print(phi_N)
    """"""""""""""""""""""""

    # Reescalamos posiciones y velocidades
    r_N /= d_SM
    vr_N /= d_SM

    # Momentos en órbita terrestre (reescalados)
    p_r_N = vr_N
    p_phi_N = omega_N * r_N**2

    return r_N, phi_N, p_r_N, p_phi_N, delta_v1, v_ini, v_p


# Función que determina la posición angular inicial de Marte (ventana de lanzamiento) y parámetros
# de la órbita de Hohmann. 
# Recibe: -
# Devuelve: phi_ini_M, a_hoh, duración de la transferencia
@njit
def hohmann():

    # Parámetros de la órbita de Hohmann
    a_hoh = (d_ST + d_SM)/2                      # Semieje mayor
    T_hoh = 2*pi*(a_hoh**3/(G*M_S))**0.5         # Periodo

    """"""""""""""""""""""""
    print("Semieje mayor Hohmann (m):")
    print(a_hoh)
    print("Semiperiodo Hohmann (s)")
    print(T_hoh/2)
    """"""""""""""""""""""""

    # Ventana de lanzamiento
    phi_ini_M = pi - 0.5 * T_hoh * omegas[1]

    """"""""""""""""""""""""
    print("Phi_ini_M (rad):")
    print(phi_ini_M)
    """"""""""""""""""""""""

    return phi_ini_M, a_hoh, T_hoh*0.5


# Función que calcula la velocidad para situarse en órbita marciana.
# Recibe: d_NM, r_M, phi_M, r, phi, p_r, p_phi
# Devuelve: p_r_def, p_phi_def, delta_v2, v_N, v_def
@njit
def startOrbit(d_NM, r_M, phi_M, r, phi, p_r, p_phi):
    
    # Velocidad lineal de Marte con respecto al Sol en cartesianas (reescalada)
    vx_M, vy_M = polToCartVel(r_M, phi_M, 0, omegas[1])

    # Velocidad lineal inicial de Nave con respecto al Sol en cartesianas (reescalada)
    vr_N = p_r                          # La velocidad radial reescalada es p_r reescalado
    omega_N = p_phi / (r*r)

    vx_N, vy_N = polToCartVel(r, phi, vr_N, omega_N)

    # Posición cart. Nave con respecto a Marte (reescalada)
    x_N, y_N = polToCart(r, phi)        # con respecto al Sol
    x_M, y_M = polToCart(r_M, phi_M)    # con respecto al Sol

    x_NM, y_NM = (x_N-x_M), (y_N-y_M)

    # Posición pol. Nave con respecto a Marte (reescalada)
    r_NM, phi_NM = cartToPol(x_NM, y_NM)

    # Vel. ang. Nave con respecto a Marte tal que hay MCU
    omega_NM = (G*M_M / (d_NM**3))**0.5

    # Velocidad lineal cart. con respecto a Marte tal que hay MCU (reescalada)
    vx_NM, vy_NM = polToCartVel(r_NM, phi_NM, 0., omega_NM)

    # Velocidad lineal cart. con respecto al Sol tal que hay MCU con respecto a Marte y se orbita el Sol
    # (Sumarle a la velocidad para orbitar Marte con respecto a Marte, 
    # la velocidad de Marte con respecto al Sol)
    vx_def, vy_def = (vx_NM + vx_M), (vy_NM + vy_M)

    # Pasamos a polares y obtenemos p_r y p_phi nuevos
    vr_def, omega_def = cartToPolVel(x_N, y_N, vx_def, vy_def)

    p_r_def = vr_def
    p_phi_def = omega_def * (r*r)

    # Calculamos el delta_v2 para poner en órbita (desescalado)
    v_N = (vx_N*vx_N + vy_N*vy_N)**0.5 * d_SM
    v_def = (vx_def*vx_def + vy_def*vy_def)**0.5 * d_SM
    delta_v2 = ((vx_def - vx_N)**2 + (vy_def - vy_N)**2)**0.5 * d_SM

    return p_r_def, p_phi_def, delta_v2, v_N, v_def


# Función que calcula el coeficiente k_i_j para cada variable j
# Recibe: r, phi, p_r, p_phi, instante, paso, omegas, delta, mus
# Devuelve: paso*array j de coeficientes k_i_j de Runge-Kutta
@njit
def evalEquation(r, phi, p_r, p_phi, t, step, omegas, delta, mus):

    k = np.zeros(4)

    # Cálculos previos
    phases = np.array([phi - omegas[0]*t, phi - (omegas[1]*t + phi_M_ini)]) # Marte NO parte del eje X, como la Tierra
    r2 = r*r
    r3 = r2*r
    r3_primes = np.array([((d_ST/d_SM)**2 + r2 - 2*r*d_ST/d_SM*np.cos(phases[0]))**(3/2), (1 + r2 - 2*r*np.cos(phases[1]))**(3/2)])

    # Coeficientes
    k[0] = p_r
    k[1] = p_phi/r2
    k[2] = p_phi*p_phi/r3 - delta*(1/r2 + mus[0]/r3_primes[0] * (r - d_ST/d_SM*np.cos(phases[0])) + mus[1]/r3_primes[1] * (r - np.cos(phases[1])))
    k[3] = -delta*(mus[0]*r/r3_primes[0] * d_ST/d_SM * np.sin(phases[0]) + mus[1]*r/r3_primes[1] * np.sin(phases[1]))

    return step*k


# Función que computa la simulación RK4 del sistema
# Recibe: paso, duración de la simulación, phi_M_ini, radio para orbitar Marte, factor de corrección
# Devuelve: r, phi, r_T, phi_T, r_M, phi_M, delta_v1, v_ini, v_p, delta_v2, v_N, v_def, 
#           instante de inicio de la órbita alrededor de Marte, p_r
@njit
def rk4(step, duration, phi_M_ini, r_orb_M, correction):

    # Parámetros
    delta = G*M_S/(d_SM**3)
    mus = np.array([M_T/M_S, M_M/M_S])

    iters = int(duration//step)

    # Variables 
    r, phi = np.zeros(iters+1), np.zeros(iters+1)
    p_r, p_phi = np.zeros(iters+1), np.zeros(iters+1) 

    k1, k2, k3, k4 = np.zeros(4), np.zeros(4), np.zeros(4), np.zeros(4)

    # Introducimos condiciones iniciales con la velocidad de transferencia
    r[0], phi[0], p_r[0], p_phi[0], delta_v1, v_ini, v_p = condIni(r_orb_T, correction)

    # Calculamos trayectorias Tierra y Marte por adelantado
    r_T, phi_T = planetTraj(d_ST, 0., omegas[0], step, duration)
    r_M, phi_M = planetTraj(d_SM, phi_M_ini, omegas[1], step, duration)    

    # Simulación en el tiempo
    delta_v2 = 0
    in_orbit = False
    for i in range(iters):

        t = i*step  # Instante actual

        # Evaluamos k1
        k1 = evalEquation(r[i], phi[i], p_r[i], p_phi[i], t, step, omegas, delta, mus)

        # Evaluamos k2
        taux = t + step*0.5
        a = r[i] + k1[0]*0.5
        b = phi[i] + k1[1]*0.5
        c = p_r[i] + k1[2]*0.5
        d = p_phi[i] + k1[3]*0.5

        k2 = evalEquation(a, b, c, d, taux, step, omegas, delta, mus)

        # Evaluamos k3
        a = r[i] + k2[0]*0.5
        b = phi[i] + k2[1]*0.5
        c = p_r[i] + k2[2]*0.5
        d = p_phi[i] + k2[3]*0.5

        k3 = evalEquation(a, b, c, d, taux, step, omegas, delta, mus)

        # Evaluamos k4
        taux = t + step
        a = r[i] + k3[0]
        b = phi[i] + k3[1]
        c = p_r[i] + k3[2]
        d = p_phi[i] + k3[3]

        k4 = evalEquation(a, b, c, d, taux, step, omegas, delta, mus)

        # Actualizamos el valor de las variables solución
        r[i+1] = r[i] + (k1[0] + 2.*k2[0] + 2.*k3[0] + k4[0])/6.               
        phi[i+1] = phi[i] + (k1[1] + 2.*k2[1] + 2.*k3[1] + k4[1])/6. 
        p_r[i+1] = p_r[i] + (k1[2] + 2.*k2[2] + 2.*k3[2] + k4[2])/6. 
        p_phi[i+1] = p_phi[i] + (k1[3] + 2.*k2[3] + 2.*k3[3] + k4[3])/6. 


        # Comprobamos si la Nave está lo suficientemente cerca de Marte
        if in_orbit == False:
            
            # Distancia nave-Marte desescalada
            d_NM = (r[i+1]*r[i+1] + r_M[i+1]*r_M[i+1] - 2*r[i+1]*r_M[i+1]*np.cos(phi[i+1] - phi_M[i+1]))**0.5 * d_SM

            # Si la nave está a una distancia fija de Marte, cambiamos su velocidad instantáneamente a la necesaria
            # para orbitar Marte circularmente y orbitar el Sol. También calculamos delta_v2 (desescalado)
            if d_NM < r_orb_M:
                
                p_r[i+1], p_phi[i+1], delta_v2, v_N, v_def = startOrbit(d_NM, r_M[i+1], phi_M[i+1], r[i+1], phi[i+1], p_r[i+1], p_phi[i+1])
                in_orbit = True
                in_orbit_time = t

                """"""""""""""""""""""""
                print("Delta v2 (m/s):")
                print(delta_v2)
                print("Radio de orbita marciana (m): ")
                print(d_NM)
                print("Posicion angular de Marte en el acercamiento (rad):")
                print(phi_M[i+1])
                print("Duracion del viaje (s): ")
                print(t)
                """"""""""""""""""""""""             


    return r, phi, r_T, phi_T, r_M, phi_M, delta_v1, v_ini, v_p, delta_v2, v_N, v_def, in_orbit_time, p_r


# Función que computa la trayectoria orbital de un planeta asumiendo MCU
# Recibe: R, phi_ini, omega, paso, duraión de la simulación
# Devuelve: r, phi
@njit
def planetTraj(R, phi_ini, omega, step, duration):

    iters = int(duration//step)

    # MCU reescalado
    r = np.full(iters, R / d_SM)
    
    phi = np.zeros(iters)
    for i in range(iters):
        t = i*step
        phi[i] = t
    phi *= omega
    phi += phi_ini

    return r, phi


# Función que realiza el cambio de coordenadas polares a cartesianas
# Recibe: r, phi
# Devuelve: x, y
@njit
def polToCart(r, phi):

    x = r*np.cos(phi)
    y = r*np.sin(phi)

    return x, y


# Función que realiza el cambio de velocidades polares a cartesianas
# Recibe: r, phi, vr, omega
# Devuelve: vx, vy
@njit
def polToCartVel(r, phi, vr, omega):

    vx = vr*np.cos(phi) - r*omega*np.sin(phi)
    vy = vr*np.sin(phi) + r*omega*np.cos(phi)

    return vx, vy


# Función que realiza el cambio de coordenadas cartesianas a polares
# Recibe: x, y
# Devuelve: r, phi
@njit
def cartToPol(x, y):

    r = (x*x + y*y)**0.5
    phi = np.arctan2(y, x)

    return r, phi


# Función que realiza el cambio de velocidades cartesianas a polares
# Recibe: x, y, vx, vy
# Devuelve: vr, omega
@njit
def cartToPolVel(x, y, vx, vy):

    vr = (x*vx + y*vy)/(x*x + y*y)**0.5
    omega = (x*vy - vx*y)/(x*x + y*y)

    return vr, omega


"-------------------------------------MAIN------------------------------------------"

# Condiciones iniciales Nave
r_orb_T = R_T + 30000                       # [m] (Cambiar en cada caso)
r_orb_M = 0.1 * (M_M/M_S)**(2/5) * d_SM     # [m] (Fracción de la SOI de Marte) (Cambiar en cada caso)

""""""""""""""""""""""""
print("Radio de la orbita de aparcamiento terrestre (m):")
print(r_orb_T)
""""""""""""""""""""""""

# Condiciones iniciales Tierra y Marte
phi_T_ini = 0
phi_M_ini, a_hoh, T_hoh = hohmann()

# Condiciones de simulación
step = 15               # (Ajustable)
duration = T_hoh*2      # [s] ~ 18 meses
correction = 1.0125     # (Cambiar en cada caso)

# Simulación del movimiento de nave, Tierra y Marte
r, phi, r_T, phi_T, r_M, phi_M, delta_v1, v_ini, v_p, delta_v2, v_N, v_def, in_orbit_time, p_r = rk4(step, duration, phi_M_ini, r_orb_M, correction)
x, y = polToCart(r, phi)

# Trayectoria de la Tierra en cartesianas
x_T, y_T = polToCart(r_T, phi_T)

# Trayectoria de Marte en cartesianas
x_M, y_M = polToCart(r_M, phi_M)

# Energía específica (por unidad de masa) total empleada por el motor en la simulación
delta_E_esp = (np.abs(v_p**2 - v_ini**2) + np.abs(v_def**2 - v_N**2))/2

""""""""""""""""""""""""
print("Delta E esp (J/kg):")
print(delta_E_esp)
""""""""""""""""""""""""

# Escribimos los resultados en formato animacion.py
iters = int(duration//step)

fileName = "resultados.dat"
with open(fileName, "w") as file:

    mostrar = 0
    for i in range(iters):

        # Escribimos las posiciones en el fichero (1 de cada 50)
        if mostrar == 5000:
            file.write(f"{0}, {0}\n") # Sol
            file.write(f"{x_T[i]},  {y_T[i]}\n") # Tierra
            file.write(f"{x_M[i]},  {y_M[i]}\n") # Marte   
            file.write(f"{x[i]},  {y[i]}\n") # Nave

            file.write("\n")

            mostrar = 0
        
        mostrar += 1


# Escribimos las distancias nave-Marte en un fichero aparte para comprobar que ambos quedan ligados.

fileName = "in_orbit_dist.dat"
with open(fileName, "w") as file:

    dist = np.sqrt((x[:-1] - x_M)**2 + (y[:-1] - y_M)**2)

    mostrar = 0
    for i in range(iters):

        if mostrar == 50000:
            file.write(f"{i} {dist[i]}")
            file.write("\n")
            mostrar = 0
        mostrar += 1

# Si escribimos en el fichero los "p_r" en el tiempo (en lugar de "dist"), podemos comprobar que el
# movimiento es oscilatorio (MCU)

fileName = "in_orbit_momentum.dat"
with open(fileName, "w") as file:

    mostrar = 0
    for i in range(iters):

        if i*step >= in_orbit_time:
            if mostrar == 20000:
                file.write(f"{i} {p_r[i]}")
                file.write("\n")
                mostrar = 0
            mostrar += 1


print("Los archivos han sido creados con exito en la localizacion del script.")
