"//////////////////////////////////////////////////////////////"
# Programa que grafica el momento radial de la nave una vez
# asentada en la órbita marciana.
# Autor: Luis Gil Martín
# Universidad de Granada
"//////////////////////////////////////////////////////////////"

import matplotlib.pyplot as plt
import numpy as np

file = "in_orbit_momentum.dat"

fig, ax = plt.subplots()
with open(file, "r") as file:

    # Leo el archivo
    lines=file.readlines()

    t = [float(line.split()[0]) for line in lines]
    p = [float(line.split()[1]) for line in lines]
    t=np.array(t)/10**6
    p=np.array(p)/10**(-9)

    t, p = zip(*sorted(zip(t, p)))

    # Graficamos
    ax.plot(t, p, marker = "^", ms = 3, lw = 0.8)

    # Fijamos para la simulación actual (9) cuál fue el paso temporal y el momento de llegada a Marte (cambiar en cada caso)
    h = 15  # [s]
    duration = 22603065 # [s]

    ax.axvline(x=duration/(h*10**6), linestyle="dashed", color = "blue", label="Llegada a Marte")

    ax.set_ylabel(r"Momento radial normalizado $\tilde{p_r} ~(\times 10^{-9})$")
    ax.set_xlabel(r"Tiempo normalizado $t/h ~(\times 10^6)$")

    ax.set_ylim(-5,5)

    ax.legend()

    plt.title("Momento radial de la nave en órbita marciana", fontname = "Arial Rounded MT Bold")
    plt.show()