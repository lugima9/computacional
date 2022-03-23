import numpy as np

def leerCondIni(nombreFichero):
     # Abrimos el fichero en modo lectura
     f = open(nombreFichero, "r")

     datos = f.read()

     f.close()

     return datos

datos = leerCondIni("datos_planetas.txt")

print(datos)

