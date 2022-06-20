#Se importa libreria Biopython, que permite leer secuencias geneticas en formato fasta
""" Este Fracmento de codigo se ejecuta solo si es en google colaboratory
try:
    import google.colab
    !pip install biopython
except ImportError:
    pass
"""

#Se importa las librerias numpy y pandas para operar matrices que se obtendran en el algoritmo, adémas de las funciones necesarias de Biopython para leer el archivo tipo fasta
import numpy as np
import pandas as pd
from Bio import Seq
from Bio.Seq import MutableSeq 
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import math
import sys

#Se lee el archivo que contiene las secuencias a analizar y se guardan en la lista 'record' 
#Este algoritmo compara la primera secuencia del archivo con la demás
record = list(SeqIO.parse(sys.argv[1], "fasta"))
"""Solo se ejecuta si es en google colaboratory, eliminando la linea anterior a este mensaje
record = list(SeqIO.parse("/content/prueba4.fasta", "fasta"))"""
def matriz_puntuacion(record):
  ##Guarda la matriz de puntuación de la secuencia con mayor similitud
  matrix_max = 0
  ##Guarda en un diccionario el mayor coeficiente obtenido en la matriz de puntuación para cada secuencia comparada
  coef_simil={}
  ##Guarda la posición del mayor coeficiente en la matriz de puntuación de la secuencia con mayor similitud
  coef_max_posicion = [0,0]
  ##Guarda el mayor coeficiente obtenido en la matriz de puntuación de la secuencia con mayor similitud
  coef_max = 0 
  ##Guarda posición de la secuencia en la lista de secuencias 'record'
  sec_pos = 0

  #para acceder a cada secuencia se hace por medio de un ciclo:
  n = 0
  for r in range(len(record)-1):
      ##Guarda longitud de la secuencia con la que se comparan las demás
      fila = len(record[0]) + 1
      ##Guarda longitud de la secuencia que se comparará con la primera
      columna = len(record[r+1]) + 1
      ##se crea matriz de puntuación con todas sus entradas en cero
      matrix = np.zeros(shape=(fila, columna), dtype=int)
      ##Guarda la posición del mayor coeficiente en la matriz de puntuación
      max_posicion = [0,0]
      ##Guarda el mayor coeficiente obtenido en la matriz de puntuación
      val_max = 0 

      #Este ciclo anidado se encarga de comparar la primera secuencia con la 'r' secuencia por medio del  Algoritmo de Smith–Waterman
      ##i representará el i-ésimo nucleotido de la primera secuencia de ARN 
      for i in range(1,fila):
          ##j representará el j-ésimo nucleotido de la 'r' secuencia de ARN a comparar
          for j in range(1,columna):
              ##Se compara el nucleotido i con cada nucleotido j de la 'r' secuencia y se le asigna una puntuación de acuedo la similitud entre ellos
              if record[0][i-1] == record[r+1][j-1]:###buscar otro metodo(numpy)
                  coincidencia = 2
              else:
                  coincidencia = -1
              ##se calcula la relación de recurrencia que indica el algoritmo de Smith–Waterman
              diagonal=matrix[i-1][j-1] + coincidencia
              arriba = matrix[i-1][j] + (-1)
              izquierda = matrix[i][j-1] + (-1)
              ##se elige la puntuación de similitud con la relación de recurrencia y se guarda en la matriz creada para tal fin (matrix)
              matrix[i][j] = max(diagonal, arriba, izquierda,0)
              ##Se guarda el mayor coeficiente de puntuación y su posición en la matriz de puntuación
              if val_max <= matrix[i][j]:
                  val_max = matrix[i][j]
                  max_posicion[0]=i
                  max_posicion[1]=j
      n +=1
      ###fin del ciclo anidado###
      ##Se guarda el mayor coeficiente de puntuación para cada secuencia en un diccionario, donde el key es la identificacón de la secuencia
      id = record[r+1].id
      coef_simil.setdefault(id,val_max)
      ##Se guarda el con mayor coeficiente de puntuación entre todas las secuencias con su respectiva matriz y posición del coeficiente en esta ultima
      if coef_max <= val_max:
          coef_max = val_max
          coef_max_posicion = max_posicion
          matrix_max = matrix
          sec_pos = n
      ###fin del ciclo de comparación###
  return [coef_max_posicion, sec_pos, matrix_max, coef_simil]

def matriz_recorrido(record, matrix_max, sec_pos):
    ##se crea matriz que indica el camino de similitud para la secuencia con mayor similitud respecto a la primera secuencia
    matrix2ruta = np.zeros(shape=(len(matrix_max),len(matrix_max[0]) ), dtype=int)
    ##se recorre la matriz de maxima similitud siguendo los parametros descritos en el algoritmo de Smith–Waterman para crear la nueva matriz que indica el recorrido
    for i in range(1,len(matrix_max)):
        for j in range(1,len(matrix_max[0])):
            ##se calcula nuevamente la puntuación de coincidencia (solo para la matriz de maxima similitud)
            if record[0][i-1] == record[sec_pos][j-1]:
                coincidencia = 2
            else:
                coincidencia = -1
            ##se calcula la relación de recurrencia que indica el algoritmo de Smith–Waterman
            diagonal=matrix_max[i-1][j-1] + coincidencia
            arriba = matrix_max[i-1][j] + (-1)
            izquierda = matrix_max[i][j-1] + (-1)
            ##se llena la matriz de recorrido
            if matrix_max[i][j] == 0 :
                matrix2ruta[i][j]=0
            elif matrix_max[i][j] == diagonal:
                matrix2ruta[i][j]=1
            elif matrix_max[i][j] == arriba:
                matrix2ruta[i][j]=2
            elif matrix_max[i][j] == izquierda:
                matrix2ruta[i][j]=3
    ###fin del ciclo para llenado de matriz de recorrido###
    return matrix2ruta

def resultados(record, matrix2ruta, coef_max_posicion, coef_simil, sec_pos):
    ##se 'renombra' la lista que contiene la posición del mayor coeficiente de puntuación
    max_posicion1= coef_max_posicion.copy()
    ##se crea el string que guardará la secuencia principal de comparación (al revés, debido a que, por el algoritmo, el camino siempre compara de atras hacia adelante las secuencias)            
    principal=""
    ##se crea el string que guardará la secuencia que fue comparada con la principal y tuvo mayor similitud (al revés, debido a que, por el algoritmo, el camino siempre compara de atras hacia adelante las secuencias)
    comparacion=""
    ##se llenan los strings de acuerdo a la matriz de recorrido odenandolos para que queden alineados
    while matrix2ruta[max_posicion1[0]][max_posicion1[1]] != 0:
        if matrix2ruta[max_posicion1[0]][max_posicion1[1]]== 1:
            principal=principal+record[0][max_posicion1[0]-1]
            comparacion=comparacion+record[sec_pos][max_posicion1[1]-1]
            max_posicion1[0]=max_posicion1[0]-1
            max_posicion1[1]=max_posicion1[1]-1
            
        elif matrix2ruta[max_posicion1[0]][max_posicion1[1]]== 2:
            principal=principal + record[0][max_posicion1[0]-1]
            comparacion=comparacion + '-'
            max_posicion1[0]=max_posicion1[0]-1
            max_posicion1[1]=max_posicion1[1]
            
        elif matrix2ruta[max_posicion1[0]][max_posicion1[1]]== 3:
            principal=principal + '-'
            comparacion=comparacion + record[sec_pos][max_posicion1[1]-1]
            max_posicion1[0]=max_posicion1[0]
            max_posicion1[1]=max_posicion1[1]-1
    ##se crea el string que guardará los simbolos correspondiente al alineamiento entre las secuencias
    alin=""  
    ##Se crea contador que permite saber el total de nucleotidos que se alinean entre ambas secuencias 
    favor=0
    ##Se invierte los strings que indican las secuencias comparadas (al revés)
    comparacion=comparacion[::-1]
    principal=principal[::-1]

    ##se llena la linea que indica el alineamiento entre secuencias comparando las secuencias por medio de un ciclo
    for i in range(len(comparacion)):
      if principal[i]==comparacion[i]:
        ##'|'indica que en ambas cadenas de ARN se encuentra el mismo nucleotido en la misma posición
        alin=alin+"|"
        favor = favor + 1
      else:
        ##' ' indica que no hay similitud en la i-ésima posición de las secuencias
        alin=alin+" "

    ##se convierte el diccionario con los coeficientes de similitud de cada secuencia en un DataFrame
    dataf = pd.DataFrame([[key, coef_simil[key]] for key in coef_simil.keys()], columns=['Id de la secuencia', 'Coef. de similitud'])
    ##Se calcula el porcentaje de similitus entre las secuencias
    P = (favor/len(record[0]))*100
    per = round(P,3)
    ###dividir el alineamiento de la secuenciencia#####
    M=math.ceil(len(principal)/100)
    ##se crea y se escribe en un archivo de texto toda la información analizada en el algoritmo
    doc = open('comparacion_similitud.txt','w')
    doc.write(f'La secuencia principal, con la que se compararon las demás: {record[0].id}')
    doc.write(f'\n')
    doc.write(f'La secuencia con mayor similitud a la secuencia principar fue {record[sec_pos].id} con un')
    doc.write(f'\n')
    doc.write(f'porcentaje de similitud del {per}%')
    doc.write(f'\n')
    doc.write(f'\n')
    doc.write(f'A continuación se muestra el alineamiento entre ambas secuencias:')
    doc.write(f'\n')
    for i in range(M):
        doc.write(f'{principal[100*i:100*(i+1)]}')
        if i == M-1:
          doc.write(f'----{record[0].id}')
        doc.write('\n')
        doc.write(f'{alin[100*i:100*(i+1)]}')
        doc.write('\n')
        doc.write(f'{comparacion[100*i:100*(i+1)]}')
        if i == M-1:
          doc.write(f'----{record[sec_pos].id}')
        doc.write('\n')
    doc.write(f'\n')
    doc.write(f'\n')
    doc.write(f'Además se presenta la siguiente tabla que relaciona cada una de las secuencias comparadas con su')
    doc.write(f'\n')
    doc.write(f'coeficiente de similitud, que corresponde a la puntuación maxima de similitud que tuvo cada secuencia')
    doc.write(f'\n')
    doc.write(f'respecto a la secuencia principal')
    doc.write(f'\n')
    doc.write(f'\n')
    doc.write(f'{dataf}')
    doc.close()

punt = matriz_puntuacion(record)
matr = matriz_recorrido(record, punt[2], punt[1])
resultados(record, matr, punt[0], punt[3], punt[1])