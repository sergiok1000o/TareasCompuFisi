#- coding: latin-1 -*-
import sys, string, os
import numpy as np
libro = sys.argv[1];
letras = [""]*1000
repeticion = [0.00]*1000
f = [0.00]*20
r = [0.00]*20
k=0;
i=0
l=1
total=0;
prueba=0
def introducir_letra(let):
		k=0
		cambiarletra = 0
		for a in letras:
			if(letras[k] == "" and cambiarletra == 0):
				letras[k] = let
				repeticion[k] = 1
				cambiarletra = 1
			k=k+1
def sumar(reng):
	encontro = 0;
	for a in reng:
		k=0
		for b in letras:
			if(a == b):
				repeticion[k] = repeticion[k]+1
				encontro = 1
			k=k+1
		if(encontro == 0):
			introducir_letra(a)
		encontro=0
libro1 = open(libro, 'r')
for x in libro1:
       sumar(x)
arreglo_final = [repeticion,letras]

from operator import itemgetter
inv = zip(*arreglo_final)
inv.sort(reverse = 1,key=itemgetter(0))
arreglo_final = zip(*inv)
nombrearchivo = sys.argv[1]
final = open("frecuencia_"+nombrearchivo,"w")
final.write("El siguiente archivo cuenta cuantas veces se repite cada letra en el texto ingresado\n\n       Caracter                           Frecuencia\n")
for a in repeticion:
	o= arreglo_final[0][i]
	p= arreglo_final[1][i]
	if(o != "0.0" and  p != " " and  p != "\n" and  p != "." and  p != "," and  p != ":" and  p != ";" and  p != "!" and  p != "?" and  p != "¿" and p != ":" and  p != "(" and  p != ")" and p != "±" and p != "_" and p != "³" and p != "" and p != "/" and p != "=" and p != "»" and p != "Œ" and p != "#" and p != "'" and p != "" and p != "" and p != "" and p != "Â" and p != "¡" and p != "	"):
		total=total+o	
	i = i+1
for a in repeticion:
	o= arreglo_final[0][k]
	o= (o/total)*100
	prueba=prueba+o
	p= arreglo_final[1][k]	
	o = str(o)
	letra_repeticion = "         " + p +"                                 "+ o + " %\n"
	if(o != "0.0" and  p != " " and  p != "\n" and  p != "." and  p != "," and  p != ":" and  p != ";" and  p != "!" and  p != "?" and  p != "¿" and p != ":" and  p != "(" and  p != ")" and p != "±" and p != "_" and p != "³" and p != "" and p != "/" and p != "=" and p != "»" and p != "Œ" and p != "#" and p != "'" and p != "" and p != "" and p != "" and p != "Â" and p != "¡" and p != "	"):
				final.write(letra_repeticion)
	k = k+1



#desde aqui comienza el quiz :

for a in repeticion:
	if(l != 20):
		o= arreglo_final[0][l]
		o= (o/total)
		o=np.log(o)
		f[l-1]=o	
		r[l-1]=np.log(l)
		l=l+1
print f
print r
z=np.polyfit(r,f,1)
final.write("\nEste libro tiene un k de valor\n")

z = str(z[0])

final.write(z)
