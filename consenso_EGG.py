#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:18:55 2020

@author: Silvia Talavera Marcos

@usage: consenso_EGG.py input_folder output_folder (outputname) (percentage)

@parameters:
input_folder . Folder containing all annotated genomes for the target PCG; they
               must include "annotations" in their filename
output_folder . Where the consensus annotation file is going to be saved
outputname . Output file name (.tsv extension will be added).
percentage . Percentage of files in which an annotation must be present to be
             considered core; 0.80 (80%) by default
"""

import os
import sys
#from carveme.reconstruction.eggnog import load_eggnog_data
from utils_core_functions import split_and_expand, load_eggnog_data
import datetime, time

start = time.time()

# =============================================================================
# Leer todos los modelos
# =============================================================================
# INPUT: consenso_EGG.py input_folder output_folder (outputname) (percentage)
# Only reads files that contain "annotations" in its name
# The path specified for the input_folder must be the realpath

wd = sys.argv[1].rstrip("/")

if len(sys.argv) > 2:
    outputdir = sys.argv[2].rstrip("/") # si el usuario da el outputdir
    
    if len(sys.argv) > 3:
        outputname = sys.argv[3].rstrip("/") # si el usuario da el nombre del output
        if len(sys.argv) > 4: 
            perc = float(sys.argv[4]) # si el usuario da el porcentaje para filtrar
        else:
            perc = 0.80
    else:
        outputname = wd.split("/")[-1]
else:
    outputdir = wd

models = [] # guardamos aquí los nombres de los archivos

for filename in os.listdir(wd):
    if "annotations" in filename: # ignore seed orthologs files
        models.append(wd+"/"+filename)

if len(models) == 0:
    quit("The input file is empty.")
else:
    num_of_models = len(models)


# =============================================================================
# Crear un almacén de todas las reacciones a partir de un modelo cualquiera
# =============================================================================
all_reactions = load_eggnog_data(models[0],drop_unused_cols=False)

# Agrupamos por reacción para evitar repeticiones:
all_reactions = all_reactions.sort_values(by='score', ascending=False) \
                             .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
# Índice para acceder eficientemente a la fila de cada reacción:
reac2idx = {reac:idx for idx,reac in enumerate(all_reactions["BiGG_gene"])}
last_index = len(all_reactions)-1

# =============================================================================
#  Voy abriendo los demás y contando las apariciones de cada reacción
# =============================================================================

# Inicializamos el diccionario de conteo con las reacciones del primer modelo
conteo = {reac:1 for reac in all_reactions["BiGG_gene"]}
for model in models[1:]:
    # abrimos el modelo
    # FIX june 2023; load_eggnog_data might parse a different format depending
    #     on the version. just in case I load here (line 22) a function that
    #     should work
    open_model = load_eggnog_data(model,drop_unused_cols=False)
    open_model = open_model.sort_values(by='score', ascending=False) \
                          .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])

    for i, df_row in open_model.iterrows():
        reaction = df_row["BiGG_gene"]
        
        if reaction not in conteo.keys():
            conteo[reaction] = 1         # la contamos
            all_reactions = all_reactions.append(df_row,ignore_index=True) # la añadimos al almacén
            last_index += 1              # y anotamos su índice
            reac2idx[reaction] = last_index
            
        
        else:
            conteo[reaction] += 1        # la contamos
            # puntuación nueva ponderada (el resultado es la media para esa reacción):
            old_score = all_reactions.at[reac2idx[reaction],"score"]
            new_score = ( old_score *(conteo[reaction]-1) + df_row["score"] ) /conteo[reaction]
            all_reactions.at[reac2idx[reaction],"score"] = new_score # actualizar
    # cerramos el modelo (lo quitamos de memoria)
    del(open_model)
    
    
# =============================================================================
# Filtro y elimino las reacciones que están en menos de un perc% de modelos:
# =============================================================================
for r in reac2idx.keys():
    if conteo[r] < int(perc*num_of_models):
        all_reactions = all_reactions.drop(reac2idx[r],axis=0)
# =============================================================================
# Guardo el archivo de anotaciones consenso
# =============================================================================
if not os.path.exists(outputdir):
    os.mkdir(outputdir)
f = open(outputdir+"/"+outputname+".tsv","a+")

f.write("# emapper version: emapper-2.0.1 emapper DB: 2.0\n")
f.write("# consensus annotation file for "+outputname+"\n")
f.write("# time: "+str(datetime.datetime.now())+"\n")
f.write("#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	best_tax_level	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	taxonomic scope	eggNOG OGs	best eggNOG OG	COG Functional cat.	eggNOG free text desc.\n")
all_reactions.to_csv(f,sep=" ",header=False, index=False)

f.close()

end = time.time()
print("Running time: ",end - start)