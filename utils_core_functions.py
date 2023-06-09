#!/usr/bin/env python3

import pandas as pd

## functions
# https://raw.githubusercontent.com/silvtal/carveme/7e832e8974fbe6075eb69055a0f1a84f0b9e2314/carveme/reconstruction/eggnog.py
def split_and_expand(df, col, sep):
    split_col = df[col].str.split(sep).apply(pd.Series, 1).stack()
    split_col.index = split_col.index.droplevel(-1)
    split_col.name = col
    df = df.drop(col, axis=1).join(split_col)
    df.reset_index(drop=True, inplace=True)
    return df


def load_eggnog_data(filename, drop_unannotated=True, drop_unused_cols=True, sep=' '):
    """ Load and parse eggnog results for new eggnog-mapper version.

    Args:
        filename (str): input file
        drop_unannotated (bool): remove entries without BiGG annotation (default: True)
        drop_unused_cols (bool): remove columns not used for model carving (default: True)

    Returns:
        pandas.DataFrame: eggnog data

    """
    columns = ['query_gene', 'seed_eggNOG_ortholog', 'evalue', 'score', 'best_tax_level',
                   'predicted_gene_name', 'GO_terms', 'EC', 'KEGG_ko','KEGG_pathways', 
                   'KEGG_Module', 'KEGG_Reaction', 'KEGGrclass',
                   'BRITE*','KEGG_TC','CAZy','BiGG_gene']

    data = pd.read_csv(filename, comment='#', sep=sep, usecols=range(17), names=columns)

    if drop_unannotated:
        data.dropna(subset=['BiGG_gene'], inplace=True)

    if drop_unused_cols:
        data = data[['query_gene', 'BiGG_gene', 'score']]

    data = split_and_expand(data, 'BiGG_gene', ',')

    return data

def consensus_eggnog(wd, outputdir, outputname="", perc=0.8, sep="\t"):
    """
    Created on Fri Nov 20 11:18:55 2020
    
    @author: Silvia Talavera Marcos
    
    @parameters:
    input_folder . Folder containing all annotated genomes for the target PCG;
                   they must include "annotations" in their filename
    output_folder . Where the consensus annotation file is going to be saved
    outputname . Output file name (.tsv extension will be added).
    perc . Percentage of files in which an annotation must be present to be
                 considered core; 0.80 (80%) by default
                 
    @original (script call): for NODE in $(ls $INPUTDIR/models); do ./consenso_EGG.py annotations_2022_12_14/$NODE"_annots" annotations_2022_12_14/ pan$NODE 0.9; done

    """
    
    import os
    import sys
    #from carveme.reconstruction.eggnog import load_eggnog_data
    #from utils_core_functions import load_eggnog_data
    import datetime, time
    
    start = time.time()
    
    # =============================================================================
    # Leer todos los modelos
    # =============================================================================
    # INPUT: consenso_EGG.py input_folder output_folder (outputname) (percentage)
    # Only reads files that contain "annotations" in its name
    # The path specified for the input_folder must be the realpath
    
    wd = wd.rstrip("/")
    outputdir = outputdir.rstrip("/")
    if outputname == "":
        outputname = wd.split("/")[-1]
    perc = float(perc)
    
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
    all_reactions = load_eggnog_data(models[0],drop_unused_cols=False, sep="\t")
    
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
        open_model = load_eggnog_data(model,drop_unused_cols=False, sep="\t")
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
        os.makedirs(outputdir)
    f = open(outputdir+"/"+outputname+".tsv","a+")
    
    f.write("# consensus annotation file for "+outputname+"\n")
    f.write("# time: "+str(datetime.datetime.now())+"\n")
    f.write("#query_name	seed_eggNOG_ortholog	seed_ortholog_evalue	seed_ortholog_score	best_tax_level	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	taxonomic scope	eggNOG OGs	best eggNOG OG	COG Functional cat.	eggNOG free text desc.\n")
    all_reactions.to_csv(f,sep="\t",header=False, index=False)
    
    f.close()
    
    end = time.time()
    print("Running time: ",end - start)