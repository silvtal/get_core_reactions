#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@description:
    This script generates a list of core functions for a given phylogenetic
    core group (PCG). This list is considered that PCG's pangenome, defined
    as annotations present in at least 90% of available genomes for each PCG.
    The input for this tool is a path containing .tsv annotation files
    generated with eggNOG-mapper, for which a consensus annotation file will
    be generated. The core reactions can be retrieved in KEGG, EC and in BiGG
    format. A file containing descriptions corresponding to EC and KEGG entries
    is also generated.
"""

# Setup =======================================================================

## modules
import os
from utils_core_functions import split_and_expand, load_eggnog_data, consensus_eggnog
from carveme.reconstruction.scoring import *
import pandas as pd
import csv
import Bio.KEGG.REST

nodes = ["Node13821", "Node27828", "Node28853", "Node28866", "Node35562"]

for node in nodes:
    
    # variables ===============================================================
    ## assume wkdir == Example_Goldford
    annotations_folder = "input_annotations/" + node
    perc = 0.90
    pangenome_folder = "output_core_pangenomes"
    if not os.path.exists(pangenome_folder):
        os.makedirs(pangenome_folder)
        
    gprs_bigg = "/home/silvia/bigg_gprs.csv.gz"
    ko_list = "../ko"
    enzyme_dat = "../EC"
    cog_list = "../cog"
    
    # output ==================================================================
    output_folder = pangenome_folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    descriptions_output = output_folder+"/"+node+'_descriptions.txt'
    bigg_reac_output    = output_folder+"/"+node+'_bigg_reac.txt'
    kegg_reac_output    = output_folder+"/"+node+'_kegg_reac.txt'
    kegg_ko_output      = output_folder+"/"+node+'_kegg_ko.txt'
    cog_output = output_folder+"/"+node+'_cogs.txt'
    
    # Obtain consensus annotations ============================================
    consensus_eggnog(annotations_folder, pangenome_folder, outputname=node, perc=perc)
    consensus_file = pangenome_folder + "/" + node + ".tsv" # input from previous line
    
    # Initial parsing =========================================================
    # Parse consensus annotations
    annotation = load_eggnog_data(consensus_file, drop_unused_cols=False, sep="\t")
    
    # take the annotations and filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
    	              .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
    
    # drop duplicates (keep first==best score!)
    bef = len(gene2gene)
    gene2gene = gene2gene.drop_duplicates("query_gene", keep="first")
    print("Removed "+str(bef - len(gene2gene))+"/"+str(bef)+" duplicated entries ('query_gene') from the annotations table")
    
    # Parse BiGG database
    gprs = pd.read_csv(gprs_bigg)
    gprs['BiGG_gene'] = gprs.apply(lambda row: f"{row['model']}.{row['gene'][2:]}", axis=1)
    
    # save available BiGG_Reactions ===========================================
        # This way we can link each annotation to a BiGG reaction
        # en gprs; aquí fusiono columnas en base a su col en común, bigg_gene
        # 'left' significa que me quedo solo con las filas presentes en gene2gene
        # (annots) y añado lo que pueda de la tabla de BiGG. "inner" == solo para las
        # que pueda aportar algo de la tabla de bigg.
        # -> the length of the table can grow: some annotations might have multiple
        # reactions
    df     = pd.merge(gene2gene, gprs, how='left')
    # df=gene_scores_in_gprs.drop_duplicates("query_gene")
    # print("Further removed "+str(bef - len(df))+"/"+str(bef)+" duplicated entries after merging with BiGG GPRs table")
        # Duplicates might actually be two reactions encoded in the same annotation.
        # Like R_MCTP1App and R_MCTP2App
    with open(bigg_reac_output, 'w') as the_file:
        for r in df["reaction"]:
            if not pd.isna(r):
                the_file.write(r+"\n")
    print("Saved BiGG reactions to "+bigg_reac_output+"\n")

    
    
    # save available ECs ======================================================
    my_IDs=list(df[df["EC"].notnull()]["EC"])
    newList=[]
    for i in my_IDs:
        newList.append(i.split(',')[0])
    my_IDs = newList
    
    print("Retrieved ECs: "+str(len(my_IDs))+" (genes with no ECs: "+str(len(gene2gene[gene2gene["EC"].isnull()]["EC"]))+")")
    
    # save available EC descriptions in human language ========================
    import Bio.ExPASy.Enzyme
    handle = open(enzyme_dat)
    descriptions = {record.get("ID"):record for record in Bio.ExPASy.Enzyme.parse(handle)}
    
    with open(descriptions_output, 'w') as the_file:
        for id_ in my_IDs:
            the_file.write(descriptions[id_].get("DE")+"\n")
    print("Saved reaction descriptions to "+descriptions_output+"\n")
    
    # save some more descriptions from KOs
    # http://rest.kegg.jp/list/ko
    with open(ko_list, mode='r') as infile:
        reader = csv.reader(infile,delimiter="\t")
        mydict = {rows[0]:rows[1] for rows in reader}
    
    with open(descriptions_output, 'a') as the_file:
        for keggkos in list(gene2gene[gene2gene["EC"].isnull()]["KEGG_ko"]):
            try:
                for keggko in keggkos.split(","):
                    the_file.write(mydict[keggko]+"\n")
            except:
                print("WARNING: keggkos.split(',') failed. This is keggkos: "+str(keggkos))
    
    print("Saved additional reaction descriptions to "+descriptions_output+"\n")
    
    
    # save all possible kegg reactions ========================================
    with open(kegg_reac_output, 'w') as the_file:
        for rs in df["KEGG_Reaction"]:
            if not pd.isna(rs):
        	    for r in rs.split(","):
        	        ## reactions
        	        the_file.write(r+"\n")
        	        ## compounds
        ##                kegg_compounds = Bio.KEGG.REST.kegg_link("compound", r)
        ##                [the_file.write(c.split("\tcpd:")[-1].strip()+"\n") for c in set(kegg_compounds)]
    
    print("Saved available KEGG_Reactions to "+kegg_reac_output)
        	    
    
    # save KOs ================================================================
    with open(kegg_ko_output, 'w') as the_file:
    #        for rs in list(gene2gene[gene2gene["KEGG_Reaction"].isnull()]["KEGG_ko"]):
        for rs in list(gene2gene[gene2gene["KEGG_ko"].notnull()]["KEGG_ko"]):
            for r in rs.split(","):
                the_file.write(r.split("ko:")[1]+"\n")
            
    print("Saved available KEGG_ko to "+kegg_ko_output)
    
    
    # COG annotations =========================================================
    cog = pd.read_csv(cog_list,
                      sep="\t",
                      header=None)
    cog.columns=["COG","code","description","name","pathway","no","PDB"]
    cog["name"]=[str(i).lower() for i in cog["name"]]
    cog.index=cog.name
    
    ## Assign COGs and save
    with open(cog_output, 'w') as the_file:
        with_code={}
        no_code={}
        for k in list(gene2gene['query_gene'][1:]):
            try: # si ese gen tiene entrada en esta DB de COG
                gene_name=gene2gene[gene2gene['query_gene'] == k]['predicted_gene_name']
                print(gene_name)
                fetched=cog.loc[gene_name.iloc[0].lower()]
                ## if we have synonyms
                if type(fetched) == pd.pandas.core.frame.DataFrame:
                    for i in range(len(fetched)):
                        with_code[k]=gene_name.iloc[0]
                        the_file.write('\t'.join([str(i) for i in list(fetched.iloc[i])])+"\t"+k+"\n")
                else:
                    with_code[k]=gene_name
                    the_file.write('\t'.join([str(i) for i in list(fetched)])+"\t"+k+"\n")
            except:
                no_code[k]=gene_name
    
    ## La mayoría están en cog:
    print("Number of entries present in COG db")
    print(len(with_code.keys()))#329
    print("\nNot present")
    print(len(no_code.keys()))#134
    
