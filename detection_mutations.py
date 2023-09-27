# from io import TextIOWrapper
import os 
import re
import csv
from Bio import Align
from typing import Literal
from Utils import extract_number
from Bio.Seq import Seq

# from Bio.Align import Alignment
# import numpy as np


# TODO : automatiser la génération d'un fasta avec la ref et les sequences étudiées à partir d'un fichier dont je n'ai pas le nom.
#====================================================    Resume    ====================================================#

# TODO!: faire en sorte d'extraire la variable duplication du nom des fastas G afin de supprimer la variable globale "dulication"
#================================================= Variables Globales =================================================#
chemin = os.getcwd()
disque = chemin[0]


# path pour la vérification/création du chemin de dossiers ci-dessous
path = f"./../10-test_positions/results"

# path pour accéder aux dossiers fasta A et B
path_A = f"./../3-Fasta_Sequences_Prots/Fasta_A"
path_B = f"./../3-Fasta_Sequences_Prots/Fasta_B"


# variable globale à changer qui permet de déterminer si les séquences de références vont comporter la duplication au niveau de la protéine G
duplication = True


# path pour les sequences de référence des types A et B
path_ref = f"./../references_phylogenie/ref_combined_insertion.fasta"
path_ref_without_duplication = f"./../references_phylogenie/RSV_ref.fasta"




# Directory:
# resultat_mutations,result_proteine_type_A,result_proteine_type_B,result_final_A,result_final_B
resultat_mutations = f"./../10-test_positions/results_mutations"
result_proteine_type_A = f"./../10-test_positions/result_proteines/results_protéines_type_A"
result_proteine_type_B = f"./../10-test_positions/result_proteines/results_protéines_type_B"
result_final_A = f"./../10-test_positions/result_proteines_final_A"
result_final_B = f"./../10-test_positions/result_proteines_final_B"
csv_final = f"./../10.1-dossier_xlsx_result_mutations"



# path à supprimer :
# path_a_delete = [result_proteine_type_A,result_proteine_type_B]




directories = [resultat_mutations,result_proteine_type_A,result_proteine_type_B,result_final_A,result_final_B,csv_final]

#================================================= fonctions =================================================#

# ouvre le fichier Path_A
#si le fichier finir par un G, 
# alors ouvre tous les dossiers 1 à 1 
# recherche et print le numéro du barcode
# pour chacune des lignes du fichier: 
# si c'est la première ligne, alors l'id est récupéré comme clé et stockée dans une liste de clés 
# si c'est la seconde ligne alors elle est enregistrée comme séquence dans une liste 

#TODO! vérifier si des bouts de code sont inutiles !!!!!!!

def extraction_regex(inputed_string: str, var_str : str, expression_regex: str)-> str:

    return var_str

def extraction_id(file_name: str,file_id: str)-> str:

    file_id_pattern_regex = r'result_([^_]+)'
    # mutation_boolean_regex = r'True'

    file_id_match = re.search(file_id_pattern_regex, file_name)

    if file_id_match:
        file_id = file_id_match.group(1)
    else:
        file_id = "Fasle"        
    return file_id

def extraction_boolean_from_file_name(file_name: str,boolean_presence: str)-> str:

    if "True" in file_name:
        boolean_presence = "True"
    else:
        boolean_presence = "False"

    return boolean_presence


def append_sequences_prot(path_A:str, path_B:str)-> tuple[list, list, list, list, list, list, list, list, dict[str,bool], dict[str,bool]]:
    
    clefs_F_A = []
    clefs_F_B = []
    clefs_G_A = []
    clefs_G_B = []
    duplication_G_A = {}
    duplication_G_B = {}
    sequences_bc_prot_G_A = []
    sequences_bc_prot_G_B = []
    sequences_bc_prot_F_A = []
    sequences_bc_prot_F_B = []

    # chacun des pour les dossiers dans 3-Fasta_Sequences_Prots
    for dossiers in os.listdir(path_A):

        # pour les dossiers se terminant par un G (soit: ceux contenant la protéine G)
        if dossiers.endswith("G"):

            folder_path = os.path.join(path_A, dossiers)

            for fichier in os.listdir(folder_path):

                file_path = os.path.join(folder_path, fichier)  # Chemin du fichier

                with open(file_path, 'r') as file:

                    presence_duplication = re.search("True",file_path)
                    if presence_duplication:

                        file_id=""
                        file_id=extraction_id(file_path,file_id)
                       
                        duplication_G_A[f"{file_id}"]=True
                    else:
                        duplication_G_A[f"{file_id}"]=False
                    
                    i = 0
                    for lines in file:
                        if i == 0:
                            lines = lines.rstrip("\n")
                            clefs_G_A.append(f"{lines[1:-2]}")

                        elif i == 1:
                            sequences_bc_prot_G_A.append(Seq(lines))
                        
                        i += 1
                    
        elif dossiers.endswith("F"):
            folder_path = os.path.join(path_A, dossiers)  # Chemin du sous-dossier F

            for fichier in os.listdir(folder_path):
                file_path = os.path.join(folder_path, fichier)  # Chemin du fichier

                with open(file_path, 'r') as file:
                    
                    i = 0

                    for lines in file:
                        if i == 0:
                            lines = lines.rstrip("\n")
                            clefs_F_A.append(f"{lines[1:-2]}")
                             
                        elif i == 1:
                            sequences_bc_prot_F_A.append(Seq(lines))
                        
                        i += 1


    for dossiers in os.listdir(path_B):

        if dossiers.endswith("G"):
            folder_path = os.path.join(path_B, dossiers)  # Chemin du sous-dossier G

            for fichier in os.listdir(folder_path):
                file_path = os.path.join(folder_path, fichier)  # Chemin du fichier

                with open(file_path, 'r') as file:
                    presence_duplication = re.search("True",file_path)

                    file_id=""
                    file_id=extraction_id(file_path,file_id)

                    if presence_duplication:
                        duplication_G_B[f"{file_id}"]=True
                        
                    else:
                        duplication_G_B[f"{file_id}"]=False
                    # correspondance = re.search(r"bc(\d+)", file_path)
                    # if correspondance:
                    #     expression = correspondance.group(0)+"_B"
                    i = 0

                    for lines in file:
                        if i == 0:
                            lines = lines.rstrip("\n")
                            clefs_G_B.append(f"{lines[1:-2]}")

                        elif i == 1:
                            sequences_bc_prot_G_B.append(Seq(lines))
                        
                        i += 1

        elif dossiers.endswith("F"):
                folder_path = os.path.join(path_B, dossiers)  # Chemin du sous-dossier F
                for fichier in os.listdir(folder_path):
                    file_path = os.path.join(folder_path, fichier)  # Chemin du fichier
                    with open(file_path, 'r') as file:

                        i = 0

                        for lines in file:
                            if i == 0:
                                lines = lines.rstrip("\n")
                                clefs_F_B.append(f"{file_path[62:66]}")

                            elif i == 1:
                                sequences_bc_prot_F_B.append(Seq(lines))
                            
                            i+=1
    
    return clefs_G_A, clefs_G_B, clefs_F_A, clefs_F_B, sequences_bc_prot_G_A, sequences_bc_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, duplication_G_A, duplication_G_B 


def extraction_sequence_reference(path_ref_with_duplication:str, path_ref_without_duplication:str)-> tuple[Seq,Seq,Seq,Seq]:


    with open(path_ref_with_duplication,"r") as ref_file_with_duplication:   
        

        table1 = []
        sequence = {}
        stockage = ""
        stockage_sequence = ""
        marqueurA = False
        marqueurB = False
        i = 0
        
        for elem in ref_file_with_duplication:
            i += 1
            if marqueurA and not elem.startswith(">typeB"):
                
                ligne = (elem.lower()).rstrip("\n")

                stockage_sequence += ligne            
                
            elif elem.startswith(">typeA_"):
                
                stockage = elem.rstrip("\n")
                
                marqueurA = True
                

            elif marqueurB == True and not elem.startswith(">typeB"): 
                
                ligne = (elem.lower()).rstrip("\n")
                stockage_sequence += ligne 

            
            elif elem.startswith(">typeB"):
                
                sequence[stockage] = stockage_sequence
                table1.append(sequence)
                sequence = {}
                stockage_sequence = ""
                stockage = elem.rstrip("\n")


                marqueurA = False
                marqueurB = True 


        sequence[stockage] = stockage_sequence
        table1.append(sequence)
        

        trousseau = []
        for element in table1:
            trousseau.append(next(iter(element.keys())))
        
        
    ref_A = Seq((table1[0])[trousseau[0]])
    ref_B = Seq((table1[1])[trousseau[1]])
    


    with open(f"{path_ref_without_duplication}","r") as ref_file_without_duplication:
        
        table = []
        sequence = {}
        stockage = ""
        stockage_sequence = ""
        marqueurA = False
        marqueurB = False
        i = 0
        
        for elem in ref_file_without_duplication:
            i += 1

            if marqueurA and not elem.startswith(">typeB"):
                
                ligne = (elem.lower()).rstrip("\n")

                stockage_sequence += ligne            
                
            elif elem.startswith(">typeA_"):
                
                stockage = elem.rstrip("\n")
                
                marqueurA = True
                

            elif marqueurB == True and not elem.startswith(">typeB"): 
                
                ligne = (elem.lower()).rstrip("\n")
                stockage_sequence += ligne 

            
            elif elem.startswith(">typeB"):
                
                sequence[stockage] = stockage_sequence
                table.append(sequence)
                sequence = {}
                stockage_sequence = ""
                stockage = elem.rstrip("\n")


                marqueurA = False
                marqueurB = True 

        sequence[stockage] = stockage_sequence
        table.append(sequence)
        trousseau = []

        for element in table:
            trousseau.append(next(iter(element.keys())))
        table.append(table)

        
    ref_A_without_duplication = Seq((table[0])[trousseau[0]])
    ref_B_without_duplication = Seq((table[1])[trousseau[1]])

    return ref_A, ref_B, ref_A_without_duplication , ref_B_without_duplication

# TODO! détécter automatiquemetn si la dulication est présente ou non et supprimer cette variable hardcodée
def extraction_proteines_ref(ref_A: Seq, ref_B: Seq, ref_A_without_duplication:Seq , ref_B_without_duplication:Seq)-> tuple[list[Seq], list[Seq], list[Seq], list[Seq]]:
    # , duplication_G_A:dict[bool, bool], duplication_G_B:dict[bool, bool]

    #TODO!: A
    proteines_A_without_duplication = []
    proteine_A_dupli = []
    duplication_A = False
    
    # coordonnées protéine G chez virus type A
    
    start_prot_G_type_A = 4658
    end_prot_G_type_A = 5555

    start_prot_F_type_A = 5631
    fin_prot_F_type_A = 7356

    longueur_duplication_A = 0
    stockage_sequence_G_A_reference = ref_A_without_duplication[int(start_prot_G_type_A): int(end_prot_G_type_A)+longueur_duplication_A]


    # protéine G de A
    proteines_A_without_duplication.append(stockage_sequence_G_A_reference) #]4659;5554] alors nombre = multiple de 3
    
    # protéine F de F
    # proteines_A.append(ref_A[start_prot_F_type_A+longueur_duplication_A:fin_prot_F_type_A+longueur_duplication_A]) # ]5633;7355] alors multiple de 3
    proteines_A_without_duplication.append(ref_A[start_prot_F_type_A:fin_prot_F_type_A])


    
    # duplication
    longueur_duplication_A = 72
    stockage_sequence_G_A_reference_dupli = ref_A[int(start_prot_G_type_A): int(end_prot_G_type_A)+longueur_duplication_A]

    
    # protéine G de A
    proteine_A_dupli.append(stockage_sequence_G_A_reference_dupli) #]4659;5554] alors nombre = multiple de 3
    
    # protéine F de F
    # proteine_A_dupli.append(ref_A[start_prot_F_type_A+longueur_duplication_A:fin_prot_F_type_A+longueur_duplication_A]) # ]5633;7355] alors multiple de 3

    proteine_A_dupli.append(ref_A[start_prot_F_type_A+longueur_duplication_A:fin_prot_F_type_A+longueur_duplication_A])
    
    proteines_B_without_duplication = []
    proteine_B_dupli = []
    duplication_B = False

    start_prot_G_type_B = 4687
    fin_prot_G_type_B = 5566

    start_prot_F_type_B = 5663
    fin_prot_F_type_B = 7388

    

    # for elem in duplication_G_B:
    #     if duplication_G_B[elem]:
    #         duplication_B = True 

    # if duplication_B:
        
    longueur_duplication_B = 60

    stockage_sequence_G_B_reference_dupli = ref_B[int(start_prot_G_type_B): int(fin_prot_G_type_B)+longueur_duplication_B]

    # protéine G avec la duplication
    proteine_B_dupli.append(stockage_sequence_G_B_reference_dupli) #]4659;5554] alors nombre = multiple de 3
    # protéine F 
    proteine_B_dupli.append(ref_B[start_prot_F_type_B + longueur_duplication_B:fin_prot_F_type_B + longueur_duplication_B]) # ]5633;7355] alors multiple de 3


    longueur_duplication_B = 0
    

    stockage_sequence_G_B_reference = ref_B_without_duplication[int(start_prot_G_type_B): int(fin_prot_G_type_B)+longueur_duplication_B]
    
    # protéine G avec la duplication
    proteines_B_without_duplication.append(stockage_sequence_G_B_reference) #]4659;5554] alors nombre = multiple de 3
    # protéine F 
    proteines_B_without_duplication.append(ref_B[start_prot_F_type_B + longueur_duplication_B:fin_prot_F_type_B + longueur_duplication_B]) # ]5633;7355] alors multiple de 3

    #* ça c'est good c'est ce que l'on cherche à input (enfin la séquence semble bonne, rien à voir avec la longueur?)
    return proteines_A_without_duplication, proteine_A_dupli, proteines_B_without_duplication, proteine_B_dupli

#TODO*: extraction_proteines_ref le problème n'a pas l'air de venir d'ici car saltation{proteine_A_dupli[1].translate()} donne bien le résultat recherché
#TODO*: traduction() la fonction n'a pas l'air de poser problème, car la sélection du résultat de la vérification du dessus est en surbrillance quand je sélectionne celui là
#TODO*: génération alignement()
    #*: premier problème détecté, les variables n'étaient pas énoncées dans le bon ordre.
    #* : tout à l'aire de fonctionner correctement jusqu'avant la fonction alignement.
    #*: alignement()
        #*: les input d'alignement on l'air d'être corrects.
        #* le result d'alignement() est correct
#todo: génération résult():
    #todo: extraction_txt()

#TODO!: ##########################################################################################################
#TODO!: ##########################################################################################################
#TODO!: ##########################################################################################################
#TODO!: ## problème avec la séquence dupliquée qui ne correspond #################################################
#TODO!: ###################toujours pas avec la séquence étudiée mais uniquement sur les txt v2.##################
#TODO!: ##########################################################################################################
#TODO!: ##########################################################################################################
#TODO!: ##########################################################################################################


"""
def transcription(
    sequences_prot_G_A: list,
    sequences_prot_G_B: list, 
    sequences_bc_prot_F_A: list, 
    sequences_bc_prot_F_B: list, 
    proteines_A:    list,
    proteines_B:    list,)-> tuple[list,list,list,list,list,list]:
    
    #clefs_G_A, clefs_G_B, clefs_F_A, clefs_F_B
    
    for i in range(len(sequences_prot_G_A)):
        sequences_prot_G_A[i] = sequences_prot_G_A[i].transcribe()

    for i in range(len(sequences_prot_G_B)):
        
        sequences_prot_G_B[i] = sequences_prot_G_B[i].transcribe()
        
    for i in range(len(sequences_bc_prot_F_A)):
        sequences_bc_prot_F_A[i] = sequences_bc_prot_F_A[i].transcribe()
        
    for i in range(len(sequences_bc_prot_F_B)):
        sequences_bc_prot_F_B[i] = sequences_bc_prot_F_B[i].transcribe()

    # references
    for i in range(len(proteines_A)):
        proteines_A[i] = proteines_A[i].transcribe()
        
    for i in range(len(proteines_B)):
        proteines_B[i] = proteines_B[i].transcribe()

    return sequences_prot_G_A, sequences_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, proteines_A, proteines_B
"""

def traduction_sequence(i,sequences_prot_X_Y):
    sequence_traduite = sequences_prot_X_Y[i].translate()
    return sequence_traduite


def traduction(
    proteines_A:list,
    proteines_B:list,
    sequences_bc_prot_G_A:list,
    sequences_bc_prot_G_B:list,
    sequences_bc_prot_F_A:list,
    sequences_bc_prot_F_B:list,
    proteine_A_dupli:list,
    proteine_B_dupli:list
    

    )-> tuple[list,list,list,list,list,list,list,list]:

    Table_sequence_proteines_G_A = []
    Table_sequence_proteines_G_B = [] 
    Table_sequence_proteines_F_A = []
    Table_sequence_proteines_F_B = []

    for i in range(len(sequences_bc_prot_G_A)):
        sequence_traduite = traduction_sequence(i,sequences_bc_prot_G_A)
        Table_sequence_proteines_G_A.append(sequence_traduite)
        
    for i in range(len(sequences_bc_prot_G_B)):
        sequence_traduite = traduction_sequence(i,sequences_bc_prot_G_B)
        Table_sequence_proteines_G_B.append(sequence_traduite)
        
    for i in range(len(sequences_bc_prot_F_A)):
        sequence_traduite = traduction_sequence(i,sequences_bc_prot_F_A)
        Table_sequence_proteines_F_A.append(sequence_traduite)
        
    for i in range(len(sequences_bc_prot_F_B)):
        sequence_traduite = traduction_sequence(i,sequences_bc_prot_F_B)
        Table_sequence_proteines_F_B.append(sequence_traduite)
        
    for i in range(len(proteines_A)):
        proteines_A[i] = proteines_A[i].translate()
        
    for i in range(len(proteines_B)):
        proteines_B[i] = proteines_B[i].translate()

    for i in range(len(proteine_A_dupli)):
        proteine_A_dupli[i] = proteine_A_dupli[i].translate()
        

    for i in range(len(proteine_B_dupli)):
        proteine_B_dupli[i] = proteine_B_dupli[i].translate()
        
    return Table_sequence_proteines_G_A, Table_sequence_proteines_G_B, Table_sequence_proteines_F_A, Table_sequence_proteines_F_B, proteines_A, proteines_B, proteine_A_dupli, proteine_B_dupli


def analyse_sequence_A_A(i, proteine_X:Seq, proteine_X_ref:Seq)-> list:

    results = []
    
    y = 0
    
    for y in range(max(len(proteine_X_ref),len(proteine_X))):
        if proteine_X[y] != proteine_X_ref[y]:
            
            results.append([f'{proteine_X_ref[y]}{y+1}{proteine_X[y]}'])


    return results


def extract_number(string)-> int:
    # Use regular expression to find a number in the string
    match = re.search(r'\d+', string)
    
    if match:
        # Extract the matched number as a string
        number_str = match.group()
        
        # Convert the string number to an integer or float if needed
        number = int(number_str)  # Use int() for integer
        # number = float(number_str)  # Use float() for floating-point number
        return number
    else :
        return 0


def extract_number_and_letter(element):
    number = ""
    letter = ""
    for char in element:
        if char.isdigit():
            number += char
        elif char.isalpha():
            letter += char
    return int(number), letter

"""
def liste_mutations_prot(proteines_A,sequences_bc_prot_G_A):
                
    mutations = []
    output_mutations = " |  0- 99|"
    # Parcourir chaque position de l'alignement
    for i in range(len(proteines_A)):
        target_aa = proteines_A[i]
        query_aa = sequences_bc_prot_G_A[i]

        # Vérifier si les acides aminés sont différents (mutation)
        if target_aa != query_aa:
            mutation = f"{target_aa}{i+1}{query_aa}"
            if mutation[-1] != "X":
                mutations.append(mutation)

    x = 100
    y = 0
    while y < len(mutations):
        if extract_number(mutations[y]) < x:
            if len(mutations[y]) == 3:
                output_mutations += f"   {mutations[y]}"
            elif len(mutations[y]) == 4:
                output_mutations += f"  {mutations[y]}"
            else:
                output_mutations += f" {mutations[y]}"
            y+=1
        else:
            output_mutations += f"\n |{x}-{x+99}|"
            x += 100
    return output_mutations,mutations
"""
#TODO! : destruction par l'atome des mutation comprennant un X
def alignement(target:str, query:str) :
    aligner = Align.PairwiseAligner()
    
    aligner.mode = 'global'  # Mode d'alignement global
    aligner.match_score = 2  # Score pour une correspondance
    aligner.mismatch_score = -1  # Score pour une non-correspondance
    aligner.open_gap_score = -2  # Score pour l'ouverture d'un espace
    aligner.extend_gap_score = -0.5  # Score pour l'extension d'un espace


    alignments = aligner.align(target,query)
    best_alignment = alignments[0]
    differences = []
    
    for i in range(len(best_alignment.target)):
        target_char = best_alignment.target[i]
        query_char = best_alignment.query[i]

        if target_char != "-" and query_char != "-" and target_char != query_char:
            differences.append((f"{target_char}{i}{query_char}"))

    
    return best_alignment


def generation_alignement_txtv1(
    proteines_A:list,
    proteine_A_dupli:list,
    proteines_B:list,
    proteine_B_dupli:list,
    sequences_bc_prot_G_A:list,
    sequences_bc_prot_G_B:list, 
    sequences_bc_prot_F_A:list, 
    sequences_bc_prot_F_B:list,
    clefs_G_A:list[str], 
    clefs_G_B:list[str], 
    clefs_F_A:list[str], 
    clefs_F_B:list[str],
    duplication_G_A:dict[str,bool], # dictionnaire comprennant les booleans de présences ou non de la duplication au niveau de la protéine G avec les ids pour clés. 
    duplication_G_B:dict[str,bool], #idem
    ):
    
    # f"{clefs_G_A[y]},{protein},{subunit},{mutation},Fasle,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab}\n"

    for y in range(len(sequences_bc_prot_G_A)):
        csv_ref = open(f"./../references_phylogenie/Mutations_RSVA.csv","r")
        
        with open(f"./../10-test_positions/result_proteines/results_protéines_type_A/result_{clefs_G_A[y].replace('>', '_')}.txt","a") as sortie_A:
            
            if duplication_G_A[clefs_G_A[y]] == False:
                sortie_A.write(f"Alignement de la protéine G de {clefs_G_A[y]} avec la protéine G de référence\n\n")
                alignement_proteine = alignement(proteines_A[0],sequences_bc_prot_G_A[y])
                sortie_A.write(f"{alignement_proteine}")

                # mutations_G_A_output,mutations_G_A = liste_mutations_prot(proteines_A[0],sequences_bc_prot_G_A[(y)])
                sortie_A.write(f"Alignement de la protéine F de {clefs_G_A[y]} avec la protéine F de référence\n\n")
                alignement_proteine = alignement(proteines_A[1],sequences_bc_prot_F_A[(y)])
                sortie_A.write(f"{alignement_proteine}")

            else:

                sortie_A.write(f"Alignement de la protéine G de {clefs_G_A[y]} avec la protéine G de référence\n\n")
                alignement_proteine = alignement(proteine_A_dupli[0],sequences_bc_prot_G_A[y])
                sortie_A.write(f"{alignement_proteine}")

                sortie_A.write(f"Alignement de la protéine F de {clefs_G_A[y]} avec la protéine F de référence\n\n")
                alignement_proteine = alignement(proteine_A_dupli[1],sequences_bc_prot_F_A[y])
                sortie_A.write(f"{alignement_proteine}")


            

            # sortie_A.write(f"Alignement de la protéine G du {clefs_G_A[y]} avec la protéine G de référence\n\n")
            # sortie_A.write(f"{alignement_proteine}")#"\n{differences_G_A}\n \n")

            # alignement_proteine,differences_F_A = alignement(proteines_A[1],sequences_bc_prot_F_A[(y)])
            # mutations_F_A_output,mutations_F_A = liste_mutations_prot(proteines_A[1],sequences_bc_prot_F_A[(y)])

            # sortie_A.write(f"Alignement de la protéine F du {clefs_G_A[y]} avec la protéine F de référence\n\n")
            # sortie_A.write(f"{alignement_proteine}")#\n{differences_F_A}\n \n")
 
    for z in range(len(sequences_bc_prot_G_B)):
        with open(f"./../10-test_positions/result_proteines/results_protéines_type_B/result_{clefs_G_B[z].replace('>', '_',)}.txt","a") as sortie_B:
    
            if duplication_G_B[clefs_G_B[z]] == False:
                sortie_B.write(f"Alignement de la protéine G de {clefs_G_B[z]} avec la protéine G de référence\n\n\n")
                alignement_proteine= alignement(proteines_B[0],sequences_bc_prot_G_B[z])
                sortie_B.write(f"{alignement_proteine}")

                # mutations_G_A_output,mutations_G_A = liste_mutations_prot(proteines_A[0],sequences_bc_prot_G_A[(y)])
                sortie_B.write(f"Alignement de la protéine F de {clefs_G_B[z]} avec la protéine F de référence\n\n")
                alignement_proteine= alignement(proteines_B[1],sequences_bc_prot_F_B[(z)])
                sortie_B.write(f"{alignement_proteine}")
            else:

                sortie_B.write(f"Alignement de la protéine G de {clefs_G_B[z]} avec la protéine G de référence\n\n")
                alignement_proteine= alignement(proteine_B_dupli[0],sequences_bc_prot_G_B[z])
                sortie_B.write(f"{alignement_proteine}")

                sortie_B.write(f"Alignement de la protéine F de {clefs_G_B[z]} avec la protéine F de référence\n\n")
                alignement_proteine= alignement(proteine_B_dupli[1],sequences_bc_prot_F_B[(z)])
                sortie_B.write(f"{alignement_proteine}")
                    
    return


def changement_pipe_en_hashtag(sequence_target:str,sequence_query:str,new_sequence:str,sequence_symboles:str,listing_mutation:list,query_number:int) -> tuple[str,list]:

    for z in range(len(sequence_target)):

        if sequence_query[z]!="X":
            new_sequence = new_sequence + sequence_symboles[z]
        else:
            new_sequence += "#"
            
        if sequence_query[z]!="X" and sequence_target[z] != sequence_query[z] and sequence_query[z] != "-" and sequence_target[z] != "-":
            listing_mutation.append(f"{sequence_target[z]}{z+query_number+1}{sequence_query[z]}")
            
    return new_sequence, listing_mutation
       
       
def add_spaces(number)->str:
                    
    spaced_number = str(number)
    
    if len(spaced_number) == 1:
        spaced_number = f"  {spaced_number}"
    elif len(spaced_number) == 2:
        spaced_number = f" {spaced_number}"
        
        
    return spaced_number


def selection_mutations_pour_csv(liste_mutation_à_ajouter_au_csv):
            
            resistance_palivizumab_boolean = False
            resistance_nirsevimab_boolean = False
            resistance_suptavumab_boolean = False
            
            type_resistance_palivizumab = ""
            type_resistance_nirsevimab = ""
            type_resistance_suptavumab = ""
            
            resistance_palivizumab = []
            resistance_nirsevimab = []
            resistance_suptavumab = []
            
            # résistant "R"
            for elem in liste_mutation_à_ajouter_au_csv:
                if elem["palivizumab"] =="R":
                    resistance_palivizumab.append(elem["mutation"]) 
                if elem["nirsevimab"] == "R":
                    resistance_nirsevimab.append(elem["mutation"])
                if elem["suptavumab"] == "R":
                    resistance_suptavumab.append(elem["mutation"])
                    
            if  len(resistance_palivizumab) > 0:
                resistance_palivizumab_boolean = True 
                type_resistance_palivizumab = "R"
                
            if  len(resistance_nirsevimab) > 0:
                resistance_nirsevimab_boolean = True
                type_resistance_nirsevimab = "R"
                
            if  len(resistance_suptavumab) > 0:
                resistance_suptavumab_boolean = True
                type_resistance_suptavumab = "R"
                       
            if resistance_palivizumab_boolean == False or \
               resistance_nirsevimab_boolean == False or \
               resistance_suptavumab_boolean == False :
                   
                # intermédiaire "I"
                for elem in liste_mutation_à_ajouter_au_csv:
                    if elem["palivizumab"] =="I":
                        resistance_palivizumab.append(elem["mutation"]) 
                    if elem["nirsevimab"] == "I":
                        resistance_nirsevimab.append(elem["mutation"])
                    if elem["suptavumab"] == "I":
                        resistance_suptavumab.append(elem["mutation"])
                
                if  len(resistance_palivizumab) > 0:
                    resistance_palivizumab_boolean = True 
                    type_resistance_palivizumab = "I"
                
                if  len(resistance_nirsevimab) > 0:
                    resistance_nirsevimab_boolean = True
                    type_resistance_nirsevimab = "I"
                    
                if  len(resistance_suptavumab) > 0:
                    resistance_suptavumab_boolean = True
                    type_resistance_suptavumab = "I"        
                
                # sensible "S"
                for elem in liste_mutation_à_ajouter_au_csv:
                    if elem["palivizumab"] =="S":
                        resistance_palivizumab.append(elem["mutation"])
                    if elem["nirsevimab"] == "S":
                        resistance_nirsevimab.append(elem["mutation"])
                    if elem["suptavumab"] == "S":
                        resistance_suptavumab.append(elem["mutation"])
                        
                if  len(resistance_palivizumab) > 0 and resistance_palivizumab_boolean == False:
                    resistance_palivizumab_boolean = True 
                    type_resistance_palivizumab = "S"
                else:
                    type_resistance_palivizumab = "Undetermined"
                    
                if  len(resistance_nirsevimab) > 0 and resistance_nirsevimab_boolean == False:
                    resistance_nirsevimab_boolean = True
                    type_resistance_nirsevimab = "S"
                else:
                    type_resistance_nirsevimab = "Undetermined"
                    
                if  len(resistance_suptavumab) > 0 and resistance_suptavumab_boolean == False:
                    resistance_suptavumab_boolean = True
                    type_resistance_suptavumab = "S"
                else:
                   type_resistance_suptavumab = "Undetermined" 
                
            
            return resistance_palivizumab,type_resistance_palivizumab,resistance_nirsevimab,type_resistance_nirsevimab,resistance_suptavumab,type_resistance_suptavumab


def concat(resistance_palivizumab):
            result = ""
            for elem in resistance_palivizumab:
                result += f" {elem} |"
            return result
        

def extraction_txt(path_txt:str, fichier:str, csv_posits, csv_result, csv_final_write,barcode:str,ligne_barcode:list)->tuple[list,list]:
    
    
    nom_fichier = os.path.join(path_txt, fichier) #path à changer
    # ouverture du fichier en "read"       
    with open(nom_fichier, "r") as file:

        lignes = file.readlines()
        # ouverture en écriture : 
        final_file = open(f"{nom_fichier[0:22]}/result_proteines_final_{nom_fichier[63]}/{barcode}.txt","a")
        
        # ouverture en lecture :
        csv_ref = open(f"./../references_phylogenie/Mutations_RSV{nom_fichier[63]}.csv","r")
        #découpage du fichier csv en lignes
        content =csv.reader(csv_ref)
        
        mutation_number = 0
        #! protéine G
        i = 0
        formated_mutations_listing = ""
        listing_mutation_returned_G = []
        listing_csv_mutations = []
        liste_mutation_à_ajouter_au_csv_G = []
        liste_mutation_à_ajouter_au_csv_F = []
        #lecture et extraction du fichier txt 
        for y in range(6):
            
            new_sequence = ""

            listing_mutation = []
            new_lines = ""
            
            #génération du fichier txt v2
            # génération de la ligne "titre" pour la protéine G
            if i == 0:
                final_file.write(f"Alignement de la protéine G de {barcode} avec la protéine G de référence\n\n")
                i+=2
                
            if i != 0 :
                # extraction des variables par caractères par 3 lignes 
                # extraction de des éléments de la ligne "target" i.e. de la séquence de référence 

                target = lignes[i]
                
                target = target.split(" ")
                target = [elem for elem in target if elem]
                
                
                target_number = target[1]
                sequence_target = target[2]
                
                # extraction des variables par caractères par 3 lignes 
                # extraction de des éléments de la ligne "symboles" i.e. des symboles correspondants aux "correspondances" ou non entre la seq ref et la query 
                symboles = lignes[i+1]
                
                symboles = symboles.split(" ")
                symboles = [elem for elem in symboles if elem]
                
                number_symboles = symboles[0]
                sequence_symboles = symboles[1]

                # extraction des variables par caractères par 3 lignes 
                # extraction de des éléments de la ligne "query" i.e. de la séquence provenant du génome étudié
                query = lignes[i+2] 
                query = query.split(" ") 
                query = [elem for elem in query if elem]
                
                query_number = query[1]
                sequence_query = query[2]
                
                # spécification pour les 6 rangés de lignes pour la mise en forme du document en supprimant les sauts de lignes.
                if y <6:
                    sequence_target = sequence_target[:-2]
                    sequence_symboles = sequence_symboles[:-2]
                    sequence_query = sequence_query[:-2]
                

                # fonctions d'ajout d'espaces avant les nombres de "bases" pour la mise en forme et l'alignement
                target_number = add_spaces(target_number)
                number_symboles = add_spaces(number_symboles)
                query_number = add_spaces(query_number)
                # fonction permettant de changer les | pour les # quand il y a un non match et que le caractère de la query est un X i.e. un a.a. non déterminé car trou au niveau des reads.
                new_sequence, listing_mutation = changement_pipe_en_hashtag(sequence_target,sequence_query,new_sequence,sequence_symboles,listing_mutation,int(query_number))
                
                #permet la génération d'un ""tableau"" comprenant les différentes mutations. (petit destructuring qui fait joli)
                formated_mutations_listing = f"{formated_mutations_listing}\n{query_number} {listing_mutation}"
                listing_mutation_returned_G = listing_mutation_returned_G+[*listing_mutation]
                # ajout dans le format.txt
                final_file.write(f"target    {target_number} {sequence_target}\n")
                final_file.write(f"          {number_symboles} {new_sequence}\n")
                
                # pour la dernière rangée de lignes, ajouts de la liste des mutations et de sauts de lignes pour le format.
                if y == 5:
                    new_lines = f"\n{formated_mutations_listing}\n\n"
                final_file.write(f"query     {query_number} {sequence_query}{new_lines}\n")
            # le +4 ne sort pas de nulpart, mais il est là pour la lecture du fichier txt. (ps : on est pas au UNO)
            i+=4

            
        #! : protéine F
        # litéralement le même process pour la protéine F
        i = 27 
        new_lines = ""   
        formated_mutations_listing = ""
        listing_mutation_returned_F = []

        for y in range(10):

            new_sequence = ""
            listing_mutation = []
            
            target = lignes[i]
            target = target.split(" ")
            target = [elem for elem in target if elem]
            
            target_number = target[1]
            
            sequence_target = target[2]
            

            symboles = lignes[i+1]
            symboles = symboles.split(" ")
            symboles = [elem for elem in symboles if elem]
            
            number_symboles = symboles[0]
            sequence_symboles = symboles[1]
            

            query = lignes[i+2] 
            query = query.split(" ") 
            query = [elem for elem in query if elem]

            query_number = query[1]
            sequence_query = query[2]
            

            if y <10:
                sequence_target = sequence_target[:-2]
                sequence_symboles = sequence_symboles[:-2]
                sequence_query = sequence_query[:-2]
            

            target_number = add_spaces(target_number)
            number_symboles = add_spaces(number_symboles)
            query_number = add_spaces(query_number)
        

            new_sequence, listing_mutation = changement_pipe_en_hashtag(sequence_target,sequence_query,new_sequence,sequence_symboles,listing_mutation,int(query_number))
            
            
            formated_mutations_listing = f"{formated_mutations_listing}\n{query_number} {listing_mutation}"
            listing_mutation_returned_F = listing_mutation_returned_F + [*listing_mutation]
            
            final_file.write(f"target    {target_number} {sequence_target}\n")
            final_file.write(f"          {number_symboles} {new_sequence}\n")

            if y == 9:
                new_lines = f"\n{formated_mutations_listing}\n\n"
            final_file.write(f"query     {query_number} {sequence_query}{new_lines}\n")

            i += 4

        for lignes in content:
            protein = lignes[0]
            subunit = lignes[1]
            mutation = lignes[2]
            """
            # RSV_A = lignes[3]
            # RSV_B = lignes[4]
            # Comment = lignes[5]
            # Impact_on_fitness = lignes[6]
            """
            palivizumab = lignes[7]
            FC_IC50_palivizumab = lignes[8]
            nirsevimab = lignes[9]
            FC_IC50_nirsevimab = lignes[10]
            suptavumab = lignes[11]
            FC_IC50_suptavumab = lignes[12]
            

            
            if protein == "G":
                    
                presence_verification_list = []
                unpredicted_verification_list = []
                
                # compte le nombre de + dans la string mutation si il y en a plusieurs, alors il y a plusieurs mutations (no shit Sherlock)
                if mutation.count("+") == "0":
                    
                    # la variable unpredicted est set a Fasle
                    
                    
                    presence = [index for (index, item) in enumerate(listing_mutation_returned_G) if item == mutation]
                    
                    
                    if presence:
                        
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        csv_posits.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        mutation_number+=1
                        liste_mutation_à_ajouter_au_csv_G.append({"mutation":mutation,"palivizumab":palivizumab,"nirsevimab":nirsevimab,"suptavumab":suptavumab})
                    else : 
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},False,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        
                elif mutation.count("+") == "1":
                    
                    liste_mutations = mutation.split("+")
                    
                    for elem in liste_mutations:
                        
                        presence = [index for (index, item) in enumerate(listing_mutation_returned_G) if item == elem]
                        if presence:
                            presence_verification_list.append(True)
                            mutation_number += 1    
                        
                            
                    if len(presence_verification_list) == len(liste_mutations):
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        csv_posits.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        mutation_number += 1 
                        liste_mutation_à_ajouter_au_csv_G.append({"mutation":mutation,"palivizumab":palivizumab,"nirsevimab":nirsevimab,"suptavumab":suptavumab})
                    else:
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},False,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                
            elif protein == "F":
                
                presence_verification_list = []
                
                if mutation.count("+") == 0:
                    
                    presence = [index for (index, item) in enumerate(listing_mutation_returned_F) if item == mutation]
                    

                    
                    if presence:
                        
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        csv_posits.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        mutation_number+=1
                        liste_mutation_à_ajouter_au_csv_F.append({"mutation":mutation,"palivizumab":palivizumab,"nirsevimab":nirsevimab,"suptavumab":suptavumab})
                    else : 
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},False,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                    
                elif mutation.count("+") >= 1:
                    
                    liste_mutations = mutation.split("+")

                    for elem in liste_mutations:
                        
                        presence = [index for (index, item) in enumerate(listing_mutation_returned_F) if item == elem]
                        if presence:
                            presence_verification_list.append(True)
                            mutation_number+=1
                            
                            
                        
                    if len(presence_verification_list) == len(liste_mutations):
                        
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        csv_posits.write(f"{barcode},{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")
                        
                        liste_mutation_à_ajouter_au_csv_F.append({"mutation":mutation,"palivizumab":palivizumab,"nirsevimab":nirsevimab,"suptavumab":suptavumab, "mutation_numbers":len(liste_mutations)})
                    else:
                        csv_result.write(f"{barcode},{protein},{subunit},{mutation},False,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]}\n")    
        
        if mutation_number == 0:
            # (f"{barcode}, Couverture_génome_à_10X(en%), exploitable(à_70%), Type, Couverture_protéine_G(en%_à10X), Couverture_protéine_F(en%_à10X), Mediane_de_couverture_en_reads, Présence_dulication_protéine_G,{protein},{subunit},{mutation},True,{palivizumab},{FC_IC50_palivizumab},{nirsevimab},{FC_IC50_nirsevimab},{suptavumab},{FC_IC50_suptavumab[:-1]} ,abscence de mutations")
            csv_final_write.write(f"{ligne_barcode[0]},{ligne_barcode[1]},{ligne_barcode[2]},{ligne_barcode[3]},{ligne_barcode[4]},{ligne_barcode[5]},{ligne_barcode[6]},{ligne_barcode[7]},##########,abscence_mutations,Undetermined,abscence_mutations,Undetermined,abscence_mutations,Undetermined\n")
        
        elif mutation_number > 0:
            #! c'est pas finit !!!
            csv_final_write.write(f"{ligne_barcode[0]},{ligne_barcode[1]},{ligne_barcode[2]},{ligne_barcode[3]},{ligne_barcode[4]},{ligne_barcode[5]},{ligne_barcode[6]},{ligne_barcode[7]},")
            #TODO! rajouter la selection des mutations en fonction des résistances
            #TODO! générer une liste des mutations et des résitances.
            
           
            resistance_palivizumab_G,type_resistance_palivizumab_G,resistance_nirsevimab_G,type_resistance_nirsevimab_G,resistance_suptavumab_G,type_resistance_suptavumab_G = selection_mutations_pour_csv(liste_mutation_à_ajouter_au_csv_G)   
            resistance_palivizumab_F,type_resistance_palivizumab_F,resistance_nirsevimab_F,type_resistance_nirsevimab_F,resistance_suptavumab_F,type_resistance_suptavumab_F = selection_mutations_pour_csv(liste_mutation_à_ajouter_au_csv_F)
                
            resistance_palivizumab_F = concat(resistance_palivizumab_F)
            resistance_nirsevimab_F = concat(resistance_nirsevimab_F)
            resistance_suptavumab_F = concat(resistance_suptavumab_F)
            csv_final_write.write(f"F,{resistance_palivizumab_F},{type_resistance_palivizumab_F},{resistance_nirsevimab_F},{type_resistance_nirsevimab_F},{resistance_suptavumab_F},{type_resistance_suptavumab_F} \n")
        final_file.close()
        csv_ref.close()

    return listing_mutation_returned_G,listing_mutation_returned_F
    
#TODO: changer le tableau de récap de présence.
#TODO? : vérifier si il n'est pas possible de faire tourner cette fonction.

def path_creation(directories:list[str]):
    
    

    for directory in directories:
        if not os.path.exists(directory):
            os.makedirs(directory)


def generation_results():

    liste_lignes_csv_v1 = []
    csv_result = open(f"./../10-test_positions/Boolean_table_of_the_presence_of_referenced_mutations_in_the_dataset.csv","a")
    csv_result.write("Id,Protein,Subunit,Mutation,Présence_mutation,palivizumab,FC.IC50.palivizumab,nirsevimab,FC.IC50.nirsevimab,suptavumab,FC.IC50.suptavumab\n")

    csv_posits = open(f"./../10-test_positions/results_mutations/table_of_presence_of_mutations_involving_resistance.csv","a")
    csv_posits.write("Id,Protein,Subunit,Mutation,Présence_mutation,palivizumab,FC.IC50.palivizumab,nirsevimab,FC.IC50.nirsevimab,suptavumab,FC.IC50.suptavumab\n")

    csv_final_write = open(f"./../10.1-dossier_xlsx_result_mutations/resume_mutations.csv","a")

    csv_recap_lecture = open(f"./../2.1-dossier_xlsx_result_pileup/resume_pileup.csv","r")
    content = csv.reader(csv_recap_lecture)
    

    for lignes in content:
        liste_lignes_csv_v1.append(lignes)


    concatenation = f""
    
    for categories in liste_lignes_csv_v1[0]:
        concatenation +=  f"{categories},"
    csv_final_write.write(f"{concatenation}Protein, Mutation_resistance_palivizumab, palivizumab, Mutation_resistance_nirsevimab, nirsevimab, Mutation_resistance_suptavumab, suptavumab \n")


    if len(os.listdir(result_proteine_type_A)) != 0:

        for fichier_A in os.listdir(result_proteine_type_A): 

            for i in range(len(clefs_G_A)):
                if extract_number(fichier_A) == extract_number(clefs_G_A[i]):

                    for y in range(len(liste_lignes_csv_v1)):

                        
                        if extract_number(liste_lignes_csv_v1[y][0]) == extract_number(fichier_A):
                            ligne_barcode = liste_lignes_csv_v1[y]


                            if ligne_barcode != None :
                                barcode = clefs_G_A[i]
                                extraction_txt(result_proteine_type_A,fichier_A,csv_posits,csv_result,csv_final_write,barcode,ligne_barcode)
                                #  TODO! à modifier: v
                                

    if len(os.listdir(result_proteine_type_B)) != 0:

        for fichier_B in os.listdir(result_proteine_type_B):

            for i in range(len(clefs_G_B)):
                if extract_number(fichier_B) == extract_number(clefs_G_B[i]):

                    for y in range(len(liste_lignes_csv_v1)):
                        if extract_number(liste_lignes_csv_v1[y][0]) == extract_number(fichier_B):
                            ligne_barcode = liste_lignes_csv_v1[y]
                            if ligne_barcode != None :
                                barcode = clefs_G_B[i]
                                extraction_txt(result_proteine_type_B,fichier_B,csv_posits,csv_result,csv_final_write,barcode,ligne_barcode)

    csv_result.close()
    csv_posits.close()
    csv_recap_lecture.close()
    csv_final_write.close()

#================================================= Code =================================================#


clefs_G_A, clefs_G_B, clefs_F_A, clefs_F_B, sequences_bc_prot_G_A, sequences_bc_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, duplication_G_A, duplication_G_B = append_sequences_prot(path_A, path_B)


ref_A, ref_B, ref_A_without_duplication, ref_B_without_duplication = extraction_sequence_reference(path_ref, path_ref_without_duplication)

proteines_A, proteine_A_dupli, proteines_B, proteine_B_dupli = extraction_proteines_ref(ref_A, ref_B, ref_A_without_duplication , ref_B_without_duplication)
# , duplication_G_A, duplication_G_B

"""
sequences_prot_G_A,sequences_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, proteines_A, proteines_B = transcription(sequences_prot_G_A,
    sequences_prot_G_B,
    sequences_bc_prot_F_A,
    sequences_bc_prot_F_B,
    proteines_A,
    proteines_B)
"""
sequences_bc_prot_G_A, sequences_bc_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, proteines_A, proteines_B, proteine_A_dupli, proteine_B_dupli = traduction(proteines_A, proteines_B, sequences_bc_prot_G_A, sequences_bc_prot_G_B, sequences_bc_prot_F_A, sequences_bc_prot_F_B, proteine_A_dupli, proteine_B_dupli)

path_creation(directories)

generation_alignement_txtv1(proteines_A, proteine_A_dupli, proteines_B, proteine_B_dupli, sequences_bc_prot_G_A, sequences_bc_prot_G_B, sequences_bc_prot_F_A , sequences_bc_prot_F_B, clefs_G_A, clefs_G_B, clefs_F_A, clefs_F_B, duplication_G_A, duplication_G_B) 




generation_results()
#TODO?: same que le TODO bleu dessus.


#TODO! : chaines de Markov 
