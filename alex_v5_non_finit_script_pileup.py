import itertools
import os
from numpy import median,array,random
from typing import Literal
import re 

# protéine G:
# TYPE_A : 4659...5555
# TYPE_B : 4688..5566

# protéine F:
# TYPE_A : 5632...7356
# TYPE_B : 5664...7388  
#============================================================ Variables Globales ============================================================#

chemin = os.getcwd()
disque = chemin[0]


nom_fichier_csv  = r"resume_pileup.csv"

# script_path = os.path.abspath(__file__)
# parent_dir = os.path.dirname(script_path)
# parent_dir = os.path.dirname(parent_dir)

script_path = os.path.abspath(__file__)
parent_dir = os.path.dirname(os.path.dirname(script_path))


# path_windows = f"{disque}:\\ANALYSE_RSV\\1-fastq\\pileup"
# nouveau path 
path_windows = f"{parent_dir}/1-fastq/pileup"



# path_windows_test =f"{disque}:\\ANALYSE_RSV\\env_test"

path_linux = "/media/virologie/MyPassport/ANALYSE_RSV"
path_linux_test = "/media/virologie/MyPassport/ANALYSE_RSV/env_test"

extension = ".pileup"

# [5459; 5530]
# 72 pour le type A 

# [5468; 5527]
# 60 pour le type B

G_protein_type_A = [4659,5555]
F_protein_type_A = [5632,7356]
duplication_type_A = 72
location_duplication_A = [5459, 5530]

G_protein_type_B = [4688,5566]
F_protein_type_B = [5664,7388]
duplication_type_B = 60
location_duplication_B = [5468, 5527]
#================================================== fonctions  ==================================================#
def supression_correction_G(resultG, sequence_consensus : str, type: Literal['A', 'B'], location_duplication_A : list, location_duplication_B: list, decalage_pre_genes: int):
    value = ""
    if type == "A":
        location_Duplication_G = location_duplication_A
    
        if len(resultG) == 2: 

            last_modif = 0
            sorted_resultG = sorted(resultG, key=lambda x: x[0])

            for elem in resultG:

                if sorted_resultG[elem][0] in [G_protein_type_A[0], location_Duplication_G[0]]:
                    sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif

                elif sorted_resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                    sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif + (sorted_resultG[elem][0] - location_Duplication_G[0])

                elif sorted_resultG[elem][0] in [location_Duplication_G[1], G_protein_type_A[1]]:
                    sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + 72 + last_modif

                    if (sorted_resultG[elem][1][0]) == "-":
                        value = "+2"
                        number = extract_number(sorted_resultG[elem][1])
                        # on ajoute N  
                        if number != 0:
                            i = 0
                        # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                            replacements = ""
                            while i < number:
                                replacements = [(sorted_resultG[elem][0],"N")]
                                i += 1
                            last_modif = number
                            sequence_consensus = replace_characters(sequence_consensus, replacements)             
                    elif (sorted_resultG[elem][1][0]) == "+":

                        number = extract_number(sorted_resultG[elem][1])
                        
                        sequence_consensus = replace_range(sequence_consensus, sorted_resultG[elem][0][1], sorted_resultG[elem][0][1 + len(str(number))], "N")
                        value = "-2"
            return sequence_consensus, value
        else: 
            # mais il n'y a qu'un seul élément 
            for elem in resultG:

                if resultG[elem][0] in [G_protein_type_A[0], location_Duplication_G[0]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes

                elif resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + (resultG[elem][0] - location_Duplication_G[0])

                elif resultG[elem][0] in [location_Duplication_G[1], G_protein_type_A[1]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + 72

                    if (resultG[elem][1][0]) == "-":
                        number = extract_number(resultG[elem][1])
                        # on ajoute N  
                        if number != 0:
                            i = 0

                        # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                            replacements = ""
                            while i < number:
                                replacements = [(resultG[elem][0],"N")]
                                i += 1
                            value = "+" + str(number)
                            sequence_consensus = replace_characters(sequence_consensus, replacements)
                            
                    elif (resultG[elem][1][0]) == "+":

                        number = extract_number(resultG[elem][1])
                        
                        sequence_consensus = replace_range(sequence_consensus, resultG[elem][0][1], resultG[elem][0][1 + len(str(number))], "N")
                        value = "-" + str(number)

    elif type == "B":
        location_Duplication_G = location_duplication_B
        if len(resultG) == 2: 

            last_modif = 0
            
            sorted_resultG = sorted(resultG, key=lambda x: x[0])
            for elem in resultG:

                    if sorted_resultG[elem][0] in [G_protein_type_B[0], location_Duplication_G[0]]:
                        sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif

                    elif sorted_resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                        sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif + (sorted_resultG[elem][0] - location_Duplication_G[0])

                    elif sorted_resultG[elem][0] in [location_Duplication_G[1], G_protein_type_B[1]]:
                        sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + 60 + last_modif

                        if (sorted_resultG[elem][1][0]) == "-":
                            value = "+2"
                            number = extract_number(sorted_resultG[elem][1])
                            # on ajoute N  
                            if number != 0:
                                i = 0

                            # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                                replacements = ""
                                while i < number:
                                    replacements = [(sorted_resultG[elem][0],"N")]
                                    i += 1
                                last_modif = number
                                sequence_consensus = replace_characters(sequence_consensus, replacements)
                
                                
                                
                        elif (sorted_resultG[elem][1][0]) == "+":

                            number = extract_number(sorted_resultG[elem][1])
                            
                            sequence_consensus = replace_range(sequence_consensus, sorted_resultG[elem][0][1], sorted_resultG[elem][0][1 + len(str(number))], "N")
                            value = "-2"
            return sequence_consensus, value
        else: 
            # mais il n'y a qu'un seul élément 
            for elem in resultG:

                if resultG[elem][0] in [G_protein_type_B[0], location_Duplication_G[0]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes

                elif resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + (resultG[elem][0] - location_Duplication_G[0])

                elif resultG[elem][0] in [location_Duplication_G[1], G_protein_type_B[1]]:
                    resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + 60

                    if (resultG[elem][1][0]) == "-":
                        number = extract_number(resultG[elem][1])
                        # on ajoute N  
                        if number != 0:
                            i = 0

                        # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                            replacements = ""
                            while i < number:
                                replacements = [(resultG[elem][0],"N")]
                                i += 1
                            value = "+" + str(number)
                            sequence_consensus = replace_characters(sequence_consensus, replacements)
                            
                            
                    elif (resultG[elem][1][0]) == "+":

                        number = extract_number(resultG[elem][1])
                        
                        sequence_consensus = replace_range(sequence_consensus, resultG[elem][0][1], resultG[elem][0][1 + len(str(number))], "N")
                        value = "-" + str(number)

                return sequence_consensus, value
    return str(sequence_consensus), value

#TODO! finir de corriger cette version F
def supression_correction_F(resultF, sequence_consensus : str, type: Literal['A', 'B'], location_duplication_A : list, location_duplication_B: list, decalage_pre_genes: int , decalage_proteine_G: int, decalage_inter_proteine_G_F: int):
    """
        if resultF[elem][0] in [G_protein_type_A[0], location_Duplication_G[0]]:
            resultF[elem][0] = resultF[elem][0] + decalage_pre_genes

        elif resultF[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
            resultF[elem][0] = resultF[elem][0] + decalage_pre_genes + (resultF[elem][0] - location_Duplication_G[0])

        elif resultF[elem][0] in [location_Duplication_G[1], G_protein_type_A[1]]:
            resultF[elem][0] = resultF[elem][0] + decalage_pre_genes + 72
        """
    if len(resultF) == 2: 

        last_modif = 0
        
        sorted_resultF = sorted(resultF, key=lambda x: x[0])

        for elem in resultF:
            # TODO! peut-êre rajouter ici une spécification de l'intervalle de F
            sorted_resultF[elem][0] = sorted_resultF[elem][0] + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + last_modif

            if (sorted_resultF[elem][1][0]) == "-":
                number = extract_number(sorted_resultF[elem][1])
                # on ajoute N  
                if number != 0:
                    i = 0

                # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                    replacements = ""
                    while i < number:
                        replacements = [(sorted_resultF[elem][0],"N")]
                        i += 1
                    last_modif = number
                    sequence_consensus = replace_characters(sequence_consensus, replacements)

            elif (sorted_resultF[elem][1][0]) == "+":

                number = extract_number(sorted_resultF[elem][1])
                
                sequence_consensus = replace_range(sequence_consensus, sorted_resultF[elem][0][1], sorted_resultF[elem][0][1 + len(str(number))], "N")
        return sequence_consensus
    else: 
    # mais il n'y a qu'un seul élément 
        for elem in resultF:
            resultF[elem][0] = resultF[elem][0] + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F 

            if (resultF[elem][1][0]) == "-":
                    number = extract_number(resultF[elem][1])
                    # on ajoute N  
                    if number != 0:
                        i = 0

                    # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                        replacements = ""
                        while i < number:
                            replacements = [(resultF[elem][0],"N")]
                            i += 1
                        sequence_consensus = replace_characters(sequence_consensus, replacements)
                        
            elif (resultF[elem][1][0]) == "+":

                number = extract_number(resultF[elem][1])

                sequence_consensus = replace_range(sequence_consensus, resultF[elem][0][1], resultF[elem][0][1 + len(str(number))], "N")

    return str(sequence_consensus)


def correction_insersion_deletion_G(resultG, sequence_consensus: str, type: Literal['A', 'B'], location_duplication_A: list, location_duplication_B: list, decalage_pre_genes:int):

    value = ""

    if type == "A":
        location_Duplication_G = location_duplication_A
    elif type == "B":
        location_Duplication_G = location_duplication_B

    if len(resultG) == 2: 

        last_modif = 0
        
        sorted_resultG = sorted(resultG, key=lambda x: x[0])

        for elem in resultG:

            if sorted_resultG[elem][0] in [G_protein_type_A[0], location_Duplication_G[0]]:
                sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif

            elif sorted_resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + last_modif + (sorted_resultG[elem][0] - location_Duplication_G[0])

            elif sorted_resultG[elem][0] in [location_Duplication_G[1], G_protein_type_A[1]]:
                sorted_resultG[elem][0] = sorted_resultG[elem][0] + decalage_pre_genes + 72 + last_modif

                if (sorted_resultG[elem][1][0]) == "-":
                    value = "+2"
                    number = extract_number(sorted_resultG[elem][1])
                    # on ajoute N  
                    if number != 0:
                        i = 0

                    # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                        replacements = ""
                        while i < number:
                            replacements = [(sorted_resultG[elem][0],"N")]
                            i += 1
                        last_modif = number
                        sequence_consensus = replace_characters(sequence_consensus, replacements)

            elif (sorted_resultG[elem][1][0]) == "+":

                number = extract_number(sorted_resultG[elem][1])
                
                sequence_consensus = replace_range(sequence_consensus, sorted_resultG[elem][0][1], sorted_resultG[elem][0][1 + len(str(number))], "N")
                value = "-2"

        return str(sequence_consensus), value

    else: 
        # mais il n'y a qu'un seul élément 
        for elem in resultG:

            if resultG[elem][0] in [G_protein_type_A[0], location_Duplication_G[0]]:
                resultG[elem][0] = resultG[elem][0] + decalage_pre_genes

            elif resultG[elem][0] in [location_Duplication_G[0],location_Duplication_G[1]]:
                resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + (resultG[elem][0] - location_Duplication_G[0])

            elif resultG[elem][0] in [location_Duplication_G[1], G_protein_type_A[1]]:
                resultG[elem][0] = resultG[elem][0] + decalage_pre_genes + 72

                if (resultG[elem][1][0]) == "-":
                    number = extract_number(resultG[elem][1])
                    # on ajoute N  
                    if number != 0:
                        i = 0

                    # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                        replacements = ""
                        while i < number:
                            replacements = [(resultG[elem][0],"N")]
                            i += 1
                        value = "+" + str(number)
                        sequence_consensus = replace_characters(sequence_consensus, replacements)

                elif (resultG[elem][1][0]) == "+":

                    number = extract_number(resultG[elem][1])
                    
                    sequence_consensus = replace_range(sequence_consensus, resultG[elem][0][1], resultG[elem][0][1 + len(str(number))], "N")
                    value = "-" + str(number)
        return str(sequence_consensus), value

# TODO! corriger cette fonction
def correction_insersion_deletion_F(resultF, sequence_consensus: str, decalage_pre_genes:int, decalage_proteine_G: int, decalage_inter_proteine_G_F: int):


        if len(resultF) == 2: 

            last_modif = 0
            
            sorted_resultF = sorted(resultF, key=lambda x: x[0])

            for elem in resultF:

                sorted_resultF[elem][0] = sorted_resultF[elem][0] + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + last_modif

                if (sorted_resultF[elem][1][0]) == "-":
                    number = extract_number(sorted_resultF[elem][1])
                    # on ajoute N  
                    if number != 0:
                        i = 0

                    # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                        replacements = ""
                        while i < number:
                            replacements = [(sorted_resultF[elem][0],"N")]
                            i += 1
                        last_modif = number
                        sequence_consensus = replace_characters(sequence_consensus, replacements)

                elif (sorted_resultF[elem][1][0]) == "+":

                    number = extract_number(sorted_resultF[elem][1])
                    
                    sequence_consensus = replace_range(sequence_consensus, sorted_resultF[elem][0][1], sorted_resultF[elem][0][1 + len(str(number))], "N")

            return str(sequence_consensus)

        else: 
            # mais il n'y a qu'un seul élément 
            for elem in resultF:

                resultF[elem][0] = resultF[elem][0] + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F 

                if (resultF[elem][1][0]) == "-":
                    number = extract_number(resultF[elem][1])
                    # on ajoute N  
                    if number != 0:
                        i = 0

                    # traiter le fait d'ajouter "N" charactères si sorted_resultG[elem][1] est une délétion
                        replacements = ""
                        while i < number:
                            replacements = [(resultF[elem][0],"N")]
                            i += 1
                        sequence_consensus = replace_characters(sequence_consensus, replacements)

                elif (resultF[elem][1][0]) == "+":

                    number = extract_number(resultF[elem][1])
                    
                    sequence_consensus = replace_range(sequence_consensus, resultF[elem][0][1], resultF[elem][0][1 + len(str(number))], "N")
            return str(sequence_consensus)

# fonction permettant d'extraire un nombre d'une string
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


def replace_range(string, start, end, replacement):
    return string[:start] + replacement + string[end:]


def replace_characters(string, replacements):
    # Convertir la chaîne en une liste mutable
    chars = list(string)
    
    # Effectuer les remplacements
    for position, replacement in replacements:
        chars[position] = replacement
    
    # Reconvertir la liste en une chaîne de caractères
    new_string = ''.join(chars)
    
    return new_string


def generations_combinaisons(data,target_sum,valeur_min,valeur_max,min_max : Literal["min","max"]):
    

    combinations = []

    for i in range(1, len(data) + 1):
        combinations.extend(itertools.combinations(data, i))

    # ajoute les combinaisons d'une seule possibilités.
    combinations.extend([(item,) for item in data])

    # Applique une logique pour filtrer les combinaisons

    

    if min_max == "min":
        filtered_combinations = [
        comb for comb in combinations if sum(int(item[1][0+len(str(extract_number(item[1])))]) for item in comb) == target_sum \
        and min(data[i][2] for i, item in enumerate(comb)) \
        and all(valeur_min <= item[0] <= valeur_max for item in comb)
        ]
    elif min_max == "max":
        filtered_combinations = [
    comb for comb in combinations if sum(int(item[1][0+len(str(extract_number(item[1])))]) for item in comb) == target_sum \
    and max(data[i][2] for i, item in enumerate(comb)) \
    and all(valeur_min <= item[0] <= valeur_max for item in comb)
    ]
    # Affiche les combinaisons filtrées
    for comb in filtered_combinations:
        # return une liste vide si rien n'est bon
        return comb


def proportion_inser_del_minoritaires(
    depth,
    hashmap,
    base_number,
    stockage_hasmap_inser_del_non_majoritaire,
    liste_hashmap_inser_del_minoritaire_filtree:list
):
    liste_stockage = []
    proportions = [(key, value / int(depth)) for key, value in hashmap.items() if key.startswith(("+", "-"))]

    if proportions:

        key_with_highest_proportion = max(proportions, key=lambda x: x[1])[0]

        stockage_hasmap_inser_del_non_majoritaire = next(x[1] for x in proportions if x[0] == key_with_highest_proportion)

        # filtre les hashmaps d'insertions/délétions non majoritaires, avec une proportion de présence supérieur à 30: 
        if round(stockage_hasmap_inser_del_non_majoritaire*100,2) >= 30 :
            if base_number is not None and key_with_highest_proportion is not None and round(stockage_hasmap_inser_del_non_majoritaire*100,2) is not  None:
                liste_stockage.append(base_number)
                liste_stockage.append(key_with_highest_proportion)
                liste_stockage.append(round(stockage_hasmap_inser_del_non_majoritaire*100,2))
                # print(f"key_with_highest_proportion :{base_number} {key_with_highest_proportion} {round(stockage_hasmap_inser_del_non_majoritaire*100,2)}")
                if liste_stockage:
                    liste_hashmap_inser_del_minoritaire_filtree.append(liste_stockage)                     
    return liste_hashmap_inser_del_minoritaire_filtree


def calcul_proportion_hashmap(hashmap):

# calcul proportions  - / + / , / autres afin de déterminer l'élément le plus présent dans le hashmap
    total_insersions = sum(value for key, value in hashmap.items() if key.startswith("+"))
    proportion_plus = (total_insersions / sum(hashmap.values())) * 100

    total_deletions = sum(value for key, value in hashmap.items() if key.startswith("-"))
    proportion_moins = (total_deletions / sum(hashmap.values())) * 100
    
    total_mutations_A = sum(value for key, value in hashmap.items() if key.startswith("A"))
    proportion_mutations_A = (total_mutations_A / sum(hashmap.values())) * 100

    total_mutations_T = sum(value for key, value in hashmap.items() if key.startswith("T"))
    proportion_mutations_T = (total_mutations_T / sum(hashmap.values())) * 100

    total_mutations_G = sum(value for key, value in hashmap.items() if key.startswith("G"))
    proportion_mutations_G = (total_mutations_G / sum(hashmap.values())) * 100

    total_mutations_C = sum(value for key, value in hashmap.items() if key.startswith("C"))
    proportion_mutations_C = (total_mutations_C / sum(hashmap.values())) * 100

    total_base_standard = sum(value for key, value in hashmap.items() if key.startswith("."))
    proportion_base_standard = (total_base_standard / sum(hashmap.values())) * 100

    return proportion_plus , proportion_moins, proportion_base_standard, proportion_mutations_A, proportion_mutations_T, proportion_mutations_G, proportion_mutations_C


def determination_hashmap_modification_proteineG_sequence(
    decalage_proteine_G,
    decalage_pre_genes,
    duplication_type_A,
    G_protein_type_A,
    liste_hashmap_inser_del_minoritaire_filtree,
    stockage_decalage_proteine_G
    ):
    
    longueur_G = (G_protein_type_A[1] + duplication_type_A + decalage_pre_genes + decalage_proteine_G) - (G_protein_type_A[0] + decalage_pre_genes) 
    
    i = 0
    valeur_min = 0
    valeur_max = 0 
    resultG = []
    modifG = False
    valeur_correction_multiple_3 = multiple_de_3_plus_proche(longueur_G)
    # return forcément +1/-2 ou -1/+2

    
    valeur_min = (G_protein_type_A[0] + decalage_pre_genes)
    # print(f"valeur_min{valeur_min}")
    valeur_max = (G_protein_type_A[1] + duplication_type_A + decalage_pre_genes + decalage_proteine_G)
    # print(f"valeur_max{valeur_max}")


    if valeur_correction_multiple_3 == "-1": 
    
        resultG =  generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max,min_max="max")

        if resultG == []:
            valeur_max = "+2"
            resultG = generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max,min_max="max")

        if resultG != []:
            modifG = True 

        # si pas d'inser/deletions à rajouter, alors on enlève l'insersion/del la moins probable.
        else:
            resultG = generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"+1",min_max="min")
            if generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"+1",min_max="min") == []:

                resultG = generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"-2",min_max="min")




    elif valeur_correction_multiple_3 == "+1":
        
        resultG =  generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max,min_max="max")

        if resultG == []:
            valeur_max = "-2"
            resultG = generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max,min_max="max")
        if resultG != []:
            modifG = True 
        else:
            resultG = generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"+1",min_max="min")
            if generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"+1",min_max="min") == []:

                resultG = generations_combinaisons(stockage_decalage_proteine_G,valeur_correction_multiple_3,valeur_min,"-2",min_max="min")    



    """
        #     resultG = [lst for lst in liste_hashmap_inser_del_minoritaire_filtree if (valeur_min <= lst[0] <= valeur_max) and (lst[1][:2] == "+1")]    
        #     if resultG == []:
        #         resultG = [lst for lst in liste_hashmap_inser_del_minoritaire_filtree if (valeur_min <= lst[0] <= valeur_max) and (lst[1][:2] == "-2")] 


        # if valeur_correction_multiple_3 == "+1":
        #     resultG = [lst for lst in liste_hashmap_inser_del_minoritaire_filtree if (valeur_min <= lst[0] <= valeur_max) and (lst[1][:2] == "-1")]    
        #     if resultG == []:
        #         resultG = [lst for lst in liste_hashmap_inser_del_minoritaire_filtree if (valeur_min <= lst[0] <= valeur_max) and (lst[1][:2] == "+2")]   

        # if resultG != []:
        #     if len(resultG) > 1:
        #         resultG = max(resultG, key=lambda x: x[2])
        #     longueur_G += int(resultG[0][2][:2])
        #     decalage_proteine_G += int(resultG[0][2][:2])
        #     modifG = True
        
    
#si sortie de while, et pas de solution, alors suppression de la modification ( à faire )
# ou si liste vide
    """
    return resultG, modifG


def determination_hashmap_modification_proteineF_sequence(
    decalage_pre_genes,
    decalage_proteine_G,
    decalage_inter_proteine_G_F,
    decalage_proteine_F,
    duplication_type_A_B,
    F_protein_type_B,
    liste_hashmap_inser_del_minoritaire_filtree,
    stockage_decalage_proteine_F
    ):

    longueur_F = (F_protein_type_B[0]+ duplication_type_A_B + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F) \
    - (F_protein_type_B[0] + duplication_type_A_B + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F)
    
    i = 0
    valeur_min = 0
    valeur_max = 0 
    resultF = []
    modifF = False
    valeur_correction_multiple_3 = multiple_de_3_plus_proche(longueur_F)
    # return forcément +1/-2 ou -1/+2
    
    valeur_min = (F_protein_type_B[0] + duplication_type_A_B + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F)
    valeur_max = (F_protein_type_B[0]+ duplication_type_A_B + decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F)


    if valeur_correction_multiple_3 == "-1": 

        resultF =  generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_max,valeur_max, min_max="max")

        if resultF == []:
            valeur_max = "+2"
            resultF = generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_max,valeur_max, min_max="max")

        if resultF != []:
            modifF = True

        else:
            resultF = generations_combinaisons(stockage_decalage_proteine_F,valeur_correction_multiple_3,valeur_min,"+1",min_max="min")
            if generations_combinaisons(stockage_decalage_proteine_F,valeur_correction_multiple_3,valeur_min,"+1",min_max="min") == []:

                resultF = generations_combinaisons(stockage_decalage_proteine_F,valeur_correction_multiple_3,valeur_min,"-2",min_max="min")


    elif valeur_correction_multiple_3 == "+1":
        
        resultF =  generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max, min_max="min")

        if resultF == []:
            valeur_max = "-2"
            resultF = generations_combinaisons(liste_hashmap_inser_del_minoritaire_filtree,valeur_correction_multiple_3,valeur_min,valeur_max, min_max="min")


    if resultF != []:
        modifF = True


    return resultF, modifF


def multiple_de_3_plus_proche(nombre)-> Literal['+1', '-1']:
    quotient = nombre // 3  # Division entière pour obtenir le quotient
    multiple_inf = quotient * 3  # Multiple de 3 inférieur ou égal à nombre
    multiple_sup = (quotient + 1) * 3  # Multiple de 3 supérieur à nombre
    if abs(nombre - multiple_inf) < abs(nombre - multiple_sup):
        return "+1"
    else:
        return "-1"


def base_consensus(
    sequence:str,
    base_ref:str,
    skip,
    depth:str,
    base_number:int,
    stockage_hasmap_inser_del_non_majoritaire,
    stockage_inser_del_majoritaire,
    stockage_valeurs_propres_inser_del_majoritaires,
    liste_hashmap_inser_del_minoritaire_filtree 
):

    # homogénéiation de la séquence en mettant les caractères en caps   
    sequence = sequence.upper()

    # détermination et stockage de la longueur de la sequence pour une optimisation du code et de la génération de la liste + calculs induits.
    longueur = len(sequence)
    # variables 
    lstdel = []
    i = 0
    compteur_moins = 0
    compteur_plus = 0  
    compteur_mutation = 0

    # liste de type de caractères trouvable dans les reads pour une base données, avec leurs nombre correspondants.
    hashmap = {}

    # génération de la liste contenant les éléments de la string 
    # liste qui permettra de déterminer quel est l'élément majoritaire 
    while i < longueur: 

        if sequence[i] == "^":

            if not "." in hashmap:

                hashmap["."]= 1

            else:

                hashmap["."] += 1

            lstdel.append(".")
            # car 
            i += 2  
        elif sequence[i] == "-":

            # [+3ATT]
            # sequence[i] : + 
            # sequence[i + 1] : n  (n = 3       [+3ATT][i, n , , ,i + n + 1])
            # sequence[i+1] + 2 : n + 2 
            y = 1
            if (sequence[i + 1]).isdigit():
                number = sequence[i + 1]

                while (sequence[i +1 + y]).isdigit():
                    
                    number += sequence[i+1+y]
                    y += 1

                if (sequence[i+1]).isdigit():

                    lstdel.append(sequence[ i : i + len(str(number)) + int(number)+1 ])

                    if not (sequence[i : i + len(str(number)) + int(number)+1 ]) in hashmap:
                        hashmap[sequence[i : i + len(str(number)) + int(number)+1 ]] = 1
                    else:
                        hashmap[sequence[i : i + len(str(number)) + int(number)+1 ]] += 1

                    # modification du compteur de bases dans le reads 
                    i += len(sequence[ i : i + len(str(number)) + int(number) + 1 ])
                    
                    compteur_moins += 1
                else: 
                    lstdel.append("")   

        elif sequence[i] == "+":

            # ajout de l'élément dans la liste
            y = 1

            if (sequence[i + 1]).isdigit():

                number = sequence[i + 1]
            
                while (sequence[i + 1 + y]).isdigit():
                    
                    number += sequence[i+1+y]

                    y += 1

                if (sequence[i+1]).isdigit():
                    
                    lstdel.append(sequence[i : i + len(str(number)) + int(number)+1 ])

                    if not (sequence[i : i + len(str(number)) + int(number)+1 ]) in hashmap:

                        hashmap[sequence[i : i + len(str(number)) + int(number)+1 ]] = 1

                    else:

                        hashmap[sequence[i : i + len(str(number)) + int(number)+1 ]] += 1

                    # incrémentation de i pour compenser la taille de l'insersion
                    i += len(sequence[i : i + len(str(number)) + int(number)+1])
                    # compteur déterminant le nombre d'élements + dans la liste 
                    compteur_plus += 1 

                else: 
                    
                    lstdel.append("")


            # ajout du reste des élément 
        elif sequence[i] in ["A", "T", "G", "C"]:
            
            # Si i n'est pas dans le dico, alors append 
            if not sequence[i] in hashmap:

                hashmap[sequence[i]]= 1

            else:

                hashmap[sequence[i]]+=1

            lstdel.append(sequence[i])

            # bientôt useless
            compteur_mutation +=1
            i += 1
        elif sequence[i] in [",", ".", "*", "$"]:
            
            if not "." in hashmap:

                hashmap["."] = 1

            else:
                hashmap["."] += 1

            lstdel.append(".")
            i += 1
        else:

            if not sequence[i] in hashmap:

                hashmap[sequence[i]]= 1

            else:

                hashmap[sequence[i]] += 1
            

            lstdel.append(sequence[i])
            # sans compteur car sortie par défaut     

            i += 1

    
    
    # suppression des points/ virgules ayant été comptés en trop par la fonction
    for elem in hashmap:
        if elem.startswith("+"):
            key = elem
            hashmap["."] = hashmap["."] - hashmap[key]

    for elem in hashmap:
        if elem.startswith("-"):
            key = elem
            hashmap["."] = hashmap["."] - hashmap[key]
    
    proportion_plus, proportion_moins, proportion_base_standard, proportion_mutations_A, proportion_mutations_T, proportion_mutations_G, proportion_mutations_C = calcul_proportion_hashmap(hashmap)

    # (-)
    if (proportion_moins >= 50) \
    or (proportion_moins == proportion_mutations_A and (proportion_moins != 0 and proportion_moins > 35)) \
    or (proportion_moins == proportion_mutations_T and (proportion_moins != 0 and proportion_moins > 35)) \
    or (proportion_moins == proportion_mutations_G and (proportion_moins != 0 and proportion_moins > 35)) \
    or (proportion_moins == proportion_mutations_C and (proportion_moins != 0 and proportion_moins > 35)) \
    or (proportion_moins == proportion_base_standard and proportion_moins != 0 and proportion_moins >= 35):
    
        
        stockage_inser_del_majoritaire[base_number] = hashmap
        # list_determination_skip = [item[1] for item in lstdel if item.startswith("-")]


        # skip = max(list_determination_skip, key = list_determination_skip.count)


        proportions = [(key, value / int(depth)) for key, value in hashmap.items() if key.startswith("-")]

        key_with_highest_proportion = max(proportions, key=lambda x: x[1])[0]
        

        stockage_valeurs_propres_inser_del_majoritaires[base_number] = {"base_number": base_number, "clé":key_with_highest_proportion, "value": stockage_inser_del_majoritaire[base_number][key_with_highest_proportion]/int(depth) }

        if len(key_with_highest_proportion) >= 104:
            skip = len(key_with_highest_proportion) - 4
        elif len(key_with_highest_proportion) >= 13:
            skip = len(key_with_highest_proportion) - 3      
        elif len(key_with_highest_proportion) < 12:
            skip = len(key_with_highest_proportion) - 2


        most_common = base_ref
        
        return most_common, str(skip), stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires,liste_hashmap_inser_del_minoritaire_filtree
    
    
    #(+)
    if (proportion_plus >= 50) \
    or (proportion_plus == proportion_mutations_A and (proportion_plus != 0 and proportion_plus > 35)) \
    or (proportion_plus == proportion_mutations_T and (proportion_plus != 0 and proportion_plus > 35)) \
    or (proportion_plus == proportion_mutations_G and (proportion_plus != 0 and proportion_plus > 35)) \
    or (proportion_plus == proportion_mutations_C and (proportion_plus != 0 and proportion_plus > 35)) \
    or (proportion_plus == proportion_base_standard and proportion_plus != 0 and proportion_plus >= 35):

        # most_common = base_ref + max(list_determination_skip, key = list_determination_skip.count)
        
        stockage_inser_del_majoritaire[base_number] = hashmap

        proportions = [(key, value / int(depth)) for key, value in hashmap.items() if key.startswith(("+"))]

        key_with_highest_proportion = max(proportions, key=lambda x: x[1])[0]
        
        stockage_valeurs_propres_inser_del_majoritaires[base_number] = {"base_number":base_number, "clé":key_with_highest_proportion, "value": stockage_inser_del_majoritaire[base_number][key_with_highest_proportion]/int(depth) }
        y = 1
        most_common = ""
        if len(key_with_highest_proportion) >= 104:
            most_common = base_ref + key_with_highest_proportion[4:]
        elif len(key_with_highest_proportion) >= 13:
            most_common = base_ref + key_with_highest_proportion[3:]        
        elif len(key_with_highest_proportion) < 12:
            most_common = base_ref + key_with_highest_proportion[2:]
        
        skip = "0"
        return most_common, str(skip), stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires,liste_hashmap_inser_del_minoritaire_filtree
    
    # (mutations):
    elif proportion_mutations_A >= 50\
    or proportion_mutations_T >= 50 \
    or proportion_mutations_C >= 50\
    or proportion_mutations_G >= 50:
        most_common = ""
        # A seulement >= 50 
        if proportion_mutations_A >= 50:
            most_common = "A"  
        elif proportion_mutations_T >= 50:
            most_common = "T"
        elif proportion_mutations_C >= 50:
            most_common = "C"
        elif proportion_mutations_G >= 50:
            most_common = "G"
            
        elif proportion_mutations_A >= 50 and proportion_mutations_T >= 50 and proportion_mutations_A == proportion_mutations_T:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "A"
            elif aleatoire == 0:
                most_common = "T"             
        elif proportion_mutations_A >= 50 and proportion_mutations_C >= 50 and proportion_mutations_A == proportion_mutations_C:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "A"
            elif aleatoire == 0:
                most_common = "C"
        elif proportion_mutations_A >= 50 and proportion_mutations_G >= 50 and proportion_mutations_A == proportion_mutations_G:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "A"
            elif aleatoire == 0:
                most_common = "G"
        elif proportion_mutations_T >= 50 and proportion_mutations_C >= 50 and proportion_mutations_T == proportion_mutations_C:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "T"
            elif aleatoire == 0:
                most_common = "C"
        elif proportion_mutations_T >= 50 and proportion_mutations_G >= 50 and proportion_mutations_T == proportion_mutations_G:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "T"
            elif aleatoire == 0:
                most_common = "G" 
        elif proportion_mutations_C >= 50 and proportion_mutations_G >= 50 and proportion_mutations_C == proportion_mutations_G:
            aleatoire = random.randint(0,1)
            if aleatoire == 1:
                most_common = "C"
            elif aleatoire == 0:
                most_common = "G"
        
        # (A et B >= 50) et A > B 
        elif proportion_mutations_A >= 50 and proportion_mutations_T >= 50 and proportion_mutations_A > proportion_mutations_T:
            most_common = "A"
        elif proportion_mutations_A >= 50 and proportion_mutations_C >= 50 and proportion_mutations_A > proportion_mutations_C:
            most_common = "A"
        elif proportion_mutations_A >= 50 and proportion_mutations_G >= 50 and proportion_mutations_A > proportion_mutations_G:
            most_common = "A"
        elif proportion_mutations_T >= 50 and proportion_mutations_C >= 50 and proportion_mutations_T > proportion_mutations_C:
            most_common = "T"
        elif proportion_mutations_T >= 50 and proportion_mutations_G >= 50 and proportion_mutations_T > proportion_mutations_G:
            most_common = "T" 
        elif proportion_mutations_T >= 50 and proportion_mutations_A >= 50 and proportion_mutations_T > proportion_mutations_A:
            most_common = "T" 
        elif proportion_mutations_C >= 50 and proportion_mutations_G >= 50 and proportion_mutations_C > proportion_mutations_G:
            most_common = "C"
        elif proportion_mutations_C >= 50 and proportion_mutations_T >= 50 and proportion_mutations_C > proportion_mutations_T:
            most_common = "C"
        elif proportion_mutations_C >= 50 and proportion_mutations_A >= 50 and proportion_mutations_C > proportion_mutations_A:
            most_common = "C"
        elif proportion_mutations_G >= 50 and proportion_mutations_A >= 50 and proportion_mutations_G > proportion_mutations_A:
            most_common = "G"
        elif proportion_mutations_G >= 50 and proportion_mutations_T >= 50 and proportion_mutations_G > proportion_mutations_T:
            most_common = "G"
        elif proportion_mutations_G >= 50 and proportion_mutations_C >= 50 and proportion_mutations_G > proportion_mutations_C:
            most_common = "G"

        skip = "0"
        liste_hashmap_inser_del_minoritaire_filtree=(proportion_inser_del_minoritaires(depth,hashmap,base_number,stockage_hasmap_inser_del_non_majoritaire,liste_hashmap_inser_del_minoritaire_filtree))


        return most_common, str(skip), stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires,liste_hashmap_inser_del_minoritaire_filtree
    
    #(base référence)
    elif proportion_base_standard > 50:
        skip = "0"
        
        
        liste_hashmap_inser_del_minoritaire_filtree=(proportion_inser_del_minoritaires(depth,hashmap,base_number,stockage_hasmap_inser_del_non_majoritaire,liste_hashmap_inser_del_minoritaire_filtree))


        return base_ref, skip, stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires,liste_hashmap_inser_del_minoritaire_filtree
    
    # (condition normalement jamais utilisée si le pileup est bien fait)
    else:
        
    
        skip = "0"
        liste_hashmap_inser_del_minoritaire_filtree=(proportion_inser_del_minoritaires(depth, hashmap, base_number, stockage_hasmap_inser_del_non_majoritaire,liste_hashmap_inser_del_minoritaire_filtree))


        return base_ref, skip, stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires, liste_hashmap_inser_del_minoritaire_filtree


def calcul_percent(nom_fichier):
	with open (nom_fichier,"r") as file :
            i = 0
            
            nb_bases_couvertes_A = 0
            nb_bases_couvertes_B = 0
            ligne_type_A = 0
            ligne_type_B = 0
            

            for nb_ligne in file:   
                
                i +=1

                ligne = nb_ligne.split()

                depth = ligne[3]
                lignee_ref = ligne[0]
                

                if lignee_ref.startswith("typeA"):

                    ligne_type_A += 1
                    

                    if int(depth) >= 10:
                        
                        nb_bases_couvertes_A += 1 

                
                else:
                    # variable de comptage (obtention du total) pour le calcul du % de couverture 
                    ligne_type_B += 1 

                    # comptage du nombre de bases couvertes par 10 reads pour la lignée A 
                    if int(depth) >= 10:
                        
                        nb_bases_couvertes_B += 1
            
            if ligne_type_A != 0 and ligne_type_B != 0:
                
                pourcentage_couverture_A = (nb_bases_couvertes_A*100)/(ligne_type_A)

                pourcentage_couverture_B = (nb_bases_couvertes_B*100)/(ligne_type_B)

                if pourcentage_couverture_A >= 70 and pourcentage_couverture_A > pourcentage_couverture_B:

                    type = "A"
                    exploitable = True
                    return type, pourcentage_couverture_A, ligne_type_A, ligne_type_B, exploitable

                elif pourcentage_couverture_B >= 70 and pourcentage_couverture_B > pourcentage_couverture_A:

                    type = "B"
                    exploitable = True

                    return type, pourcentage_couverture_B, ligne_type_A, ligne_type_B, exploitable
            
                elif pourcentage_couverture_A > pourcentage_couverture_B: 
                    
                    type = "A"
                    exploitable = False
                    
                    return type, pourcentage_couverture_A, ligne_type_A, ligne_type_B, exploitable 
                

                elif pourcentage_couverture_B > pourcentage_couverture_A:

                    type = "B"
                    exploitable = False

                    return type, pourcentage_couverture_B, ligne_type_A, ligne_type_B, exploitable 
                else:
                    type = pourcentage_couverture_B = ligne_type_A = ligne_type_B = exploitable = "/"
                    return type, pourcentage_couverture_B, ligne_type_A, ligne_type_B, exploitable
            else:

                exploitable = False
                type, pourcentage = "/","/"
                return type, pourcentage, ligne_type_A, ligne_type_B, exploitable


def correction_sequence(
    sequence_consensus,
    stockage_valeurs_propres_inser_del_majoritaires,
    type: Literal['A', 'B'],
    G_protein_type_A,
    F_protein_type_A,
    duplication_type_A,
    G_protein_type_B,
    F_protein_type_B,
    duplication_type_B,
    Duplication,
    liste_hashmap_inser_del_minoritaire_filtree
    ):

    valeur_decalage_pre_genes = 0
    decalage_proteine_G = 0
    stockage_decalage_proteine_G = []
    decalage_inter_proteine_G_F = 0
    stockage_decalage_inter_proteine_G_F = []
    decalage_proteine_F = 0
    stockage_decalage_proteine_F = []
    modulation = 0
    modifG = False
    resultG = []
    modifF = False
    resultF = []

    for elem in stockage_valeurs_propres_inser_del_majoritaires : 
        # renvoie la clé 
        

        # (-)
        if stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0] == "-":
            y = 1
            nombre = ""
            while stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0+y].isdigit():
                
                nombre += stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0+y]
                y += 1

            if type == "A":

                if not Duplication:
                    duplication_type_A = 0 

                # zone pré protéine G
                if int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < G_protein_type_A[0] + valeur_decalage_pre_genes:
                    
                    valeur_decalage_pre_genes -= int(nombre)
                    print(f"decalage_pre_genes {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                # zone intra_protéine G
                elif  G_protein_type_A[0] + valeur_decalage_pre_genes \
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= G_protein_type_A[1]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_proteine_G.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_G -= int(nombre)
                    print(f"decalage_proteine_G {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                #zone inter protéine G et F 
                elif G_protein_type_A[1]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G \
                < int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < F_protein_type_A[0]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_inter_proteine_G_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                                stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                                stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_inter_proteine_G_F -= int(nombre) 
                    print(f"decalage_inter_proteine_G_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


                # zone intra_protéine F 
                elif F_protein_type_A[0] + duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F\
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= F_protein_type_A[0]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F:
                    stockage_decalage_proteine_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_F -=int(nombre)
                    print(f"decalage_proteine_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


            if type == "B":

                if not Duplication:
                    duplication_type_B = 0 

                # zone pré protéine G
                if int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < G_protein_type_B[0] + valeur_decalage_pre_genes:
                    
                    valeur_decalage_pre_genes -= int(nombre)
                    print(f"decalage_pre_genes {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                # zone intra_protéine G
                elif  G_protein_type_B[0] + valeur_decalage_pre_genes \
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= G_protein_type_B[1]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_proteine_G.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_G -= int(nombre)
                    print(f"decalage_proteine_G {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                #zone inter protéine G et F 
                elif G_protein_type_B[1]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G \
                < int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < F_protein_type_B[0]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_inter_proteine_G_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_inter_proteine_G_F -= int(nombre) 
                    print(f"decalage_inter_proteine_G_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


                # zone intra_protéine F 
                elif F_protein_type_B[0] + duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F\
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= F_protein_type_B[0]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F:
                    stockage_decalage_proteine_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_F -=int(nombre)
                    print(f"decalage_proteine_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")



            modulation -= int(nombre)
        
        # (+)
        if stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0] == "+":
            y = 1
            nombre = ""
            while stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0+y].isdigit():
                
                nombre += stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"][0+y]
                y += 1

            if type == "A":

                if not Duplication:
                    duplication_type_A = 0 

                # zone pré protéine G
                if int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < G_protein_type_A[0] + valeur_decalage_pre_genes:
                    
                    valeur_decalage_pre_genes += int(nombre)
                    print(f"decalage_pre_genes {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                # zone intra_protéine G
                elif  G_protein_type_A[0] + valeur_decalage_pre_genes \
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= G_protein_type_A[1]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_proteine_G.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_G += int(nombre)
                    print(f"decalage_proteine_G {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                #zone inter protéine G et F 
                elif G_protein_type_A[1]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G \
                < int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < F_protein_type_A[0]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_inter_proteine_G_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_inter_proteine_G_F += int(nombre) 
                    print(f"decalage_inter_proteine_G_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


                # zone intra_protéine F 
                elif F_protein_type_A[0] + duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F\
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= F_protein_type_A[0]+ duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F:
                    
                    decalage_proteine_F +=int(nombre)
                    print(f"decalage_proteine_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


            if type == "B":

                if not Duplication:
                    duplication_type_B = 0 

                # zone pré protéine G
                if int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < G_protein_type_B[0] + valeur_decalage_pre_genes:
                    
                    valeur_decalage_pre_genes += int(nombre)
                    print(f"decalage_pre_genes {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                # zone intra_protéine G
                elif  G_protein_type_B[0] + valeur_decalage_pre_genes \
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= G_protein_type_B[1]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_proteine_G.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_G += int(nombre)
                    print(f"decalage_proteine_G {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

                #zone inter protéine G et F 
                elif G_protein_type_B[1]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G \
                < int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                < F_protein_type_B[0]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G:
                    stockage_decalage_inter_proteine_G_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_inter_proteine_G_F += int(nombre) 
                    print(f"decalage_inter_proteine_G_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")


                # zone intra_protéine F 
                elif F_protein_type_B[0] + duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F\
                <= int(stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number']) \
                <= F_protein_type_B[0]+ duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F:
                    stockage_decalage_proteine_F.append([stockage_valeurs_propres_inser_del_majoritaires[elem]['base_number'],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]["clé"],\
                                                        stockage_valeurs_propres_inser_del_majoritaires[elem]['value']])
                    decalage_proteine_F +=int(nombre)
                    print(f"decalage_proteine_F {stockage_valeurs_propres_inser_del_majoritaires[elem]}")

            modulation += int(nombre)

    print(f"stockage_decalage_proteine_G{stockage_decalage_proteine_G}")
    print(f"stockage_decalage_inter_proteine_G_F{stockage_decalage_inter_proteine_G_F}")
    print(f"{stockage_decalage_proteine_F}")
    print(f"{stockage_decalage_proteine_F}")
    print(f"liste_hashmap_inser_del_minoritaire_filtree{liste_hashmap_inser_del_minoritaire_filtree}")

    if type == "A":
        
        longueur_G = ((G_protein_type_A[1] + duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G) - (G_protein_type_A[0] + valeur_decalage_pre_genes)) + 1
        if longueur_G % 3 !=0:

            resultG,modifG = determination_hashmap_modification_proteineG_sequence(
                decalage_proteine_G,
                valeur_decalage_pre_genes,
                duplication_type_A,
                G_protein_type_A,
                liste_hashmap_inser_del_minoritaire_filtree,
                stockage_decalage_proteine_G)
        
        elif longueur_G % 3 == 0: resultG = [0]

        if modifG == False and resultG != [0]:
            
            sequence_consensus, decalage_post_correction_G = supression_correction_G(resultG, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes)
        
        # si il existe une possibilité de correction :
        elif modifG == True:
            # traiter insersion / délétion 
            sequence_consensus, decalage_post_correction_G = correction_insersion_deletion_G(resultG, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes)

        # si pas besoin de correction
        
        else:   
            decalage_post_correction_G = "+0"
            # ajout du décalage de la correction de la protéine G au décalage globale

        if decalage_post_correction_G[0] == "+":
            decalage = extract_number(decalage_post_correction_G)
            decalage_proteine_G += decalage
            
        else:
            decalage = extract_number(decalage_post_correction_G)
            decalage_proteine_G -= decalage


        longueur_F = ((F_protein_type_A[1] + duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F) \
        - (F_protein_type_A[0] + duplication_type_A + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F))+1
    
        if longueur_F %3 != 0:
            resultF, modifF = determination_hashmap_modification_proteineF_sequence(
            valeur_decalage_pre_genes,
            decalage_proteine_G,
            decalage_inter_proteine_G_F,
            decalage_proteine_F,
            duplication_type_A,
            F_protein_type_A,
            liste_hashmap_inser_del_minoritaire_filtree,
            stockage_decalage_proteine_F)

        
        if longueur_F % 3 == 0:
            resultF = [0]

        # pas d'inser/del à rajouter donc on corrige le cadre de lecture en supprimant l'insersion/délétion d'origine.
        if modifF == False and resultF != [0]:
            
            sequence_consensus = supression_correction_F(resultF, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes, decalage_proteine_G, decalage_inter_proteine_G_F)
        
        # si il existe une possibilité de correction :
        elif modifF == True:
            
            # traiter insersion / délétion 
            sequence_consensus = correction_insersion_deletion_F(resultF, sequence_consensus, valeur_decalage_pre_genes, decalage_proteine_G,decalage_inter_proteine_G_F)


    if type == "B":

        longueur_G = (G_protein_type_B[1] + duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G) - (G_protein_type_B[0] + valeur_decalage_pre_genes)+1
    
        if longueur_G % 3 !=0:

            resultG,modifG = determination_hashmap_modification_proteineG_sequence(
                decalage_proteine_G,
                valeur_decalage_pre_genes,
                duplication_type_B,
                G_protein_type_B,
                liste_hashmap_inser_del_minoritaire_filtree,
                stockage_decalage_proteine_G)
        
        if longueur_G % 3 == 0:
            resultG = [0]



# TODO puis mettre les calculs et corrections de F après.

        # pas d'inser/del à rajouter donc on corrige le cadre de lecture en supprimant l'insersion/délétion d'origine.
        if modifG == False and resultG != [0]:
            
            sequence_consensus, decalage_post_correction_G = supression_correction_G(resultG, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes)
        
        # si il existe une possibilité de correction :
        elif modifG == True:
            
            # traiter insersion / délétion 
            sequence_consensus, decalage_post_correction_G = correction_insersion_deletion_G(resultG, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes)

        # si pas besoin de correction
        else:    
            decalage_post_correction_G = "+0"

# ajout du décalage de la correction de la protéine G au décalage globale
        if decalage_post_correction_G[0] == "+":
            decalage = extract_number(decalage_post_correction_G)
            decalage_proteine_G += decalage
            
        else:
            decalage = extract_number(decalage_post_correction_G)
            decalage_proteine_G -= decalage

        
        longueur_F = ((F_protein_type_B[1] + duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F + decalage_proteine_F) \
        - (F_protein_type_B[0] + duplication_type_B + valeur_decalage_pre_genes + decalage_proteine_G + decalage_inter_proteine_G_F)) + 1

        if longueur_F %3 != 0:
            resultF, modifF = determination_hashmap_modification_proteineF_sequence(
            valeur_decalage_pre_genes,
            decalage_proteine_G,
            decalage_inter_proteine_G_F,
            decalage_proteine_F,
            duplication_type_B,
            F_protein_type_B,
            liste_hashmap_inser_del_minoritaire_filtree,
            stockage_decalage_proteine_F)
        
        if longueur_F % 3 == 0:
            resultF = [0]

        # pas d'inser/del à rajouter donc on corrige le cadre de lecture en supprimant l'insersion/délétion d'origine.
        if modifF == False and resultF != [0]:
            
            sequence_consensus = supression_correction_F(resultF, sequence_consensus, type, location_duplication_A, location_duplication_B, valeur_decalage_pre_genes, decalage_proteine_G, decalage_inter_proteine_G_F)
        
        # si il existe une possibilité de correction :
        elif modifF == True:
            
            # traiter insersion / délétion 
            sequence_consensus = correction_insersion_deletion_F(resultF, sequence_consensus, valeur_decalage_pre_genes, decalage_proteine_G,decalage_inter_proteine_G_F)



    return sequence_consensus


def determination_sequence_consensus(
    type,
    nom_fichier,
    lignes_type_A,
    lignes_type_B,
    content,
    G_protein_type_A,
    F_protein_type_A,
    duplication_type_A,
    G_protein_type_B,
    F_protein_type_B,
    duplication_type_B
):

    sequence_consensus = ""
    i = 0
    base_G = 0
    base_F = 0
    base_G_dupli = 0
    base_F_dupli = 0

    skip = "0"   
    liste_median = []
    stockage_hasmap_inser_del_non_majoritaire = {}
    stockage_inser_del_majoritaire = {}
    stockage_valeurs_propres_inser_del_majoritaires = {}
    liste_hashmap_inser_del_minoritaire_filtree = []

    # si typage déterminé => A
    if type == "A":
        
        for lignes in range(0, lignes_type_A):

            ligne = content[lignes].split() 

            lignee_ref, base_number, base_ref, depth, sequence = ligne[0], ligne[1], ligne[2], ligne[3], ligne[4]
            liste_median.append(int(depth)) 


        # prot G TYPE_A : 4659...5555
            if i >= G_protein_type_A[0] and i < G_protein_type_A[1] and int(depth) > 9: 

                base_G += 1

            if i >= G_protein_type_A[0] and i < (G_protein_type_A[1] + duplication_type_A) and int(depth) > 9: 

                base_G_dupli += 1

        #prot F TYPE_A : 5632...7356
            if i >= F_protein_type_A[0] and i < F_protein_type_A[1] and int(depth) > 9:
                
                base_F += 1

            if i >= (F_protein_type_A[0] + duplication_type_A) and i < (F_protein_type_A[1] + duplication_type_A) and int(depth) > 9:

                base_F_dupli += 1

            if lignee_ref.startswith("typeA"):
                
                if int(depth) >= 10:
                    
                        
                    if skip == "0":

                        result_base, skip, stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires, liste_hashmap_inser_del_minoritaire_filtree = base_consensus(str(sequence), base_ref, skip, depth, int(base_number), stockage_hasmap_inser_del_non_majoritaire,stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires, liste_hashmap_inser_del_minoritaire_filtree)
                        if result_base == "." :
                            print(f"il y a un point base : {base_number}")
                        sequence_consensus += str(result_base)

                    else : 
                        sequence_consensus += ""
                        skip = str(int(skip) - 1)
                        
                else:
                    sequence_consensus += "N"
            i += 1

    # si typage déterminé => B
    if type == "B":

        for lignes in range(lignes_type_A+1, lignes_type_A+1+lignes_type_B):

            # compteur de bases
            i += 1

            ligne = content[lignes-1].split()   
            lignee_ref, base_number, base_ref, depth, sequence = ligne[0], ligne[1], ligne[2], ligne[3], ligne[4]

            liste_median.append(int(depth)) 


            #protG TYPE_B : 4688..5566
            if i >= G_protein_type_B[0] and i < G_protein_type_B[1] and int(depth) > 9:

                base_G += 1 
            if i >= G_protein_type_B[0] and i < (G_protein_type_B[1] + duplication_type_B) and int(depth) > 9:

                base_G_dupli += 1 
            
            #prot F TYPE_B : 5664...7388
            if i >= F_protein_type_B[0] and i < F_protein_type_B[1] and int(depth) > 9: 
                
                base_F += 1
            if i >= ( F_protein_type_B[0] + duplication_type_B ) and i < ( F_protein_type_B[1] + duplication_type_B ) and int(depth) > 9: 

                base_F_dupli += 1

            if lignee_ref.startswith("typeB"):

                if int(depth) >= 10:
                    
                    if str(skip) == "0":

                        result_base , skip, stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires, liste_hashmap_inser_del_minoritaire_filtree = base_consensus(str(sequence), base_ref,skip, depth,int(base_number),stockage_hasmap_inser_del_non_majoritaire, stockage_inser_del_majoritaire, stockage_valeurs_propres_inser_del_majoritaires,liste_hashmap_inser_del_minoritaire_filtree)
                        sequence_consensus += str(result_base)
                    else: 
                        sequence_consensus += ""
                        skip = str(int(skip) - 1)
                else:
                    sequence_consensus += "N"


    """
    if type == "A":
        pourcent_cover_G = (base_G*100) / (5555 - 4659)
        pourcent_cover_F = (base_F*100) / (7356 - 5632)

    elif type == "B":
        
        pourcent_cover_G = (base_G*100) / (5566 - 4688)
        pourcent_cover_F = (base_F*100) / (7388 - 5664)

    else:
        pourcent_cover_G = "/"
        pourcent_cover_F = "/"
    """
    
    liste_median = median(liste_median) #TODO??????? peut-être à corriger 
    
    return sequence_consensus, liste_median, base_G, base_F, base_G_dupli, base_F_dupli, stockage_valeurs_propres_inser_del_majoritaires,  liste_hashmap_inser_del_minoritaire_filtree


def creation_dossiers(disque: str,nom_fichier_csv: str):
    if not os.path.exists(f"./../2.1-dossier_xlsx_result_pileup"):
        os.makedirs(f"./../2.1-dossier_xlsx_result_pileup")
    sortie_csv = open(f"./../2.1-dossier_xlsx_result_pileup/{nom_fichier_csv}","a")
    
    if not os.path.exists(f"./../2-Dossier_results_FASTA"):
        os.makedirs(f"./../2-Dossier_results_FASTA")
    return sortie_csv
#============================================================ script ============================================================#



sortie_csv = creation_dossiers(disque, nom_fichier_csv)



sortie_csv.write("barcode,Couverture_génome_à_10X(en%),exploitable(à_70%),Type,Couverture_protéine_G(en%_à10X),Couverture_protéine_F(en%_à10X),Mediane_de_couverture_en_reads,Présence_dulication_protéine_G\n")

# import les fichiers du path choisit et en fait une liste dont les éléments sont les déterminants de la boucle 
for fichier in os.listdir(path_windows): # path à changer 
	
    # si le fichier comporte l'extension .pileup, alors il est traité 
    if fichier.endswith(extension):

            # import du fichier de l'occurence de la boucle 
            nom_fichier = os.path.join(path_windows, fichier) #path à changer
            print(nom_fichier)
            # ouverture du fichier en "read"        
            file = open (nom_fichier,"r")
            content = file.readlines()

            # fonction permettant le traitement des fichiers pileups en lisant le fichier de A à Z et en sortant un type A/B,
            # qui permet déviter de refaire des opérations pour l'entièreté du fichier pour les fonction d'après.
            # De plus, il sort le pourcentage de couverture du génome ainsi que si le fichier permet d'exploiter le génome consensus qui sera déterminée plus tard,
            # cela en fonction du pourcentage de couverture du génome (seuil à 70%)
            # Pour finir, il sort le numéro des lignes A et lignes B qui seront réutilisée pour des calculs ultérieurement, et éviter des traitement de données inutiles.

            type, pourcent, lignes_type_A, lignes_type_B, exploitable = calcul_percent(nom_fichier)
        
            # si le génome du virus est considéré comme ayant un génome de Type A ou B        
            if type == "A" or type == "B":
                
                # cette fonction permet de déterminer en lisant la partie A ou la partie B du fichier (A ou B étant déterminés précédemment dans le programme) la séquence consensus, 
                # ainsi que le pourcentage de couverture des protéines G et F, mais aussi de la médiane de couverture du génome en reads. 
                sequence_consensus, mediane, base_G, base_F, base_G_dupli, base_F_dupli, stockage_valeurs_propres_inser_del_majoritaires, liste_hashmap_inser_del_minoritaire_filtree = determination_sequence_consensus(
                        type,
                        nom_fichier,
                        lignes_type_A,
                        lignes_type_B,
                        content,
                        G_protein_type_A,
                        F_protein_type_A,
                        duplication_type_A,
                        G_protein_type_B,
                        F_protein_type_B,
                        duplication_type_B
                        )
                
                if type == "A":
                    longueur_genome_reference_typeA_sans_duplication = 15191

                    # duplication contient un Boolean
                    Duplication = len(sequence_consensus) > longueur_genome_reference_typeA_sans_duplication + 50

                    # TODO!!!! changer les valeurs hardcodées et la méthode de calcul.
                    pourcent_cover_G = (base_G*100) / (5555 - 4659)
                    pourcent_cover_F = (base_F*100) / (7356 - 5632)

                    if Duplication:
                        pourcent_cover_G = (base_G_dupli*100) / (( 5555 + 72 ) - 4659 )
                        pourcent_cover_F = (base_F_dupli*100) / (( 7356 + 72 ) - ( 5632 + 72 ))

                        sequence_consensus = correction_sequence(
                            sequence_consensus,
                            stockage_valeurs_propres_inser_del_majoritaires,
                            type,
                            G_protein_type_A,
                            F_protein_type_A,
                            duplication_type_A,
                            G_protein_type_B,
                            F_protein_type_B,
                            duplication_type_B,
                            Duplication,
                            liste_hashmap_inser_del_minoritaire_filtree
                            )
                        
                    else:
                        duplication_type_A = 0
                        duplication_type_B = 0

                        sequence_consensus = correction_sequence(
                            sequence_consensus,
                            stockage_valeurs_propres_inser_del_majoritaires,
                            type,
                            G_protein_type_A,
                            F_protein_type_A,
                            duplication_type_A,
                            G_protein_type_B,
                            F_protein_type_B,
                            duplication_type_B,
                            Duplication,
                            liste_hashmap_inser_del_minoritaire_filtree
                            )
                        
                    # "barcode,Couverture_génome_à_10X(en%),exploitable(à_70%),Type,Couverture_protéine_G(en%_à10X),Couverture_protéine_F(en%_à10X),Mediane_de_couverture_en_reads,Présence_dulication_protéine_G\n"
                    if exploitable:
                        
                        sortie_csv.write(f"{nom_fichier[-9:-7]},{round(float(pourcent),2)},{exploitable},A,{round(float(pourcent_cover_G),2)},{round(float(pourcent_cover_F),2)},{mediane},{Duplication}\n")

                        with open(f"./../2-Dossier_results_FASTA/result_bc{nom_fichier[-9:-7]}_A.fasta","w") as sortie_fastaA:

                            sortie_fastaA.write(f">{nom_fichier[-9:-7]}_A\n")
                            sortie_fastaA.write(f"{sequence_consensus}")
                    else:

                        sortie_csv.write(f"{nom_fichier[-9:-7]},{round(float(pourcent),2)},{exploitable},A,{round(float(pourcent_cover_G),2)},{round(float(pourcent_cover_F),2)},{mediane},{Duplication}\n")
                        
                elif type =="B":
                    
                    longueur_genome_reference_typeB_sans_duplication = 15225

                    Duplication = len(sequence_consensus) > longueur_genome_reference_typeB_sans_duplication + 45


                    #pourcent_cover_G = (base_G*100) / (5566 - 4688)
                    #pourcent_cover_F = (base_F*100) / (7388 - 5664)
                    pourcent_cover_G = (base_G*100) / (5566 - 4688)
                    pourcent_cover_F = (base_F*100) / (7388 - 5664)

                    if Duplication:
                        pourcent_cover_G = (base_G_dupli*100) / ((5566 + 60) - 4688 )
                        pourcent_cover_F = (base_F_dupli*100) / ((7388 + 60) - ( 5664 + 60 ))
                    
                        sequence_consensus = correction_sequence(
                            sequence_consensus,
                            stockage_valeurs_propres_inser_del_majoritaires,
                            type,
                            G_protein_type_A,
                            F_protein_type_A,
                            duplication_type_A,
                            G_protein_type_B,
                            F_protein_type_B,
                            duplication_type_B,
                            Duplication,
                            liste_hashmap_inser_del_minoritaire_filtree
                            )
                    else:
                        duplication_type_A = 0
                        duplication_type_B = 0

                        sequence_consensus = correction_sequence(
                            sequence_consensus,
                            stockage_valeurs_propres_inser_del_majoritaires,
                            type,
                            G_protein_type_A,
                            F_protein_type_A,
                            duplication_type_A,
                            G_protein_type_B,
                            F_protein_type_B,
                            duplication_type_B,
                            Duplication,
                            liste_hashmap_inser_del_minoritaire_filtree
                            )
                    if exploitable:

                        sortie_csv.write(f"{nom_fichier[-9:-7]},{round(float(pourcent),2)},{exploitable},B,{round(float(pourcent_cover_G),2)},{round(float(pourcent_cover_F),2)},{mediane},{Duplication}\n")

                        with open(f"./../2-Dossier_results_FASTA/result_bc{nom_fichier[-9:-7]}_B.fasta","w") as sortie_fastaB:

                            sortie_fastaB.write(f">bc{nom_fichier[-9:-7]}_B\n")
                            sortie_fastaB.write(f"{sequence_consensus}")
                    else: 
                        

                        sortie_csv.write(f"{nom_fichier[-9:-7]},{round(float(pourcent),2)},{exploitable},B,{round(float(pourcent_cover_G),2)},{round(float(pourcent_cover_F),2)},{mediane},{Duplication}\n")
            else:  
                sortie_csv.write(f"{nom_fichier[-9:-7]},/,/,/,/,/,/,/\n")

            file.close()

sortie_csv.close()                 

    
    # TODO* : optimiser le code pour qu'il tourne plus vite (HAHAHAHAHAHAHAHHHAHAHHHAHAHAHHAAHAHHAAHAHHAHAHAHAAHHAHA que je suis drôle.)
    
    # TODO! : voir si il y a moyen de paramétrer le csv avant de l'ouvrir.
    # TODO!! : voir si il n'y a pas moyen d'utiliser les parties de code sur les protéines pour générer des fastas des prots pour les prots G et F


    #TODO! corriger les valeurs des coordonnées pour le RSV 
    #TODO! passer totalement à la technique d'hashmap







#     si pas de couverture alors N +/- 30

# insertiosn de plus de longuer 