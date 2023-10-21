import re
import sys 
import os 


ext_variable = sys.argv[1]

distance_A  = "distance_A_result.txt"
distance_B = "distance_B_result.txt"
sous_lignee_A = "sous_lignées_A.txt"
sous_lignee_B = "sous_lignées_B.txt"
ref_fasta_A = "ref_A.fasta"
ref_fasta_B = "ref_B.fasta"
dossier = f"./../9-genetic_distances/attribution_lignée"


def import_variable(ext_variable:str)->tuple[str, str, str]:

    file_name=""
    sous_lignee = ""
    ref_fasta=""
    if ext_variable == "A":
        file_name = distance_A
        sous_lignee = sous_lignee_A
        ref_fasta = ref_fasta_A

    elif ext_variable == "B":
        file_name = distance_B
        sous_lignee = sous_lignee_B
        ref_fasta = ref_fasta_B

    return file_name, sous_lignee, ref_fasta



def main():

    file_name,sous_lignee, ref_fasta = import_variable(ext_variable)

    if os.path.exists(dossier):
        print(f"Le dossier '{dossier}' existe.")
    else:
        print(f"Le dossier '{dossier}' n'existe pas.")
        os.mkdir(f"{dossier}")
        print(f"le dossier {dossier} a été créé")

    with open(f"./../9-genetic_distances/result_txt/{file_name}","r") as open_file:
    
        with open (f"./../9-genetic_distances/attribution_lignée/{sous_lignee}","w") as write_file:
            write_file.write(f"Num ID sous_lignée sous_lignée_apparentée\n")
            
            
            # parcour des lignes du fichier distanc, pour chaque ligne (correspondant à un id spécifique on regarde le match dans le fichier de ref: si match, alors )
            i = 0 
            for nb_ligne in open_file:
                i += 1
                ligne = nb_ligne.split()
                
                id_name = id_distance = ligne[3]
                

                # si les id sont ceux des sequences étudiées, on prend leurs
                if id_distance.endswith("_A") or id_distance.endswith("_B"):
                    
                    

                    id_distance = ligne[5]
                    id_distance = id_distance[:-1] # pour enlever les ":" à la fin de l'id

                # on prend le premier id du fichier distance et on cherche dans la str du second fichier l'occurence si id finit par _A ou _B alors on prend l'id en position +2 avec un espace séparateur
                # ensuite on récupère la ligne dans le fichier et on extrait tous les caractères sauf le premier ">" jusqu'à l'underscore .
                # Distance genetique entre KY654518.1 et KU950692.1: 0.0026060425999999996
                # Distance genetique entre 14_A et KY654518.1: 0.0081626275
                
                with open(f"./../references_phylogenie/{ref_fasta}","r") as open_file_ref:
                    
                        for nombre_ligne in open_file_ref:
                            
                            
                            
                            match = re.search(f"{id_distance[:-2]}", f"{nombre_ligne}")
                            if match :
                                if id_name.endswith("_A") or id_name.endswith("_B"):
                                    write_file.write(f"{i} {id_name} {nombre_ligne[1:].split('_')[0]} {id_distance[:-2]}\n")
                                    print(f"{i} {id_distance[:-2]} {nombre_ligne[1:].split('_')[0]} {id_name}")   


                                else:
                                    write_file.write(f"{i} {id_distance[:-2]} {nombre_ligne[1:].split('_')[0]} {'/'}\n")
                                    print(f"{i} {id_distance[:-2]} {nombre_ligne[1:].split('_')[0]} {'/'}")

                                    
                        
main()


















