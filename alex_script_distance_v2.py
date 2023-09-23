from ete3 import Tree
import os
# c'est BIGGz
directory = "./../8-IQTree"
dossier = f"./../9-Distances/result_txt"
# Lire le fichier .treefile
if os.path.exists(dossier):
    print(f"Le dossier '{dossier}' existe.")
else:
    print(f"Le dossier '{dossier}' n'existe pas.")
    os.mkdir(f"{dossier}")
    print(f"le dossier {dossier} a été créé")

for filename in os.listdir(directory):
    if filename.endswith(".treefile"):
        
        tree = Tree(f"./../8-IQTree/{filename}") #mafft_result_B_GB.fasta-gb.treefile
        distances_min=[]
        # Récupérer les noms des feuilles (sorties)
        leaf_names = tree.get_leaf_names()

        # détermination de la lettre 
        parts = filename.split("_")
        
        viral_type = parts[2]

        with open(f"./../9-Distances/result_txt/distance_{viral_type}_result.txt","w") as file:
            # Parcourir toutes les sorties
            for i, seq1 in enumerate(leaf_names):
                print("Distances génétiques pour la sortie:", seq1)
                
                list_temporaire = []
                liste_seq_tempo = []
                
                for j, seq2 in enumerate(leaf_names):
                    if i != j:  # Exclure la sortie courante
                        if seq1.endswith("_A") or seq1.endswith("_B"):
                            if seq2.endswith("_A") or seq2.endswith("_B"):
                                continue  # Ignorer les séquences se terminant par "_A"
                        
                        dist = tree.get_distance(seq1, seq2)
                        liste_seq_tempo.append(seq2)
                        list_temporaire.append(dist)
                
                distance_min_index = list_temporaire.index(min(list_temporaire))
                print(f"Distance génétique entre {seq1} et {liste_seq_tempo[distance_min_index]}: {list_temporaire[distance_min_index]}\n")
                file.write(f"Distance genetique entre {seq1} et {liste_seq_tempo[distance_min_index]}: {list_temporaire[distance_min_index]}\n")