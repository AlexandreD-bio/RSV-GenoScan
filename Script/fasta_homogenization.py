from Bio import SeqIO
import os 

"""
# variables de la fonction fasta_fusion
 

# path_sortie_fusion_A = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\fusion_result_fasta_A.fasta" #TODO: fusion des génomes étudiés
# path_sortie_fusion_B = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\fusion_result_fasta_B.fasta" #TODO: fusion des génomes étudiés




# # variables de la fonction clean_ref
# path_input_A = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\sequences_ref_A_pre_mafft.fasta" #TODO*: génomes ref 
# path_input_B = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\sequences_ref_B_pre_mafft.fasta" #TODO*: génomes ref 

# path_output_A = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\cleaned_sequences_ref_A_pre_mafft.fasta" #TODO*: génomes ref sans les protéines
# path_output_B = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\cleaned_sequences_ref_B_pre_mafft.fasta" #TODO*: génomes ref sans les protéines



# # variables de la fonction fasta_fusion
#                         #TODO*: génomes ref sans les protéines                                                                     
# files_to_combine_A = [f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\cleaned_sequences_ref_A_pre_mafft.fasta", path_sortie_fusion_A] #TODO: fusion des génomes étudiés
# files_to_combine_B = [f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\cleaned_sequences_ref_B_pre_mafft.fasta", path_sortie_fusion_B] #TODO: fusion des génomes étudiés
#                         #TODO*: génomes ref sans les protéines        

# output_A = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\combined_files_A.fasta" #TODO? ref sans pro + fusion génomes étudiés 
# output_B = f"{disque}:\\ANALYSE_RSV\\6-Mafft\\cleaned\\combined_files_B.fasta" #TODO? ref sans pro + fusion génomes étudiés 

# New_references_A_path = f"./../5-genbank/New_references_A"

"""
#============================================================ Fonctions ============================================================#


# vérifie si un dossier n'est pas vide
def is_directory_not_empty(directory_path)-> bool:
    # Obtenir la liste des fichiers et dossiers dans le répertoire
    files = os.listdir(directory_path)
    
    # Vérifier si la liste est vide ou non
    if len(files) > 0:
        return True
    else:
        return False

# pour chaque fichier dans le dossier input, si le fichier finit par X.fasta avec X:(A ou B), et mettant sous format sequence de Biopython, les séquences de plus de 15Kb et les ajoutant au fichier output.
def fasta_fusion(input: str, output_file, fasta: str):

    fichiers = os.listdir(input)
    print(output_file)
    with  open(output_file,"a") as output:

        for fichier in fichiers:

            if fichier.endswith(fasta):
                    
                for seq_record in SeqIO.parse(f"{input}\\{fichier}","fasta"):
                        
                    if len(seq_record.seq) > 15000:

                        output.write(f">{seq_record.id}\n{seq_record.seq}\n")

def fasta_fusion_add_sequences(input:str, output_file):

    fichiers = os.listdir(input)

    with  open(output_file,"a") as output:

        for fichier in fichiers:

            if fichier.endswith(".fasta"):
                    
                for seq_record in SeqIO.parse(f"{input}\\{fichier}","fasta"):
                        
                    if len(seq_record.seq) > 15000:

                        output.write(f">{seq_record.id}\n{seq_record.seq}\n")

# fonction prennant un fichier fasta comptenant des sequences de référence, crée un fichier fasta comprenantn uniquement les sequences faisant plus de 15kbases (Génome complet)
def clean_ref(path_input: str, path_output):

    for seq_record in SeqIO.parse(path_input,"fasta"):
        
        if len(seq_record.seq) > 15000:
            print(len(seq_record.seq))
            print(seq_record)
            path_output.write(f">{seq_record.id}\n{seq_record.seq}\n")

# fusion en 1 fichier des sequences de références, et des séquences étudiées.
# def fusion_par_l_atome(files_to_combine: list, output: str):

# #     vvv output vvv
#     with open(output, "w") as output_file:

#     #                       vvv input vvv    
#         for file_name in files_to_combine:


#             with open(file_name, "r") as input_file:

#                 output_file.write(input_file.read())


#============================================================ script ============================================================#
""" # todo! : à remettre si le script fonctionne pas 
    # TODO changer toute la logique du script en suprimant les fonctions
    
def main():
    sortie_file_A = open(f"{chemin_output_structure_dépendant_A}/file_all_genomes_A.fasta","w") # TODO? outputA
    sortie_file_B = open(f"{chemin_output_structure_dépendant_B}/file_all_genomes_B.fasta","w") # TODO? outputB

    fasta_fusion(studied_genomes_path, sortie_file_A, fasta_A) # ajout des séquences des génomes étudiés (A)
    fasta_fusion(studied_genomes_path, sortie_file_B, fasta_B) # ajout des séquences des génomes étudiés (B)

    if is_directory_not_empty(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_A):
        fasta_fusion_add_sequences(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_A, sortie_file_A) # ajout de référence rajouté par l'utilisateur (A)
    if is_directory_not_empty(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_B):  
        fasta_fusion_add_sequences(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_A, sortie_file_B) # ajout de référence rajouté par l'utilisateur (B)


    clean_ref(path_genomes_references_logiciel_A, sortie_file_A)
    clean_ref(path_genomes_references_logiciel_B, sortie_file_B)
    sortie_file_A.close()
    sortie_file_B.close()
main()

# TODO  : fasta par défaut dans le logicel (protéines à supprimer)
# TODO! : fasta étudiés 
# TODO* : fasta à rajouter manuellement (protéines à supprimer)
# TODO? : générer l'output
"""

def fasta_homogenization_main():


    #============================== Variables Globales
    chemin = os.getcwd()
    disque = chemin[0]

    # permet d'obtenir le répertoire parent du parent du script
    script_path = os.path.abspath(__file__)
    parent_dir = os.path.dirname(os.path.dirname(script_path))

    #TODO!Corriger les path qui dépendront de la structure finale du projet et qui devront-être automatique.

    studied_genomes_path = f"{parent_dir}/2-FASTA_result_folder"  #TODO: génomes étudiés

    fasta_A = "A.fasta"
    fasta_B = "B.fasta"

    chemin_output_structure_dépendant_A = f"./../5-genbank/result" 
    chemin_output_structure_dépendant_B = f"./../5-genbank/result"

    # chemin_genomes_input_par_user_A = f"{aefzfez}"
    # chemin_genomes_input_par_user_B = f"{aefzfez}"

    directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_A = f"./../5-genbank/New_References_A"
    directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_B = f"./../5-genbank/New_References_B"

    path_genomes_references_logiciel_A = f"{parent_dir}/references_phylogenie/sequences_ref_A_pre_mafft.fasta"
    path_genomes_references_logiciel_B = f"{parent_dir}/references_phylogenie/sequences_ref_B_pre_mafft.fasta"

    # test A: 
    sortie_file_A = open(f"{chemin_output_structure_dépendant_A}/file_all_genomes_A.fasta","w") # TODO? outputA
    # fasta_clean()

    fichiers = os.listdir(f"./../2-FASTA_result_folder")
    for fichier in fichiers:
        
        if fichier.endswith("A.fasta"):
            
            for seq_record in SeqIO.parse(f"./../2-FASTA_result_folder/{fichier}","fasta"):
                    
                if len(seq_record.seq) > 15000:

                    sortie_file_A.write(f">{seq_record.id}\n{seq_record.seq}\n")

    # if is_directory_not_empty(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_A):
    #     fichiers = os.listdir(f"./../5-genbank/New_References_A")
    #     for fichier in fichiers:

    #         if fichier.endswith(".fasta"):
                    
    #             for seq_record in SeqIO.parse(f"./../5-genbank/New_References_A/{fichier}","fasta"):
                        
    #                 if len(seq_record.seq) > 15000:

    #                     sortie_file_A.write(f">{seq_record.id}\n{seq_record.seq}\n")

    for seq_record in SeqIO.parse(f"./../references_phylogenie/sequences_ref_A_pre_mafft.fasta","fasta"):
            
        if len(seq_record.seq) > 15000:
            print(len(seq_record.seq))
            print(seq_record)
            sortie_file_A.write(f">{seq_record.id}\n{seq_record.seq}\n")

    sortie_file_A.close()

    #test B
    sortie_file_B = open(f"{chemin_output_structure_dépendant_B}/file_all_genomes_B.fasta","w") # fasta_clean()

    fichiers = os.listdir(f"./../2-FASTA_result_folder")
    for fichier in fichiers:

        if fichier.endswith("B.fasta"):
                
            for seq_record in SeqIO.parse(f"./../2-FASTA_result_folder/{fichier}","fasta"):
                    
                if len(seq_record.seq) > 15000:

                    sortie_file_B.write(f">{seq_record.id}\n{seq_record.seq}\n")

    # if is_directory_not_empty(directory_du_dossier_contenant_les_fasta_input_par_user_en_plus_B):
    #     fichiers = os.listdir(f"{parent_dir}/5-genbank/New_References_B")
    #     for fichier in fichiers:

    #         if fichier.endswith(".fasta"):
                    
    #             for seq_record in SeqIO.parse(f"{parent_dir}/5-genbank/New_References_B/{fichier}","fasta"):
                        
    #                 if len(seq_record.seq) > 15000:

    #                     sortie_file_B.write(f">{seq_record.id}\n{seq_record.seq}\n")

    for seq_record in SeqIO.parse(f"{parent_dir}/references_phylogenie/sequences_ref_B_pre_mafft.fasta","fasta"):
            
        if len(seq_record.seq) > 15000:
            print(len(seq_record.seq))
            print(seq_record)
            sortie_file_B.write(f">{seq_record.id}\n{seq_record.seq}\n")
            
    sortie_file_B.close()
fasta_homogenization_main()