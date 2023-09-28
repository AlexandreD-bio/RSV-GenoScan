import os,sys,subprocess


# var_ext_cores = sys.argv[1]
# var_fastq_type = str(sys.argv[2])
var_ext_cores=4
var_fastq_type = "A"
dico_fichier_ilumina={}




dossier_fastq = os.listdir(f"./../1-fastq/fastq.gz")
liste_fichier_ilumina=[]

liste_even=[]

for fichier in dossier_fastq:
    if var_fastq_type == "2":
        if fichier.endswith(".fastq.gz"):
            liste_fichier_ilumina.append(fichier.replace('.fastq.gz', ''))
            dico_fichier_ilumina[str(fichier.replace('.fastq.gz', ''))]=[]
            dico_fichier_ilumina[str(fichier.replace('.fastq.gz', ''))].append(f"{fichier.replace('.fastq.gz', '')[:-1]}")
    else:
        if fichier.endswith(".fastq.gz"):
            print(fichier)
            commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta ./../1-fastq/fastq.gz/{fichier} > ./../1-fastq/sam/{fichier.replace('.fastq.gz', '')}.sam"
            try:
                result=subprocess.check_output(commande, shell=True, universal_newlines=True)
            except subprocess.CalledProcessError as e:
                print(f"Erreur lors de l'exécution de la commande Bash : {e}")
# print(liste_fichier_ilumina)
# print(dico_fichier_ilumina)
# for i in range(len(liste_fichier_ilumina)):

for elem in liste_fichier_ilumina:
    fichier1_key = elem
    fichier1_value = dico_fichier_ilumina[elem]

    # print(elem)

    # print(f"fichier value : {fichier1_value}")

    del(dico_fichier_ilumina[f"{elem}"])

    # print(dico_fichier_ilumina)
    for pages in dico_fichier_ilumina:
        # print(f"pages {pages}")
        if fichier1_value == dico_fichier_ilumina[pages]:
            # print(f"c'est Biggz {fichier1_key} et {pages}" )
            liste_even.append(fichier1_key)
            liste_even.append(pages)
            
            commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta ./../1-fastq/fastq.gz/{fichier1_key}.fastq.gz ./../1-fastq/fastq.gz/{pages}.fastq.gz > ./../1-fastq/sam/{fichier1_value[0].replace('[','').replace(']','').replace('_','')}.sam"
            try:
                result=subprocess.check_output(commande, shell=True, universal_newlines=True)
            except subprocess.CalledProcessError as e:
                print(f"Erreur lors de l'exécution de la commande Bash : {e}")





                
            

    














