import os,sys,subprocess


var_ext_cores = sys.argv[1]
var_fastq_type = str(sys.argv[2])

dico_fichier_ilumina={}




dossier_fastq = os.listdir(f"./../1-fastq/fastq.gz")
# dossier_fastq = os.listdir(f"./../1-fastq/flux")

# print(dossier_fastq)
number_ = 0
files_count = 0
convention=False

# détection de la convention de nomenclature
for file in dossier_fastq:
    
    files_count+=1
    if file.count("_") == 4:
        number_ += 1
    

# les fichiers illumina sont nommés selon la convention
if number_ == files_count:

    files={}
    id_list = []
    files_names = []

    for file in dossier_fastq:

        dico_convention_name = {}
        dico_convention_name["id"]=f"{file.split('_')[0]}"
        dico_convention_name["sample_index"]=f"{file.split('_')[1]}"
        dico_convention_name["lane"]=f"{file.split('_')[2]}"
        dico_convention_name["Read_number"]=f"{file.split('_')[3]}"
        dico_convention_name["sequencing_iteration"]=f"{file.split('_')[4]}"
        id_list.append(dico_convention_name)
        files_names.append(file.split('_')[0])


    if var_fastq_type == "2":
        print("2.2")
        
        alredy_used_item = []

        # print(set(files_names))
        unique_files_names=set(files_names)

        double=0
        for elem in unique_files_names:

            name_stockage = []
            
            # print(f"elem des ids : {elem}")
            for names in id_list :
                # print(f"names : {names['id']}")
                if names['id'] in alredy_used_item and double ==2 :
                    double=0
                    alredy_used_item.append(names["id"])
                    continue
                if names["id"] == elem:
                    double +=1
                    
                    
                    name_stockage.append(names)
                    if len(name_stockage) ==2:
                        
                        file1=f"./../1-fastq/fastq.gz/{name_stockage[0]['id']}_{name_stockage[0]['sample_index']}_{name_stockage[0]['lane']}_{name_stockage[0]['Read_number']}_{name_stockage[0]['sequencing_iteration']}"
                        file2=f"./../1-fastq/fastq.gz/{name_stockage[1]['id']}_{name_stockage[1]['sample_index']}_{name_stockage[1]['lane']}_{name_stockage[1]['Read_number']}_{name_stockage[1]['sequencing_iteration']}"

                        commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta {file1} {file2} > ./../1-fastq/sam/{name_stockage[0]['id']}.sam"
                        


                        try:
                            result=subprocess.check_output(commande, shell=True, universal_newlines=True)
                        except subprocess.CalledProcessError as e:
                            print(f"Error executing Bash command : {e}")
                        name_stockage =[]
    
    # Pour le non paired-end
    else:
        print("2.1")
        for file in dossier_fastq:
            commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta ./../1-fastq/fastq.gz/{file} > ./../1-fastq/sam/{file.split('_')[0]}.sam"
            try:
                result=subprocess.check_output(commande, shell=True, universal_newlines=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing Bash command : {e}")
# les fichiers illumina ne sont pas nommés selon la convention 
else:
    id_list=[]
    files_names=[]
    for files in dossier_fastq:
        
        dico_non_convention_name={}
        dico_non_convention_name['id']=files.split('_')[0]
        dico_non_convention_name['Read']=files.split('_')[1]
        id_list.append(dico_non_convention_name)
        
        
        files_names.append(files.split('_')[0])

    if var_fastq_type == "2":
        
        alredy_used_item = []
        
        unique_files_names=set(files_names)
        
        double=0
        for elem in unique_files_names:

            name_stockage = []
            
            # print(f"elem des ids : {elem}")
            for names in id_list :
                # print(f"names : {names['id']}")
                if names['id'] in alredy_used_item and double ==2 :
                    double=0
                    alredy_used_item.append(names["id"])
                    continue
                if names["id"] == elem:
                    double +=1
                    
                    
                    name_stockage.append(names)
                    if len(name_stockage) ==2:
                        
                    
                        file1=f"./../1-fastq/fastq.gz/{name_stockage[0]['id']}_{name_stockage[0]['Read']}"
                        file2=f"./../1-fastq/fastq.gz/{name_stockage[1]['id']}_{name_stockage[1]['Read']}"

                        commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta {file1} {file2} > ./../1-fastq/sam/{name_stockage[0]['id']}.sam"
                        print(commande)


                        try:
                            result=subprocess.check_output(commande, shell=True, universal_newlines=True)
                        except subprocess.CalledProcessError as e:
                            print(f"Error executing Bash command : {e}")
                        name_stockage =[]
    else:
        for file in dossier_fastq:
            
            commande = f"bwa mem -t{var_ext_cores} ./../references_phylogenie/ref_combined_insertion.fasta ./../1-fastq/fastq.gz/{file} > ./../1-fastq/sam/{file.split('_')[0]}.sam"
            try:
                result=subprocess.check_output(commande, shell=True, universal_newlines=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing Bash command : {e}")
 