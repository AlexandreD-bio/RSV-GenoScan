#! /bin/bash

#Script Alexandre
#Creation data: 02/05/2023
#Last modification: 02/05/2023

#####-----variables-----#####

# Chemin de l'arborescence d'organisation
dossier0=./../_Input_Folder_
dossier1=./../1-fastq
dossier2="$dossier1/fastq.gz"
dossier3="$dossier1/sam"
dossier4="$dossier1/bam"
dossier5="$dossier1/sorted.bam"
dossier6="$dossier1/pileup"
dossier7=./../5-genbank
dossier8="$dossier7/references"


#détection de minimap
PATH_minimap=$(find / -name minimap2 2>/dev/null)
#Adding MINIMAP2 to the PATH environment
export PATH="$PATH_minimap:$PATH"

#Location of the Reference genome
REF_GEN=./../references_phylogenie/ref_combined_insertion.fasta

#FILES contains the fastq file for each sample
extension=.fastq.gz
FILES_FASTQ=./../_Input_Folder_/*$extension
FILES_FASTQ_POST_ORGANISATION=$dossier2/*$extension

# Valeur utilisée pour le multithreading



file_name_searched="references_mafft"
file_name_searched_A="sequences_ref_A_pre_mafft.fasta"
file_name_searched_B="sequences_ref_B_pre_mafft.fasta"




######-----######-----######-----######-----######-----fonctions-----######-----######-----######-----######-----###### 

verification_conditions() {
    echo "Drop your (.fastq.gz) files in the '_Input_Folder_'"
    read -p "After deposit, enter 'Y' to continue: " choice
    
    while [ "$choice" != "Y" ] || [ ! -n "$(ls -A $dossier0)" ] || [ ! "$(echo $dossier/*$extension)" ]; do 
            
        if [ "$choice" != "Y" ]; then 
            echo "You entered the wrong character" 
        elif [ ! -n "$(ls -A $dossier0)" ]; then
            echo "The 'Input_Folder' is empty"
        elif [ ! "$(echo $dossier/*$extension)" ]; then 
            echo "There are no (.fastq.gz) files in the 'Input_Folder'"
        fi
        
        echo "Drop your (.fastq.gz) files in the '_Input_Folder_'"
        read -p "After deposit, enter 'Y' to continue: " choice
    done
}    

######-----######-----######-----######-----######-----SCRIPT-----######-----######-----######-----######-----###### 

coeur_machine=$(nproc)
read -p "In order to optimize execution speed, this program uses multithreading whenever possible.
This machine uses $coeur_machine cores.
How many core(s) do you want to use ? [1-$coeur_machine] : " cores 
    
while ! [[ "$cores" =~ ^[0-9]+$ ]] || (("$cores > $coeur_machine")); do
    if ! [[ "$cores" =~ ^[0-9]+$ ]];
        then

        read -p "This input is not a number.
Please enter a number.
In order to optimize execution speed, this program uses multithreading whenever possible.
This machine uses $coeur_machine cores
How many core(s) do you want to use ? : " cores 
        
    elif (("$cores > $coeur_machine"));
        then
        read -p "This input is greater than the number of cores on this machine.
In order to optimize execution speed, this program uses multithreading whenever possible.
This machine uses $coeur_machine cores
How many core(s) do you want to use ? : " cores 
    fi

done 
echo "$cores"


# génération du dossier d'input
if [ -d "$dossier0" ]; then
    echo "Le répertoire $dossier0 existe déjà."
else
    mkdir "$dossier0"
    echo "Le répertoire $dossier0 a été créé avec succès."
fi

# choix et test de type ilumina ou
read -p "which type of data do you want to treat ?
    Ilumina (1) | nanopore (2) -> : " data_type

while [ "$data_type" != "1" ] && [ "$data_type" != "2" ]; do 
    read -p "which type of data do you want to treat ?
    Ilumina (1) | nanopore (2) -> : " data_type
done  

verification_conditions

######-----######-----######-----######-----######-----génération_dossiers-----######-----######-----######-----######-----######

# generation of the ./../1-fastq folder 
if [ -d "$dossier1" ]; then
    echo "Le répertoire $dossier1 existe déjà."
else
    mkdir "$dossier1"
    echo "Le répertoire $dossier1 a été créé avec succès."
fi

# generation of the ./../1-fastq/fastq.gz folder 
if [ -d "$dossier2" ]; then
    echo "Le répertoire $dossier2 existe déjà."
else
    mkdir "$dossier2"
    echo "Le répertoire $dossier2 a été créé avec succès."
fi

# generation of the ./../1-fastq/sam folder
if [ -d "$dossier3" ]; then
    echo "Le répertoire $dossier3 existe déjà."
else
    mkdir "$dossier3"
    echo "Le répertoire $dossier3 a été créé avec succès."
fi

# generation of the ./../1-fastq/bam folder
if [ -d "$dossier4" ]; then
    echo "Le répertoire $dossier4 existe déjà."
else
    mkdir "$dossier4"
    echo "Le répertoire $dossier4 a été créé avec succès."
fi

# generation of the ./../1-fastq/sorted.bam folder
if [ -d "$dossier5" ]; then
    echo "Le répertoire $dossier5 existe déjà."
else
    mkdir "$dossier5"
    echo "Le répertoire $dossier5 a été créé avec succès."
fi

# generation of the ./../1-fastq/pileup folder
if [ -d "$dossier6" ]; then
    echo "Le répertoire $dossier6 existe déjà."
else
    mkdir "$dossier6"
    echo "Le répertoire $dossier6 a été créé avec succès."
fi

# generation of the ./../5-genbank folder 
if [ -d "$dossier7" ]; then
    echo "Le répertoire $dossier7 existe déjà."
else
    mkdir "$dossier7"
    echo "Le répertoire $dossier7 a été créé avec succès."
fi




# Ilumina 
if [ "$data_type" == "1" ]; then
    echo "You have chosen to process ilumina data"
    # TODO: vérifier si le format de fichier est bien fastq.gz

    Déplacer les fichiers dans les bons répertoires "1-fastq/fastq.gz" 
    for file in `ls $FILES_FASTQ`; do

        mv "$file" "$dossier2"

    done
    
    # préparation des données à l'analyse:
    python3 analyse_bwa_ilumina.py $cores
    echo "BWA PROCESSED"
        

# Nanopore
elif [ "$data_type" == "2" ]; then

    echo "You have chosen to process Nanopore data"
    # Déplacer les fichiers dans les bons répertoires "1-fastq/fastq.gz" 
    for file in `ls $FILES_FASTQ`; do

        mv "$file" "$dossier2"

    done

    # préparation des données à l'analyse:
    for files in `ls $FILES_FASTQ_POST_ORGANISATION`; do 
        
        current_name=$(basename $files .fastq.gz)
        
        echo "$current_name"
        minimap2 -ax map-ont -t $cores $REF_GEN $files > "$dossier3/$current_name.sam"
        echo "MINIMAP2 PROCESSED $files"   
        
    done

fi



for files in `ls $dossier3`; do 
        
    current_name=$(basename $files .sam)
    
    echo "$current_name"
      
    samtools view -b -S -@ $cores "$dossier3/$current_name.sam" > "$dossier4/$current_name.bam"
    echo "SAM to BAM PROCESSED $current_name"
    #rm $current_name.sam

    samtools sort -@ $cores "$dossier4/$current_name.bam" -o "$dossier5/$current_name.sorted.bam"
    echo "SORTING PROCESSED $current_name"
    #rm $current_name.bam

    samtools index "$dossier5/$current_name.sorted.bam"
    echo "INDEXATION PROCESSED $current_name"

    samtools mpileup -a -B "$dossier5/$current_name.sorted.bam" -f $REF_GEN -o "$dossier6/$current_name.pileup"
    echo "PILEUP PROCESSED $current_name"

done


#TODO: analyse des pileups
echo "pileups analysis ..."
python3 alex_v5_non_finit_script_pileup.py
echo "ANALYSE PROCESSED"

#TODO: génération des fastas des protéines 
echo "protein fasta generation ..."
python3 alex_script_génération_fasta_F_G.py
echo "PROTEIN FASTA GENERATION PROCESSED"


read -p "Voulez-vous lancer la génération des graphiques ? (Y/n) " choix2

# TODO: génération des graphiques
while [ "$choix2" != "Y" ] && [ "$choix2" != "n" ];
    do 
        echo "L'input est incorrect." 
        read -p "Voulez-vous lancer la génération des graphiques  ? (Y/n) " choix2
    done  

if [ "$choix2" == "n" ]; then
    echo "Vous avez choisi de ne pas lancer la génération des graphiques"

else 

    dossier9=./../4-graphiques
    if [ -d "$dossier9" ]; then
    echo "Le répertoire $dossier9 existe déjà."
    else
        mkdir "$dossier9"
        echo "Le répertoire $dossier9 a été créé avec succès."
    fi
    
    echo "graph generation ..."
    python3 alex_script_graph_alex.py
    echo "GRAPHICS GENERATION PROCCESSED"
    
fi


export $cores

bash alex_script2.sh

