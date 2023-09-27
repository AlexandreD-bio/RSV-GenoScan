#! /bin/bash

#Script Alexandre
#Creation data: 02/05/2023
#Last modification: 02/05/2023

#####-----variables-----#####

# path variables
folder0=./../_Input_Folder_
folder1=./../1-fastq
folder2="$folder1/fastq.gz"
folder3="$folder1/sam"
folder4="$folder1/bam"
folder5="$folder1/sorted.bam"
folder6="$folder1/pileup"
folder7=./../5-genbank
folder8="$folder7/references"
folder9=./../4-graphiques

#détection de minimap
PATH_minimap=$(find / -name minimap2 2>/dev/null)
#Adding MINIMAP2 to the PATH environment
export PATH="$PATH_minimap:$PATH"

#Location of the Reference genome
REF_GEN=./../references_phylogenie/ref_combined_insertion.fasta

#FILES contains the fastq file for each sample
extension=.fastq.gz
FILES_FASTQ=./../_Input_Folder_/*$extension
FILES_FASTQ_POST_ORGANISATION=$folder2/*$extension

# Valeur utilisée pour le multithreading

file_name_searched="references_mafft"
file_name_searched_A="sequences_ref_A_pre_mafft.fasta"
file_name_searched_B="sequences_ref_B_pre_mafft.fasta"

######-----######-----######-----######-----######-----fonctions-----######-----######-----######-----######-----###### 

verification_conditions() {
    echo "Drop your (.fastq.gz) files in the '_Input_Folder_'"
    read -p "After deposit, enter 'Y' to continue: " choice
    
    while [ "$choice" != "Y" ] || [ ! -n "$(ls -A $folder0)" ] || [ ! "$(echo $dossier/*$extension)" ]; do 
            
        if [ "$choice" != "Y" ]; then 
            echo "You entered the wrong character" 
        elif [ ! -n "$(ls -A $folder0)" ]; then
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

# input file generation
if [ -d "$folder0" ]; then
    echo "The $folder0 directory already exists."
else
    mkdir "$folder0"
    echo "The $folder0 directory has been successfully created."
fi

# ilumina or nanopore selection and testing
read -p "which type of data do you want to treat ?
    Ilumina (1) | nanopore (2) -> : " data_type

while [ "$data_type" != "1" ] && [ "$data_type" != "2" ]; do 
    read -p "which type of data do you want to treat ?
    Ilumina (1) | nanopore (2) -> : " data_type
done  

verification_conditions

######-----######-----######-----######-----######-----génération_dossiers-----######-----######-----######-----######-----######

# generation of the ./../1-fastq folder 
if [ -d "$folder1" ]; then
    echo "The $folder1 directory already exists."
    
else
    mkdir "$folder1"
    echo "The $folder1 directory has been successfully created."
fi

# generation of the ./../1-fastq/fastq.gz folder 
if [ -d "$folder2" ]; then
    echo "The $folder2 directory already exists."
    rm -r $folder2/*
else
    mkdir "$folder2"
    echo "The $folder2 directory has been successfully created."
fi

# generation of the ./../1-fastq/sam folder
if [ -d "$folder3" ]; then
    echo "The $folder3 directory already exists."
    rm -r $folder3/*
else
    mkdir "$folder3"
    echo "The $folder3 directory has been successfully created."
fi

# generation of the ./../1-fastq/bam folder
if [ -d "$folder4" ]; then
    echo "The $folder4 directory already exists."
    rm -r $folder4/*
else
    mkdir "$folder4"
    echo "The $folder4 directory has been successfully created."
fi

# generation of the ./../1-fastq/sorted.bam folder
if [ -d "$folder5" ]; then
    echo "The $folder5 directory already exists."
    rm -r $folder5/*
else
    mkdir "$folder5"
    echo "The $folder5 directory has been successfully created."
fi

# generation of the ./../1-fastq/pileup folder
if [ -d "$folder6" ]; then
    echo "The $folder6 directory already exists."
    rm -r $folder6/*
else
    mkdir "$folder6"
    echo "The $folder6 directory has been successfully created."
fi

# generation of the ./../5-genbank folder 
if [ -d "$folder7" ]; then
    echo "The $folder7 directory already exists."
    rm -r $folder7/*
else
    mkdir "$folder7"
    echo "The $folder7 directory has been successfully created."
fi




# Ilumina 
if [ "$data_type" == "1" ]; then
    echo "You have chosen to process ilumina data"

    # Move files to the correct "1-fastq/fastq.gz" directories
    for file in `ls $FILES_FASTQ`; do

        mv "$file" "$folder2"

    done
    
    # Preparing data for analysis:
    python3 analyse_bwa_ilumina.py $cores
    echo "BWA PROCESSED"
        

# Nanopore
elif [ "$data_type" == "2" ]; then

    echo "You have chosen to process Nanopore data"
    # Move files to the correct "1-fastq/fastq.gz" directories 
    for file in `ls $FILES_FASTQ`; do

        mv "$file" "$folder2"

    done

    # preparing data for analysis:
    for files in `ls $FILES_FASTQ_POST_ORGANISATION`; do 
        
        current_name=$(basename $files .fastq.gz)
        
        echo "$current_name"
        minimap2 -ax map-ont -t $cores $REF_GEN $files > "$folder3/$current_name.sam"
        echo "MINIMAP2 PROCESSED $files"   
        
    done

fi



for files in `ls $folder3`; do 
        
    current_name=$(basename $files .sam)
    
    echo "$current_name"
      
    samtools view -b -S -@ $cores "$folder3/$current_name.sam" > "$folder4/$current_name.bam"
    echo "SAM to BAM PROCESSED $current_name"
    #rm $current_name.sam

    samtools sort -@ $cores "$folder4/$current_name.bam" -o "$folder5/$current_name.sorted.bam"
    echo "SORTING PROCESSED $current_name"
    #rm $current_name.bam

    samtools index "$folder5/$current_name.sorted.bam"
    echo "INDEXATION PROCESSED $current_name"

    samtools mpileup -a -B "$folder5/$current_name.sorted.bam" -f $REF_GEN -o "$folder6/$current_name.pileup"
    echo "PILEUP PROCESSED $current_name"

done


# Pileup analysis
echo "pileups analysis ..."
python3 alex_v5_non_finit_script_pileup.py
echo "ANALYSE PROCESSED"

# Generating protein fastas
echo "protein fasta generation ..."
python3 alex_script_génération_fasta_F_G.py
echo "PROTEIN FASTA GENERATION PROCESSED"


read -p "Do you want to start the graph generation ? (Y/n) " choix2

# Graphics generation
while [ "$choix2" != "Y" ] && [ "$choix2" != "n" ];
    do 
        echo "Input is incorrect." 
        read -p "Do you want to start the graph generation ? (Y/n) " choix2
    done  

if [ "$choix2" == "n" ]; then
    echo "You have chosen not to run graphics generation"
else 

    
    if [ -d "$dossier9" ]; then
        echo "The $dossier9 directory already exists."
        rm -r $dossier9/*
    else
        mkdir "$dossier9"
        echo "The $dossier9 directory has been successfully created."
    fi
    
    echo "graph generation ..."
    python3 alex_script_graph_alex.py
    echo "GRAPHICS GENERATION PROCCESSED"
    
fi


export $cores

bash alex_script2.sh

