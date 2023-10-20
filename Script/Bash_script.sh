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
folder9=./../4-Graphs
folder10=./../2.1-dossier_xlsx_result_pileup
folder11=./../2-Dossier_results_FASTA
folder12=./../3-Fasta_Sequences_Prots
folder13="./../references_phylogenie"

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

illumina_indexing_check(){
    amb="False"
    ann="False"
    bwt="False"
    pac="False"
    sa="False"

    for files in  "$folder13"/*; do
        echo "coucou $files"
        echo "salut  $folder13/ref_combined_insertion.fasta.amb"
        if [ "$files" == "$folder13/ref_combined_insertion.fasta.amb" ];then 
            echo "TRUE"
            amb="True"
        fi

        if [ "$files" == "$folder13/ref_combined_insertion.fasta.ann" ];then 
            echo "TRUE"
            ann='True'
        fi
        if [ "$files" == "$folder13/ref_combined_insertion.fasta.bwt" ];then
            echo "TRUE"
            bwt='True'
        fi
        if [ "$files" == "$folder13/ref_combined_insertion.fasta.pac" ];then
            echo "TRUE"
            pac='True'
        fi
        if [ "$files" == "$folder13/ref_combined_insertion.fasta.sa" ];then 
            echo "TRUE"
            sa='True'
        fi
    done
    echo "amb = $amb"
    echo "ann = $ann"
    echo "bwt = $bwt"
    echo "pac = $pac"
    echo "sa = $sa"
    
    if  [ "$amb" = "False" ] || [ "$ann" = "False" ] || [ "$bwt" = "False" ] || [ "$pac" = "False" ] || [ "$sa" = "False" ]; then  
        echo "start indexing references"
        bwa index ./../references_phylogenie/ref_combined_insertion.fasta 
        echo "INDEXING REFERENCES COMPLETED"
    fi
}
######-----######-----######-----######-----######-----SCRIPT-----######-----######-----######-----######-----###### 

coeur_machine=$(nproc)
read -p "In order to optimize execution speed, this program uses multithreading whenever possible.
This machine uses $coeur_machine cores.
How many core(s) do you want to use ? [1-$coeur_machine] : " cores

while ! [[ "$cores" =~ ^[0-9]+$ ]] || (("$cores > $coeur_machine")); do
    if ! [[ "$cores" =~ ^[0-9]+$ ]]; then
        read -p "This input is not a number.
Please enter a number.
In order to optimize execution speed, this program uses multithreading whenever possible.
This machine uses $coeur_machine cores
How many core(s) do you want to use ? : " cores
    elif (("$cores > $coeur_machine")); then
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
    rm -r $folder2/* 2>/dev/null
else
    mkdir "$folder2"
    echo "The $folder2 directory has been successfully created."
fi

# generation of the ./../1-fastq/sam folder
if [ -d "$folder3" ]; then
    echo "The $folder3 directory already exists."
    rm -r $folder3/* 2>/dev/null
else
    mkdir "$folder3"
    echo "The $folder3 directory has been successfully created."
fi

# generation of the ./../1-fastq/bam folder
if [ -d "$folder4" ]; then
    echo "The $folder4 directory already exists."
    rm -r $folder4/* 2>/dev/null
else
    mkdir "$folder4"
    echo "The $folder4 directory has been successfully created."
fi

# generation of the ./../1-fastq/sorted.bam folder
if [ -d "$folder5" ]; then
    echo "The $folder5 directory already exists."
    rm -r $folder5/* 2>/dev/null
else
    mkdir "$folder5"
    echo "The $folder5 directory has been successfully created."
fi

# generation of the ./../1-fastq/pileup folder
if [ -d "$folder6" ]; then
    echo -e "The $folder6 directory already exists. \n"
    rm -r $folder6/* 2>/dev/null
else
    mkdir "$folder6"
    echo "The $folder6 directory has been successfully created."
fi

# Ilumina 
if [ "$data_type" == "1" ]; then
    echo "You have chosen to process ilumina data"

    # Move files to the correct "1-fastq/fastq.gz" directories
    for file in `ls $FILES_FASTQ`; do

        mv "$file" "$folder2"

    done
    # faire ça
    # ilumina or nanopore selection and testing
    read -p "which type of fastq do you want to treat ?
    paired (1) | paired-end (2) -> : " pair_type


    while [ "$data_type" != "1" ] && [ "$data_type" != "2" ]; do 
        read -p "which type of fastq do you want to treat ?
        paired (1) | paired-end (2) -> : " pair_type
    done 

    # Preparing data for analysis:
    illumina_indexing_check
    

    python3 analyse_bwa_ilumina.py $cores $pair_type
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

    samtools sort -@ $cores "$folder4/$current_name.bam" -o "$folder5/$current_name.sorted.bam"
    echo "SORTING PROCESSED $current_name"

    samtools index "$folder5/$current_name.sorted.bam"
    echo "INDEXATION PROCESSED $current_name"

    samtools mpileup -a -B "$folder5/$current_name.sorted.bam" -f $REF_GEN -o "$folder6/$current_name.pileup"
    echo "PILEUP PROCESSED $current_name"

done


# Pileup analysis
echo "pileups analysis ..."

if [ -d "$folder10" ]; then
    echo "The $folder10 directory already exists."
    rm -r $folder10/* 2>/dev/null
fi

if [ -d "$folder11" ]; then
    echo "The $folder11 directory already exists."
    rm -r "$folder11"
    mkdir "$folder11"
else
    mkdir "$folder11"
    echo "The $folder11 directory has been successfully created."
fi


python3 pileup_analysis.py
echo "ANALYSE PROCESSED"

# Generating protein fastas
if [ -d "$folder12" ]; then
    echo "The $folder12 directory already exists."
    rm -r $folder12/* 2>/dev/null
fi

echo "protein fasta generation ..."
python3 fasta_protein_generation.py
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

    
    if [ -d "$folder9" ]; then
        echo "The $folder9 directory already exists."
        rm -r $folder9/* 2>/dev/null
    else
        mkdir "$folder9"
        echo "The $folder9 directory has been successfully created."
    fi
    
    echo "graph generation ..."
    python3 graphics_generation.py
    echo "GRAPHICS GENERATION PROCCESSED"
    
fi

#######################################################################################################
##########################################  Partie2  ##################################################
#######################################################################################################
PATH_Gblocks=$(find / -name Gblocks 2>/dev/null)
export PATH="$PATH_Gblocks:$PATH"

dossier1=./../5-genbank

dossier4_1=$dossier1/result
dossier4="./../6-Mafft"
dossier5="./../7-Gblocks"
dossier6="./../8-IQTree"
dossier7="./../9-Distances"
dossier8="./../10-test_positions"
dossier9="./../10.1-dossier_xlsx_result_mutations"

if [ -z "$cores" ]; then
    Coeur=2
else
    Coeur=$cores
fi

########## Script ##########

if [ -d "$dossier1" ]; then
    echo "The $dossier1 directory already exists."
    rm -r $dossier1/*
else
    mkdir "$dossier1"
    echo "The $dossier1 directory has been successfully created."
fi

# Phylogeny generation


    # Sequence preparation

echo "preparing sequences ..."

if [ -d "$dossier4" ]; then
    echo "The $dossier4 directory already exists."
    rm -r $dossier4/*
else
    mkdir "$dossier4"
    echo "The $dossier4 directory has been successfully created."
fi

if [ -d "$dossier4_1" ]; then
    echo "The $dossier4_1 directory already exists."
    rm -r $dossier4_1/*
else
    mkdir "$dossier4_1"
    echo "The $dossier4_1 directory has been successfully created."
fi

python3 fasta_homogenization.py   
echo "PREPARATION PROCESSED" 


echo "multiple alignement ..."
#TODO?: mafft
mafft --auto ./../5-genbank/result/file_all_genomes_A.fasta > ./../6-Mafft/mafft_result_A.fasta 
mafft --auto ./../5-genbank/result/file_all_genomes_B.fasta > ./../6-Mafft/mafft_result_B.fasta 

echo "MULTIPLE ALIGNEMENT PROCESSED"

#TODO?: gblocks

echo "Cleaning sequences ..."

# création du directory 7-Gblocks
if [ -d "$dossier5" ]; then
    echo "The $dossier5 directory already exists."
    rm -r $dossier5/*
else
    mkdir "$dossier5"
    echo "The $dossier5 directory has been successfully created."
fi

export PATH="./Gblocks_0.91b:$PATH"
Gblocks ./../6-Mafft/mafft_result_A.fasta -t=d -b1=./../7-Gblocks/mafft_result_A_gblocks.fasta
Gblocks ./../6-Mafft/mafft_result_B.fasta -t=d -b1=./../7-Gblocks/mafft_result_B_gblocks.fasta
echo "CLEANING PROCESSED"

# IQTREE

echo "trees creation ..."
if [ -d "$dossier6" ]; then
    echo "The $dossier6 already exists."
    rm -r $dossier6/*
else
    mkdir "$dossier6"
    echo "The $dossier6 directory has been successfully created."
fi

mv ./../6-Mafft/mafft_result_A.fasta-gb ./../7-Gblocks
mv ./../6-Mafft/mafft_result_B.fasta-gb ./../7-Gblocks
mv ./../6-Mafft/mafft_result_A.fasta-gb.htm ./../7-Gblocks
mv ./../6-Mafft/mafft_result_B.fasta-gb.htm ./../7-Gblocks

iqtree -s ./../7-Gblocks/mafft_result_A.fasta-gb -m GTR+I+G -nt $Coeur -pre $dossier6/mafft_result_A_iqtree
iqtree -s ./../7-Gblocks/mafft_result_B.fasta-gb -m GTR+I+G -nt $Coeur -pre $dossier6/mafft_result_B_iqtree
echo "TREES CREADTED"

#TODO: Détermination distances génétiques 

echo "genetic distance determination ..." 
if [ -d "$dossier7" ]; then
    echo "The $dossier7 already exists."
    rm -r $dossier7/*
else
    mkdir "$dossier7"
    echo "The $dossier7 directory has been successfully created."
fi

python3 genetic_distances.py
echo "GENETIC DISTANCES DETERMINED"


echo "assignment in progress ..."


python3 lineage_determination.py "A"  
python3 lineage_determination.py "B" 
echo "ASSIGNMENT COMPLETED"

#TODO: mutations 
echo "mutation detection ..."

if [ -d "$dossier8" ]; then
    echo "The $dossier8 already exists."
    rm -r $dossier8/*

fi

if [ -d "$dossier9" ]; then
    echo "The $dossier9 already exists."
    rm -r $dossier9/*

fi

python3 mutation_detection.py
echo "DETECTION COMPLETED"

echo "ANALYSE COMPLETED !"







