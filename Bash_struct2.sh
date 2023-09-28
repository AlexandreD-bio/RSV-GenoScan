#! /bin/bash

########## Variables ##########
source script1.sh
 
#Adding Gblocks to the PATH environment
    #détection de minimap
PATH_Gblocks=$(find / -name Gblocks 2>/dev/null)
export PATH="$PATH_Gblocks:$PATH"

dossier1=./../5-genbank
dossier2=$dossier1/New_References_A
dossier3=$dossier1/New_References_B
dossier4_1=$dossier1/result
dossier4="./../6-Mafft"
dossier5="./../7-Gblocks"
dossier6="./../8-IQTree"
dossier7="./../9-Distances"
dossier8="./../10-test_positions"
dossier9="./../10.1-dossier_xlsx_result_mutations"

if [ -z "$cores" ]; then
    Coeur=4
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

# Generation of a "New_references_A" folder
if [ -d "$dossier2" ]; then
    echo "The $dossier2 directory already exists."
else
    mkdir "$dossier2"
    echo "The $dossier2 directory has been successfully created."
fi

# Generation of a "New_references_B" folder
if [ -d "$dossier3" ]; then
    echo "The $dossier3 directory already exists."
else
    mkdir "$dossier3"
    echo "The $dossier3 directory has been successfully created."
fi

# Phylogeny generation
echo "Would you like to integrate new referenced genomes into the phylogeny?"
echo "If so, please drag and drop your fasta files into the folders corresponding to the genome type (A or B)." 
read -p "Type (Y) after deposit or not to continue : " choix1

# Wrong input
while [ "$choix1" != "Y" ] 
    do 
        echo "The input is incorrect." 
        echo "Would you like to integrate new referenced genomes into the phylogeny?"
        echo "If so, please drag and drop your fasta files into the folders corresponding to the genome type (A or B)." 
        read -p "Type (Y) after deposit or not to continue : " choix1
    done  

# Correct input
if [ "$choix1" == "Y" ]; then

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

fi