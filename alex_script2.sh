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

if [ -z "$cores" ]; then
    Coeur=4
else
    Coeur=$cores
fi


########## Script ##########

if [ -d "$dossier1" ]; then
    echo "Le répertoire $dossier1 existe déjà."
else
    mkdir "$dossier1"
    echo "Le répertoire $dossier1 a été créé avec succès."
fi

# génération d'un dossier New_references_A
if [ -d "$dossier2" ]; then
    echo "Le répertoire $dossier2 existe déjà."
else
    mkdir "$dossier2"
    echo "Le répertoire $dossier2 a été créé avec succès."
fi

# génération d'un dossier New_references_B
if [ -d "$dossier3" ]; then
    echo "Le répertoire $dossier3 existe déjà."
else
    mkdir "$dossier3"
    echo "Le répertoire $dossier3 a été créé avec succès."
fi


    #TODO: génération phylogénie
echo "Voulez-vous intégrer de nouveaux génomes référencés à la phylogénie ?"
echo "si oui veuillez glissez vos fichiers fasta dans les dossiers respectivement au type des génomes (A ou B)" 
read -p "Tapez (Y) après dépot ou non pour continuer : " choix1

# mauvais input
while [ "$choix1" != "Y" ] 
    do 
        echo "L'input est incorrect." 
        echo "Voulez-vous intégrer de nouveaux génomes référencés à la phylogénie ?"
        echo "si oui veuillez glissez vos fichiers fasta dans les dossiers respectivement au type des génomes (A ou B)" 
        read -p "Tapez (Y) après dépot ou non pour continuer : " choix1
    done  
# bon input
if [ "$choix1" == "Y" ]; then

    # TODO?: préparation des séquences

    echo "preparing sequences ..."
    

    

    echo "multiple alignement ..."

    if [ -d "$dossier4" ]; then
        echo "Le répertoire $dossier4 existe déjà."
    else
        mkdir "$dossier4"
        echo "Le répertoire $dossier4 a été créé avec succès."
    fi

    if [ -d "$dossier4_1" ]; then
        echo "Le répertoire $dossier4_1 existe déjà."
    else
        mkdir "$dossier4_1"
        echo "Le répertoire $dossier4_1 a été créé avec succès."
    fi



    
    python3 fasta_clean.py   
    echo "PREPARATION PROCESSED" 
    #TODO?: mafft
    mafft --auto ./../5-genbank/result/file_all_genomes_A.fasta > ./../6-Mafft/mafft_result_A.fasta 
    mafft --auto ./../5-genbank/result/file_all_genomes_B.fasta > ./../6-Mafft/mafft_result_B.fasta 

    echo "MULTIPLE ALIGNEMENT PROCESSED"

    #TODO?: gblocks

    echo "Cleaning sequences ..."

    # création du directory 7-Gblocks
    if [ -d "$dossier5" ]; then
        echo "Le répertoire $dossier5 existe déjà."
    else
        mkdir "$dossier5"
        echo "Le répertoire $dossier5 a été créé avec succès."
    fi

    export PATH="./Gblocks_0.91b:$PATH"
    Gblocks ./../6-Mafft/mafft_result_A.fasta -t=d -b1=./../7-Gblocks/mafft_result_A_gblocks.fasta
    Gblocks ./../6-Mafft/mafft_result_B.fasta -t=d -b1=./../7-Gblocks/mafft_result_B_gblocks.fasta
    echo "CLEANING PROCESSED"

    #TODO?: IQTREE

    echo "trees creation ..."
    if [ -d "$dossier6" ]; then
        echo "Le répertoire $dossier6 existe déjà."
    else
        mkdir "$dossier6"
        echo "Le répertoire $dossier6 a été créé avec succès."
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
        echo "Le répertoire $dossier7 existe déjà."
    else
        mkdir "$dossier7"
        echo "Le répertoire $dossier7 a été créé avec succès."
    fi

    
    python3 alex_script_distance_v2.py
    echo "GENETIC DISTANCES DETERMINED"
    


    echo "assignment in progress ..."
    if [ -d "$dossier9" ]; then
        echo "Le répertoire $dossier9 existe déjà."
    else
        mkdir "$dossier9"
        echo "Le répertoire $dossier9 a été créé avec succès."
    fi

    python3 script_lignee_post_distance.py "A" #TODO? normalement bon ? 
    python3 script_lignee_post_distance.py "B" #TODO? normalement bon ?
    echo "ASSIGNMENT COMPLETED"

    #TODO: mutations 
    echo "mutation detection ..."
    python3 detection_mutations.py
    echo "DETECTION COMPLETED"

    echo "ANALYSE COMPLETED !"

fi