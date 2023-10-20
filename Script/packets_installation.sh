#! /bin/bash

#Script Alexandre
#Creation data: 11/08/2023
#Last modification: 11/08/2023

read -p "
In order to use this automated RSV genome analysis pipeline,
several (python)/software packages are required.
The purpose of this script is to facilitate the use of the pipeline by automatically installing packets/software.
would you like to proceed with this installation ? (Y/n) " choix_installation

if [ "$choix_installation" != "Y" ]; then 

    echo "You have chosen not to install packets"
    exit 1
else  
    # vérification de si numpy est déjà installé, sinon intallation
    if python3 -c "import biopython" &> /dev/null; then
        echo "biopython is already installed."
    else
        echo "biopython is not installed."
        echo "Installation in progress..."
        pip3 install biopython
        echo "biopython is installed"
    fi

    if python3 -c "import numpy" &> /dev/null; then
        echo "numpy is already installed."
    else
        echo "numpy is not installed."
        echo "Installation in progress..."
        pip3 install numpy
        echo "numpy is installed"
    fi

    if python3 -c "import bokeh" &> /dev/null; then
        echo "bokeh is already installed."
    else
        echo "bokeh is not installed."
        echo "Installation in progress..."
        pip3 install bokeh
        echo "bokeh is installed"
    fi

    if python3 -c "import ete3" &> /dev/null; then
        echo "ete3 is already installed."
    else
        echo "ete3 is not installed."
        echo "Installation in progress..."
        pip3 install ete3
        echo "ete3 is installed"
    fi

    if python3 -c "import selenium" &> /dev/null; then
        echo "selenium is already installed."
    else
        echo "selenium is not installed."
        echo "Installation in progress..."
        pip3 install selenium
        echo "selenium is installed"
    fi

    if minimap2 --version &> /dev/null; then
        echo "minimap2 is already installed."
    else
        echo "minimap2 is not installed."
        echo "Installing..."
        sudo apt update
        sudo apt install minimap2
        echo "minimap2 is installed."
    fi

    if samtools --version &> /dev/null; then
        echo "samtools is already installed."
    else
        echo "samtools is not installed."
        echo "Installation in progress..."
        sudo apt update
        sudo apt install samtools
        echo "samtools is installed."
    fi

    if mafft --version &> /dev/null; then
    echo "mafft is already installed."
    else
        echo "mafft is not installed."
        echo "Installation in progress..."
        sudo apt update
        sudo apt install mafft 
        echo "mafft is installed."
    fi

    if iqtree2 --version &> /dev/null; then
    echo "IQTree is already installed."
    else
        echo "IQTree is not installed."
        echo "Installation in progress..."
        sudo apt update
        sudo apt install iqtree 
        echo "IQTree is installed."
    fi

    if bwa --version &> /dev/null; then
    echo "bwa is already installed."
    else
        echo "bwa is not installed."
        echo "Installation in progress..."
        sudo apt update
        sudo apt install bwa 
        echo "bwa is installed."
    fi
fi
