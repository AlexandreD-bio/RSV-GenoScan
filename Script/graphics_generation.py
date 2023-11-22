from bokeh.plotting import figure, show
from bokeh.models import Legend, LegendItem, BoxAnnotation, Label
from bokeh.io import export_png
import os
import re

#======================================== Fonctions ========================================#

def graphiqueB(ligne_B, num_ligne_virus, total_reads, nom_fichier, file_id,result_folder):

# création de la figure
    p = figure(
        title=f"Number of reads for each base of the virus genome (type B) corresponding to {file_id}",
        x_axis_label='Genome Position',
        x_range=(0, ligne_B),
        y_axis_label='number of Reads',
        y_axis_type="log",
        width=2100,
        height=800,
        toolbar_location=None
    )

    # trace de la courbe se faisant en fonction de 2 listes (ici num_ligne_virus, et total reads)
    p.line(
        num_ligne_virus,
        total_reads,
        line_width=1,
    )

    # trace du palier de 10 reads (avec les coordonnées de sont point de début, et de fin)
    p.segment(
        x0=0,
        x1=ligne_B,
        y0=10,
        y1=10,
        line_color="red",
        line_width=2,
    )

    # ajout des boîtes
    Proteine_G = BoxAnnotation(left=4688, right=5566, fill_color='green', fill_alpha=0.2, level="underlay")
    Proteine_F = BoxAnnotation(left=5664, right=7388, fill_color='orange', fill_alpha=0.2, level="underlay")
    p.add_layout(Proteine_G)
    p.add_layout(Proteine_F)

    # paramètres de taille
    p.xaxis.axis_label_text_font_size = '25px'
    p.yaxis.axis_label_text_font_size = '25px'
    p.title.text_font_size = "20px" # type: ignore
    p.xaxis.major_label_text_font_size = '20px'
    p.yaxis.major_label_text_font_size = '20px'

    # ajout de la légende
    legend = Legend(
         
        title="Caption:",
        items=[
            ("Number of Reads", [p.renderers[0]]),# type: ignore
            ("steps of 10 Reads", [p.renderers[1]]),# type: ignore
        ],
        label_text_font_size="16pt",
        title_text_font_size="16pt"  
    )
    
    # génération d'un encadré pour la protéine G
    legende_Prot_G = Label(
         
        x=640,
        y=20,
        x_units='screen',
        y_units='screen',
        text='G protein',
        background_fill_color='yellowgreen',
        background_fill_alpha=0.5   
    )
    
    legende_Prot_F = Label(
        
        x=825,
        y=20,
        x_units='screen',
        y_units='screen',
        text='F Protein',
        background_fill_color='orange',
        background_fill_alpha=0.5,    
    )
    
    p.add_layout(legende_Prot_G)
    p.add_layout(legende_Prot_F)
    p.add_layout(legend, 'below')


    show(p) # type: ignore 

    export_png(p, filename=f"{result_folder}/4-Graphs/Graph_{file_id}_B.png") # type: ignore

    return p 


def graphiqueA(ligne_A, num_ligne_virus, total_reads, nom_fichier, file_id,result_folder:str):
    p = figure(
        title=f"Number of reads for each base of the virus genome (type A) corresponding to {file_id}",
        x_axis_label='Genome Position',
        x_range=(0, ligne_A),
        y_axis_label='Number of Reads',
        y_axis_type="log",
        width=2100,
        height=800,
        toolbar_location=None
    )

    # trace de la ligne
    p.line(
        num_ligne_virus,
        total_reads,
        line_width=1,
    )

    # trace du palier de 10 reads
    p.segment(
        x0=0,
        x1=ligne_A,
        y0=10,
        y1=10,
        line_color="red",
        line_width=2,
    )

    # ajout des boîtes
    Proteine_G = BoxAnnotation(left=4659, right=5555, fill_color='green', fill_alpha=0.2, level="underlay")
    Proteine_F = BoxAnnotation(left=5632, right=7356, fill_color='orange', fill_alpha=0.2, level="underlay")
    p.add_layout(Proteine_G)
    p.add_layout(Proteine_F)

    # paramètres de taille
    p.xaxis.axis_label_text_font_size = '25px'
    p.yaxis.axis_label_text_font_size = '25px'
    p.title.text_font_size = "20px" # type: ignore
    p.xaxis.major_label_text_font_size = '20px'
    p.yaxis.major_label_text_font_size = '20px'

    # ajout de la légende
    legend = Legend(
        title="Légende:",
        items=[
            ("Number of Reads", [p.renderers[0]]), # type: ignore
            ("Palier de 10 reads", [p.renderers[1]]),# type: ignore
        ],
        label_text_font_size="16pt",
        title_text_font_size="16pt"     
    )
    
    legende_Prot_G = Label(x=641, y=20, x_units='screen', y_units='screen', text='G Protein',
        background_fill_color='yellowgreen', background_fill_alpha=0.5)
    
    legende_Prot_F = Label(x=828, y=20, x_units='screen', y_units='screen', text='F Protein',
     background_fill_color='orange', background_fill_alpha=0.5)
    
    p.add_layout(legende_Prot_G)
    p.add_layout(legende_Prot_F)
    p.add_layout(legend, 'below')

    show(p) # type: ignore

    export_png(p, filename=f"{result_folder}/4-Graphs/Graph_{file_id}_A.png") # type: ignore

    return p 


def calcul_percent(nom_fichier, file_id: str,result_folder:str):
	with open (nom_fichier,"r") as file :
            i = 0
                            
            nb_bases_couvertes_A = 0
            nb_bases_couvertes_B = 0

            ligne_A = 0
            ligne_B = 0

            total_reads = []
            num_ligne_virus = []

            for nb_ligne in file:
                                    
                i +=1
                          
                ligne = nb_ligne.split()

                lignee_ref = ligne[0]
                total_reads.append(ligne[3])
                num_ligne_virus.append(ligne[1])


                # prot G TYPE_A : 4659...5555
                
                if lignee_ref.startswith("typeA") :
                    ligne_A += 1
                    
                     
                    
                    if int(ligne[3]) >= 10:
                        
                        nb_bases_couvertes_A += 1 
               
                else:
                    
                    ligne_B += 1
                    
                    # comptage du nombre de bases couvertes par 10 reads pour la lignée A 
                    if int(ligne[3]) >= 10:
                        
                        nb_bases_couvertes_B += 1

            if ligne_A != 0 and ligne_B != 0:
                
                pourcentage_couverture_A = (nb_bases_couvertes_A*100)/(ligne_A)
                print(f"percentage coverage_A {pourcentage_couverture_A}")

                pourcentage_couverture_B = (nb_bases_couvertes_B*100)/(ligne_B)
                print(f"percentage coverage_B {pourcentage_couverture_B}")
                #type A
                if  pourcentage_couverture_A > pourcentage_couverture_B:
                    #x
                    num_ligne_virus = num_ligne_virus[0:ligne_A]
                    #y
                    total_reads = total_reads[0:ligne_A]

                    p = graphiqueA(ligne_A, num_ligne_virus, total_reads, nom_fichier, file_id,result_folder)

                # type B
                elif pourcentage_couverture_B > pourcentage_couverture_A:

                    total_reads = total_reads[ligne_A+1:ligne_A+1+ligne_B]
                    num_ligne_virus = num_ligne_virus[ligne_A+1:ligne_A+1+ligne_B]
                    
                        # Plot
                    p = graphiqueB(ligne_B, num_ligne_virus, total_reads, nom_fichier, file_id,result_folder)


def extraction_id(file_name: str, file_id: str)-> str:

    file_id_pattern_regex = r'pileup/([^/]+)\.pileup'

    file_id_match = re.search(file_id_pattern_regex, file_name)

    if file_id_match:
        file_id = file_id_match.group(1)
        
    

    return file_id


#======================================== Code ========================================#
def graphics_main():
    #==================== Variables Globales

    chemin = os.getcwd()
    disque = chemin[0]
    
    result_folder = "./.."
    nom_fichier_txt  = r"resume_pileup.txt"

    script_path = os.path.abspath(__file__)
    parent_dir = os.path.dirname(script_path)
    parent_dir = os.path.dirname(parent_dir)

    path_windows_input = f"{result_folder}/1-fastq/pileup"

    path_windows_output = f"{result_folder}/4-Graphs"

    extension = ".pileup"

    for fichier in os.listdir(path_windows_input):
        if fichier.endswith(extension):
            
                nom_fichier = os.path.join(path_windows_input, fichier)
                file_id = ""
                file_id = extraction_id(nom_fichier, file_id)
                
                file = open (nom_fichier,"r")
                content = file.readlines()

                calcul_percent(nom_fichier, file_id,result_folder)             

graphics_main()