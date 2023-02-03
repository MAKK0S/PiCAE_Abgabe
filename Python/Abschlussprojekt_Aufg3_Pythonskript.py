# -*- coding: utf-8 -*-

from abaqus import *
from abaqus import *
from abaqusConstants import *
from odbAccess import *
from abaqusConstants import *
from material import *
from odbMaterial import *
from odbSection import *

import os
from datetime import datetime

date_prefix = datetime.today().strftime("%Y%m%d")

# Schaumdichte ?ndern
def setFoamDensity (FoamDensity,nameofmodel):
    Reference = [135,1000.0,0.4,0.4,0.4,10] # [Referenzdichte 135 kg/m^3,E-Modul, nu, k_m, g_m, tau_m]

    if FoamDensity == 45 or 90 or 135 or 180:
        E = Reference[1]*FoamDensity*FoamDensity/Reference[0]/Reference[0]
    else:
        print('Wrong density: Density has to be 45,90,135 or 180')

    PronyParameters = (E,Reference[2],Reference[3],Reference[4],Reference[5]) # PronyParameters muss tuple sein

    mdb.models[nameofmodel].materials['Foam'].UserMaterial().setValues(mechanicalConstants=PronyParameters) # PronyParameters muss tuple sein


# Schaumdicke ändern
#def setFoamThickness (FoamThickness,nameofmodel):
#    wholeModelRegion=regionToolset.Region[faces=wholeModel]
#    #wholeModelSet=mdb.models[nameofmodel].parts['b-pillar'].Set['wholeModel']
#    mdb.models[nameofmodel].parts['b-pillar'].compositeLayups['CompositeLayup-1'].CompositePly(thickness=10.0,region=wholeModelRegion,material='Foam',plyName='ply-2',orientationType='ANGLE_0',thicknessType='SPECIFY_THICKNESS')

#    #mdb.models[nameofmodel].parts['b-pillar'].compositeLayups['CompositeLayup-1'].plies[1]

# Volumenberechnung

# Preisberechnung
# volumen [mm^3], ThicknessFoam [mm], densityFoam[kg/mm^3]
def cost(volume,ThicknessFoam,densityFoam):
    cost = 0.0
    massFoam = volume*ThicknessFoam/5.0*densityFoam # 5.0 ist die Gesamtdicke in mm
    massCFK = volume*(5.0-ThicknessFoam)/5.0*0.00000016 # Dichte von CFK 1.6*10^-6 kg/mm^3
    costFoam = massFoam*250.0
    costCFK = massCFK*100.0
    cost = costFoam + costCFK
    return cost

# Bestimmen der Leichtbaukennzahl

def LBK (mass,defromation,cost):

    return

# Automatisierte Vernetzung und Joberstellung:
def autojob(nameofmodel, nameofpart, MeshSize): #Eingabe: ('Modellname','Partname',Meshgröße')
    root_model = mdb.models[nameofmodel] #zuweisen des Modells
    
    #hier beginnt der groessen check
    listofCoordinates= root_model.parts[nameofpart].faces.getBoundingBox() # Bestimmen der Bauteilgrenzen in x,y,z (Quader um Bauteil)
    #listofCoordinates format: {'high' [0.0,1.0,2.0], 'low' [-5.0,-10.251,-5]}
    
    runner = 0
    dimesion = [0,0,0]
    while (runner<=2): # Bestimmt aus Grenzen die Kantenl?ngen des Quaders
        dimesion[runner]= listofCoordinates['high'][runner]- listofCoordinates['low'][runner]
        runner +=1
    #dimension= [5,255, 10]   
    small = min(dimesion) # bestimmen der kleinsten Kantenl?nge des Quaders
    
    if (small> MeshSize):
        root_model.parts[nameofpart].seedPart(size= MeshSize) #abaqus interne befehle zum beseeden eines Teils
        root_model.parts[nameofpart].generateMesh() # Mesh wird erstellt
        newjob = mdb.Job(name ="{}_NetzKonvStudie_it".format(date_prefix), model=root_model,) # automatisches erstellen eines Jobs
        newjob.submit() # job mit vorbestimmter name weitergeleitet
    else:
        print("Mesh size is too big for part dimensions") # wenn angaben fuer mesh Kantenlaengen zu gross Fehler

def evaluatemaxDeformation(input_odb):
	data_i=input_odb.steps["Step-1"].frames[1].fieldOutputs['U']
	number_nodes=len(input_odb.steps['Step-1'].frames[1].fieldOutputs['U'].values)
	maximales_u = 0.0

	for ii in range(number_nodes):
		u=data_i.values[ii].data
		
		if maximales_u < abs(max(u)): # suche nach hoechstem Wert (Betrag) waehrend die forschleife laeuft
			maximales_u = abs(max(u))
			nodeID=data_i.values[ii].nodeLabel # mitspeichern der lokation

	with open('maximales_u.txt', 'a') as f:
		f.write(str(maximales_u) +" an Node:" + str(nodeID) + " ")# Sicherung der Ergebnisse
	return(maximales_u)



# main-Methode
def main():


    if __name__=="__main__":
        main()