import re
import math
import Bio
import Bio.PDB
from Bio.PDB import PDBParser
fh=open("Hydrophobicity.txt", "r")
j=0
hydro=[[0 for x in range(21)] for y in range(21)]
for line in fh.readlines():
#    print(line)
    data=line.rstrip().split()
    print(data)
    for i in range(1,len(data)):
        
        if(len(data[i])>0) and j>0:
            hydro[i][j]=10-int(data[i])
    j=j+1
fh.close()
#print(hydro)
AA_list={'A':1,'R':2,'N':3,'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20}
def dist_hydro(A,B):
    return hydro[AA_list[A]][AA_list[B]]
NameDict={}
NameDict['ALA']='A'
NameDict['ARG']='R'
NameDict['ASN']='N'
NameDict['ASP']='D'
NameDict['CYS']='C'
NameDict['GLU']='E'
NameDict['GLN']='Q'
NameDict['GLY']='G'
NameDict['HIS']='H'
NameDict['ILE']='I'
NameDict['LEU']='L'
NameDict['LYS']='K'
NameDict['MET']='M'
NameDict['PHE']='F'
NameDict['PRO']='P'
NameDict['SER']='S'
NameDict['THR']='T'
NameDict['TRP']='W'
NameDict['TYR']='Y'
NameDict['VAL']='V'
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure('1rnv', '1rnv.pdb')
total_dist=0
tsp_list=[]
for model in structure:
    for chain in model:
        residues1=chain.get_residues()
        for residue1 in residues1:
            if Bio.PDB.is_aa(residue1):
#        print(residue1.get_resname())
                dist2=100
                for atom1 in residue1:
                    coords1=atom1.get_coord()
            residues2=chain.get_residues()
            for residue2 in residues2:
                dist=100
                if Bio.PDB.is_aa(residue2) and residue1 !=residue2 and residue1.get_full_id()[3][1]!=(residue2.get_full_id()[3][1]+1) and residue1.get_full_id()[3][1]!=(residue2.get_full_id()[3][1]-1) and residue2 not in tsp_list:
#            if residue2.get_segid()<125 and residue1 != residue2:
                    for atom2 in residue2:
                        coords2=atom2.get_coord()
                        dist1=math.sqrt((coords1[0]-coords2[0])*(coords1[0]-coords2[0])+(coords1[1]-coords2[1])*(coords1[1]-coords2[1])+(coords1[2]-coords2[2])*(coords1[2]-coords2[2]))
#                        print(dist1)
#                        if Bio.PDB.is_aa(residue1) and Bio.PDB.is_aa(residue2):
#                            dist1=hydro[AA_list[NameDict[residue1.get_resname()]]][AA_list[NameDict[residue2.get_resname()]]]
                        if dist1<dist:
                            dist=dist1

                    dist3=dist
                    residueA=residue2
                    if dist3<dist2:
                        dist2=dist3
                        residueB=residueA
#            print(residue1)
#            print(residueB)
            if(Bio.PDB.is_aa(residue1) and Bio.PDB.is_aa(residueB)):
                tsp_list.append(residueB)
                print(NameDict[residue1.get_resname()]+" "+NameDict[residueB.get_resname()]+" "+str(dist2)+" "+str(hydro[AA_list[NameDict[residue1.get_resname()]]][AA_list[NameDict[residueB.get_resname()]]]/10)+" "+str(residue1.get_full_id()[3][1])+" "+str(residueB.get_full_id()[3][1]))
                total_dist=total_dist+hydro[AA_list[NameDict[residue1.get_resname()]]][AA_list[NameDict[residueB.get_resname()]]]

print("Total distance:"+str(total_dist))
