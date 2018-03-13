import math
import Bio
import Bio.PDB
from Bio.PDB import PDBParser
parser=PDBParser(PERMISSIVE=1)
structure=parser.get_structure('1rnv', '1rnv.pdb')
for model in structure:
    residues1=model.get_residues()
    for residue1 in residues1:
        if Bio.PDB.is_aa(residue1):
#        print(residue1.get_resname())
            dist2=100
            for atom1 in residue1:
                coords1=atom1.get_coord()
            residues2=model.get_residues()
            for residue2 in residues2:
                dist=100
                if Bio.PDB.is_aa(residue2) and residue1 !=residue2:
#            if residue2.get_segid()<125 and residue1 != residue2:
                    for atom2 in residue2:
                        coords2=atom2.get_coord()
                        dist1=math.sqrt((coords1[0]-coords2[0])*(coords1[0]-coords2[0])+(coords1[1]-coords2[1])*(coords1[1]-coords2[1])+(coords1[2]-coords2[2])*(coords1[2]-coords2[2]))
#                        print(dist1)
                        if dist1<dist:
                            dist=dist1
#                    print(str(residue1)+" "+str(residue2)+" "+str(dist)+"\n")
                    dist3=dist
                    residueA=residue2
                    if dist3<dist2:
                        dist2=dist3
                        residueB=residueA
            print(str(residue1)+" "+str(residueB)+" "+str(dist2))
