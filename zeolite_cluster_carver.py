from ase.io import vasp
import ase
import numpy as np
from math import log10, floor
from itertools import chain
import sys

def round_to_sig(x, sig):
    if abs(x)>0:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    else:
        return 0.000


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def union(lst1, lst2):
    final_list = list(set(lst1) | set(lst2))
    return final_list


# read in CONTCAR file
atoms_obj=vasp.read_vasp(str(sys.argv[1]))

# parameters
# radius of sphere around atoms
sphere_radius=4
# atoms to include
atom_id_list=[108,109,110,111,112,113,114,115,116]
# id of center atom around which coordinates are built
center_id=108




# get the set of all atoms within a x angstrom range from the atoms in atom_id_list
set_tot=[]
vec_tot=[]
for atom_id in atom_id_list:
    set1=[]
    vec1=[]
    for i in range(0,len(atoms_obj)):
        if atoms_obj.get_distance(atom_id, i, mic=True, vector=False)<sphere_radius:
            set1.append(i)
            vec1.append(atoms_obj.get_distance(atom_id, i, mic=True, vector=True))


    set_tot.append(set1)
    vec_tot.append(vec1)

combined=list(set(chain.from_iterable(set_tot)))


# find all atoms in combined that share a common atom
add_to_list = []
combined_bond = []
for i in range(0,len(atoms_obj)):
    if i not in combined:
        temp_id=[]
        for atom_id in combined:
            if atoms_obj.get_distance(i, atom_id, mic=True, vector=False)<2: 
                temp_id.append(atom_id)
        if len(temp_id)>1:
            combined_bond.append(temp_id)
            add_to_list.append(i)

# add coordinates centered on all atoms in list wihtin sphere radius
atom_coords=[]
for atom_id in atom_id_list:
    temp_coords=[]
    for atom_comb in combined:
        coords = atoms_obj.get_distance(atom_id, atom_comb, mic=True, vector=True)+atoms_obj.get_distance(center_id, atom_id, mic=True, vector=True)

        dist_list=[np.linalg.norm(coords-atoms_obj.get_positions()[this]+atoms_obj.get_positions()[center_id]) for this in atom_id_list]

        dist_check = [dist for dist in dist_list if dist < sphere_radius ]
        
        if len(dist_check) > 0:
            temp_coords.append([atoms_obj[atom_comb].symbol, round_to_sig(coords[0], 3), round_to_sig(coords[1], 3), round_to_sig(coords[2], 3)])
            
    if len(temp_coords) > 0:
        atom_coords.append(temp_coords)

combined_atom_coords = [item for sublist in atom_coords for item in sublist]


# remove repeated coordinates
new_k = []
for elem in combined_atom_coords:
    if elem not in new_k:
        new_k.append(elem)

with open('carved_cluster.xyz','w') as f:
    for coords in new_k:
        f.write('   ' + coords[0]+'   '+str(round_to_sig(coords[1], 3))+'     ' + str(round_to_sig(coords[2], 3))+'     '+str(round_to_sig(coords[3], 3))+'\n')

    for i, connected_atom_id in enumerate(add_to_list):
        coords = atoms_obj.get_distance(combined_bond[i][0], connected_atom_id, mic=True, vector=True)+atoms_obj.get_distance(center_id, combined_bond[i][0], mic=True, vector=True)
        f.write('   '+atoms_obj[connected_atom_id].symbol +'   '+str(round_to_sig(coords[0], 3))+'     ' + str(round_to_sig(coords[1], 3))+'     '+str(round_to_sig(coords[2], 3))+'\n')

