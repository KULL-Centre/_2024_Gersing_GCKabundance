### script to update .gro file with Na ion information
### NA_wt coordinates should be updated depending on the wt crystal structure used
### (so NA_wt differs for closed and superopen crystal structure files) 

# open .gro file made from pdb without HETATMs
with open('processed_original.gro', 'r') as f:
    lines = f.readlines()
    
# get total number of atoms and add 1    
atoms_total = int(lines[1][1:-1]) + 1
atoms_total_str = f' {atoms_total}\n'

# update total number of atoms
lines[1] = atoms_total_str

# get last line of file
lines_last = lines[-1]

# add ion coordinates copied from tmp.gro
NA_wt = f'  600NA      NA {atoms_total}   3.687   1.059   4.647\n'
lines[-1] = NA_wt

# append original last line
lines.append(lines_last)

# write new file
f = open("processed.gro", "a")
f.writelines(lines)
f.close()

# update .top file with Na ion information
# also add CA position restraint information

with open('topol_original.top') as f:
    lines = f.readlines()

lines[-21] = '\n; Include Position restraint file (CA atoms) \n#ifdef POSRES-CA \n#include "ca_posre.itp" \n#endif \n \n'
    
lines.append('NA                  1\n')

f = open("topol.top", "a")
f.writelines(lines)
f.close()