import pandas as pd
from ase.io import read

from density_of_states import *


data = {}

## 1) Set file path
mos2_DOSCAR = "../../Data/GeometryRelax_Dimer_DOS/3L_MoS2/3L_MoS2/SCF/DOS/DOSCAR"
mos2_sv_DOSCAR = "../../Data/GeometryRelax_Dimer_DOS/3L_MoS2/3L_MoS2_sv/SCF/DOS/DOSCAR"

mos2_POSCAR = "../../Data/GeometryRelax_Dimer_DOS/3L_MoS2/3L_MoS2/SCF/DOS/POSCAR"
mos2_sv_POSCAR = "../../Data/GeometryRelax_Dimer_DOS/3L_MoS2/3L_MoS2_sv/SCF/DOS/POSCAR"

## 2) Extract from 3L-MoS2 DOS
atoms = read(mos2_POSCAR)
pos = atoms.get_positions()
sheet_names = ['total', 'top', 'middle', 'bottom']
atoms_indices_list = [list(np.arange(0, len(pos), 1)),
                   [i for i in range(len(atoms)) if pos[i][2] > 22.],
                   [i for i in range(len(atoms)) if 15. < pos[i][2] < 22.],
                   [i for i in range(len(atoms)) if pos[i][2] < 15.]
                   ]
for i in range(len(atoms_indices_list)):
    atoms_indices = atoms_indices_list[i]
    doscar, energy_range, fermi_ene = plot_dos_graph(filename = mos2_DOSCAR,
                                                    atoms_list = mos2_POSCAR,
                                                    particular_atoms = None,
                                                    particular_orbital = None,
                                                    spin = True,
                                                    total_dos = False,
                                                    atoms_indices = atoms_indices,
                                                    plot = False)
    p_dos = np.zeros((2, 3000))
    for j in range(len(doscar)-1):
        p_dos = np.add(p_dos, doscar[j])
    p_dos = p_dos[0] + p_dos[1]
    data[f'{sheet_names[i]}'] = {'E-fermi': energy_range-fermi_ene, 'DOS': p_dos}

## 3) Extract from 3L-MoS2-SV DOS
atoms = read(mos2_sv_POSCAR)
pos = atoms.get_positions()
sheet_names = ['sv-total', 'sv-top', 'sv-middle', 'sv-bottom']
atoms_indices_list = [list(np.arange(0, len(pos), 1)),
                   [i for i in range(len(atoms)) if pos[i][2] > 22.],
                   [i for i in range(len(atoms)) if 15. < pos[i][2] < 22.],
                   [i for i in range(len(atoms)) if pos[i][2] < 15.]
                   ]
for i in range(len(atoms_indices_list)):
    atoms_indices = atoms_indices_list[i]
    doscar, energy_range, fermi_ene = plot_dos_graph(filename = mos2_sv_DOSCAR,
                                                    atoms_list = mos2_sv_POSCAR,
                                                    particular_atoms = None,
                                                    particular_orbital = None,
                                                    spin = True,
                                                    total_dos = False,
                                                    atoms_indices = atoms_indices,
                                                    plot = False)
    p_dos = np.zeros((2, 3000))
    for j in range(len(doscar)-1):
        p_dos = np.add(p_dos, doscar[j])
    p_dos = p_dos[0] + p_dos[1]
    data[f'{sheet_names[i]}'] = {'E-fermi': energy_range-fermi_ene, 'DOS': p_dos}

## 3) Save as ".xlsx" format
with pd.ExcelWriter('density_of_states.xlsx', engine='openpyxl') as wt:
    for sheet_name, df_data in data.items():
        df = pd.DataFrame(df_data)
        df.to_excel(wt, sheet_name=sheet_name, index=False)

