import sys
import os
import zipfile
import linecache
import numpy as np
import matplotlib.pyplot as plt
from ase.io import read


def orbit_dos(orbit, dos, nedos, spin):    
    if spin:
        if orbit == 's':
            orbit = [[1], [2]]
        elif orbit == 'p':
            orbit = [[3, 5, 7], [4, 6, 8]]
        elif orbit == 'd':    
            orbit = [[9, 11, 13, 15, 17], [10, 12, 14, 16, 18]]
        elif orbit == 'tot':
            orbit = [[1], [2]]
        orbit_up = np.array(nedos*[0])
        orbit_down = np.array(nedos*[0])
        for i in range(len(dos)):
            for j in range(len(orbit[0])):
                orbit_up = np.add(orbit_up, np.array(dos[i][:,orbit[0][j]]))
        for i in range(len(dos)):
            for j in range(len(orbit[1])):
                orbit_down = np.add(orbit_down, np.array(dos[i][:,orbit[1][j]]))
        orbit_dos = np.array([orbit_up, orbit_down])
    else:
        if orbit == 's':
            orbit = [1]
        elif orbit == 'p':
            orbit = [2, 3, 4]
        elif orbit == 'd':    
            orbit = [5, 6, 7, 8, 9]
        elif orbit == 'tot':
            orbit = [1]
        orbit_tot = np.array(nedos*[0])
        for i in range(len(dos)):
            for j in range(len(orbit)):
                orbit_tot = np.add(orbit_tot, np.array(dos[i][:,orbit[j]]))
        orbit_dos = orbit_tot
    return orbit_dos


def plot_dos_graph(filename, atoms_list, ene_range=None,
                   particular_atoms=None, particular_orbital=None, 
                   spin=True, dos_range=None, total_dos=False, 
                   spin_sum=False, plot=True, atoms_indices=None):
    
    if type(atoms_list) == str:        
        atoms = read(atoms_list)
        POSCAR = open(atoms_list)
        poscar = POSCAR.readlines()
        symbols = poscar[5].split()
        num_of_atoms_list = poscar[6].split()
        print('\n', symbols, num_of_atoms_list)
        atoms_list = [symbols[i] for i in range(len(symbols)) \
                                 for _ in range(int(num_of_atoms_list[i]))]
    
    default_orbit_list = ['s','p','d']
    if total_dos:
        dos_list = ['tot']
    else:
        dos_list = ['tot']+atoms_list

    file = open(filename)
    line_5 = file.readlines().copy()[5].split()
    EMAX = float(line_5[0])
    EMIN = float(line_5[1])
    nedos = int(line_5[2])
    fermi_ene = float(line_5[3])
    print("\n EMAX", "\tEMIN", "\tNEDOS", "\tE_FERMI")
    print(f' {EMAX:.4g}\t{EMIN:.4g}\t{nedos:.4g}\t{fermi_ene:.4g}')
    
    for i in range(len(dos_list)):
        dos_list[i] = np.loadtxt(filename, skiprows=i*nedos+i+6, max_rows=nedos)
    
    if particular_atoms is not None:
        element_list = particular_atoms
    else:
        element_list = []
        for i in atoms_list:
            if i not in element_list:
                element_list += [i]
    
    if particular_orbital is not None:
        orbit_list = particular_orbital
    else:
        orbit_list = default_orbit_list
    
    doscar = []
    label = []
    element_list = element_list.copy()
    if atoms_indices != None:
        print("\n PARTIAL DOS INFO.")

    if not total_dos:
        for i in element_list: 
            if atoms_indices == None:
                num = atoms_list.index(i)+1 
                max_num = atoms_list.count(i)
                dos_ls = dos_list[num:num+max_num]
            else:
                num = atoms_list.index(i)+1 
                max_num = atoms_list.count(i)
                index_list = list(np.arange(num, num+max_num, 1))
                dos_ls = [dos_list[atoms_indices[j]+1] for j in range(len(atoms_indices)) \
                          if atoms_indices[j]+1 in index_list]
                print(f' {i}:\t{len(dos_ls)}')
            for j in orbit_list:
                label.append(f'{i}-{j}')
                doscar.append(orbit_dos(orbit = j, 
                                        dos = dos_ls, 
                                        nedos = nedos,
                                        spin = spin))
        label.append(f'Total DOS')
        doscar.append(orbit_dos(orbit = 'tot', 
                                dos = dos_list[0:1], 
                                nedos = nedos,
                                spin = spin))
    else:
        doscar.append(orbit_dos(orbit = 'tot', 
                                dos = dos_list, 
                                nedos = nedos,
                                spin = spin))

    print(f'\n {label}\n')
    if plot:
        plt.figure(figsize=(10,7))
        for i in range(len(doscar)):
            if i < len(doscar)-1:
                if spin:
                    if spin_sum:
                        dos = doscar[i][0] + doscar[i][1]
                        plt.plot(dos_list[0][:,0]-fermi_ene, dos,linewidth=1.2, label=label[i])
                    else:
                        plt.plot(dos_list[0][:,0]-fermi_ene, doscar[i][0],linewidth=1.2, label=label[i])
                        plt.plot(dos_list[0][:,0]-fermi_ene, -doscar[i][1],linewidth=1.2, label=label[i])
                else:
                    plt.plot(dos_list[0][:,0]-fermi_ene,doscar[i],linewidth=1.2,label=label[i])
            else:
                if spin:
                    if spin_sum:
                        dos = doscar[i][0] + doscar[i][1]
                        plt.plot(dos_list[0][:,0]-fermi_ene, dos, 'black', linewidth=1.2, label='Total DOS')
                    else:
                        plt.plot(dos_list[0][:,0]-fermi_ene, doscar[i][0],'black',linewidth=1.2,label='Tot/spin_up')
                        plt.plot(dos_list[0][:,0]-fermi_ene, -doscar[i][1],'black',linewidth=1.2,label='Tot/spin_down')  
                else:
                    plt.plot(dos_list[0][:,0]-fermi_ene, doscar[i],'black',linewidth=1.2,label='Total DOS')

            if dos_range is not None:
                plt.ylim(dos_range)  
            if ene_range is not None:
                plt.xlim(ene_range)
        plt.axvline(x=0, color='blue', alpha=0.9, ls='--', label=r'$E_{Fermi}$')
        plt.xlabel(r'$E_{DFT}$ - $E_{Fermi}$ (eV)', fontsize=13)
        plt.ylabel('Density of States', fontsize=13)
        plt.legend(loc='upper right', fontsize=12)
        plt.show() 

    energy_range = dos_list[0][:,0]
    return doscar, energy_range, fermi_ene


def unzip_file(zip_path):
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        file_name = zip_ref.namelist()[0]
        extracted_path = zip_ref.extract(file_name, 
                                         os.path.dirname(zip_path))
    return extracted_path


def get_band_gap(energy_range, dos, fermi_ene):
    dos_dat = dos[0][0] + dos[0][1]
    energy_r = energy_range.copy()
    energy_r -= fermi_ene
    band_tag = False
    for i in range(len(dos_dat)):
        if not band_tag and energy_r[i] >= 0:
            ene_0 = energy_r[i]
            band_tag = True
        if band_tag == True and dos_dat[i] != 0 and energy_r[i] - ene_0 > 0.5:
            band_gap = energy_r[i] - ene_0
            print(f' Bandgap: {band_gap:.5g} eV')
            break
    return band_gap
    
