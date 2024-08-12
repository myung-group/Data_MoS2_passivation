from density_of_states import *

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
            print(band_gap)
            break
    return band_gap

spin = True
spin_sum = True
# atoms_list = 72*['Au'] + 40*['S'] + 20*['W']

#atoms_list = ['W']+2*['S']
# particular_atoms = ['W']
# particular_orbital = ['p']
particular_atoms = None
particular_orbital = None
doscar, energy_range, fermi_ene = plot_dos_graph(filename = 'DOSCAR',
                                                 atoms_list = 'POSCAR',
                                                 particular_atoms = particular_atoms,
                                                 particular_orbital = particular_orbital,
                                                 spin = spin,
                                                 ene_range = None,
                                                 dos_range = None,
                                                 total_dos = True,
                                                 spin_sum = True,
                                                 plot = False)


sv_doscar, sv_energy_range, sv_fermi_ene = plot_dos_graph(filename = 'DOSCAR-sv',
                                                          atoms_list = 'POSCAR-sv',
                                                          particular_atoms = particular_atoms,
                                                          particular_orbital = particular_orbital,
                                                          spin = spin,
                                                          ene_range = None,
                                                          dos_range = None,
                                                          total_dos = True,
                                                          spin_sum = True,
                                                          plot = False)

get_band_gap(energy_range, doscar, fermi_ene)
get_band_gap(sv_energy_range, sv_doscar, sv_fermi_ene)

lw = 1.1
plt.figure(figsize=(14, 8))
if not spin_sum:
    plt.plot(energy_range-fermi_ene, dos, 'blue', linewidth=lw, label=r'1L-MoS$_{2}$')
    plt.plot(sv_energy_range-sv_fermi_ene, sv_dos, 'red', linewidth=lw, label=r'1L-MoS$_{2}$-S$_{v}$')
else: 
    plt.plot(energy_range-fermi_ene, doscar[0][0], 'blue', linewidth=lw, label=r'3L-MoS$_{2}$ (spin up)', zorder=4)
    plt.plot(energy_range-fermi_ene, -doscar[0][1], 'darkblue', linewidth=lw, label=r'3L-MoS$_{2}$ (spin down)', zorder=3)
    plt.plot(sv_energy_range-sv_fermi_ene, sv_doscar[0][0], 'red', linewidth=lw, label=r'3L-MoS$_{2}$-SV (spin up)', zorder=2)
    plt.plot(sv_energy_range-sv_fermi_ene, -sv_doscar[0][1], 'darkred', linewidth=lw, label=r'3L-MoS$_{2}$-SV (spin down)', zorder=1)
plt.axvline(x=0, color='grey', alpha=0.9, ls='--', lw=lw, label=r'$E_{Fermi}$')
plt.axhline(y=0, color='black', alpha=0.9, ls='-', lw=lw)
plt.xlabel(r'$E_{DFT}$ - $E_{Fermi}$ (eV)', fontsize=17)
plt.ylabel('Density of States (a.u.)', fontsize=17)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.xlim(-3, 3)
#plt.legend(loc='upper right', fontsize=14)
plt.show()

