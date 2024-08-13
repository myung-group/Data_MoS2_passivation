from ase import Atoms
from ase.io import read, write, Trajectory
from ase.calculators.singlepoint import SinglePointCalculator

import numpy as np
import sys


## 1) Set filename
atoms_file = sys.argv[1] # any atoms object type or POSCAR...; can recognize atoms indices
filename = sys.argv[2]   # out.log of geometry optimization

## 2) Get Atoms object
atoms = read(atoms_file)

## 3) Get geometry optimization trajectory
FILE = open(filename)
file = FILE.readlines()
traj = Trajectory(f'{atoms_file[-4:]}.traj', 'w')

for i in range(len(file)): 
    if 'Geometry' in file[i].split():
        pos = []
        for j in range(len(atoms)):
            line = file[i+3+j].split()
            pos.append([float(line[1]), float(line[2]), float(line[3])])
        relaxed_atoms = Atoms(symbols=atoms.symbols, positions=np.array(pos))
            
    if "converged" in file[i].split() and "energy" in file[i].split() and "SCF" in file[i].split():
        energy = float(file[i].split()[-1])
        calc = SinglePointCalculator(atoms=relaxed_atoms, energy=energy)
        relaxed_atoms.calc = calc
        traj.write(relaxed_atoms)
traj.close()
