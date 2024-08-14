from pyscf import gto, scf
from pyscf.tools import cubegen, molden
from pyscf.geomopt.geometric_solver import optimize

from ase.io import Trajectory, read
import numpy as np
import sys


def ase_to_pyscf(atoms, basis):
    pos = atoms.get_positions()
    atoms_and_coords = f'{atoms[0].symbol} {pos[0][0]} {pos[0][1]} {pos[0][2]};'
    for i in range(len(atoms)-1):
        atoms_and_coords += f' {atoms[i+1].symbol} {pos[i+1][0]} {pos[i+1][1]} {pos[i+1][2]};'
    mol = gto.M(atom=atoms_and_coords[:-1], basis=basis)
    return mol


## 1) Set filename
fout_name = 'dipole_moment.out'
atoms_file = sys.argv[1]

## 2) Set Mole object
atoms = read(atoms_file)
mol = ase_to_pyscf(atoms, 'augccpvtz')

## 3) Geometry optimize 
mf = scf.RKS(mol, 'b3lyp')
mol_eq = optimize(mf)

## 4) SCF calculation
mf = scf.RKS(mol_eq, 'b3lyp')
mf.kernel()

## 5) Dipole moment calculation
fout = open(fout_name, 'a', 1)
dm = mf.dip_moment()
print(f'dipole moment of {atoms_file}: {dm}\n {np.sqrt(dm[0]**2+dm[1]**2+dm[2]**2)}', file=fout)
fout.close()

## 6) Generate cube files (density and mep)
cubegen.density(mol, f'{atoms_file[-4:]}_density.cube', mf.make_rdm1(), nx=200, ny=200, nz=200)
cubegen.mep(mol, f'{atoms_file[-4:]}_pot.cube', mf.make_rdm1(), nx=200, ny=200, nz=200)

