from pyscf import gto, scf
from pyscf.tools import cubegen
from pyscf.geomopt.geometric_solver import optimize

import numpy as np
import sys


## 1) Set filename
atoms_file = sys.argv[1]

## 2) Set Mole object
mol = gto.M(atom=atoms_file, basis='augccpvtz')

## 3) Geometry optimize 
mf = scf.RKS(mol, 'b3lyp')
mol_eq = optimize(mf)

## 4) SCF calculation
mf = scf.RKS(mol_eq, 'b3lyp')
mf.kernel()

## 5) Dipole moment calculation
dm = mf.dip_moment()

## 6) Generate cube files (density and mep)
cubegen.density(mol, f'{atoms_file[:4]}_density.cube', mf.make_rdm1(), nx=200, ny=200, nz=200)
cubegen.mep(mol, f'{atoms_file[:4]}_pot.cube', mf.make_rdm1(), nx=200, ny=200, nz=200)

