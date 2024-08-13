for i in 4ABT 4MBT 4FBT PFBT; do

cd $i/
python ../dipole_moment.py POSCAR-$i > log-$i.out
python ../pyscf_traj.py POSCAR-$i log-$i.out
cd ../

done
