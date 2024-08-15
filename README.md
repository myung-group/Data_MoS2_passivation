This repository accompanies the paper " **Back-End-of-Line-Compatible Passivation of Sulfur Vacancies in $MoS_2$ Transistors with Benzenethils** " by Hacksoon Jung, Mingyu Kim, Gi Beom Sim, Youngwoo Lee, Sumin Hong, Hyeonho Gu, Jaehyun Lee, Donghyeop Lee, Sanghyun Lee, Taoyu Zou, Kibum Kang, Chang Woo Myung, Yong-Young Noh, and Jimin Kwon.

* ```Data/``` contains inputs and outputs of 1) Geometry relaxation, 2) Density of states (DOS), and 3) Improved dimer calculation of MoS2 passivation by the substituted thiolbenzene (TB) molecules (PFBT, 4-ABT, 4-MBT, 4-FBT). Also, contains inputs and outputs of 4) Dipole moment, 5) Electron density, and 6) Molecular electrostatic potential calculations of the TB molecules.
  * To calculate dipole moment, electron density, and molecular electrostatic potential of TB molecules, you can use ```Data/DipoleMoment_MEP/*/dipole_moment.py```
  * For example,
    ```
    python dipole_moment.py POSCAR-PFBT > PFBT.out
    ```
  * Some of the outputs were uploaded as compressed (.zip) files because the data is huge. (DOSCAR, density.cube, and mep.cube). 

* ```Figures/``` contains figures of the paper that obtained from ```Data/```
* ```Scripts/``` contains Jupyter notebook and Python scripts that can reproduce figures.
  * The DOS of $MoS_2$ and defective $MoS_2$ can be plotted using Jupyter notebook code ```Scripts/Density_of_states/Plot_DOS.ipynb```
  * Before plotting the DOS, the compressed files ```/Data/*/SCF/DOS/doscar.zip``` should be extracted as "DOSCAR".
