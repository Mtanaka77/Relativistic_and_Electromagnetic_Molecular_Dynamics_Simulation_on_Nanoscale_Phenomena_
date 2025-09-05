## Relativistic and Electromagnetic Molecular Dynamics Simulation for Nanoscale Phenomena ##

For "Open Internet Access by Molecular Dynamics Simulations", a couple of 
various codes are shown in https://github.com/Mtanaka77/. The dynamics 
simulation codes are, 
(1) "Relativistic and Electromagnetic Molecular Dynamics Simulation 
for Nanoscale Phenomena", 
(2) "Macro-particle Simulation for Magnetic Reconnection", (3) "SIESTA-4.1 on Vector and 
Parallel Clusters", (4) "Water and Hydrate by Molecular Dynamics TIP5P Simulation", 
and (5) " AlmaLinux by MD and Siesta". 

This page is discussed on the relativistic and nanoscale molecular 
dynamics simulations, Computer Physics Communications (2019, Ref. 1). 
Updated files of a non-parallelized code @cnt3em_3pCa.f03 for 
Maxwell equations, and a parallelized code @a_cnt3-3p7Ca.f03 
for the Maxwell equations with the z direction are given in this directory. 
They are compiled by parameters and configuration files which are recently 
updated in Sep. to Nov., 2024. 
Each simulation code is suitable for the non-paralleled or parallelized 
Maxwell equations both with parallelized particles, which really 
depends on different architectures. 


### Molecular Dynamics Simulation: CGS Units and Necessary Files ###

The present molecular dynamics simulation code is implemented for 
relativistic and electromagnetic fields in three dimensions. 
It is applied to nanoscale particle phenomena such as nanotube accelerators. 
In the code, Maxwell equations are solved and momentum equations of 
relativistic particles are advanced in time. 
Four physical CGS units are used in this code: a_unit= 1.00d-08 cm, 
t_unit= 1.00d-15 sec, electron mass m_unit= 0.9110d-27 g and its charge 
e_unit= 4.8032d-10 esu. The mass of hydrogen, for example, is 1.6726d-24 g.

One needs files in the simulation that is,
1) @cnt3em_3pCa.f03 or @a_cnt3-3p7Ca.f03: Molecular dynamics simulation code, 
which is the non-paralleled or parallelized simulation code,
2) param_em3p7_Ca.h (non-parallel) or param_3p7_Ca.h (parallel): 
Common parameters of this simulation, 
3) Cntemp_config.STARTC: configuring parameters, 
4) p_config_ss.xyz_D150 and P135 of pellet atoms (ions and electrons): 
H, C and Au ions and electrons. 
The program is written in Fortran 2003/Fortran 2008 (the format is written 
in the same write line) and MPI of Ver.3 for parallelization.

The description of each subroutine and important comment lines of 
@cnt3em_3pCa.f03 and/or @a_cnt3-3p7Ca.f03, and the post-processed program 
@3ddisppC.f03, etc. (to be shown later), are written with comments 
and remarks of the simulation. 
Initial 80-100 lines of the file @cnt3em_3pCa.f03 or @a_cnt3-3p7Ca.f03 
are devoted to the title, references, summary of subroutines and 
remarks of the simulation code. 
The major subroutine /moldyn/ consists of the followings:
(i) the magnetic field is advanced, 
(ii) current density is calculated, and separation of the longitudinal
electric field is made using the Poisson equation,
(iii) the transverse electric field is advanced, 
(iv) the forces are calculated for the lontudinal electric field, and 
(v) positions and momenta of particles are advanced toward the next time step.
The step (iv) is most time-consuming than other steps..


### Gauss's Law, Courant Condition and Realistic Mass Simulation ###

It is noted here that the Ampere's law becomes inaccurate in 
the Cartesian coordinate space, (1/c)\partial{\textbf E}/\partial t=
rot{\textbf B} -(4\pi/c){\textbf J}, residing from the longitudinal electric field.
Thus, the Gauss's law, div{\textbf E}=4\pi q, is utilized to separate the 
transverse current on the total one (thelongitudinal and transverse currents), 
J_transvese= J -(J*EL)*EL/|ELl^2.
Thus, the Gauss's law must be solved in five-step interval for the 
discrete coordinate space (Refs. 1, 2 and 3).

On the other hand, the relativistic formulae of velocity and momeutum
/vec{v}= \vec{p}/(sqrt(m^2 +(px^2 +py^2 +pz^2)/c^2)) is valid
in the nanoscale cases (Ref. 1 and Ref. 2).
We stress on the use of real mass of hydrogen, carbon and gold atoms, 
and electrons.
Also, all the explicit simulation code must satisfy the Courant condition,
that is, Dx(length)/Dt(time step) > c, the speed of light.
Otherwise, simulation steps overflow quite shortly.


### Simulation on Nanoscales ###

A simulation of the nanotube accelerator is set up by putting pellets of H, C 
and Au atoms and associated electrons at null velocity. 
Electromagnetic monochromatic waves at the wavelength 800 nm are travelling from 
the negative direction toward the origin and then go out to the positive direction. 
The pellet includes a cold pellet of hydrogen, carbon and gold, and electrons.
The pellets at the origin are irradiated by these waves and are ejected to ion 
perpenducular and electron parallel directions toward an open space. 
The fine particle ejection in 3-D and x-y plots are shown newly in cntemp.77Ca-2024a.pdf. 

In the second case, electromagnetic circulary-polarized waves are injected periodically.
Thirdly, complex double pellets are excited by the electromagnetic waves.  

The final energies for laser intensity 10^22 W/cm^2 are around 30-40 MeV in 20-40 fs, 
and they are discussed in the latter half of the Computer Physics Commun. in 2019 (Ref. 1).
Also, it is shown by animation movies at my homepage, https://physique.isc.chubu.ac.jp/,
and/or http://www1.m4.mediacat.ne.jp/dphysique/.


### Execution Scripts ###

(1) Linux (PGI Fortran): MPI and FFTW; configure, make, make install for installation..

  %  mpich-4.0.2: ./configure --prefix=/opt/pgi/mpich-4.0.2 2>&1 | tee conf.txt

  % fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/pgi/fftw3

Compilation by mpif90: 

The script 'mpif90 @cnt3em-3pCa.f03' needs param_em3p7_Ca.h, Cntemp_config.STARTC, p_config_ss.xyz_P135 and p_config_ss.xyz_D150.

  % mpif90 -byteswapio -mcmodel=medium -fast @cnt3em3_3pCa.f03 -I/opt/pgi/fftw3/include -L/opt/pgi/fftw3/lib -lfftw3

(2) Linux (gcc gfortran): MPI and FFTW; configure, make, and make install for installation.

  %  mpich-4.0.2: ./configure --prefix=/opt/mpich-4 2>&1 | tee conf.txt

  % fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/fftw3

  % mpif90 -mcmodel=medium -fPIC @a_cnt3-3p7Ca.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3

  The parallelization simulation code is utilized, for example. 
Two different fortrans, PGI and gfotran, are incompatible with processors of the same name but different styles.
One may want to use about 50 processors at least, where the z-direction is exactly divided by parallelization. 

  % Execution by mpiexec (with lots of co-processors): Only for a test case. % mpiexec -n 6 a.out &
 
### Post-processing Simulation Analysis ###

To analyze simulation results, this program provides the post-processing tool. 
They are named @3ddisppC.f03, @3dfdisp.f03, and @3dfv2C.f03. 
The velocity distributions in parallel and perpendicular directions, @3dfdispC.f03, are plotted 
in sequential times of ions and electrons. The @3ddisppC.f03 program is time sequential plots 
of H, C, Au ions and electrons from side and top views with energy histories as well at the end of the run. 

The @3dfv2C.f03 plot makes a kind of animation (100 plots, for example) which tells different distributions 
of parallel and perpendicular directions of axis.
These graphic outputs by PDF files are shown on the PC screen, either cntemp.77Cfb.pdf or 
cntemp.77Csa.pdf. 

### Parallelization of Electromagnetic Fields ###

The heavy load of many particles is generally divided on parallel processors which is easily coded. 
On the other hand, the electromagnetic fields are parallelized for one dimension 
(the z-direction) where the long axis of pellets is open to eject heavy ions in that direction.
Small segments 'do n=nz1(ipar),nz2(ipar)' in the z-direction of each processor node 
are written in the individual electric and magnetic fields of the Maxwell equations.   
Note: Parellization of the Maxwell equation depends on the machines. The Fujitsu by FX100
and maybe the current processos are excellent on parallelization of the Maxwell equations.

Finally, the machine run time depends on the physical time and cpu's architecture. 
For the elapsed time of parallel simulation, it executed in 3.2 sec/step for 
52 ranks and 4.0 10^5 particles (equal to 1 fs for elapsed 1.8 hours) by 
Fujitsu FX100 Supercomputer (2019).

### References: ###

1. M. Tanaka and M. Murakami, Relativistic and electromagnetic molecular dynamics simulations for a carbon-gold nanotube accelerator, Computer Physics Communications, 241, 56-63 (2019).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Tanaka, Macro-EM particle simulation method and a study of collisionless magnetic reconnection, Comput.Phys.Commun., 87, 117-138 (1995).

4. M. Murakami and M. Tanaka, Generation of high-quality mega-electron volt proton beams with intense-laser-driven nanotube accelerator, Applied Phys. Letters, 102, 163101 (2013).

