
## Description

 This repository contains the scripts to prepare the b12 antibody from
PDB ID 1hzh.pdb for molecular dynamics (MD) simulation using the
CHARMM36m energy function.

## Summary  

The scripts download the required PDB file, process it to extract
coordinates, connectivity, etc. Minimization, solvation and ionization
is performed using fully open source software, and the resulting files
are saved in the local 'struc' directory. Modeling of the missing loop
in the HC2 antibody chain was performed with Truman's CRIMM code,
available freely from the Charlie Brooks Lab @MSU. The antibody is
solvated in a thin shell of water solvent.  More solvent can be added in
,e.g., VMD if simulations in periodic cells are desired.


## Required software:
 Name:	|	Purpose:	|	Source:
 ---------------|----------------|----------------------------------
 GNU Octave|	Structure preparation|	Linux distribution/package
 VMD (v. >=1.9.3)	|	Structure preparation|	(UIUC/TCBG)
 NAMD2	|	Coordinate minimization|	(UIUC/TCBG)
 CRIMM	|	Model missing residues|	(Brooks group @ MSU)
 Python	|	To run CRIMM|		Linux distribution
 wget, gunzip|	To download & uncompress|	Linux distribution


## How to Run

The steps to create the structure topology file (PSF) coordinate file
(PDB) and combined parameter files (PAR/TOP) are collected in the bash
script prep.

The procedure has been tested on an Arch Linux workstation, where

`./prep >&prep.out`

was used to create the simulation files.

The complete output from prep upon successful execition is
provided for reference in the tarball b12-prepare-output.txz

An example temperature-accelerated MD simulation (TAMD) using an earlier
version of the b12 structure files performed with the OpenMM MD code can
be found at the link below:

https://drive.google.com/drive/folders/1v4ktKVx4JepigmvBuUEyX7T0H4CI1GuS?usp=sharing
