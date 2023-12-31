# Keywords
- amber simulations
- remove excess solvent/excess water
- truncate simulation box
- change box-dimensions/box-boundaries

# Preliminaries
Python packages to install 
- MDAnalysis
- Ambertools

# truncate_amberbox
Often in biophysical simulations a lot of excess water doesn't need to be simulated after complex formation, but deleting all of it would destroy the nice water structure around eventual macromolecules. This python script keeps a user specified amount of solvent molecules around the complex and reports the approx. distance to the outer solv. shell.

The truncation process only keeps the specified amount of solvent molecules which are in closest vicinity to the macromolecules of interest referred to as 'complex'. Solvent molecules (does not have to be water) that are further out are removed, resulting in a spherical solvation shell around the complex. Of course spherical solvation shells are not a thing, therefore tleap is used from the python program to create an octahedral simulation box. The forcefield that is used by tleap is hardcoded, therefore it can only be adjusted by manually changing the corresponding lines in truncate_simbox.py.  

In the examples/ are dummy output files from an amber simulation run that can be truncated interactively with the jupyter notebook that can be found also in examples/. The jupyter notebook contains the same code as the python script, with the difference that the notebook does not parse user input. Apart of the simulation output and the notebook, all files that would be generated from a successful execution are also already available and would be updated upon executing every jupyter cell without adjustment.

Usage: `./truncate_simbox.py -p path/to/prmtop -x path/to/mdcrd -N 7000 --csel 'protein or nucleic or resname G5 or resname C3' -c 0.082 -s tip3p`

The system that is then outputted is somewhat non-standard and needs careful treatment when equilibrating it with pmemd.cuda or sander. What seemed to work is to skip the NVT equilibration and instead run an NPT with a decreased time step. Thus: EM -> NPTSLOW -> NPT -> PRE -> PROD. The timestep can be set to around 0.1 fs for NPTSLOW.
# To-Do
Make the molecule identification use Residue Connections and not Atom bonds.
