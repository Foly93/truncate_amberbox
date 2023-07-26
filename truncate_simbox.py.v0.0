#!/usr/bin/env python

import warnings
import subprocess
import numpy as np
import MDAnalysis as mda
from optparse import OptionParser
from MDAnalysis.transformations import center_in_box

def main():
    # ignore annoying warnings sparked by the <atoms>.write(stru.file)
    warnings.filterwarnings("ignore", message="Found no information for attr:")
    warnings.filterwarnings("ignore", message="Found missing chainIDs.")

    # initialise Optionparser from optparse as parser
    parser = OptionParser("Usage: ./truncate_simbox.py.v0.0 -p path/to/prmtop -x path/to/mdcrd4 -N 7000 --csel 'protein or nucleic or resname G5 or resname C3' -c 0.082")
    # add option to the parser, see help arguments for explanations
    parser.add_option("-p","--parmin", dest="top", type="string", default="prmtop",
        help="Input PRMTOP file.")
    parser.add_option("-x","--trajin", dest="trj", type="string", default="mdcrd.nc",
        help="Input trajectory.")
    parser.add_option("-N", "--nsolv", dest="N_slvnt", type="int", default=5000,
        help="Number of solvent residues to Keep from initial simbox.")
    parser.add_option("--csel", dest="cmplx_sel", type="string", default="protein or nucleic",
        help="selection string with MDAnalysis' syntax to select the macro molecules of interest usually a complex.")
    parser.add_option("--ssel", dest="slvnt_sel", type="string", default="resname WAT",
        help="selection string with MDAnalysis' syntax to select the system's solvent.")
    parser.add_option("-c","--conc", dest="conc", type="float", default=0.150,
        help="concentration of NaCl in the final simulation system in mol/L.")
    parser.add_option("-o", "--parmout", dest="out", type="string", default="prmtop_out",
        help="name of the output prmtop file for amber.")
    parser.add_option("-r", "--rstout", dest="rst", type="string", default="restrt_out",
        help="name of the output coords file for amber.")
    parser.add_option("-q", "--quiet", action="store_false", dest="verbose", default=True,
        help="Do not communicate any progress to the user.")
    # parse the arguments into options and args
    (options,args) = parser.parse_args()

    # instantiate universe from user input; hardcoded formats, yikes.
    u = mda.Universe(options.top, options.trj, format='NCDF', topology_format='PRMTOP')
    # declare atom groups for the complex and solvent from user input
    cmplx = u.select_atoms('({}) and not element H'.format(options.cmplx_sel))
    slvnt = u.select_atoms(options.slvnt_sel)

    # read the last frame of the trajectory
    for ts in u.trajectory[-1:]:
        # calculate box dimensions the center coords
        dim = ts.triclinic_dimensions
        box_cog = np.sum(dim, axis=0) / 2

        # translate complex atom 1 to the box COG to avoid a split complex, unwrap the complex and calc the COM
        u.atoms.translate(box_cog - cmplx[:1].center_of_mass())
        cmplx.unwrap(compound='fragments')
        cmplx_com = cmplx.center_of_mass()

        # translate the Complex COM to the box COG and wrap the residues into the box
        u.atoms.translate(box_cog - cmplx_com)
        slvnt.wrap(compound='residues')

        # calculate the solvent COMs and find closest N solvent molecules to the complex COM
        slvnt_xyz = slvnt.center_of_mass(compound='residues')
        close = np.argsort(np.linalg.norm(slvnt_xyz - cmplx.center_of_mass()[None,:], axis=1))[:options.N_slvnt]
        close_slvnt = u.select_atoms(options.slvnt_sel).residues[close].atoms

        # combine cmplx atom group and closest solvennt atom group and write it to PDB for tleap input
        trunc = cmplx + close_slvnt
        trunc.write("trunc.pdb")
        break

    # print user message if not supressed with -q and run pdb4amber to create a tleap ready PDB file
    if options.verbose: print("pdb4amber creates auxiliary PDB of truncated simbox and other files...")
    solv_sys = subprocess.run(['pdb4amber', '-i', 'trunc.pdb','-o','trunc_a.pdb'],
                               stderr=subprocess.STDOUT,
                               stdout=open('log.amb2pdb','w'),
                               check=True,
                               text=True)

    # write the necessary leap file via python, necessary since user input needs to be included
    with open('solv_box.leap','w') as f:
        f.write('source leaprc.RNA.OL3\n')
        f.write('source leaprc.water.tip3p\n')
        f.write('source leaprc.protein.ff14SB\n')
    
        f.write('system = loadpdb trunc_a.pdb\n')

        f.write('solvateoct system TIP3PBOX 0.0\n')
        f.write('addions system Na+ 0\n')
        f.write('addions system Cl- 0\n')
    
        f.write('savepdb system noions.pdb\n')
        f.write('saveAmberParm system noions.prmtop noions.rst7\n')
        f.write('quit')

    # print user message if not supressed by -q and run tleap to create solvated dummy system
    if options.verbose: print("(1/2) tleap solvates truncated simbox, output redirected into log.solv...")
    solv_sys = subprocess.run(['tleap', '-f', 'solv_box.leap'],
                               stderr=subprocess.STDOUT,
                               stdout=open('log.solv','w'),
                               check=True,
                               text=True)

    # instantiate universe of the solvated dummy system and calculate the number of Ions from user inputted conc
    ### maybe recall noions into dummy
    solvated = mda.Universe('noions.prmtop','noions.rst7',format='RESTRT')
    nW = solvated.select_atoms('resname WAT and element O').n_atoms
    nI = np.round(nW * 1 / (( 1 / options.conc / 0.018) + 2)).astype(int)

    # use python to write a tleap input file, necessary since user input needs to be included
    with open('genions.leap','w') as f:
        f.write('source leaprc.RNA.OL3\n')
        f.write('source leaprc.water.tip3p\n')
        f.write('source leaprc.protein.ff14SB\n')
    
        f.write('system = loadpdb trunc_a.pdb\n')
    
        f.write('solvateoct system TIP3PBOX 0.0\n')
        f.write('addions system Na+ 0\n')
        f.write('addions system Cl- 0\n')
        f.write('addions system Na+ {} Cl- {}\n'.format(nI,nI))
    
        f.write('savepdb system system.pdb\n')
        f.write('saveamberparm system {} {}\n'.format(options.out, options.rst))
        f.write('quit')
    
    # print user message if not supressed via -q and use tleap to create the final system with correct ion conc
    if options.verbose: print("(2/2) tleap generates amber files for truncated simbox, output redirected into log.ions...")
    ions_sys = subprocess.run(['tleap','-f','genions.leap'], 
                               stderr=subprocess.STDOUT,
                               stdout=open('log.ions','w'),
                               check=True,
                               text=True)

    # instantiate the new truncated system from the tleap output files and declare atom groups
    trunc_sys = mda.Universe(options.out, options.rst, format='RESTRT',topology_format='PRMTOP')
    trunc_slvnt = trunc_sys.select_atoms(options.slvnt_sel)
    trunc_cmplx = trunc_sys.select_atoms(options.cmplx_sel)
    # calculate the comdists betweenn solvent and complex and find the 100 furthest solvent molecules
    com_dists = trunc_slvnt.center_of_mass(compound='residues') - trunc_cmplx.center_of_mass()[None,:]
    far = np.argsort(np.linalg.norm(com_dists, axis=1))[-100:]
    far_slvnt = trunc_slvnt.residues[far].atoms

    # calculate the positional difference from the far solvents and the complex then the distance stored as matrix 
    surf_diffs = far_slvnt.center_of_mass(compound='residues')[:,None,:] - trunc_cmplx.positions[None,:,:]
    surf_dists = np.linalg.norm(surf_diffs.reshape(-1,3), axis=1).reshape(100,-1)
    # find the minimum value basically making it the surface distance between complex and outer solvent layer
    surf_min = surf_dists.min()
    # find the minimmum index and thus the minimum solvent residue ID
    min_idx = np.unravel_index(surf_dists.argmin(), surf_dists.shape)[0]
    min_resid = far_slvnt[min_idx].resid

    # inform the user about the surface-solv-layer distance and the water resid that this distance corresponds to
    print("The approximate distance to the outer solvation layer is {:.2f} Angstrom measured as the minimum distance of the complex surface to the 100 outermost water molecules. The distance corresponds to the water with resid {}.".format(surf_min,min_resid))
    

if __name__ == '__main__':
    main()
