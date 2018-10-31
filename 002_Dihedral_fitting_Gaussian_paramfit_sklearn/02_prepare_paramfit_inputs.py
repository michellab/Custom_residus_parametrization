d=' \
############################ \
Prepare Paramfit inputs input from gaussians input files and the mol2 file used to produce those \
!!!! It HAS TO BE THE SAME MOL2 FILE OR AT LEAST ATOMS HAVE TO BE IN THE SAME ORDER! !!!!! \
Call process_gaussians_files which produces  \
        the stationary structures from gaussians,\
        -energy files \
        - listof files (succesful gaussians calculations) \
this script produce: \
    - the initial prmtop with gaff2 parameters, \
            KNOWN BUG :  you will have to remove  atoms types in the initial mol2 file...  \
    - mdcrd files \
    - corresponding formated paramfit input files \
############################\
Author : Marie Bluntzer \
############################\
'

from  process_gaussians_files  import extract_stationary_structures
from  get_all_dihedral import find_dihedrals
import MDtraj as md
import shutil


def get_gaff2_parmcheck(input):

     os.system('parmchk2 -f mol2 -i %s -o gaff2.frcmod  -s  gaff2' %(input))

def make_traj(prefix):
    list_pdb=open(prefix +'-listfiles')
    t = md.load(list_pdb.readline(), top='input.prmtop')
    for filename in list_pdb.readline():
        md.join(filename, check_topology=True, discard_overlapping_frames=False)
        t.save_mdcrd(prefix+ '-traj.mdcrd' )
        return t.n_frames

def pepare_paramfit_job_files (prefix, nframes) :
    ## This fonction prepares
    # - fit_K.in for paramfit : find constant to fit Kcal to Hartree
    # fit_P.in for torsion

    temp_fitK= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'fit_K_A.in' , 'r')
    temp_fitP= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'fit.in' , 'r')
    new_fitK= open(prefix + '-K-fit.in','w')
    new_fitP= open(prefix + '-P-fit.in','w')
    for line in temp_fitK.readline():
        new_fitK.writelines(line.replace('$PREFIX', prefix))
        new_fitK.writelines(line.replace('$NSTRUCT',  nframes))
    temp_fitK.close()
    for line in temp_fitP.readline():
        new_fitP.writelines(line.replace('$PREFIX', prefix))
        new_fitP.writelines(line.replace('$NSTRUCT',  nframes))
        new_fitP.writelines(line.replace('$PARAMFILE', prefix + '-P-fit.in' ))
    temp_fitP.close()


def pepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type) :
    shutil.copy(os.path.realpath(sys.argv[0]))+'param.in' ,'parameter_'+prefix +'.in' )

    temp_fit=open('parameter_'+prefix +'.in' , 'a'):
    index=torsion_names.index(prefix)
        for types in all_dihedrals_type[index]:
            for mult in range(type[4])
                temp_fit.writelines('	%s %s %s %s	%s	1	0	0' %(type[0], type[1], type[2],type[3], mult  ) )
    temp_fit.close()



if __name__==__main__:

    all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name ,torsion_names , dihedrals_heavy_index= find_dihedrals(inputmol2)
    for p in torsion_names:
        prefix = p + '_'
        if os.isfile(prefix'*') :
            extract_stationary_structures(prefix)
            nframes= make_traj(prefix)
            pepare_paramfit_job_files (prefix, nframes)
            pepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type)

    input='PSI_-70.mol2'
