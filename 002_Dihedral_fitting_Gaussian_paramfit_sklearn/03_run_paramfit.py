import MDtraj as md
import sklearn.metrics as sk1
import parmed
import 02_prepare_paramfit_inputs
import energy-decomposition
from  enenergy-decomposition import singlepoint

def get_converted_QM_data(prefix):
    if not os.isfile(prefix+'energy')  :sys.exit(prefix+'energy   file not found'  )
    qm = open(prefix+'-energy','r')
	QM_data=[ float(i)*627.509 for i in qm.readline() ]
    minimum = min(QM_data)
    converted_QM= [val - minimum for val in QM_data]
    return converted_QM

def single_point_calculations(topol, prefix):

	singlepointenergy=singlepoint(topol, prefix+'-traj.mdcrd')
    return singlepointenergy

	#return K
def paramfit_fit_K(prefix):
    os.system('paramfit -i %s-K-fit.in -p  %s -c %s  -q %s-energy  -c %s-traj.mdcrd' >paramrun.out' %(prefix, topol,prefix ,prefix))
    file=open(paramrun.out,'r')
    for line in file.realine():
        if '*K' in line :
            print (line , line.split())
            K=float(line.split()[2] )
            break
    file.close()
    return K

def run_paramfit(topolin, prefix,K,mode,optimizations, max_gen, gen_to_conv , gen_simplex, gen_without_simplex, mutation_rate , parent_percent, search_space ):
    # K value of K
    # mode = BOTH , GENETIC, SIMPLEX
    # mutation_rate

    old_fitP= open(prefix + '-P-fit.in','r')
    new_fitP= open(prefix + '-P-fit-run.in','w')
    for line in temp_fitP.readline():
        new_fitP.writelines(line.replace('$K', K))
        new_fitP.writelines(line.replace('$ALGORITHM', mode))
        new_fitP.writelines(line.replace('$OPT',optimizations ))
        new_fitP.writelines(line.replace('$GEN_CONV',gen_to_conv ))
        new_fitP.writelines(line.replace('$MAX_GEN',max_gen ))
        new_fitP.writelines(line.replace('$GEN_SIMPLEX',gen_simplex ))
        new_fitP.writelines(line.replace('$GEN_WT_SIMPLEX', gen_without_simplex))
        new_fitP.writelines(line.replace('$MUT_RATE',mutation_rate))
        new_fitP.writelines(line.replace('$PARENT',parent_percent ))
        new_fitP.writelines(line.replace('$SSPACE', search_space))
        os.system('paramfit -i %s-P-fit.in -p  %s -c %s  -q %s-energy  -c %s-traj.mdcrd' >paramrun.out' %(prefix, topol,prefix ,prefix))
        file=open(prefix+-paramfit.frcmod,'r')
        top = parmed.load_file(topol)
        for line in file.readline():
            s=line.split('-')
            if len( s) ==4:
                a=np.concatenate([s[0].remplace(' ','')],[s[1].remplace(' ','')], [s[2].remplace(' ','')],[s[3].split()]  )
                phi=a[4]
                phase= a[6]
                if phi < 0 :
                    phi=-phi
                    phase= phase + 180 %180    #phase == 0 or 180 ::::   (phi*cos(angle +180)= -phi*cos(angle))
                a=s[0]
                if [s[0],s[1],s[2],s[3]] != previousline :  top.deleteDihedral(s[0],s[1],s[2],s[3])  #delete all multiplicity for one dihedral
                top.addDihedral(s[0],s[1],s[2],s[3],  phi , a[5] , phase )

                previousline=[s[0],s[1],s[2],s[3]]

    top.writeprmtop('newtopol.prmtop')
    top.tools.writefrcmod('new.frcmod')

if __name__==__main__ :

    mode='BOTH'
    optimizations=30
    max_gen=50
    gen_to_conv=50
    gen_simplex=20
    gen_without_simplex=0
    mutation_rate=0.1
    parent_percent=0.99
    search_space=0.1000000

    solutions=[1,2,3,4,5]
    sol=[0]
    while len(solutions) != 0 :
        QM_data=get_converted_QM_data(prefix)
        singlepointenergy=single_point_calculations(topol, prefix)

    	singlepointenergyparamfit=single_point_calculations('newtopol.prmtop', prefix)

    	score_before= sk1.r2_score(QM_data, singlepointenergy)

    	score_after =  sk1.r2_score(QM_data, singlepointenergyparamfit )
        if score_after >  score_before :
            shutil.copy('newtopol.prmtop' ,'best.prmtop')
            shutil.copy('new.frcmod' ,'new.frcmod')

        else :
            sol = solutions.pop[0]  #got to the next step


        if sol==0:
            gen_simplex+=10
        if sol==1
            parent_percent+=-0.3
            gen_simplex=20
        if sol==2:
            gen_simplex+=10
        if sol==3:
            gen_without_simplex+=10
            gen_simplex=20
            mutation_rate+=0.1
        if sol==4:
            gen_without_simplex+=10
            search_space+=0.1
        if sol==5:
            optimizations+=20
