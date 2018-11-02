#! /home/marie/Utilities/sire.app/bin/python
import mdtraj as md
import sklearn.metrics as sk1
import parmed
#import 02_prepare_paramfit_inputs
from  get_all_dihedral import find_dihedrals
from  energyDecomposition import singlepoint
import sys
import os
import subprocess

def get_converted_QM_data(prefix):
    if not os.isfile(prefix+'energySP.dat')  :sys.exit(prefix+' energy file not found'  )
    qm = open(prefix+'energySP.dat','r')
    QM_data=[ float(i)*627.509 for i in qm.readline() ]
    minimum = min(QM_data)
    converted_QM= [val - minimum for val in QM_data]
    return converted_QM

def single_point_calculations(topol, prefix):
    singlepointenergy=singlepoint(topol, prefix+'-traj.mdcrd')
    return singlepointenergy

	#return K
def paramfit_fit_K(prefix,topol):
    cmd='paramfit -i %s-K-fit.in -p  %s   -c %s-traj.mdcrd  -q %s-energySP.dat   > paramrun.out' %(prefix, topol ,prefix,prefix)
    print(cmd)
    os.system(cmd)
    file=open('paramrun.out','r')
    for line in file.readlines():
        if '*K' in line :
            print (line , line.split())
            K=float(line.split()[2] )
    file.close()
    return K

def run_paramfit(topolin, prefix,K,mode,optimizations, max_gen, gen_to_conv , gen_simplex, gen_without_simplex, mutation_rate , parent_percent, search_space ):
    # K value of K
    # mode = BOTH , GENETIC, SIMPLEX
    # mutation_rate

    old_fitP= open(prefix + '-P-fit.in','r')
    new_fitP= open(prefix + '-P-fit-run.in','w')
    for line in old_fitP.readlines():
        newline=line.replace('$K', str(K)).replace('$ALGORITHM', str(mode))
        newline=newline.replace('$SSPACE', str(search_space)).replace('$PARENT', str(parent_percent) ).replace('$MUT_RATE',str(mutation_rate))
        newline=newline.replace('$GEN_WT_SIMPLEX', str(gen_without_simplex)).replace('$GEN_SIMPLEX',str(gen_simplex ))
        newline=newline.replace('$MAX_GEN',str(max_gen)).replace('$GEN_CONV',str(gen_to_conv)).replace('$OPT', str(optimizations))
        new_fitP.writelines(newline )
        print(newline)
    cmd='paramfit -i %s-P-fit-run.in -p %s -c %s-traj.mdcrd -q %s-energySP.dat ' %(prefix, topol ,prefix,prefix)
    #cmd='paramfit -i %s-P-fit-run.in -p %s -c %s-traj.mdcrd -q %s-energySP.dat' %(prefix, topol, prefix ,prefix)
    print (cmd)
    #os.system(cmd)
    subprocess.check_call(cmd, shell=True)
    file=open(prefix+'-paramfit.frcmod','r')
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

if __name__=='__main__' :
    inputfile=sys.argv[1]
    all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index=find_dihedrals(inputfile)
    topol='input.prmtop'
    for prefix in ['CHI1_']:
        K=paramfit_fit_K(prefix, topol)
        besttopol='input.prmtop'
        for prefix in torsion_names:
            prefix += '_'
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

                run_paramfit(besttopol, prefix,K,mode,optimizations, max_gen, gen_to_conv , gen_simplex, gen_without_simplex, mutation_rate , parent_percent, search_space )
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
            if sol==1:
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
