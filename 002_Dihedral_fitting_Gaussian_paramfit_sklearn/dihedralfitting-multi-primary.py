#! /usr/bin/env python
d=' \
############################ \
Prepare gaussian input from a mol2 file. \
This is an ugly but efficient script to generate conformations spaced from 10 degrees using tleap and \
manually parse the generated mol2 files to gaussian input files.  \n \n \
It is not necessarry but advised to rename your atoms in your input : CA CB C N ... as the script will try to assign PHI PSI and CHI\
atomtypes will be assigned but won\'t be used in calculations \n\n \
Known bugs : -parmchk2 doesn\'t know CX atom type (CA) so add it to to parmcheck file or change the CX to CT\
             -The molecule to be parametrized can only be called "LIG" \
\n\n                                          !!!!!!WARNING!!!!!!!!  \
                                          This script is not robust ! \
                    PLEASE PREVIEW your mol2 files before running long gaussians calculations...  \
                                        !!!!!!!!!!!!!!!!!!!!\
############################\
Author : Marie Bluntzer \
############################\
'
import os
import subprocess
import argparse
from  get_all_dihedral import find_dihedrals
from multiprocessing import Pool
import mdtraj as md
import os.path
from  process_gaussians_files  import extract_stationary_structures
from  get_all_dihedral import find_dihedrals_mod
import mdtraj as md
import shutil
import sys
import numpy as np
import numpy as np
from lmfit import Minimizer, Parameters, report_fit
import math
import mdtraj
import os
import types
import numpy
#import MDAnalysis as MDA
import parmed
import os
import sys
import itertools
import matplotlib.pyplot as plt
import sklearn.metrics
#def get_close_atoms(target_atom):

	#close_atoms = ns.search(target_atom.coord, 2)
	#return close_atoms
sire_folder='/media/marie/0f307c0d-a7bd-40fc-a5cd-5c79fde466d8/marie/Utilities/sire.app/bin/python'
lib_folder=os.path.dirname(sys.argv[0])
inputfile='newcharges.mol2'


def guess_order(frcmod , dihe_types):
    f= open(frcmod,'r')
    lines = f.readlines()
    order=[]
    amp=[]
    for typ in dihe_types :
        lastline=0
        order.append([])
        amp.append([])
        typstring='-'.join(typ)
        typ.reverse()
        typstring2='-'.join(typ)


        for i,l in enumerate(lines) :
            if typstring in l.replace(' ', '') or typstring2 in l.replace(' ', '')  :
                if i==lastline+1 :
                    lastline=i
                    order[-1].append(abs(int(float(l[17:].split()[2]))))
                    amp[-1].append(float(l[17:].split()[0]))
                else :
                    lastline=i
                    order[-1]=[abs(int(float(l[17:].split()[2])))]
                    amp[-1]=[float(l[17:].split()[0])]
    print( dihe_types)
    print(order)
    print(amp)
    return order, amp


def makedihedralzero(frcmod ,frcmodzero, dihe_types):
    f= open(frcmod,'r')
    f= open(frcmodzero,'w')
    lines = f.readlines()

    for i,l in enumerate(lines) :
        while 'IMPROPER' not in line or flag==True :
            flag=True
            for typ in dihe_types :
                typstring='-'.join(typ)
                typ.reverse()
                typstring2='-'.join(typ)
                if typstring in l.replace(' ', '') or typstring2 in l.replace(' ', '')  :
                    par= l[17:].split()
                    par[0]='0.000'
                    l=l[:17]+ ' ' +'   '.join(par)
                    print(l)
            f.writelines(l)

def check_clashes(torsion_name, conf ,cutoffnonh=0.175 ,cutoffoneh=0.140 ,cutoffoxy=0.185 ,cutofftwoh=0.150) :
    bonds=[]
    mol=md.load('%s_%s_%s.mol2' %(torsion_name,conf ,0))
    for bond in mol.topology.bonds :
                bonds.append( [bond.atom1.index,bond.atom2.index] )
    #print('%s_%s_%s.mol2' %(torsion_name,conf ,0))
    #print('cutoff cutoffnonh=%s ,cutoffoneh=%s ,cutoffoxy=%s ,cutofftwoh=%s'  %(cutoffnonh ,cutoffoneh,cutoffoxy,cutofftwoh))
    for i in range(-180,180,10) :
            pairs_non_h=[]
            pairs_one_h=[]
            pairs_two_h=[]
            pairs_one_h_one_ox=[]
            bonds=[]
            mol=md.load('%s_%s_%s.mol2' %(torsion_name,conf ,i))
            for bond in mol.topology.bonds :
                bonds.append( [bond.atom1.index,bond.atom2.index] )

            for at1 in mol.topology.atoms :
                     for at2 in mol.topology.atoms :
                        if  at1 !=at2 and [at1.index, at2.index] not in bonds and  [at2.index, at1.index] not in bonds:
                            if at1.element.atomic_number==1 and at2.element.atomic_number==1  :
                                 pairs_two_h.append([at1.index, at2.index])
                            elif ( at1.element.atomic_number==1 or at2.element.atomic_number==1    ) and  ( at1.element.atomic_number==8 or at2.element.atomic_number==8    ):
                                pairs_one_h_one_ox.append([at1.index, at2.index])
                            elif ( at1.element.atomic_number==1 or at2.element.atomic_number==1    ) :
                                pairs_one_h.append([at1.index, at2.index])
                            else :
                                pairs_non_h.append([at1.index, at2.index])

            if min(md.compute_distances(mol,  pairs_non_h)[0])< cutoffnonh :
                print('non hydrogen')
                return   False
            if min(md.compute_distances(mol,  pairs_one_h)[0])< cutoffoneh :
                print('one hydrogen')
                return   False
            if min(md.compute_distances(mol,  pairs_two_h)[0])< cutofftwoh :
                print('two hydrogen')
                return   False
            if min(md.compute_distances(mol, pairs_one_h_one_ox)[0])< cutoffoxy:
                print('oxygen' ,min(md.compute_distances(mol, pairs_one_h_one_ox)[0]) )
                return   False
            pairs_non_h.clear()
            pairs_one_h.clear()
            pairs_two_h.clear()

    return   True


def prepfiles(list_of_atm,  dihedrals_heavy_index, torsion_names ,theo ):
    global conf_per_torsion
    list_of_key={}

    inputfile=open(inputmol2_list[0], 'r')
    inputfiliecontent=inputfile.readlines()
    for l in range(len(inputfiliecontent)):
        if res_name in inputfiliecontent[l] :
            list_of_key.update({inputfiliecontent[l].split()[1]: inputfiliecontent[l].split()[0]})

    inputfile.close()
    conf_per_torsion=[]
    for bond in range(len(list_of_atm)) :
        if torsion_names[bond] != False :
            j=0

            for conf in range( n_conf ) :
                clashes_check=False
                while clashes_check==False   and  j<  len(inputmol2_list) :
                    c=conf
                    tleap_input_file=open('tleap.in', 'w')
                    tleap_input_file.writelines(['source leaprc.protein.ff14SB \n' ])
                    for angle in range(-180,+180,10) :

                        inputmol2=inputmol2_list[j]

                        tleap_input_file.writelines( [  'm = loadmol2 %s  \n' %(inputmol2) , 'relax m \n' ])
                        tleap_input_file.writelines( [ 'impose m  { 1 } {  { %s %s } }\n' %(list_of_atm[bond][0] , angle)])
                        tleap_input_file.writelines( [ 'savemol2 m %s_%s_%s.mol2 1\n' %(torsion_names[bond] , conf , angle)] )

                    tleap_input_file.writelines(['quit' ])
                    tleap_input_file.close()

                    subprocess.check_output('tleap -f tleap.in' ,shell=True  )
                    clashes_check=check_clashes(torsion_names[bond], conf )
                    if clashes_check==False :

                        os.system('rm %s_%s*mol2 ' %(torsion_names[bond] , conf ))
                        c+=-1
                    #print(bond , conf ,  clashes_check, j )
                    else : print(inputmol2)
                    j+=1
            if  j==  len(inputmol2_list)   : print('I did not manage to produce %s set of conformations for torsion %s : only %s ' %  (n_conf ,torsion_names[bond], c+1))
            else : print('I managed to produce %s set of conformations for torsion %s ' %  (n_conf ,torsion_names[bond] ))
            conf_per_torsion.append(c+1)
        else: conf_per_torsion.append(False)
    #for finp in os.listdir('./') :
    for numb in  range(len(torsion_names)):

    #for numb in  range(1):
        if  torsion_names[numb]!= False:
            for conf in range(conf_per_torsion[numb]) :
                for angle in range(-180,180,10 ):
                    finp='%s_%s_%s.mol2'  %( torsion_names[numb],conf,angle)
                    if finp[-4:]=='mol2' and finp[:3] in ['PHI', 'PSI' , 'CHI'] and os.path.isfile(finp):
                        fi=open(finp, 'r')
                        fo=open('%s' %(finp[:-4]+'gau'), 'w')
                        #fo.write('--Link1-- \n%%RWF=%s\n%%NoSave\n%%chk=%s \n#P %s/6-31G* Opt=(ModRedundant,MaxCycle=50) Freq SCF=(Conver=6,MaxCycle=62) Pop=NoMBS\n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],finp[:-5],theo, charge,multiplicity))
                        #fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G* Geom=(ModRedundant) Opt=(ModRedundant,MaxCycle=50) SCF\n  \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],theo, charge,multiplicity))
                        fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G*  Geom=(ModRedundant)  \n\nparamfit run\n  \n%s    %s \n'  %(finp[:-5],theo, charge,multiplicity))
                        #fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G* \n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],theo, charge,multiplicity))
                        while 'ATOM' not in fi.readline()  : continue
                        i=True
                        while i==True :
                                line=fi.readline().split()
                                if 'TRIPOS' not  in line[0] :
                                    #with charges
                                    lineout=line[1][0]+'-' + line[5] + '-' + line[8]  + '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                                    #without charges
                                    #lineout=line[1][0]+ '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                                    fo.write(lineout)
                                else:
                                    i=False
                                    fo.write('\n')
                        print(finp[0:4].replace('_',''))
                        #if torsion_names.index(finp[0:4].replace('_','')) not in ['CHI'] : numb=torsion_names.index(finp[0:3].replace('_',''))
                        #else : numb=torsion_names.index(finp[0:4].replace('_',''))

                        #listdih=list_of_atoms[numb]

                        listdih= dihedrals_heavy_index[numb]
                        #print(str (torsion_names.index(finp[0:4].replace('_',''))) + ' ' +str(listdih[0])+ ' ' +str(listdih[1])+ ' '+ str(listdih[2])+ ' '+str(listdih[3]))
                        #fo.write('F'+str(list_of_atoms[numb][0])+'\nF'+str(list_of_atoms[numb][1])+'\nF'+str(list_of_atoms[numb][2])+'\nF'+str(list_of_atoms[numb][3])+'\n')
                        #fo.write('F ' + str(listdih[0]+1) +' '+ str(listdih[1]+1)+ ' '+ str(listdih[2]+1)+ ' '+str(listdih[3]+1)+ ' F\n')
                        #print(dihedrals_heavy_index[numb])
                        #fo.write('D ' + str(listdih[0]+1) +' '+ str(listdih[1]+1)+ ' '+ str(listdih[2]+1)+ ' '+str(listdih[3]+1)+ ' F\n')
                        #fo.write('D  %s  %s  %s  %s '   %( str(listdih[0]+1), str(listdih[2]+1), str(listdih[3]+1)+ ' '+ str(listdih[2]+1)+ ' '+'*'+ ' =%s B\n' %angle)
                        fo.write('D ' + '*' +' '+ str(listdih[1]+1)+ ' '+ str(listdih[2]+1)+ ' '+'*'+ ' F\n')
                        fo.close()
                        fi.close

'''
    for bond in range(len(list_of_atm)) :
        for conf in range( len(inputmol2_list)) :
            for angle in range(-180,+180,10) :
                torsion_name=torsion_names[bond]
                struct=inputmol2_list[conf]
                fi=open(struct, 'r')
                fo=open('%s-%s-%s.gau'  %(torsion_name,struct[:-4] , str(angle)), 'w')
                #fo.write('--Link1-- \n%%RWF=%s\n%%NoSave\n%%chk=%s \n#P %s/6-31G* Opt=(ModRedundant,MaxCycle=50) Freq SCF=(Conver=6,MaxCycle=62) Pop=NoMBS\n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],finp[:-5],theo, charge,multiplicity))
                #fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G* Geom=(ModRedundant) Opt=(ModRedundant,MaxCycle=50) SCF\n  \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],theo, charge,multiplicity))
                fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G*  Geom=(ModRedundant)  \n\nparamfit run\n  \n%s    %s \n'  %(struct[:-5],theo, charge,multiplicity))
                #fo.write('--Link1-- \n%%chk=%s\n%%NoSave \n#SP %s/6-31G* \n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],theo, charge,multiplicity))
                while 'ATOM' not in fi.readline()  : continue
                i=True
                while i==True :
                        line=fi.readline().split()
                        if 'TRIPOS' not  in line[0] :
                            #with charges
                            lineout=line[1][0]+'-' + line[5] + '-' + line[8]  + '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                            #without charges
                            #lineout=line[1][0]+ '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                            fo.write(lineout)
                        else:
                            i=False
                            fo.write('\n')

                #listdih=list_of_atoms[numb]

                fo.write('D  %s  %s  %s  %s  =%s B\n'  %( str(listdih[0]+1), str(listdih[1]+1), str(listdih[2]+1), str(listdih[3]+1) , angle))
                fo.write('D ' + '*' +' '+ str(listdih[1]+1)+ ' '+ str(listdih[2]+1)+ ' '+'*'+ ' F\n')
                fo.close()
                fi.close
'''


def run_g09(input):
    #prepare gaussian files run gaussian jobs and read gaussian files


    if os.path.isfile('./%s.log'%input) == False or  'Normal termination' not in  subprocess.getoutput('tail -1 ./%s.log'%input)  or  'HF=' not in  subprocess.getoutput('grep \'\\HF=\' ./%s.log'%input) :
        print('Running g09 %s.gau' %input  )
        os.system('g09 %s.gau' %input  )


    #time.wait(1)

def run_allconf(torsion_names, n_cpu=15) :
    #call run_resp for all conformation using multiprocessing
    # !!!!! For a RE-RUN  DO NOT USE multiprocessing :
    # 10 instances of antechamber running at the same time seems to break it ! use the run_resp in the 'for' loop



    args=[]
    print(n_conf)
    for i,name in enumerate(torsion_names):

        if  name !=False:
            for angle in range(-180,+180,10):
                for conf in range(conf_per_torsion[i]) :
                    args.append('%s_%s_%s' %(name,conf, angle))

    p =Pool(processes=n_cpu)
    r = p.map_async(run_g09,args)
    #print(args)
    r.wait()



def add_order():
    F= open('amber.frcmod', 'r')
    N=open('oneorder.frcmod','w')
    line=F.readline()
    while line[0:4] !='DIHE':
        N.writelines(line)
        line=F.readline()

    N.writelines(line)
    previous=F.readline()
    line=F.readline()
    next=F.readline()

    while line[:6] not in ['IMPROP','NONBON' ]:
        if line[0:11]!=previous[0:11]:
                N.writelines( previous[0:11]+ '   1    0.000       180.000          -6.000     \n')
        N.writelines(previous)
        previous=line
        line=next
        next=F.readline()


    N.writelines(previous)
    N.writelines(line)
    N.writelines(F.readlines())
    N.close()


'''
def make_traj(prefix):
    list_pdb=open(prefix +'-listfiles')
    lines=list_pdb.readlines()
    if len(lines) > 0 : first=lines.pop(0).replace(' \n','')
    else :
        #sys.exit('no file in listfile file ' + prefix +   str(lines) )
        return False
    t = md.load(first, top='input.prmtop')

    for i in range(len(lines)):
        filename=lines.pop(0).replace(' \n','')
        print(filename)
        tnext= md.load(filename.replace(' \n',''), top='input.prmtop')
        t= md.join([t , tnext], check_topology=True, discard_overlapping_frames=False)
    t.save_mdcrd(prefix+ '-traj.mdcrd' )
    t.save_netcdf(prefix+ '-traj.nc' )
    return t.n_frames
'''
def make_traj(prefix, index ):
    list_pdb=open(prefix +'-listfiles')
    lines=list_pdb.readlines()
    if len(lines) > 0 : first=lines.pop(0).replace(' \n','')
    else :
        return False
    t = mdtraj.load(first, top='input.prmtop')

    table, bonds = t.topology.to_dataframe()

    for i in range(len(lines)):
        filename=lines.pop(0).replace(' \n','')
        #print(filename)
        tnext= mdtraj.load(filename.replace(' \n',''), top='input.prmtop')
        t= mdtraj.join([t , tnext], check_topology=True, discard_overlapping_frames=False)
    t.save_mdcrd(prefix+ '-traj.mdcrd' )
    t = mdtraj.load(prefix+ '-traj.mdcrd', top='input.prmtop')
    table, bonds = t.topology.to_dataframe()
    print(index)
    dih= mdtraj.compute_dihedrals(t,index)
    t.save_netcdf(prefix+ '-traj.nc' )
    return t.n_frames ,dih


def prepare_paramfit_job_files (prefix, nframes) :
    ## This fonction prepares
    # - fit_K.in for paramfit : find constant to fit Kcal to Hartree
    # fit_P.in for torsion

    temp_fitK= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'/fit_K_A.in' , 'r')
    temp_fitP= open(os.path.dirname(os.path.realpath(sys.argv[0]))+'/job_fit_1.in' , 'r')
    new_fitK= open(prefix + '-K-fit.in','w')
    new_fitP= open(prefix + '-P-fit.in','w')

    for line in temp_fitK.readlines():
        new_fitK.writelines(line.replace('$PREFIX', prefix).replace('$NSTRUCT',  str(nframes)))
    temp_fitK.close()
    for line in temp_fitP.readlines():
        new_fitP.writelines(line.replace('$PREFIX', prefix).replace('$NSTRUCT',  str(nframes)).replace('$PARAMFILE', prefix + '-P-fit.in' ))
    temp_fitP.close()


def prepare_paramfit_param_files (prefix,torsion_names, all_dihedrals_type) :
    template= os.path.dirname(os.path.realpath(sys.argv[0])) +'/param.in'
    jobfile='parameter_'+prefix +'.in'
    shutil.copy(template ,jobfile )
    temp_fit=open(jobfile, 'a')
    index=torsion_names.index(prefix.split('_')[0])
    print(prefix , all_dihedrals_type)
    added = []
    count=[]
    for types in all_dihedrals_type[index]:
        #print(types)
        #    for mult in range(types[4]):
        if types[0] + types[1] +types[2]+types[3] not in added :
        # order != mult   if order 1 3 6 / mult = 3
            count.append(0)
            added.append(types[0] + types[1] +types[2]+types[3])
            mult=0
        else :
            count[added.index(types[0] + types[1] +types[2]+types[3] )]+=1
            mult=count[added.index(types[0] + types[1] +types[2]+types[3] )]
        temp_fit.writelines('	%s %s %s %s	%s	1	0	0\n' %(types[0], types[1], types[2],types[3] ,mult) )
    temp_fit.close()

def prepare_mdgx_job_files (prefix,torsion_names, all_dihedrals_type) :
    ## This fonction prepares
    # - fit_K.in for paramfit : find constant to fit Kcal to Hartree
    # fit_P.in for torsion

    new_fit= open(prefix + '-mdgx-fit.in','w')

    new_fit.writelines('&files\n')
    new_fit.writelines(' -parm   /home/marie/Utilities/amber16/dat/leap/parm/frcmod.ff14SB\n')
    new_fit.writelines(' -o      fit.out\n')
    new_fit.writelines(' -fmod amber.frcmod\n')
    new_fit.writelines('&end\n\n')

    new_fit.writelines('&param\n')
    new_fit.writelines(' System  input.prmtop %s-traj.nc %s-energySP.dat\n' %(prefix, prefix)   )
    new_fit.writelines(' ParmOutput frcmod\n')
    new_fit.writelines(' eunits    hartree,\n')
    new_fit.writelines(' accrep    report.m\n')
    new_fit.writelines(' verbose    1,\n\n')
    new_fit.writelines(' % Angle fitting input\n')
    angle_added=[]
    dih_added=[]
    index=torsion_names.index(prefix.split('_')[0])
    for types in all_dihedrals_type[index]:

        if sorted([types[0] , types[1] ,types[2] ])  not in angle_added :
         angle_added.append(sorted([types[0], types[1] ,types[2]]))
         new_fit.writelines('  fita       %s %s %s\n'%(types[0], types[1] ,types[2]))
        if sorted([types[1] ,types[2], types[3] ])  not in angle_added :
         angle_added.append(sorted([types[1] ,types[2],types[3]]))
         new_fit.writelines('  fita       %s %s %s\n'%(types[0], types[1] ,types[2]))

    new_fit.writelines('  FitAnglEq  1,\n  arst       0.0002,\n  arstcpl    1.0,\n\n')
    new_fit.writelines(' % Torsion fitting input\n')

    for types in all_dihedrals_type[index]:

        if [types[0] , types[1] ,types[2],types[4] ]  not in dih_added :
         dih_added.append([types[0], types[1] ,types[2],types[3] ])
         new_fit.writelines('  fith       %s %s %s %s\n'%(types[0], types[1] ,types[2],types[3]))
    new_fit.writelines('  hrst       0.0002,\n&end')
    new_fit.close()
    print(angle_added)


def comp_with_derivdata3( outputparams):

        #print(angle )
        yval=np.zeros(len(outputparams[0][1][0])*len(outputparams[0][1]))
        print(yval)
        for j, iparm in enumerate(outputparams) :
                xval=[]

                for u in range(len(iparm[1])):
                    xval.extend(iparm[1][u])
                order=iparm[3]
                amp = iparm[4]

                yval  +=  np.multiply(np.sin(np.multiply(xval,order  )) , -amp * order ) # 0.5
        print( outputparams[0][1])


        return   np.array(outputparams[0][1]).flatten(), yval


# In[77]:



def comp_with_data3( outputparams):

        yval=np.zeros(len(outputparams[0][1][0])*len(outputparams[0][1]))
        for j, iparm in enumerate(outputparams) :
            xval=[]

            for u in range(len(iparm[1])):
                    xval.extend(iparm[1][u])
            order=iparm[3]
            if iparm[4] >0 :
                amp = iparm[4]
                phase=0
            else :
                amp = -iparm[4]
                phase=  np.pi
            yval  +=  np.multiply(1 + np.cos(np.multiply(xval,order  ) +phase) ,  amp ) #0.5

        return np.array(outputparams[0][1]).flatten(), yval


def fctderiv2min_multi(params, inputfct,data,weights=False):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1][0]))
    sum_amp=0
    residuals=[]

    for u in range(len(data)):
        for inputdih in inputfct:
            ampparam=inputdih[2]
            x= inputdih[1][u]


            order= inputdih[3]
            amp = params[ ampparam].value

            sum_amp+=amp

            model +=   np.multiply(np.sin(np.multiply(x,order )  ) , -  amp * order)
        model += params['K' +str(u)]

        if weights==False :
            residuals.extend(model - data[u])
        else :
            residuals.extend(np.divide(model - data, weights))
    return residuals


def fctderiv2min_multi_prim(params, inputfct,data,weights=False):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1][0]))
    sum_amp=0
    residuals=[]

    for u in range(len(data)):
        for inputdih in inputfct:
            ampparam=inputdih[2]
            x= inputdih[1][u]


            order= inputdih[3]
            amp = params[ ampparam].value

            sum_amp+=amp
            order= inputdih[3]
            amp_abs = abs(params[ ampparam].value)
            if params[ ampparam].value < 0 : phase = math.pi
            else : phase = 0
            model +=   np.multiply(   1  +  np.cos(np.multiply(x,order ) + phase  ) ,   amp_abs)
        model += params['K' +str(u)]

        if weights==False :
            residuals.extend(model - data[u])
        else :
            residuals.extend(np.divide(model - data, weights))
    return residuals


    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion

def fctderiv2min(params, inputfct,data,weights=False):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1]))
    sum_amp=0

    for inputdih in inputfct:
        ampparam=inputdih[2]
        x= inputdih[1]
        order= inputdih[3]
        amp = params[ ampparam].value

        sum_amp+=amp
        model +=   np.multiply(np.sin(np.multiply(x,order )  ) , -  amp * order)
    model += params['K']
    if weights==False :
        return model - data
    else :
        return np.divide(model - data, weights)

def fct2min(params, inputfct,data):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1]))
    sum_amp=0

    for inputdih in inputfct:
        ampparam=inputdih[2]
        x= inputdih[1]
        order= inputdih[3]
        amp = params[ ampparam].value

        sum_amp+=amp
        model +=   np.multiply(np.cos(np.multiply(x,order )  ) ,  amp )

    model += params['K']
    return model - data




def fctderiv2minimp(params, inputfct,data):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solution1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1]))
    sum_amp=0
    for inputdih in inputfct:

        ampparam=inputdih[2]
        x= inputdih[1]
        order= inputdih[3]
        amp = params[ ampparam].value

        sum_amp+=amp
        model +=   np.multiply(np.sin(np.multiply(x,order )  ) , -  amp * order)

    model += params['K'].value
    return model - data


# In[131]:


def fctderiv2min_new(params, inputfct,data):
    # define objective function: returns the array to be minimized
    # To do: define the angle better for the non heavy atoms,
    # it works now for most dihedral but not for psi as angle are not ~2pi/3
    # Solut   ion1 : mdtraj.computedihedrals  : in rad :-)
    # Solution2 : read angle harmonic equilibration values from frcmod  : in degrees --> conversion
    model = np.zeros(len(inputfct[0][1]))
    sum_amp=0
    for inputdih in inputfct:

        ampparam=inputdih[2]
        x= inputdih[1]
        order= inputdih[3]
        print ( ampparam)
        amp = params[ ampparam].value

        #amp=2
        #print(inputdih,amp ,order)
        sum_amp+=amp
        model +=   np.multiply(np.sin(np.multiply(x,order)  ) ,  - 0.5* amp * order)

    model += params['K'].value
    return model - data



def make_missing_parms_parmck2(inputfile, frcmod='gaff2.frcmod'):
    print('parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o gaff2.frcmod -a \'Y\' ')
    os.system( 'parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o gaff2.frcmod -a \'Y\' ' %(inputfile))
    tleapinput=open('temp.in','w')
    tleapinput.writelines([  'source leaprc.protein.ff14SB \n'  , 'loadAmberParams %s \n'%(frcmod), 'input = loadmol2 %s \n '  %(inputfile ) ,    'saveAmberParm  input  input.prmtop input.inpcrd\n ', 'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp.in ' )






def write_frcmod(result, torsion_name , all_dihedrals_type,all_dihedrals_type_index):
    frcmod_file=open('%s.frcmod' %(torsion_name),'w')
    frcmod_file.writelines('\nDIHE\n')
    new_frcmod=[]
    threshold=0.005
    params=[[],[]]
    for parm  in result.params:

        if  parm[0:3]== 'amp' and result.params[parm].value != 0 and parm[-3:] != 'bis' : # dependency.keys():
                #print( parm[3] , parm[-1]  ,  result.params['order'+parm[3:]].value, result.params[parm].value)
                #if math.fabs(result.params[parm].value) > threshold :
                #new_frcmod.append([all_dihedrals_type[torsion][all_dihedrals_type_index[torsion][int(parm[-1])][-1]],  result.params[parm].value ,  result.params['order'+parm[3:]].value])
                #print( [i.ljust(2) for i in parm.split('_')[-4:]])
                #new_frcmod.append(['-'.join([i.ljust(2) for i in parm.split('_')[-4:]]),  result.params[parm].value ,  result.params['order'+parm[3:]].value])
                order=result.params[parm].name[3]
                phase='  0'
                if result.params[parm].value < 0 : phase='180'
                new_frcmod.append(['-'.join([i.ljust(2) for i in parm.split('_')[-4:]]), math.fabs( result.params[parm].value ) ,  phase ,order ] )
                params[0].append(parm)
                params[1].append(result.params[parm].value)


    new_frcmod.sort(key=lambda x: x[0])

    for l in  range(len(new_frcmod)) :

        line=new_frcmod[l]
        if l != len(new_frcmod)-1 : nextdih=new_frcmod[l+1][0]
        else : nextdih=False
        if line[0] == nextdih  : s= -1
        else : s= 1
        frcmod_file.writelines('%s  1  %s   %s   %s \n' %(line[0] , str(math.fabs(line[1]))[0:6], line[2], str(s*int(line[3] )).ljust(2) ))
    print(params)
    return params




# In[135]:


def plot(data, final,x ,x2=[],outputfile='plot.png'):
    if  len(x2)==0: x2=x

    for i, angle in enumerate(x2) :
        if angle==3.14  : x2[j]=-3.14
    for i, angle in enumerate(x) :
        if angle==3.14  : x[j]=-3.14

    if x2[0]>0 :  x2[0]=-x2[0]
    '''
    if x2[-1]<0 :  x2[0]=-x2[0]
    '''

    try:
        plt.figure()
        plt.plot(x, data, 'k+')
        plt.plot(x2, final, 'r')
        #plt.plot( data2, 'b')
        plt.savefig(outputfile)
        plt.close()
    except ImportError:
        pass


# In[143]:


import lmfit

def leastsquare_with_prem(x,qm,dih_idx, dih_type, maxorder=1 , prem_data=[[],[]]):
    #n_dih=len(dih_idx)
    dih_type_index =  [ j[-1]  for j in dih_idx ]
    #print(len(dih_type_index))
    if len(dih_type_index) == 9 :
        orders = [1,2,3,4]
    if len(dih_type_index) != 9 :
        orders = [1,2,4,5,3]
        maxorder+=1
    orders, amps = guess_order('gaff2.frcmod',dih_type)


    while maxorder*len(dih_type) > 35:
        maxorder-=1
        print('warning maxorder has been reduced to : %s ' %maxorder)
    # create a set of Parameters
    params = Parameters()
    names=[]
    inputfct=[]
    ivalue=10
    print(len(dih_idx))
    for j in range (len(dih_idx)):

            typedih= dih_type[dih_type_index[j]]
            order_list=orders[dih_type_index[j]]
            amp_list=amps[dih_type_index[j]]
            name= '_'.join(typedih)

            name_inv='_'.join(reversed(typedih))

            maxval=1


            if  len(dih_type_index) < 9 :
                maxval=3
            elif typedih[0][0]=='H' and  typedih[3][0]=='H':
                maxval=0.15
            elif typedih[0][0]=='H' or  typedih[3][0]=='H':
                maxval=0.50
            minval=-maxval
            if name not in names and name_inv  not in names :
                names.append(name)

                #for i in range(1, maxorder+1):
                print(order_list)
                for k,i in enumerate(order_list):

                    if 'amp%s_%s'%(i,name) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name))]
                        params.add('amp%s_%s'%(i,name), value=ivalue ,max=maxval , min=minval)
                        print(j, x,'amp%s_%s'%(i,name), i )
                        inputfct.append([j, x[j],'amp%s_%s'%(i,name), i ])

                    elif 'amp%s_%s'%(i,name_inv) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name_inv))]
                        params.add('amp%s_%s'%(i,name_inv), value=ivalue ,max=maxval , min=minval)
                        #print ( 'amp%s_%s'%(i,name_inv))
                        inputfct.append([j, x[j],'amp%s_%s'%(i,name_inv), i ])

                    else :
                        ivalue=0
                        params.add('amp%s_%s'%(i,name), value=amp_list[k] ,max=maxval , min=minval)
                        inputfct.append([j, x[j],'amp%s_%s'%(i,name), i ])
                        names.append(name)
                        #print('value not computed last round')
                        #print('amp%s_%s'%(i,name))
                    #if name == 'CT_CJ_C_O'  and  i == 3 : sys.exit('BUG check me !')
            elif  name in names  :
                print(name)
                for i in order_list:
                    inputfct.append([j, x[j],'amp%s_%s'%(i,name), i ])
                    #print(inputfct[-1])
            elif  name_inv  in names  :
                for i in order_list:
                    inputfct.append([j, x[j],'amp%s_%s'%(i,name_inv), i ])


                    #print(inputfct[-1])
            else : sys.exit('BUG')

            '''
            impropers_propers= [[  imp[0][1], imp[0][0],  imp[0][2] ,  imp[0][3]  ]   for imp in list_improper ]
            impropers_propers.extend( [[  imp[0][0], imp[0][1],  imp[0][3] ,  imp[0][2]  ]   for imp in list_improper ] )
            if typedih in impropers_propers :
                print('IMPRPER FOUND')
                atoms_imp=[imp[0] for imp in list_improper ][impropers_propers.index(['CJ', 'C ', 'N ', 'H'])]
                name='_'.join(atoms_imp)
                inputfct.append([j, x[j],'ampIMP_%s'%(name), i ])
                params.add('ampIMP_%s'%(name), value=0 )
            '''
    params.add('K', value=0 , max=10000 , min=-10000)

    print(params)
    # do fit, here with leastsq model
    #print(names)
    method = 'L-BFGS-B'
    reduce_fcn='neglogcauchy'
    #minner = Minimizer(fctderiv2min_new, params, fcn_args=(inputfct, data))
    #result = minner.minimize(method = 'L-BFGS-B')

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'basinhopping')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'lbfgsb' ,  reduce_fcn='neglogcauchy'  )
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'leastsq'  )
    result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'differential_evolution'  )

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'ampgo')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'nelder')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cg')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'newton')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cobyla')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'tnc')
    #fprime = lambda x: optimize.approx_fprime(x[0], fctderiv2min, 0.01, args=(inputfct, qm))
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'trust-ncg' ,jac=fprime)
    #result  = lmfit.minimize( fctderiv2min,arams, args=(inputfct, qm) , method = 'trust-exact')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'dual_annealing')
    weights=[ max ( 0.95 , min( 1/math.sqrt(i+0.01) , 1.05) ) for i in data   ]
    #weights=[  1 for i in data   ]
    print(weights)
    result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm ,weights) , method = 'cg')
    return result ,inputfct



def leastsquare_with_prem_multi(X,DerivQM,dih_idx, dih_type, maxorder=1 , prem_data=[[],[]]):
    #n_dih=len(dih_idx)

    dih_type_index =  [ j[-1]  for j in dih_idx ]
    #print(len(dih_type_index))
    if len(dih_type_index) == 9 :
        orders = [1,2,3,4]
    if len(dih_type_index) != 9 :
        orders = [1,2,4,5,3]
        maxorder+=1
    orders, amps = guess_order('gaff2.frcmod',dih_type)


    while maxorder*len(dih_type) > 35:
        maxorder-=1
        print('warning maxorder has been reduced to : %s ' %maxorder)
    # create a set of Parameters
    params = Parameters()
    names=[]
    inputfct=[]
    ivalue=10
    print(len(dih_idx))

    for j in range (len(dih_idx)):

            typedih= dih_type[dih_type_index[j]]
            order_list=orders[dih_type_index[j]]
            amp_list=amps[dih_type_index[j]]
            name= '_'.join(typedih)

            name_inv='_'.join(reversed(typedih))
            maxval=1
            if  len(dih_type_index) < 9 :
                maxval=3
            elif typedih[0][0]=='H' and  typedih[3][0]=='H':
                maxval=0.15
            elif typedih[0][0]=='H' or  typedih[3][0]=='H':
                maxval=0.50
            minval=-maxval
            if name not in names and name_inv  not in names :
                names.append(name)

                #for i in range(1, maxorder+1):
                print(order_list)
                for k,i in enumerate(order_list):

                    if 'amp%s_%s'%(i,name) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name))]
                        params.add('amp%s_%s'%(i,name), value=ivalue ,max=maxval , min=minval)
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])

                    elif 'amp%s_%s'%(i,name_inv) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name_inv))]
                        params.add('amp%s_%s'%(i,name_inv), value=ivalue ,max=maxval , min=minval)
                        #print ( 'amp%s_%s'%(i,name_inv))
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name_inv), i ])

                    else :
                        ivalue=0
                        params.add('amp%s_%s'%(i,name), value=amp_list[k] ,max=maxval , min=minval)
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        names.append(name)
                        #print('value not computed last round')
                        #print('amp%s_%s'%(i,name))
                    #if name == 'CT_CJ_C_O'  and  i == 3 : sys.exit('BUG check me !')
            elif  name in names  :
                print(name)
                for i in order_list:
                    inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                    #print(inputfct[-1])
            elif  name_inv  in names  :
                for i in order_list:
                    inputfct.append([j, [ x[j] for x in X ],'amp%s_%s'%(i,name_inv), i ])


                    #print(inputfct[-1])
            else : sys.exit('BUG')



            '''
            impropers_propers= [[  imp[0][1], imp[0][0],  imp[0][2] ,  imp[0][3]  ]   for imp in list_improper ]
            impropers_propers.extend( [[  imp[0][0], imp[0][1],  imp[0][3] ,  imp[0][2]  ]   for imp in list_improper ] )
            if typedih in impropers_propers :
                print('IMPRPER FOUND')
                atoms_imp=[imp[0] for imp in list_improper ][impropers_propers.index(['CJ', 'C ', 'N ', 'H'])]
                name='_'.join(atoms_imp)
                inputfct.append([j, x[j],'ampIMP_%s'%(name), i ])
                params.add('ampIMP_%s'%(name), value=0 )
            '''
    for y in range (len(X )):
        params.add('K' +str(y), value=0 , max=10000 , min=-10000)

    print(params)
    # do fit, here with leastsq model
    #print(names)
    method = 'L-BFGS-B'
    reduce_fcn='neglogcauchy'
    #minner = Minimizer(fctderiv2min_new, params, fcn_args=(inputfct, data))
    #result = minner.minimize(method = 'L-BFGS-B')

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'basinhopping')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'lbfgsb' ,  reduce_fcn='neglogcauchy'  )
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'leastsq'  )
    #result  = lmfit.minimize( fctderiv2min_multi, params, args=(inputfct, DerivQM) , method = 'differential_evolution'  )

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'ampgo')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'nelder')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cg')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'newton')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cobyla')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'tnc')
    fprime = lambda x: optimize.approx_fprime(x[0], fctderiv2min, 0.01, args=(inputfct, qm))
    result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'trust-ncg' ,jac=fprime)
    #result  = lmfit.minimize( fctderiv2min,arams, args=(inputfct, qm) , method = 'trust-exact')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'dual_annealing')
    #weights=[ max ( 0.95 , min( 1/math.sqrt(i+0.01) , 1.05) ) for i in data   ]
    #weights=[  1 for i in data   ]
    #print(weights)
    #result  = lmfit.minimize( fctderiv2min_multi_prim, params, args=(inputfct,  QM) , method = 'cg')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm ,weights) , method = 'cg')
    return result ,inputfct





def leastsquare_with_prem_multi_data(X,QM,dih_idx, dih_type, maxorder=1 , prem_data=[[],[]]):
    #n_dih=len(dih_idx)

    dih_type_index =  [ j[-1]  for j in dih_idx ]
    #print(len(dih_type_index))
    if len(dih_type_index) == 9 :
        orders = [1,2,3,4]
    if len(dih_type_index) != 9 :
        orders = [1,2,4,5,3]
        maxorder+=1
    orders, amps = guess_order('gaff2.frcmod',dih_type)


    while maxorder*len(dih_type) > 35:
        maxorder-=1
        print('warning maxorder has been reduced to : %s ' %maxorder)
    # create a set of Parameters
    params = Parameters()
    names=[]
    inputfct=[]
    ivalue=10
    print(X)

    for j in range (len(dih_idx)):

            typedih= dih_type[dih_type_index[j]]
            order_list=orders[dih_type_index[j]]
            amp_list=amps[dih_type_index[j]]
            name= '_'.join(typedih)

            name_inv='_'.join(reversed(typedih))

            maxval=1

            if  len(dih_type_index) < 9 :
                maxval=3
            elif typedih[0][0]=='H' and  typedih[3][0]=='H':
                maxval=0.15
            elif typedih[0][0]=='H' or  typedih[3][0]=='H':
                maxval=0.50
            minval=-maxval
            if name not in names and name_inv  not in names :
                names.append(name)

                #for i in range(1, maxorder+1):
                print(order_list)
                for k,i in enumerate(order_list):

                    if 'amp%s_%s'%(i,name) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name))]
                        params.add('amp%s_%s'%(i,name), value=ivalue ,max=maxval , min=minval)
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])

                    elif 'amp%s_%s'%(i,name_inv) in prem_data[0] :
                        #print('value computed last round')
                        ivalue = prem_data[1][prem_data[0].index('amp%s_%s'%(i,name_inv))]
                        params.add('amp%s_%s'%(i,name_inv), value=ivalue ,max=maxval , min=minval)
                        #print ( 'amp%s_%s'%(i,name_inv))
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name_inv), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name_inv), i ])

                    else :
                        ivalue=0
                        params.add('amp%s_%s'%(i,name), value=amp_list[k] ,max=maxval , min=minval)
                        print([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
                        names.append(name)
                        #print('value not computed last round')
                        #print('amp%s_%s'%(i,name))
                    #if name == 'CT_CJ_C_O'  and  i == 3 : sys.exit('BUG check me !')
            elif  name in names  :
                print(X[0])
                for i in order_list:
                    inputfct.append([j,[ x[j] for x in X ],'amp%s_%s'%(i,name), i ])
            elif  name_inv  in names  :
                for i in order_list:
                    inputfct.append([j, [ x[j] for x in X ],'amp%s_%s'%(i,name_inv), i ])

                    #print(inputfct[-1])
            else : sys.exit('BUG')



            '''
            impropers_propers= [[  imp[0][1], imp[0][0],  imp[0][2] ,  imp[0][3]  ]   for imp in list_improper ]
            impropers_propers.extend( [[  imp[0][0], imp[0][1],  imp[0][3] ,  imp[0][2]  ]   for imp in list_improper ] )
            if typedih in impropers_propers :
                print('IMPRPER FOUND')
                atoms_imp=[imp[0] for imp in list_improper ][impropers_propers.index(['CJ', 'C ', 'N ', 'H'])]
                name='_'.join(atoms_imp)
                inputfct.append([j, x[j],'ampIMP_%s'%(name), i ])
                params.add('ampIMP_%s'%(name), value=0 )
            '''
    for y in range (len(X )):
        params.add('K' +str(y), value=0 , max=10000 , min=-10000)

    print(params)
    # do fit, here with leastsq model
    #print(names)
    method = 'L-BFGS-B'
    reduce_fcn='neglogcauchy'
    #minner = Minimizer(fctderiv2min_new, params, fcn_args=(inputfct, data))
    #result = minner.minimize(method = 'L-BFGS-B')

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'basinhopping')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'lbfgsb' ,  reduce_fcn='neglogcauchy'  )
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'leastsq'  )
    #result  = lmfit.minimize( fctderiv2min_multi, params, args=(inputfct, DerivQM) , method = 'differential_evolution'  )

    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'ampgo')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'nelder')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cg')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'newton')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'cobyla')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'tnc')
    #fprime = lambda x: optimize.approx_fprime(x[0], fctderiv2min, 0.01, args=(inputfct, qm))
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'trust-ncg' ,jac=fprime)
    #result  = lmfit.minimize( fctderiv2min,arams, args=(inputfct, qm) , method = 'trust-exact')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm) , method = 'dual_annealing')
    #weights=[ max ( 0.95 , min( 1/math.sqrt(i+0.01) , 1.05) ) for i in data   ]
    #weights=[  1 for i in data   ]
    #print(weights)
    result  = lmfit.minimize( fctderiv2min_multi_prim, params, args=(inputfct,  QM) , method = 'cg')
    #result  = lmfit.minimize( fctderiv2min, params, args=(inputfct, qm ,weights) , method = 'cg')
    return result ,inputfct

import parmed

def add_new_parms(inputfile,newfrcmod, MOLNAME ):

    print('parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o gaff2.frcmod -a \'Y\' ' %(inputfile))
    os.system( 'parmchk2 -i  %s  -o frcmod -f mol2 -s 2 -o gaff2.frcmod -a \'Y\' ' %(inputfile))

    tleapinput=open('temp.in','w')
    tleapinput.writelines(['source leaprc.protein.ff14SB \n' , 'source leaprc.gaff2 \n' ,'loadAmberParams gaff2.frcmod \n',  ])

    for name in newfrcmod:
        if name != False :
            tleapinput.writelines( 'loadAmberParams %s.frcmod \n' % name  )
    tleapinput.writelines( [ 'new = loadmol2  %s  \n' %(inputfile),   'saveAmberParm  new  newtopol.prmtop new.rst\n ' ] )
    tleapinput.writelines( [  'quit'])
    os.system( 'tleap -f temp.in >tleap' )
    new=parmed.load_file('newtopol.prmtop')
    parmed.tools.writeFrcmod(new,'new.frcmod').execute()

    tleapinput=open('temp2.in','w')
    tleapinput.writelines(['source leaprc.protein.ff14SB \n' , 'source leaprc.gaff2 \n' ,'loadAmberParams gaff2.frcmod \n',  ])
    for name in newfrcmod:
        if name != False :
            tleapinput.writelines( 'loadAmberParams %s.frcmod \n' % name  )

    if MOLNAME :
        tleapinput.writelines( 'loadoff %s.lib \n' %MOLNAME)
    tleapinput.writelines( [ 'new = sequence { ACE   %s ALA NME} \n '  %(MOLNAME) ,  'source leaprc.water.tip3p \n', 'savepdb new myprotein.pdb\n' ,'solvatebox new TIP3PBOX 15 \n' , 'saveAmberParm  new  newwat.prmtop newwat.rst7\n ' ] )
    tleapinput.writelines( [  'quit'])
    tleapinput.close()
    os.system( 'tleap -f temp2.in >tleap' )


def smooth (energy):
    adjust=0
    for i in range(-2,energy.shape[0]-2):      # remove outliers : 3 max otherwise the whole daset is discarded!
        rollinground= np.round(((energy[i]+ energy[i+1] +energy[i+2])/3),6)
        if abs(rollinground   -   energy[i+1]) > threeshold:
            energy[i+1] =  np.round(((energy[i]+ energy[i+2])/2),6)
            adjust+=1
    return energy,adjust



if __name__ =='__main__':
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
        parser.add_argument("-f",   nargs="+", default='*-reo.mol2', help='input mol2 (default: %(default)s)')
        parser.add_argument("-res", type=str, default='LIG', help='residue name : script will fail if wrong ')
        parser.add_argument("-charge", type=int, default=0, help='charge of the molecule  (default: %(default)s)')
        parser.add_argument("-m", type=float, default=1, help="multiplicity (default: %(default)s)")
        parser.add_argument("-theo", type=str, default='HF', help="theory possible blyp , HF ... (default: %(default)s)")
        parser.add_argument("-n_cpu", type=int, default='10', help="number of cpu alocated for gaussian if 0 the gaussian jobs won't be run and script will stop. pcopy all log files here and run the script with -gen-skip defaults %(default)s)")
        parser.add_argument("-n_conf_max", type=int, default='5', help="number of conformations for each torsion")

        parser.add_argument("-max_tor", type=int, default='5', help="number of conformations for each torsion")
        parser.add_argument("-skipgen",  action='store_true', default=False, help="Used for debug or if gaussian are run elsewhere, will skip g09 preparation")
        args = parser.parse_args()
        theo = args.theo
        inputmol2_list= args.f
        res_name=  args.res
        charge=  args.charge
        multiplicity= args.m
        n_conf= min( args.n_conf_max ,  len(inputmol2_list) )
        max_tor = args.max_tor
        inputmol2=args.f[0]

        all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name , torsion_names , dihedrals_heavy_index= find_dihedrals(inputmol2)
        list_of_torsions =  dihedral_heavy_name
        for i,name in enumerate(  torsion_names) :
                if 'CHI' in name and  int(torsion_names[i][3:]) > max_tor:
                    torsion_names[i] =False

        if args.skipgen==False :
            prepfiles(dihedral_heavy_name, dihedrals_heavy_index, torsion_names,theo)
            if args.n_cpu != 0 :
                run_allconf(torsion_names, n_cpu=args.n_cpu)
                print(conf_per_torsion)
                #sys.exit()
                #run_allconf(['PHI'], n_cpu=args.n_cpu)
            else : print('input for gaussian ready, please run manually')
        else :conf_per_torsion = np.ones(len(torsion_names))*n_conf   ## to be fixed

        all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names_full, dihedrals_heavy_index , all_dihedrals_type_index = find_dihedrals_mod('input.prmtop')

        for p in torsion_names:
            #print(p)
            if p !=False:
                energies_torsions=np.array([])
                for j in range(0, n_conf):
                    prefix = p + '_' + str(j)  + '_'
                    extract_stationary_structures(prefix, Opt=False)
                    '''
                    nframes= make_traj(prefix)
                    if  nframes :

                        prepare_paramfit_job_files (prefix, nframes)
                    '''

        # here postparamfit
        inputfile=inputmol2
        make_missing_parms_parmck2(inputfile)

        fileparamas=open('gaff2.frcmod')
        line = fileparamas.readline()
        while line[0:4] != 'IMPR' :
            line = fileparamas.readline()
        line = fileparamas.readline()
        list_improper=[]
        while line[0:1] != '\n' :
            list_improper.append([line.split('  ')[0].split('-'),line[12:].split()[:3]] )
            line = fileparamas.readline()

        paramsout=[]
        for torsion in range(len(torsion_names)):
            if  torsion_names[torsion]!=False:
                maxorder=3
                energies_torsions=np.array([])
                #test multiconf  TEST

                #print(make_traj(torsion_names[torsion]+'_' + '%s_'%0,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]]))


                n_dih=len(all_dihedrals_type_index[torsion])
                #print(n_dih)
                #test multiconf  TEST
                #data=np.loadtxt(torsion_names[torsion]+'_-energySP.dat')
                for i in range(int(conf_per_torsion[torsion ])) :
                    print(torsion_names[torsion],i,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]])
                    print(torsion_names[torsion]+'_%s_' %i,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]])
                    nframes , x = make_traj(torsion_names[torsion]+'_%s_' %i,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]])
                    print(x)

                    sdtmin=12
                    #print(torsion_names[torsion]+'_' + '%s_'%i,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]])
                    #make_traj(torsion_names[torsion]+'_' + '%s_'%i,[[dih[0],dih[1], dih[2],dih[3]] for dih  in all_dihedrals_type_index[torsion]])
                    print('%s %s/rmSP-UBE-fromqm.py  --traj %s_%s_-traj.mdcrd --onlynonbonded   --top  input.prmtop --qm  %s_%s_-energySP.dat --NBE   %s_%s-nonbondedenergySP.dat --RBE %s_%s-residual-energySP.dat'  %(sire_folder, lib_folder,torsion_names[torsion],i, torsion_names[torsion],i,torsion_names[torsion],i, torsion_names[torsion],i))
                    os.system( '%s %s/rmSP-UBE-fromqm.py  --traj %s_%s_-traj.mdcrd --onlynonbonded   --top  input.prmtop --qm  %s_%s_-energySP.dat --NBE   %s_%s-nonbondedenergySP.dat --RBE %s_%s-residual-energySP.dat'  %(sire_folder, lib_folder,torsion_names[torsion],i, torsion_names[torsion],i,torsion_names[torsion],i, torsion_names[torsion],i))
                    energy=np.loadtxt(  '%s_%s-residual-energySP.dat'  %( torsion_names[torsion],i)  )
                    threeshold=20   #kcal
                    #prefix = torsion_names[torsion] + '_' + str(i)  + '_'
                    #if np.std(energy) < sdtmin :

                    sdtmin=np.std(energy)

                    if energy.shape[0] == len(range(-180,180,10))  and max(energy) <100 :
                        energy , adjust=smooth(energy)
                        print('energysmooth')
                        if adjust < 6 :   # how many times can we cheat ?

                            if len(energies_torsions) == 0 :
                                    energies_torsions=energy
                                    data = [ energy ]
                                    deriv=np.gradient(energy)
                                    #rectify edges for periodicity:
                                    deriv[0]=(energy[1]-energy[-1])/2
                                    deriv[-1]=(energy[-2]-energy[0])/2
                                    derivdata = [deriv]
                                    X = [x.T]

                                    plt.plot(energies_torsions,range(-180,180,10) )
                                    plt.savefig('energyresidual%s'%(torsion_names[torsion ]))
                                    outputenergy=open('%s-%s-residual-energySP.dat' %(torsion_names[torsion],i) ,'w')
                                    outputenergy.writelines( [ '%.8f \n' %e for e in  energy])
                                    outputenergy.close()
                            else :
                                energies_torsions=energies_torsions + energy
                                data.append( energy )
                                X.append(x.T)
                                deriv[0]=(energy[1]-energy[-1])/2
                                deriv[-1]=(energy[-2]-energy[0])/2
                                derivdata.append(deriv)
                        else :  conf_per_torsion[torsion ]=-1
                    else : conf_per_torsion[torsion ]=-1




                # conversion hartree to kcal/mol : 1 hartree = 627.509 kcal.mol-1
                #data*= 627.509

                k=np.min(data)
                data=data-k
                k=np.min(derivdata)
                derivdata=derivdata-k
                #print(data)
                #result ,inputfct=leastsquare_with_prem_multi(X,derivdata, all_dihedrals_type_index[torsion],all_dihedrals_type[torsion], maxorder)
                if torsion_names[torsion]=='CHI2':
                    print(len(x) , len(x.T))
                    print(len(X), len(X[0]))
                    print(conf_per_torsion[torsion ])


                if conf_per_torsion[torsion ] >0 :
                    result ,inputfct=leastsquare_with_prem_multi_data(X,data, all_dihedrals_type_index[torsion],all_dihedrals_type[torsion], maxorder)
                    #result ,inputfct=leastsquare_with_prem(x.T,data, all_dihedrals_type_index[torsion],all_dihedrals_type[torsion], maxorder)
                    final = derivdata[-1] + result.residual[-len(derivdata[-1]):]

                    report_fit(result)

                    paramsout.append(write_frcmod(result,torsion_names[torsion],  all_dihedrals_type,all_dihedrals_type_index))

                    outputfct=[]
                    for entry in inputfct :
                        entry.append(result.params[entry[2]].value)
                        entry[1]=entry[1]
                        outputfct.append(entry)
                    xvalues, fitteddata=comp_with_derivdata3(outputfct)
                    #####some multi adjustement
                    derivdata=np.array(derivdata)
                    derivdata=derivdata.flatten()
                    xvalues=np.array(xvalues)
                    xvalues=xvalues.flatten()
                    residuals=np.array(result.residual)
                    residuals=residuals.flatten()
                    fitteddata=np.array(fitteddata)
                    fitteddata=fitteddata.flatten()
                    data=np.array(data)
                    data=data.flatten()
                    print( len(derivdata), len( derivdata+residuals) ,len(xvalues))
                    plot( derivdata, derivdata+residuals,xvalues, outputfile=torsion_names[torsion]+'fittingresiduals.png' )
                    fitteddata=fitteddata-np.min(fitteddata)
                    fitteddata.flatten()
                    plot(derivdata,fitteddata,  xvalues , outputfile=torsion_names[torsion]+'fitted-function.png')
                    xvalues, fitteddata=comp_with_data3(outputfct)
                    fitteddata=fitteddata-np.min(fitteddata)
                    plot(data,fitteddata,  xvalues,outputfile=torsion_names[torsion]+ 'fit-dihedral.png' )
        add_new_parms(inputfile,torsion_names,res_name )


        for torsion in range(len(torsion_names)):
            if torsion_names[torsion]!=False :
                for i in range(int(conf_per_torsion[torsion ])) :
                    print( '%s %s/rmSP-UBE-fromqm.py --traj %s_%s_-traj.mdcrd   --top  input.prmtop --qm   %s_%s_-energySP.dat --TE    %s_%s-gaff2.dat --NBE  %s_%s-residual-gaff2.dat'  %(sire_folder, lib_folder,torsion_names[torsion], i, torsion_names[torsion], i,torsion_names[torsion], i, torsion_names[torsion], i))
                    os.system( '%s %s/rmSP-UBE-fromqm.py  --traj  %s_%s_-traj.mdcrd  --top  input.prmtop --qm  %s_%s_-energySP.dat --TE    %s_%s-gaff2.dat --NBE  %s_%s-residual-gaff2.dat'  %(sire_folder, lib_folder,torsion_names[torsion], i, torsion_names[torsion], i,torsion_names[torsion], i, torsion_names[torsion], i))
                    print( '%s %s/rmSP-UBE-fromqm.py --traj  %s_%s_-traj.mdcrd   --top  newtopol.prmtop --qm   %s_%s_-energySP.dat --TE   %s_%s-fitted.dat --NBE %s_%s-residual-fitted.dat'  %(sire_folder, lib_folder,torsion_names[torsion], i, torsion_names[torsion], i,torsion_names[torsion], i, torsion_names[torsion], i))
                    os.system( '%s %s/rmSP-UBE-fromqm.py  --traj  %s_%s_-traj.mdcrd  --top  newtopol.prmtop --qm   %s_%s_-energySP.dat --TE    %s_%s-fitted.dat --NBE  %s_%s-residual-fitted.dat'  %(sire_folder, lib_folder,torsion_names[torsion], i, torsion_names[torsion], i,torsion_names[torsion], i, torsion_names[torsion], i))

                    qm= np.loadtxt('%s_%s_-energySP.dat'  %(torsion_names[torsion],i )  )
                    gaff2= np.loadtxt('%s_%s-gaff2.dat'  %(torsion_names[torsion],i )  )
                    fitted= np.loadtxt('%s_%s-fitted.dat'%(torsion_names[torsion],i )  )
                    print(torsion_names[torsion]  )
                    angle=np.arange(-180,180,10)
                    angle=angle[:qm.shape[0]]
                    plt.figure()
                    plt.plot(angle, qm*627.509-min(qm*627.509) , 'k')
                    plt.plot(angle, gaff2, 'r')
                    plt.plot(angle, fitted , 'b')
                    plt.savefig('fit-%s-%splot.png'  %(torsion_names[torsion],i))
                    plt.close()
                    mse = sklearn.metrics.mean_squared_error(qm,fitted)
                    msegaff2 = sklearn.metrics.mean_squared_error(qm,gaff2)
                    rmse = math.sqrt(mse)
                    rmsegaff2 = math.sqrt(msegaff2 )
                    print(i , rmse)
                    print(i , rmsegaff2)
