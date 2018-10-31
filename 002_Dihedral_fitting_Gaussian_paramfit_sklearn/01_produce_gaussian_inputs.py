#! /usr/bin/env python
d=' \
############################ \
Prepare gaussian input from a mol2 file. \
This is an ugly but efficient script to generate conformations spaced from 10 degrees using tleap and \
manually parse the generated mol2 files to gaussian input files.  \n \n \
It is not necessarry but advised to rename your atoms in your input : CA CB C N ... as the script will try to assign PHI PSI and CHI\
atomtypes will be assigned but won\'t be used in calculations \n\n \
Known bugs : -parmchk2 doesn\'t know CX (The bug has been reported ) atom type (CA) so add it to to parmcheck file or change the CX to CT\
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


def main(list_of_atm,  dihedrals_heavy_index):

    list_of_key={}

    tleap_input_file=open('tleap.in', 'w')
    inputfile=open(inputmol2, 'r')
    inputfiliecontent=inputfile.readlines()
    for l in range(len(inputfiliecontent)):
        if res_name in inputfiliecontent[l] :
            list_of_key.update({inputfiliecontent[l].split()[1]: inputfiliecontent[l].split()[0]})

    inputfile.close()

    tleap_input_file.writelines(['source leaprc.protein.ff14SB \n' , 'm = loadmol2 %s  \n' %(inputmol2) , 'relax m \n' ])

    for bond in range(len(list_of_atm)) :

      for angle in range(-180,+180,10) :
        tleap_input_file.writelines( [ 'impose m  { 1 } {  { %s %s } }\n' %(list_of_atm[bond][0] , angle)])
        tleap_input_file.writelines( [ 'savemol2 m %s_%s.mol2 1\n' %(torsion_names[bond] , angle)] )

    tleap_input_file.writelines(['quit' ])
    tleap_input_file.close()
    print (' tleap -f tleap.in')
    subprocess.call('tleap -f tleap.in' ,shell=True  )
    for finp in os.listdir('./') :
             if finp[-4:]=='mol2' and finp[:3] in ['PHI', 'PSI' , 'CHI']:
                fi=open(finp, 'r')
                fo=open('%s' %(finp[:-4]+'gau'), 'w')
                fo.write('--Link1-- \n%%RWF=%s\n%%NoSave\n%%chk=%s \n#P b3lyp/6-31G* Opt=(ModRedundant,MaxCycle=400) Freq SCF=(Conver=6,MaxCycle=400) Pop=NoMBS\n \nparamfit run\n  \n%s    %s \n'  %(finp[:-5],finp[:-5],charge,multiplicity))

                while 'ATOM' not in fi.readline()  : continue
                i=True
                while i==True :
                        line=fi.readline().split()
                        if 'TRIPOS' not  in line[0] :
                            #lineout=line[1][0]+'-' + line[5] + '-' + line[8]  + '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'

                            lineout=line[1][0]+ '    ' + line[2] +  '    ' + line[3] +  '    ' + line[4]  + '\n'
                            fo.write(lineout)
                        else:
                            i=False
                            fo.write('\n')
                numb=torsion_names.index(finp[0:4].replace('_',''))
                #listdih=list_of_atoms[numb]

                listdih= dihedrals_heavy_index[numb]
                #fo.write('F'+str(list_of_atoms[numb][0])+'\nF'+str(list_of_atoms[numb][1])+'\nF'+str(list_of_atoms[numb][2])+'\nF'+str(list_of_atoms[numb][3])+'\n')
                fo.write('D ' + str(listdih[0]) +' '+ str(listdih[1])+ ' '+ str(listdih[2])+ ' '+str(listdih[3])+ ' F\n')
                fo.close()
                fi.close


if __name__ =='__main__':
        parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=d, epilog=" ")
        parser.add_argument("-f", type=str, default='Mol-sm_m1-c1.mol2', help='input mol2 (default: %(default)s)')
        parser.add_argument("-res", type=str, default='LIG', help='residue name : script will fail if wrong (default: system[HMR]_REST1.XXX.prmtop)')
        parser.add_argument("-charge", type=int, default=0, help='charge of the molecule  (default: %(default)s)')
        parser.add_argument("-m", type=float, default=1, help="multiplicity (default: %(default)s)")
        args = parser.parse_args()

        inputmol2= args.f
        res_name=  ' ' + args.res
        charge=  args.charge
        multiplicity= args.m
        all_dihedrals, all_dihedrals_type, dihedrals_heavy, dihedral_heavy_name , dihedrals_heavy_centralbond_name, dihedrals_heavy_index= find_dihedrals(inputmol2)
        list_of_torsions =  dihedral_heavy_name
        #torsion_names_dict =[  'CHI1', 'CHI2' , 'CHI3' ,'CHI4', 'CHI5','CHI6', 'CHI7', 'CHI8' , 'CHI9' ,'CHI10', 'CHI11', 'CHI12','CHI13', 'CHI14', 'CHI15' , 'CHI16' ,'CHI17' , 'CHI18', 'CHI19']
        torsion_names=[]
        i=1
        for bond in dihedrals_heavy_centralbond_name:
            bond.sort()

            if ['CA' ,'N'] == bond:                 torsion_names.append('PHI')
            elif [ "C" , "CA" ] == bond :                torsion_names.append('PSI')
            else :
                torsion_names.append('CHI'+str(i) )
                i+=1



        main(dihedral_heavy_name, dihedrals_heavy_index)
