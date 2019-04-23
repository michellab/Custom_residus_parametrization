from multiprocessing import Pool
import os
import parmed
import numpy as np
import sys
import subprocess
import string
from multiprocessing import Process
import time
import copy

def run_resp(input):
    #prepare gaussian files run gaussian jobs and read gaussian files

    os.system('antechamber -i %s.pdb -fi pdb -o %s.gau  -fo gcrt -gv 2 -ge %s.gesp'  %(input, input, input))
    #Uncomment the next if statement to rerun gaussian jobs
    if os.path.isfile('./%s.log'%input) == False or  'Normal termination' not in  subprocess.getoutput('tail -1 ./%s.log'%input)  :
        os.system('g09 %s.gau' %input  )
    os.system('antechamber -i  %s.log -fi gout -o %s.mol2 -fo mol2 -c resp -at amber   -eq 2  -pf   -dr n' %(input,input) )
    #time.wait(1)

def sdf2pdb(input):
    # read the sdf (input) and convert it to multiple pdb and return the number of conformations
    tobeparsed= subprocess.getoutput('babel -isdf %s -O .pdb --split' %(input) )
    n_conf=tobeparsed.split('\n')[-2].split()[0]

    os.system('sed -i  s/HETATM/ATOM\ \ /g  *pdb' )
    #os.system('grep -v ATOM  1.pdb > model.pdb')
    return n_conf


def run_allconf(n_conf, n_cpu=10) :
    #call run_resp for all conformation using multiprocessing
    # !!!!! For a RE-RUN  DO NOT USE multiprocessing :
    # 10 instances of antechamber running at the same time seems to break it ! use the run_resp in the 'for' loop

    args=[]
    for i in range(1,n_conf+1):
      args.append('%s' %(i))
      #run_resp('%s' %(i))
    p =Pool(processes=n_cpu)
    r = p.map_async(run_resp,args)
    print(args)
    r.wait()



def return_charges(parm, i=0):
    #remove unwanted fragment (any heavy atom fwhich name finishes with a number and bounded Hydrogen)

    frag=[]
    charges=[]

    for atom in  parm.atoms:

      if atom.name[-1].isdigit() ==False and atom.element not in [1] :

        frag.append(atom.idx)
        charges.append(atom.charge)
        for atm in atom.bond_partners:

            if atm.element == 1 :

                frag.append(atm.idx)
                charges.append(atm.charge)
    ####TO DO : here check if similar H had same charge independant
    totcharge= np.sum(charges)
    #print(charges)
    newcharges=np.add( charges ,-totcharge/len(charges) )
    new_parm = recharge(newcharges,frag , parm)
    # could save here for each using
    # new_parm.save(conf%i.mol2 %i)
    return frag, newcharges


def multiple(n_conf):
    # run return_charges and average all charge for a same atom
    I=[]

    C=[]
    for file in range(1,n_conf):
        #print(file)
        ##### TO DO :check index is conserved here :
        parm = parmed.load_file('%s-reo.mol2'%file)
        index, charges =return_charges(parm)
        I.append(index)
        #print(charges)
        if C!=[]  : C= np.concatenate( (C, [charges]) , axis = 0 )
        else :  C=[copy.deepcopy(charges)]

    avg_charge= np.average(C, axis=0)
    print(avg_charge , index)
    new_parm = recharge(avg_charge,index , parm)
    parm.save('newcharges.mol2')

def recharge(newcharges,frag , parm):
    for atom in  parm.atoms:
      if atom.idx in frag :
         print(atom.name,atom.idx, frag.index(atom.idx),np.where(frag ==atom.idx ) )
         #atom.charge=newcharges[np.where(frag ==atom.idx )[0][0]]
         atom.charge=newcharges[frag.index(atom.idx)]
      else:
        atom.charge=0
        print(atom.name, atom.charge )
    return parm

def rename_atoms(n_conf, Calpha, LastAtomInSidechain=False , extin='.mol2', extout='.mol2' ):

    ABC= list(string.ascii_uppercase)
    for i in range(1,n_conf):
        #os.system('antechamber -i %s%s -fi pdb -o %sout%s  -fo pdb' %(i,extin,i,extin))
        print('loading : %s%s'%(i,extin))
        parm = parmed.load_file('%s%s'%(i,extin))
        AmberOrder=['N', 'H','C', 'O', 'CA']
        k=1
        for atom in  parm.atoms:
            if atom.name in AmberOrder :
                atom.name = 'R' +str(k)
                atom.number=55
            k+=1
        prev_order=[-1,-1,-1,-1,-1]
        current=[]
        index=1
        for atom in  parm.atoms:


            if atom.name=='CA' or atom.idx==Calpha :
                atom.name=AmberOrder[4]
                atom.type='CX'
                prev_order[4]= atom.idx
                atom.number=4
                for atm in atom.bond_partners:

                    if atm.element==7:
                      prev_order[0]= atm.idx
                      atm.name=AmberOrder[0]
                      atom.number=0
                      for at in atm.bond_partners:
                        if at.element==1 :
                          prev_order[1]= at.idx
                          at.name=AmberOrder[1]
                          at.number=1
                    elif atm.element==6:

                        if 8 in [at.element for at in atm.bond_partners] :
                            prev_order[2]= atm.idx
                            atm.name=AmberOrder[2]
                            atm.number=2
                            for at in atm.bond_partners:
                                if at.element==8 :
                                  prev_order[3]= at.idx
                                  at.name=AmberOrder[3]
                                  at.number=3
                        else:
                            current.append(atm)
                            for atom in current :
                                    if len(current) == 1 : s=''
                                    else:  s=str(len(current))
                                    atm.name='C'+s+'B'
                                    prev_order.append(atm.idx)
                                    atm.number=6


                    elif atm.element==1 :
                         atm.number=4
                         atm.name='HA'


        if -1 in prev_order:
            parm.save( '%s-debug.mol2'%i )
            sys.exit('debug: missing atoms in backbone confromation : %s ' %i)
        if len(current) ==0 : sys.exit('debug: no next atom')


        index+=1
        while len(current) >0:
                #parm.save( '%s-reo.mol2'%index )
                #print('atom : '+str(index))
                current, index , parm, prev_order=find_next(current, index, parm, prev_order)
                #print('after atom : '+str(index))
                for atom in current :
                    print(LastAtomInSidechain,atom.idx)
                    if LastAtomInSidechain and atom.idx == LastAtomInSidechain :
                        last = current.pop(current.index(atom))

                        for at in last.bond_partners:
                            if  at.element==1 :
                                Hbound.append(at.idx)

                        if len(Hbound)==1 :j=''
                        else: j=1
                        for atm in last.bond_partners:
                            if atm.idx in Hbound :
                                atm.name= 'H' +ABC[index-1] +str(j)
                                if j!='' :   j+=1
                                atm.number=index+4
                prev_order.append(index)
                #print('after atom end of loop:' +str(index) )
                print('next :' +str(len(current)))
                Hbound=[]
        new=len(prev_order)+1
        for atom in  parm.atoms:
            if atom.idx not in prev_order :
                atom.number=new
                new+=1
                #for  i in Hbound : prev_order.append(i)
        #parm.atoms.sort(key=lambda x: x.number)
        #for atom in  parm.atoms:
            #print(atom.number)

        #print(prev_order)
        parm.save( '%s-reo.mol2'%i )

def find_atom_by_name(parm, name):
    for atom in  parm.atoms:
        if atom.name==name:
            return atom
def find_atom_by_idx(parm, idx):
    for atom in  parm.atoms:
        if atom.idx==idx:
            return atom

def find_next(current_atoms, index, parm, previous):
    ABC= list(string.ascii_uppercase)
    atom_element= dict({6:'C' , 16:'S' , 8:'O', 7:'N'  })
    Hbound=[]
    Heavybound=[]
    next_atoms=[]
    for current in current_atoms:
        for at in current.bond_partners:
            #print(at.element, at.number)
            if  at.element==1 :
                Hbound.append(at.idx)
                print(Hbound)
            elif at.element!=1 and at.number==-1 or at.idx not in previous  :
                Heavybound.append(at.idx)
                #print(Heavybound)
                next_atoms.append(at)

        if len(Heavybound) == 1 : s=''
        else : s= 1 #str(len(Heavybound))
        if len(Hbound) == 1 : i=''
        else : i=1
        for atm in current.bond_partners:
            if  atm.idx in Heavybound :
                #print(index)
                atm.name=atom_element[atm.element] + str(s) +ABC[index]
                #print(atom_element[atm.element] + str(s) +ABC[index])
                if s != '' : s+=1
                #print(atm.idx)
                atm.number=index+5
            elif atm.idx in Hbound :
                atm.name= 'H' +ABC[index-1] +str(i)
                if i!='' : i+=1
                atm.number=index+4
    prev_order= previous+Hbound+Heavybound
    #if index==6 :sys.exit()
    return  next_atoms, index+1 , parm, prev_order


input=sys.argv[1]
Calpha=int(sys.argv[2])
LastAtomInSidechain=int(sys.argv[3])
Calpha += -1
LastAtomInSidechain+=-1

n_conf= sdf2pdb(input)

run_allconf(int(n_conf))
print('There is %s conformations in input file'%n_conf)
'''
#rename_atoms(int(n_conf), Calpha=5-1,LastAtomInSidechain=17-1, extin='.mol2')
rename_atoms(int(n_conf), Calpha=Calpha,LastAtomInSidechain=LastAtomInSidechain, extin='.mol2')
multiple(int(n_conf))
'''
