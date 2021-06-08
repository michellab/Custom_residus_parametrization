from rdkit import Chem
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, TorsionFingerprints
from rdkit.ML.Cluster import Butina
import argparse

def gen_conformers(mol, numConfs=100, maxAttempts=1000, pruneRmsThresh=0.1, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, enforceChirality=True):
	ids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, maxAttempts=maxAttempts, pruneRmsThresh=pruneRmsThresh, useExpTorsionAnglePrefs=useExpTorsionAnglePrefs, useBasicKnowledge=useBasicKnowledge, enforceChirality=enforceChirality, numThreads=0)
	return list(ids)




def smi2sdf(c,numconf):
    # read the sdf (input) and convert it to multiple pdb and return the number of conformations
    #tobeparsed= subprocess.getoutput('babel -m -isdf %s -O conf.pdb --split' %(input) )

    #n_conf=tobeparsed.split('\n')[-2].split()[0]

	numconf=len(gen_conformers(c, numConfs=numconf))
	w = Chem.SDWriter('out.sdf')
	for i in range(numconf):
		w.write(c ,i)
	w.close()
    #print(tobeparsed.split('\n')[-2])
    #os.system('sed -i  s/HETATM/ATOM\ \ /g  conf*pdb' )
    #os.system('grep -v ATOM  1.pdb > model.pdb')


def sdftomol2( ids, mol2templ):
	os.system('obabel  -isdf out.sdf -O out.mol2   -m'  )


	charges=open(mol2templ,'r')
	linescharges=charges.readlines()
	charges.close()
	for i in range(1,ids):

		input=open('out%s.mol2'%i,'r')
		linescoor=input.readlines()
		input.close()
		output=open('out%s.mol2'%i,'w')
		for i,line in enumerate(linescharges) :
			if len(line.split())== 9 :
				outputline=line[:20]+linescoor[i+1][18:47] +line[50:]

			else : outputline=line
			output.writelines(outputline)
		output.close()



parser = argparse.ArgumentParser(description='Produce one pdb conformer from smile')
parser.add_argument('--input', type=str, help='a smile string')
parser.add_argument('--inputpdb', type=str, help='a pdb input')
parser.add_argument('--output', type=str, default='out.pdb' , help='an a pdb output file')
parser.add_argument('--n_conf', type=int, default=5000 , help='number of confomrations to generate')
parser.add_argument('--mol2', type=str, default='newcharges_full.mol2' , help='takes atoms names and charges from this file')
args = parser.parse_args()
inputsmi=args.input
inputpdb=args.inputpdb
output=args.output
numconf=args.n_conf
if inputsmi:
	c=Chem.AddHs(Chem.MolFromSmiles( input))
if inputpdb:
	c=Chem.AddHs(Chem.MolFromPDBFile('../conf0.pdb'))
smi2sdf(c, numconf)
sdftomol2( args.n_conf,args.mol2 )
