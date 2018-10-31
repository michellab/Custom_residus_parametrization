from extract_Stationary_structures import extract_stationary_structures


input='PSI_-70.mol2'

def get_gaff2_parmcheck(input):

     os.system('parmchk2 -f mol2 -i %s -o gaff2.frcmod  -s  gaff2' %(input))


extract_stationary_structures('CHI1_')
