# Custom_residus_parametrization-

Workflow:

 001 Use Knime and Rdkit to produce inputs for R.E.D
	-inputs can be given as sdf file or drawn directly from Marvin interface
	-run multi resp script which run multiple RESP using gaussian, rename atoms sencibly and merge reslut from run, remove extra bit of residus... As input takes the CA index (compulsory ) and a list of atoms to stop the numbering. 

	
002 Production of the .lib file containing charges and atoms types 
     Python script to produce gaussians inputs files offering a sampling over each dihedral
     Paramfit Module 
     

 004 small Simulations production for ramachadran plotting


 !!! Advisable manual Step : change atomtypes and RESIDUS name in .lib and .frcmod if overlapping with AMBER atomtypes 

