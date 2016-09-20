# MD_Pipeline v 1.1

md_pipeline.py is a beta testing version of the python script to streamline and 
automise coarse grained simulations of membrane-protein systems.

--------------------------------------------------------------------------------
 script organisation
--------------------------------------------------------------------------------

1. parse_options.py
Parsing the config.txt file for input options 

2. preprocess.py
Orienting proteins for embedding in membrane is implemented 
Preprocessing PDB file & making DSSP file

3. runprocess.py
--> convert ATM to CG protein using martinize script

--> make protien box using protein_dimensions.py

--> make protein copies using genconf (gromacs)

--> make CG membrane and memb-prot system using insane script

--> make index using make_ndx (gromacs)

--> make final topology (topol.top)

--> energy minimization

--> isothermic-isobaric equilibration

--------------------------------------------------------------------------------
 usage
--------------------------------------------------------------------------------

Download all files
Copy all files to PWD
Edit the config.txt file provided
Run python md_pipeline.py

--------------------------------------------------------------------------------
 comments
--------------------------------------------------------------------------------

This verion works for both integral membrane proteins and most of the peripheral protiens.
The error with making beads for multi-chain proteins(only for some proteins) is yet to be fixed.

