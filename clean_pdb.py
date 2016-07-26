import sys, os

#Clean PDB files
with open(sys.argv[1]) as fin:
    lines = fin.readlines()


outfile = 'clean_'+sys.argv[1]
fout = open(outfile, 'w')

nmr_model = False

for line in lines:
    if line[0:5] == 'MODEL':
        nmr_model = True
    if line[0:4] == 'ATOM':
        fout.write(line)
    if line[0:6] == 'ENDMDL' and nmr_model:
        break

fout.close()

# Generate DSSP output files

do_dssp= '~/MD_Pipeline/dssp_exe '
infile  = outfile
out_dssp = infile+".dssp"
run_cmd = do_dssp + infile + " " + out_dssp
print run_cmd
os.system(run_cmd)