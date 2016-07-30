import parse_options
import runProcess



all_flags = parse_options.parseConfig2('config.txt')
pdb_file, dssp_file = runProcess.runPreProcess(all_flags)
print pdb_file, dssp_file
cg_protein, cg_topol, cg_index, nmap = runProcess.runMaritinize(all_flags, pdb_file, dssp_file)
bool_box = runProcess.checkBox(all_flags, pdb_file)


"""
pdbfile = 'test.pdb'
clean_pdb.clean_pdb(pdbfile)
clean_pdb.runDSSP(clean_pdbfile)
"""
