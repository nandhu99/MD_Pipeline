import clean_pdb
import  parse_options

martinize_flags, clean_flags = parse_options.parseconfig('config.txt')
print martinize_flags
"""
pdbfile = 'test.pdb'
clean_pdb.clean_pdb(pdbfile)
clean_pdb.runDSSP(clean_pdbfile)
"""