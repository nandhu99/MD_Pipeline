import  parse_options

martinize_flags, preprocess_flags = parse_options.parseconfig('config.txt')
print martinize_flags
print preprocess_flags

"""
pdbfile = 'test.pdb'
clean_pdb.clean_pdb(pdbfile)
clean_pdb.runDSSP(clean_pdbfile)
"""