import preprocess
import os


def runPreProcess(all_flags):
    for option in all_flags:
        if option[0] == 'preprocess_options':
            t = option[1]
            t[0] = t[0].strip()
            t[1] = t[1].strip()
            if t[0] == '-protein':
                clean_file = preprocess.clean_pdb(t[1])
                dssp_out   = preprocess.runDSSP(clean_file)
                return clean_file, dssp_out
    return False, False


def runMaritinize(all_flags, clean_pdb, dssp_file):
    mar_flags = {}
    for option in all_flags:
        if option[0] == 'martinize_options':

            for t in option[1:]:
                print t
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                if t[1] == 'False' or t[1] == 'false':
                    pass
                elif t[1] == 'True' or t[1] == 'true':
                    mar_flags[t[0]] = ''
                else:
                    mar_flags[t[0]] = t[1]

    martinize_flags = ''
    for key in mar_flags:
        martinize_flags = martinize_flags + ' ' + key + ' ' + mar_flags[key]

    cg_protein = clean_pdb[:-4]+'_CG.gro'
    cg_topol = clean_pdb[:-4]+'_CG.top'
    cg_index = clean_pdb[:-4]+'_CG.ndx'
    nmap = clean_pdb[:-4]+'_map.ndx'

    martinize_flags = martinize_flags+' -f '+clean_pdb+' -ss '+dssp_file+' -x '+ cg_protein+ ' -o '+ cg_topol+ ' -n '+cg_index+ ' -nmap '+nmap
    run_martinize = 'python martinize.py '+martinize_flags
    os.system(run_martinize)
    return cg_protein, cg_topol, cg_index, nmap

def checkBox(all_flags, clean_pdb):
    import protein_dimensions
    protein_dim = protein_dimensions.box_dimension(clean_pdb)
    print protein_dim
    sane_flags= {}
    for option in all_flags:
        if option[0] == 'insane_options':
            for t in option[1:]:
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                if t[1] == 'False' or t[1] == 'false':
                    pass
                elif t[1] == 'True' or t[1] == 'true':
                    sane_flags[t[0]] = ''
                else:
                    sane_flags[t[0]] = t[1]
    insane_dim = []
    for key in sane_flags:
        if key == '-x':
            insane_dim.append(float(sane_flags[key]))
        if key == '-y':
            insane_dim.append(float(sane_flags[key]))
        if key == '-z':
            insane_dim.append(float(sane_flags[key]))

    for i in range(0, len(insane_dim)):
        if protein_dim[i] > insane_dim[i]:
            print "Error in specified insane dimensions"
            print "Protein is bigger than box"
            print "protein dim: ", protein_dim
            print "insane dim : ", insane_dim
            return False

    return True