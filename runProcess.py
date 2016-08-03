import os

import preprocess


def runPreProcess(all_flags):
    for option in all_flags:
        if option[0] == 'preprocess_options':
            t = option[1]
            t[0] = t[0].strip()
            t[1] = t[1].strip()
            if t[0] == '-protein' or 'protein' in t[0]:
                clean_file = preprocess.clean_pdb(t[1])
                dssp_out = preprocess.runDSSP(clean_file)
                return clean_file, dssp_out
    return False, False


def runMaritinize(all_flags, clean_pdb, dssp_file):
    mar_flags = {}
    for option in all_flags:
        if option[0] == 'martinize_options':
            for t in option[1:]:
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

    cg_protein = clean_pdb[:-4] + '_CG.pdb'
    cg_topol = clean_pdb[:-4] + '_CG.top'
    cg_index = clean_pdb[:-4] + '_CG.ndx'
    nmap = clean_pdb[:-4] + '_map.ndx'

    martinize_flags = martinize_flags + ' -f ' + clean_pdb + ' -ss ' + dssp_file + ' -x ' + cg_protein + ' -o ' + cg_topol + ' -n ' + cg_index + ' -nmap ' + nmap
    run_martinize = 'python martinize.py ' + martinize_flags
    os.system(run_martinize)
    return cg_protein, cg_topol, cg_index, nmap


def multiplyProtein(all_flags, clean_pdb):
    import protein_dimensions
    protein_dim_temp = protein_dimensions.box_dimension(clean_pdb)
    x_dim = protein_dim_temp[1]
    y_dim = protein_dim_temp[2]
    z_dim = protein_dim_temp[3]
    multi_flags = {}
    for option in all_flags:
        if option[0] == 'multiProt_options':
            for t in option[1:]:
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                multi_flags[t[0]] = t[1]

    new_x_dim = []
    new_y_dim = []
    new_z_dim = []

    for key in multi_flags:
        if key == '-x_num':
            for i in range (0,len(x_dim)):
                new_x = float(x_dim[i])

    return True

    '''

    run_genconf = gmx genconf -f protein.pdb -o multiprot.gro -nbox 3 3 1
    os.system(run_genconf)
    return multi_prot
'''


def runInsane(all_flags, clean_pdb, cg_protein, cg_topol):
    import protein_dimensions
    protein_dim_temp = protein_dimensions.box_dimension(clean_pdb)
    protein_dim = protein_dim_temp[0]
    print protein_dim
    sane_flags = {}
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
        if float(abs(protein_dim[i]) / 10) > insane_dim[i]:
            print "Error in specified insane dimensions"
            print "Protein is bigger than box"
            print "protein dim: ", protein_dim
            print "insane dim : ", insane_dim
            return False

    insane_flags = ''
    lipids = []
    for key in sane_flags:
        if key == '-l' or key == '-u':
            luflags = sane_flags[key].split()
            for luflag in luflags:
                if luflag not in lipids:
                    lipids.append(luflag)
                insane_flags = insane_flags + ' ' + key + ' ' + luflag
        else:
            insane_flags = insane_flags + ' ' + key + ' ' + sane_flags[key]

    system = clean_pdb[:-4] + '_inMemb.gro'
    system_top = clean_pdb[:-4] + '_inMemb.top'
    insane_flags = insane_flags + ' -f ' + cg_protein + ' -o ' + system + ' -p ' + system_top
    run_insane = 'python insane.py ' + insane_flags
    os.system(run_insane)
    return True, system, system_top, lipids


def make_ndx(system, lipids):
    ndx_flags = ''
    for lipid in lipids:
        lipid = lipid.split(':')
        ndx_flags = ndx_flags + ' ' + lipid[0]
    system_ndx = system[:-4] + '.ndx'
    index = 'gmx make_ndx -f ' + system + ' -o ' + system_ndx + ' <<EOF\n 1 | r ' + ndx_flags + ' \n 11  &  !  r ' + ndx_flags + ' \n  q \n EOF'
    os.system(index)
    return system_ndx


def make_topology(cg_topol, system_top):
    first_lines = '#include "martini_v2.2.itp" \n#include "martini_v2.0_lipids.itp" \n#include "martini_v2.0_ions.itp"\n'

    with open(cg_topol) as fin:
        cg_lines = fin.readlines()
    new_cg_lines = []
    for i in range(0, len(cg_lines)):
        if '"martini.itp"' in cg_lines[i]:
            pass
        else:
            new_cg_lines.append(cg_lines[i])

    with open(system_top) as fin:
        sys_lines = fin.readlines()
    new_sys_lines = []
    for i in range(0, len(sys_lines)):
        if 'molecules' in sys_lines[i]:
            start = i + 3
            for j in range(start, len(sys_lines)):
                '''
                if '+' in sys_lines[j] or '-' in sys_lines[j]:
                    tline = sys_lines[j]
                    tline = tline.split()
                    charge = '+-'
                    for char in charge:
                        tline[0] = tline[0].replace(char, '')
                    newtline = tline[0]+"           "+ tline[1]+'\n'
                    new_sys_lines.append(newtline)
                else:
                '''
                new_sys_lines.append(sys_lines[j])
            break

    topology = 'topol.top'
    with open(topology, 'w') as fout:
        fout.write(first_lines)
        for line in new_cg_lines:
            fout.write(line)
        fout.write('\n')
        for line in new_sys_lines:
            fout.write(line)

    return topology


def runMinimization(all_flags, system, topology, system_ndx):
    em_flags = {}
    for option in all_flags:
        if option[0] == 'energy_min_options':
            for t in option[1:]:
                em_flags[t[0]] = t[1]

    min_flags = ''
    for key in em_flags:
        min_flags = min_flags + ' ' + key + ' = ' + em_flags[key] + '\n'

    em_mdp = 'em.mdp'
    em_lines = "### this is an energy minimization parameter file (em.mdp) ### \n"
    with open(em_mdp, 'w') as fout:
        fout.write(em_lines)
        for line in min_flags:
            fout.write(line)

    em_gro = system[:-4] + '_EM.gro'
    em_tpr = system[:-4] + '_EM.tpr'

    run_em_grompp = 'gmx grompp -f ' + em_mdp + ' -c ' + system + ' -p ' + topology + ' -n ' + system_ndx + ' -o ' + em_tpr
    os.system(run_em_grompp)
    run_em_mdrun = 'gmx mdrun -s ' + em_tpr + ' -v -deffnm ' + system[:-4] + '_EM'
    os.system(run_em_mdrun)

    return em_gro


def runEquilibration(all_flags, system, em_gro, topology, system_ndx):
    groups = []
    with open(system_ndx) as fin:
        lines = fin.readlines()

    for line in lines:
        if '[' in line:
            tline = line.strip()
            tline = tline.lstrip('[').rstrip(']')
            groups.append(tline)
    groups = groups[-2:]

    eq_flags = {}
    for option in all_flags:
        if option[0] == 'equilibration_options':
            for t in option[1:]:
                if 'grps' in t[0]:
                    eq_flags[t[0]] = groups[0] + groups[1]
                else:
                    eq_flags[t[0]] = t[1]

    equil_flags = ''
    for key in eq_flags:
        equil_flags = equil_flags + ' ' + key + ' = ' + eq_flags[key] + '\n'

    equil_mdp = 'equil.mdp'
    equil_lines = "### this is an npt equilibration parameter file (equil.mdp) ### \n"
    with open(equil_mdp, 'w') as fout:
        fout.write(equil_lines)
        for line in equil_flags:
            fout.write(line)

    equil_gro = system[:-4] + '_EQUIL.gro'
    equil_tpr = system[:-4] + '_EQUIL.tpr'

    run_equil_grompp = 'gmx grompp -f ' + equil_mdp + ' -c ' + em_gro + ' -p ' + topology + ' -n ' + system_ndx + ' -o ' + equil_tpr + ' -maxwarn 10'
    os.system(run_equil_grompp)
    run_equil_mdrun = 'gmx mdrun -s ' + equil_tpr + ' -v -deffnm ' + system[:-4] + '_EQUIL'
    os.system(run_equil_mdrun)

    return equil_gro
