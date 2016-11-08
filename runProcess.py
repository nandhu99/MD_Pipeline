import os
import sys

import preprocess
import protein_dimensions


def runPreProcess(all_flags):
    """
    This function calls the preprocess to make a clean pdb file
    It also generates a dssp file
    :param orient_protein: protein from runMemembed function
    :return: a clean pdb file & a dssp file
    """
    prot_flag = {}
    for option in all_flags:
        if option[0] == 'preprocess_options':
            for t in option[1:]:
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                if t[1] == '-protein' or 'protein':
                    prot_flag[t[0]] = t[1]
                else:
                    pass

    in_protein = prot_flag['-protein']
    chain = prot_flag['-nos_chain']

    base_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    clean_pdb = preprocess.pdbpqr(base_dir, in_protein, chain)
    return clean_pdb


def runMemembed(clean_pdb, all_flags):
    """
    This function runs memembed program to orient protein in a membrane environment.
    :param all_flags: preprocess parameters read from 'config.txt'
    :return: returns the protein file that is oriented for a membrane system.
    """

    mem_flags = {}
    for option in all_flags:
        if option[0] == 'memembed_options':
            for t in option[1:]:
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                if t[1] == 'False' or t[1] == 'false':
                    pass
                elif t[1] == 'True' or t[1] == 'true':
                    mem_flags[t[0]] = ''
                else:
                    mem_flags[t[0]] = t[1]

    memembed_flags = ''
    for key in mem_flags:
        memembed_flags = memembed_flags + ' ' + key + ' ' + mem_flags[key]

    orient_protein = clean_pdb[:-10] + '_orient.pdb'

    memembed_flags = memembed_flags + ' -o ' + orient_protein + ' ' + clean_pdb
    run_memembed = './memembed ' + memembed_flags
    os.system(run_memembed)

    processed_pdb = clean_pdb[:-10] + '_proc.pdb'
    sep_chains = preprocess.sep_chains(orient_protein, processed_pdb)
    dssp_file = preprocess.runDSSP(processed_pdb)
    return sep_chains, processed_pdb, dssp_file


def runMaritinize(all_flags, processed_pdb, dssp_file):
    """
    This function runs martinize to convert atomistic protein to coarse-grained protein
    :param all_flags: martinize parameters read from 'config.txt'
    :param clean_pdb: oriented protein with only ATOM in pdb
    :param dssp_file: dssp file for defining secondary structures
    :return: martinized/coarse-grained protein (.pdb), coarse-grained topology, coarse-grained index & mapping file
    """
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

    cg_protein = processed_pdb[:-9] + '_CG.pdb'
    cg_topol = processed_pdb[:-9] + '_CG.top'
    cg_index = processed_pdb[:-9] + '_CG.ndx'
    nmap = processed_pdb[:-9] + '_map.ndx'

    martinize_flags = martinize_flags + ' -f ' + processed_pdb + ' -ss ' + dssp_file + ' -x ' + cg_protein + ' -o ' + cg_topol + ' -n ' + cg_index + ' -nmap ' + nmap
    run_martinize = 'python martinize_TW.py ' + martinize_flags
    print run_martinize
    os.system(run_martinize)
    return cg_protein, cg_topol, cg_index, nmap


def multiplyProtein(all_flags, cg_protein):
    """
    This function makes multiple copies of the coarse-grained protein in x-y plane
    :param all_flags: multiply protein parameters from 'config.txt'
    :param cg_protein: coarse-grained protein (.pdb)
    :return: returns a (.gro) coarse-grained protein file with multiple copies of protein along x-y plane
    :return: also returns, the protein box dimensions [x, y, z] & total number of protein copies
    """
    box_protein = protein_dimensions.box_dimension(cg_protein)
    with open(box_protein) as fin:
        lines = fin.readlines()

    x_dim, y_dim, z_dim = 0, 0, 0

    for line in lines:
        if line[0:6] == 'CRYST1':
            x_dim = float(line[7:15]) / 10
            y_dim = float(line[16:24]) / 10
            z_dim = float(line[25:33]) / 10
        else:
            pass

    protein_box = [x_dim, y_dim, z_dim]

    multi_flags = {}
    for option in all_flags:
        if option[0] == 'multiProt_options':
            for t in option[1:]:
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                multi_flags[t[0]] = t[1]

    multi_prot = cg_protein[:-4] + '_copies.gro'

    run_genconf = 'gmx genconf -f ' + box_protein + ' -o ' + multi_prot + ' -nbox ' + multi_flags['-x_num'] + ' ' + \
                  multi_flags['-y_num'] + ' ' + ' 1  -dist ' + multi_flags['-d']
    os.system(run_genconf)

    total_prot = int(multi_flags['-x_num']) * int(multi_flags['-y_num'])

    return multi_prot, total_prot


def runInsane(all_flags, clean_pdb, multi_prot):
    """
    This function runs insane to make a coarse-grained membrane-protein system
    :param all_flags: insane parameters read from 'config.txt'
    :param clean_pdb: protein file name
    :param multi_prot: coarse-grained multiple protein file
    :param protein_box: protein box dimension
    :return: coarse-grained membrane-protein system, topology file, name & number of all lipids
    """
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

    with open(multi_prot) as fin:
        lines = fin.readlines()
    new_box_size = lines[-1].split()

    for i in range(0, len(insane_dim)):
        if float(new_box_size[i]) > insane_dim[i]:
            print "Protein is bigger than defined box dimensions. Using protein dimensions."
            insane_dim[i] = new_box_size[i]
    print "Box_dimensions: " + str(insane_dim)

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

    insane_flags = insane_flags + ' -f ' + multi_prot + ' -o ' + system + ' -p ' + system_top
    run_insane = 'python insane.py ' + insane_flags
    os.system(run_insane)
    return True, system, system_top, lipids


def make_ndx(system, lipids):
    """
    this function makes index file: protein_lipids & sol_ions
    :param system: membrane-protein system
    :param lipids: names & no:of lipids
    :return: makes a gromacs index file (.ndx)
    """
    ndx_flags = ''
    for lipid in lipids:
        lipid = lipid.split(':')
        ndx_flags = ndx_flags + ' ' + lipid[0]
    system_ndx = system[:-4] + '.ndx'
    index = 'gmx make_ndx -f ' + system + ' -o ' + system_ndx + ' <<EOF\n 1 | r ' + ndx_flags + ' \n 11  &  !  r ' + ndx_flags + ' \n  q \n EOF'
    os.system(index)
    return system_ndx


def make_topology(clean_pdb, cg_topol, system_top, total_prot):
    """
    This function combines the topologies from martinize and insane and makes a topology that can be used for energy minimization
    :param clean_pdb:  protein file name
    :param cg_topol: topology from martinize
    :param system_top: topology from insane
    :param total_prot: total number of protein copies
    :return: Final topology file with correct itp files and number of molecules
    """
    first_lines = '#include "martini_v2.2.itp" \n#include "martini_v2.0_lipids.itp" \n#include "martini_v2.0_ions.itp"\n'

    with open(cg_topol) as fin:
        cg_lines = fin.readlines()
    new_cg_lines = []

    for i in range(0, len(cg_lines)):
        if '"martini.itp"' in cg_lines[i]:
            pass
        elif 'number' in cg_lines[i]:
            new_cg_lines.append(cg_lines[i])
            for j in range(i + 1, len(cg_lines)):
                if 'Protein_' in cg_lines[j]:
                    prot_line = cg_lines[j].split()
                    prot_line[-1] = total_prot
                    prot_line = prot_line[0] + "   " + str(prot_line[1]) + '\n'
                    new_cg_lines.append(prot_line)
            break
        else:
            new_cg_lines.append(cg_lines[i])

    with open(system_top) as fin:
        sys_lines = fin.readlines()
    new_sys_lines = []
    for i in range(0, len(sys_lines)):
        if 'molecules' in sys_lines[i]:
            start = i + 3
            for j in range(start, len(sys_lines)):
                new_sys_lines.append(sys_lines[j])
            break

    topology = clean_pdb[:-4] + '.top'
    with open(topology, 'w') as fout:
        fout.write(first_lines)
        for line in new_cg_lines:
            fout.write(line)
        for line in new_sys_lines:
            fout.write(line)

    return topology


def runMinimization(all_flags, system, topology, system_ndx):
    """
    This function does gromacs energy minimization
    :param all_flags: energy minimization mdp parameters read from 'config.txt'
    :param system: membrane-protein system
    :param topology: final topology
    :param system_ndx: gromacs index file
    :return: energy minimized gro file for next step
    """
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

    run_em_grompp = 'gmx grompp -f ' + em_mdp + ' -c ' + system + ' -p ' + topology + ' -n ' + system_ndx + ' -o ' + em_tpr + ' -maxwarn 1'
    os.system(run_em_grompp)
    run_em_mdrun = 'gmx mdrun -ntmpi 1 -ntomp 4  -s ' + em_tpr + ' -v -deffnm ' + system[:-4] + '_EM'
    # os.system(run_em_mdrun)

    return em_gro, em_tpr


def runEquilibration(all_flags, system, em_gro, topology, system_ndx):
    """
    This function does gromacs equilibration
    :param all_flags:  equilibration mdp parameters read from 'config.txt'
    :param system: membrane-protein system
    :param em_gro: energy minimized system
    :param topology: final topology
    :param system_ndx: gromacs index file
    :return: equilibrated gro file
    """
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
    # os.system(run_equil_grompp)
    run_equil_mdrun = 'gmx mdrun -ntmpi 1 -ntomp 4  -s ' + equil_tpr + ' -v -deffnm ' + system[:-4] + '_EQUIL'
    # os.system(run_equil_mdrun)

    return equil_gro, equil_tpr
