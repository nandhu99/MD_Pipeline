import os
import os.path


def add_chain(file_in, file_out, nos_chain):
    f = open(file_in, 'r')
    lines = f.readlines()
    f.close()

    chain_id = 65  # this is the integer that corresponds to the character "A"
    max_chain = int(nos_chain) + chain_id
    pdb_chains = []

    for i in range(len(lines)):
        line = lines[i]
        if line.startswith('TER'):
            chain_id = chain_id + 1

        if chain_id < max_chain:
            chain_char = chr(chain_id)
        else:
            chain_char = ' '

        if line.startswith('ATOM'):
            line = line[0:21] + chain_char + line[22:]
            pdb_chains.append(line)

        elif line.startswith('TER') and len(line) > 21:
            line = line[0:21] + ' ' + line[22:]
            pdb_chains.append(line)

        else:
            pdb_chains.append(line)

    f = open(file_out, 'w')
    for i in range(len(pdb_chains)):
        f.write(pdb_chains[i])
    f.close()
    return file_out


def fix_missing_BB(pdb_in, pdb_out):
    """
    This function fixes the protein for missing backbone atoms
    :param pdb_in: oriented protein pdb file
    :param pdb_out: protein with all backbone atoms fixed
    :return: True
    """
    f = open(pdb_in, 'r')
    lines = f.readlines()
    f.close()

    atom_i = 0
    atoms = []
    for line in lines:
        if line.startswith('ATOM  '):
            atom_i += 1
            chain_id = line[21]
            residue_id = line[22:27].strip()
            residue_name = line[17:21].strip()
            atom_name = line[12:16].strip()
            atoms.append(
                {'atom_i': atom_i, 'chain_id': chain_id, 'residue_id': residue_id, 'residue_name': residue_name,
                 'atom_name': atom_name, 'line': line})

    pdb_residue_list = set()
    for atom in atoms:
        pdb_residue_list.add((atom['chain_id'], atom['residue_name'], atom['residue_id']))

    pdb_residues = []
    for residue in pdb_residue_list:
        pdb_residues.append(filter(
            lambda x: x['chain_id'] == residue[0] and x['residue_id'] == residue[2] and x['residue_name'] == residue[1],
            atoms))

    good_atoms = []
    for residue in pdb_residues:
        has_CA = False
        has_C = False
        has_N = False
        for atom in residue:
            if atom['atom_name'] == 'N':
                has_N = True
            if atom['atom_name'] == 'C':
                has_C = True
            if atom['atom_name'] == 'CA':
                has_CA = True
        if has_N and has_C and has_CA:
            for atom in residue:
                good_atoms.append(atom)

    good_atoms.sort(key=lambda x: x['atom_i'])

    f = open(pdb_out, 'w')
    for atom in good_atoms:
        f.write(atom['line'])
    f.close()
    return True


def pdbpqr(base_dir, pdb, chain):
    """
    This function runs pdb2pqr script to process pdb file.
    First fixes the protein backbone, then process pdb using pdbpqr,
    next changes the CD1 to CD in ILE,
    then renames atom names for gromos compatiblity
    :param base_dir: present working directory
    :param pdb: oriented protein file as input from 'config.txt'
    :return: processed clean protein pdb file
    """
    pdb_prefix = pdb[:-4]

    chain_out = pdb_prefix + "_chains.pdb"
    bb_out = pdb_prefix + "_bb.pdb"
    out_pdb = pdb_prefix + "_pqr.pdb"

    add_chain(pdb, chain_out, chain)

    fix_missing_BB(chain_out, bb_out)
    run_pdb2pqr = base_dir + '/pdb2pqr/pdb2pqr.py --nodebump --noopt --chain --ff=PARSE --ffout=PARSE -v ' + bb_out + ' ' + out_pdb
    os.system(run_pdb2pqr)
    #os.unlink(chain_out)
    #os.unlink(bb_out)

    update_line = []
    with open(out_pdb, 'r') as fin:
        for line in fin:
            if line[0:4] == 'ATOM':
                occup = "%6.2f" % (float(1)) + "%6.2f" % (float(2))
                newline = line[0:54] + occup + line[70:] + "\n"
                update_line.append(newline)
            else:
                update_line.append(line)

    clean_pdb = pdb_prefix + "_clean.pdb"

    with open(clean_pdb, 'w') as fout:
        for i in range(len(update_line)):
            fout.write(update_line[i])


    print "Written " + clean_pdb

    #os.unlink(out_pdb)
    return clean_pdb


def sep_chains(pdb_in, processed_pdb):
    f1 = open(pdb_in, 'r')
    lines = f1.readlines()
    f1.close()

    clean_pdb = []

    for line in lines:
        if line[0:4] != 'ATOM':
            pass
        elif line[0:4] == 'ATOM':
            clean_pdb.append(line)

    new_pdb = []
    for i in range(len(clean_pdb) - 1):
        if clean_pdb[i][21] == clean_pdb[i + 1][21]:
            new_pdb.append(clean_pdb[i])
        else:
            new_line = "TER " + "\n"
            new_pdb.append(new_line)

    f2 = open(processed_pdb, 'w')
    for i in range(len(new_pdb)):
        f2.write(new_pdb[i])
    f2.write("END")
    f2.close()
    print "Written " + processed_pdb
    return processed_pdb


def runDSSP(processed_pdb):
    """
    This function assigns secondary structure for the protein
    :param clean_pdb: clean protein from pdb2pqr processing
    :return: dssp file for the protein
    """
    do_dssp = './dssp  '
    infile = processed_pdb
    out_dssp = infile[:-9] + ".dssp"
    run_cmd = do_dssp + infile + " " + out_dssp
    os.system(run_cmd)
    print "Written " + out_dssp
    return out_dssp
