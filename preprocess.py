import os
import os.path
import sys


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


def fix_ILE(in_pdb, out_pdb):
    """
    This function renames the carbon CD1 to CD in ILE aminoacid for gromos compatible
    :param in_pdb: takes the pdb2pqr processed pdb file
    :param out_pdb: CD fixed protein file
    :return:
    """
    with open(in_pdb) as in_f, open(out_pdb, 'w') as out_f:
        for line in in_f:
            if line.startswith('ATOM  '):
                residue_name = line[17:21].strip()
                atom_name = line[12:16].strip()
                if residue_name == 'ILE' and atom_name == 'CD1':
                    line = line[:12] + ' CD ' + line[16:]
            out_f.write(line)
    return True


def gromos_residues(in_pdb, out_pdb):
    """
    This function fixes the atom names for amino acids
    :param in_pdb: CD fixed protein file
    :param out_pdb: protein with atom names fixed
    :return:
    """
    resis = dict(ALA=['N', 'H', 'CA', 'CB', 'C', 'O'],
                 ARG=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22',
                      'C', 'O'],
                 ARN=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'C',
                      'O'],
                 ASN=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O'],
                 ASP=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'C', 'O'],
                 ASH=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O'],
                 CYS=['N', 'H', 'CA', 'CB', 'SG', 'C', 'O'],
                 GLN=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O'],
                 GLU=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'C', 'O'],
                 GLH=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O'],
                 GLY=['N', 'H', 'CA', 'C', 'O'],
                 HIS=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'HD1', 'CD2', 'NE2', 'CE1', 'C', 'O'],
                 HID=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'HD1', 'CD2', 'NE2', 'CE1', 'C', 'O'],
                 HIP=['N', 'H', 'CA', 'CB', 'CG', 'CD2', 'ND1', 'HD1', 'CE1', 'NE2', 'HE2', 'C', 'O'],
                 HIE=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'CE1', 'CD2', 'NE2', 'HE2', 'C', 'O'],
                 ILE=['N', 'H', 'CA', 'CB', 'CG2', 'CG1', 'CD', 'C', 'O'],
                 LEU=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                 LYS=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'],
                 LYN=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'HZ1', 'HZ2', 'C', 'O'],
                 MET=['N', 'H', 'CA', 'CB', 'CG', 'SD', 'CE', 'C', 'O'],
                 MSE=['N', 'H', 'CA', 'CB', 'CG', 'SE', 'CE', 'C', 'O'],
                 PHE=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'C', 'O'],
                 PRO=['N', 'CD', 'CA', 'CB', 'CG', 'C', 'O'],
                 SER=['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O'],
                 THR=['N', 'H', 'CA', 'CB', 'OG1', 'HG1', 'CG2', 'C', 'O'],
                 TRP=['N', 'H', 'CA', 'CB', 'CG', 'CD2', 'CE2', 'CE3', 'CD1', 'NE1', 'HE1', 'CZ2', 'CZ3', 'CH2', 'C',
                      'O'],
                 TYR=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'HH', 'C', 'O'],
                 VAL=['N', 'H', 'CA', 'CB', 'CG1', 'CG2', 'C', 'O']
                 )

    atoms = []
    with open(in_pdb) as pdb:
        for line in pdb:
            if line.startswith('ATOM  '):
                chain_id = line[21]
                residue_num = int(line[22:26])
                residue_alt = line[26]
                residue_id = str(residue_num) + residue_alt
                residue_name = line[17:21].strip()
                atom_name = line[12:16].strip()
                atoms.append({'chain_id': chain_id, 'residue_num': residue_num, 'residue_alt': residue_alt,
                              'residue_id': residue_id, 'residue_name': residue_name, 'atom_name': atom_name,
                              'line': line})

    pdb_residue_list = set()
    for atom in atoms:
        pdb_residue_list.add((atom['chain_id'], atom['residue_num'], atom['residue_alt']))

    pdb_residue_list = sorted(pdb_residue_list, key=lambda x: (ord(x[0]) * 100000 * 100) + x[1] * 100 + ord(x[2]))

    pdb_residues = []
    for residue in pdb_residue_list:
        residue_id = str(residue[1]) + residue[2]
        pdb_residues.append(filter(lambda x: x['chain_id'] == residue[0] and x['residue_id'] == residue_id, atoms))

    atom_i = 1
    with open(out_pdb, 'w') as out_f:
        for residue in pdb_residues:
            residue_name = residue[0]['residue_name']
            if residue_name not in resis:
                print >> sys.stderr, "Residue %s not supported, stripped for mapping" % residue_name
                continue
            residue_atoms = resis[residue_name]
            for resi_atom in residue_atoms:
                for pdb_atom in residue:
                    if pdb_atom['atom_name'] == resi_atom:
                        atom = "%5g" % atom_i
                        out_f.write(pdb_atom['line'][0:6] + atom + pdb_atom['line'][11:])
                        atom_i += 1
                        break

    return True


def pdbpqr(base_dir, pdb):
    """
    This function runs pdb2pqr script to process pdb file.
    First fixes the protein backbone, then process pdb using pdbpqr,
    next changes the CD1 to CD in ILE,
    then renames atom names for gromos compatiblity
    :param base_dir: present working directory
    :param pdb: oriented protein file as input from 'config.txt'
    :return: processed clean protein pdb file
    """
    pdb_prefix = pdb[:-11]

    genout = pdb_prefix + "_pgen.pdb"
    pqr = pdb_prefix + "_pmin.pqr"
    gromos_pdb = pdb_prefix + "_clean.pdb"
    outmol = pdb_prefix + "_pmin.mol.pdb"
    outmol2 = pdb_prefix + "_pmin.mol2.pdb"

    fix_missing_BB(pdb, genout)
    #  with open('error_log', 'a') as err_out, open('output_log', 'a') as out_f:
    run_pdb2pqr = base_dir + '/pdb2pqr/pdb2pqr.py --nodebump --noopt --mol_charmm_pdb --chain --ff=GROMOS --ffout=GROMOS ' + genout + ' ' + pqr
    # subprocess.call([base_dir + '/pdb2pqr/pdb2pqr.py', '--nodebump', '--noopt', '--mol_charmm_pdb', '--chain', '--ff=GROMOS', '--ffout=GROMOS', genout, pqr], stdout=out_f, stderr=err_out)
    os.system(run_pdb2pqr)
    fix_ILE(outmol, outmol2)
    gromos_residues(outmol2, gromos_pdb)

    os.unlink(genout)
    os.unlink(pqr)
    os.unlink(outmol)
    os.unlink(outmol2)
    #  os.unlink('error_log')
    # os.unlink('output_log')
    return gromos_pdb


def runDSSP(clean_pdb):
    """
    This function assigns secondary structure for the protein
    :param clean_pdb: clean protein from pdb2pqr processing
    :return: dssp file for the protein
    """
    do_dssp = './dssp  '
    infile = clean_pdb
    out_dssp = infile[:-4] + ".dssp"
    run_cmd = do_dssp + infile + " " + out_dssp
    os.system(run_cmd)
    return out_dssp
