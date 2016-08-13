import os
import os.path
import subprocess
import sys


def fix_missing_BB(pdb_in, pdb_out):
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


def fix_ILE(in_pdb, out_pdb):
    with open(in_pdb) as in_f, open(out_pdb, 'w') as out_f:
        for line in in_f:
            if line.startswith('ATOM  '):
                residue_name = line[17:21].strip()
                atom_name = line[12:16].strip()
                if residue_name == 'ILE' and atom_name == 'CD1':
                    line = line[:12] + ' CD ' + line[16:]
            out_f.write(line)


def gromos_norm(in_pdb, out_pdb):
    resis = dict(ALA=['N', 'H', 'CA', 'CB', 'C', 'O'],
                 ARG=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'HH22', 'C', 'O'],
                 ARGN=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'NE', 'HE', 'CZ', 'NH1', 'HH11', 'HH12', 'NH2', 'HH21', 'C', 'O'],
                 ASN=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'ND2', 'HD21', 'HD22', 'C', 'O'],
                 ASP=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'C', 'O'],
                 ASPH=['N', 'H', 'CA', 'CB', 'CG', 'OD1', 'OD2', 'HD2', 'C', 'O'],
                 CYS=['N', 'H', 'CA', 'CB', 'SG', 'C', 'O'],
                 CYSH=['N', 'H', 'CA', 'CB', 'SG', 'HG', 'C', 'O'],
                 GLN=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'NE2', 'HE21', 'HE22', 'C', 'O'],
                 GLU=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'C', 'O'],
                 GLUH=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'HE2', 'C', 'O'],
                 GLY=['N', 'H', 'CA', 'C', 'O'],
                 HIS=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'C', 'O'],
                 HISA=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'C', 'O'],
                 HISB=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2', 'C', 'O'],
                 HISH=['N', 'H', 'CA', 'CB', 'CG', 'ND1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'NE2', 'HE2', 'C', 'O'],
                 ILE=['N', 'H', 'CA', 'CB', 'CG1', 'CG2', 'CD', 'C', 'O'],
                 LEU=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'CD2', 'C', 'O'],
                 LYSH=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'HZ1', 'HZ2', 'HZ3', 'C', 'O'],
                 LYS=['N', 'H', 'CA', 'CB', 'CG', 'CD', 'CE', 'NZ', 'HZ1', 'HZ2', 'C', 'O'],
                 MET=['N', 'H', 'CA', 'CB', 'CG', 'SD', 'CE', 'C', 'O'],
                 PHE=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'HZ', 'C', 'O'],
                 PRO=['N', 'CA', 'CB', 'CG', 'CD', 'C', 'O'],
                 SER=['N', 'H', 'CA', 'CB', 'OG', 'HG', 'C', 'O'],
                 THR=['N', 'H', 'CA', 'CB', 'OG1', 'HG1', 'CG2', 'C', 'O'],
                 TRP=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'NE1', 'HE1', 'CE2', 'CE3', 'HE3', 'CZ2', 'HZ2', 'CZ3', 'HZ3', 'CH2', 'HH2', 'C', 'O'],
                 TYR=['N', 'H', 'CA', 'CB', 'CG', 'CD1', 'HD1', 'CD2', 'HD2', 'CE1', 'HE1', 'CE2', 'HE2', 'CZ', 'OH', 'HH', 'C', 'O'],
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


def pdbpqr(base_dir, pdb):
    pdb_prefix = pdb[:-4]

    genout = pdb_prefix + "_pgen.pdb"
    pqr = pdb_prefix + "_pmin.pqr"
    nonminout = pdb_prefix + "_clean.pdb"
    outmol = pdb_prefix + "_pmin.mol.pdb"
    outmol2 = pdb_prefix + "_pmin.mol2.pdb"

    fix_missing_BB(pdb, genout)
    with open('error_log', 'a') as err_out, open('output_log', 'a') as out_f:
        subprocess.call([base_dir + '/pdb2pqr/pdb2pqr.py', '--nodebump', '--noopt', '--mol_charmm_pdb', '--chain',
                         '--ff=LIBMOL', '--ffout=LIBMOL', genout, pqr], stdout=out_f, stderr=err_out)
    fix_ILE(outmol, outmol2)
    gromos_norm(outmol2, nonminout)

    os.unlink(genout)
    os.unlink(pqr)
    os.unlink(outmol)
    os.unlink(outmol2)
    os.unlink('error_log')
    os.unlink('output_log')
    return nonminout


def runDSSP(filename):
    do_dssp = './dssp_exe '
    infile = filename
    out_dssp = infile[:-4] + ".dssp"
    run_cmd = do_dssp + infile + " " + out_dssp
    os.system(run_cmd)
    return out_dssp
