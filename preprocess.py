import os
import sys


def clean_pdb(filename):
    with open(filename) as fin:
        lines = fin.readlines()

    outfile = 'clean_' + filename
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
    return True


def runDSSP(filename):
    # Generate DSSP output files
    do_dssp = './dssp_exe '
    infile = filename
    out_dssp = infile + ".dssp"
    run_cmd = do_dssp + infile + " " + out_dssp
    os.system(run_cmd)
    return True


def main():
    clean_pdb(sys.argv[1])
    runDSSP(sys.argv[1])
    return True


if __name__ == "__main__":
    main()
