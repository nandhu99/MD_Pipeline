def box_dimension(filename):
    """
    This function centers the protein and then calculates the protein box dimensions
    :param filename: clean pdb file
    :return: protein with a box is returned
    """
    import numpy as np
    with open(filename) as fin:
        lines = fin.readlines()
        x_dim = []
        y_dim = []
        z_dim = []

        for line in lines:
            if line[0:4] == 'ATOM':
                x_dim.append(float(line[30:38]))
                y_dim.append(float(line[38:46]))
                z_dim.append(float(line[46:54]))
            else:
                pass

    x_dim.sort()
    y_dim.sort()
    z_dim.sort()
    dx = float(abs(x_dim[0] - x_dim[-1]))
    dy = float(abs(y_dim[0] - y_dim[-1]))
    dz = float(abs(z_dim[0] - z_dim[-1]))

    box_line = 'CRYST1' + '%9.3f' % (dx) + '%9.3f' % (dy) + '%9.3f' % (dz) + '  90.00  90.00  90.00 P 1           1'

    Cx = round(np.mean(x_dim), 3)
    Cy = round(np.mean(y_dim), 3)
    Cz = round(np.mean(z_dim), 3)

    box_protein = filename[:-4] + '_center.pdb'
    with open(box_protein, 'w') as fout:
        fout.write(box_line)
        fout.write('\n')
        for line in lines:
            if line[0:4] == 'ATOM':
                new_x = float(line[30:38])
                new_y = float(line[38:46])
                new_z = float(line[46:54])

                new_coors = "%8.3f" % (new_x) + "%8.3f" % (new_y) + "%8.3f" % (new_z)
                new_outline = line[0:30] + new_coors + line[54:]
                fout.write(new_outline)
        fout.write('END')

    return box_protein
