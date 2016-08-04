def box_dimension(filename):
    import numpy as np
    with open(filename) as fin:
        lines = fin.readlines()
        x_dim = []
        y_dim = []
        z_dim = []

        outfile = filename[:-4] + '_center.pdb'
        fout = open(outfile, 'w')

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
        print box_line

        Cx = round(np.mean(x_dim), 3)
        Cy = round(np.mean(y_dim), 3)
        Cz = round(np.mean(z_dim), 3)

        for line in lines:
            if line[0:4] == 'ATOM':
                new_x = float(line[30:38]) - Cx
                new_y = float(line[38:46]) - Cy
                new_z = float(line[46:54]) - Cz

                new_coors = "%8.3f" % (new_x) + "%8.3f" % (new_y) + "%8.3f" % (new_z)
                new_outline = line[0:30] + new_coors + line[54:]
                fout.write(new_outline)
    return outfile


import sys

box_dimension(sys.argv[1])
