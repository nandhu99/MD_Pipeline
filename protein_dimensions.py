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

        Cx = round(np.mean(x_dim), 3)
        Cy = round(np.mean(y_dim), 3)
        Cz = round(np.mean(z_dim), 3)

        new_x = []
        new_y = []
        new_z = []
        for i in range(0, len(x_dim)):
            x = x_dim[i] - Cx
            new_x.append(float(round(x,3)))
            y = y_dim[i] - Cy
            new_y.append(float(y))
            z = z_dim[i] - Cz
            new_z.append(float(z))

        for line in lines:







'''
        x_dim.sort()
        y_dim.sort()
        z_dim.sort()

        dx = float(abs(x_dim[0] - x_dim[-1]))
        dy = float(abs(y_dim[0] - y_dim[-1]))
        dz = float(abs(z_dim[0] - z_dim[-1]))
'''

    return True


import sys

box_dimension(sys.argv[1])
