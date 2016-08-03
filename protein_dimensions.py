def box_dimension(filename):
    with open(filename) as fin:
        lines = fin.readlines()
        x_dim = []
        y_dim = []
        z_dim = []

        for line in lines:
            x_dim.append(float(line[30:38]))
            y_dim.append(float(line[38:46]))
            z_dim.append(float(line[46:54]))

        x_dim.sort()
        y_dim.sort()
        z_dim.sort()
    return [float(abs(x_dim[0] - x_dim[-1]) / 10), float(abs(y_dim[0] - y_dim[-1]) / 10),
            float(abs(z_dim[0] - z_dim[-1]) / 10)]


