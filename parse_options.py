def parseConfig2(filename):
    with open(filename) as fin:
        lines = fin.readlines()
    options = []
    no_of_lines = len(lines)
    i = 0
    current_option = ''
    while True:
        line = lines[i]
        tline = str(line[1:-1])
        line = line.strip()
        if '_option' in tline and 'end' not in tline:
            t_type = []
            current_option = tline
            t_type.append(tline)

        if len(line) > 0:
            if '#' in line or ';' in line:
                pass
            else:
                t = line.split('=')
                t_type.append(t)

        if 'end' in tline:
            if tline == 'end_' + current_option:
                options.append(t_type)
                t_type = []
                current_option = ''
        i += 1
        if i == no_of_lines:
            break
    return options
