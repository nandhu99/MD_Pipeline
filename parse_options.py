def parseconfig(filename):
    with open(filename) as fin:
        lines = fin.readlines()
    options = {}
    for line in lines:
        line = line.rstrip()
        if '#' in line or ';' in line:
            pass
        elif len(line) == 0:
            pass
        else:
            t = line.split('=')
            t[0] = t[0].strip()
            t[1] = t[1].strip()
            if t[1] == 'False' or t[1] == 'false':
                pass
            elif t[1] == 'True' or t[1] == 'true':
                options[t[0]] = ''
            else:
                options[t[0]] = t[1]

    preprocess = ['-protein']
    option_keys = options.keys()

    martinize_flags = ''
    for key in option_keys:
        if key in options:
            pass
        else:
            martinize_flags = martinize_flags + ' ' + key + ' ' + options[key]

    preprocess_flags = ''
    for key in option_keys:
        if key in options:
            preprocess_flags = preprocess_flags + ' ' + key + ' ' + options[key]
        else:
            pass

    return preprocess_flags


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
