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
        if key in preprocess:
            pass
        else:
            martinize_flags = martinize_flags + ' ' + key + ' ' + options[key]

    preprocess_flags = ''
    for key in option_keys:
        if key in preprocess:
            preprocess_flags = preprocess_flags + ' ' + key + ' ' + options[key]
        else:
            pass

    return martinize_flags, preprocess_flags
