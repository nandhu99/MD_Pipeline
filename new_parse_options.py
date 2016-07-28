def parseconfig(filename):
    with open(filename) as fin:
        lines = fin.readlines()

    preprocess_options = {}
    martinize_options = {}
    insane_options = {}
    simulation_options = {}

    for line in lines:
        line = line.rstrip()
        if '#' in line or ';' in line:
            pass
        elif len(line) == 0:
            pass
        elif 'preprocess_options' == line:
            print line
            if '#' in line or ';' in line:
                pass
            elif len(line) == 0:
                pass
            else:
                print 'preprocess'
                t = line.split('=')
                t[0] = t[0].strip()
                t[1] = t[1].strip()
                if t[1] == 'False' or t[1] == 'false':
                    pass
                elif t[1] == 'True' or t[1] == 'true':
                    preprocess_options[t[0]] = ''
                    print preprocess_options[t[0]]
                else:
                    preprocess_options[t[0]] = t[1]
                    print preprocess_options[t[0]]
        elif 'end_preprocess_options' in line:
            print 'end_preprocess'

    preprocess_keys = preprocess_options.keys()
    print preprocess_keys

    preprocess_flags = ''
    for key in preprocess_keys:
        if key in preprocess_options:
            preprocess_flags = preprocess_flags + ' ' + key + ' ' + preprocess_options[key]
        else:
            pass

    return preprocess_flags
