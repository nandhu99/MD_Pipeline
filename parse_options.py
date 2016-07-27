import sys, string


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

    #print options

    exclude = ['-protein', '-clean_protein', '-dssp']
    option_keys = options.keys()
    martinize_flags = ''

    for key in option_keys:
        if key in exclude:
            pass
        else:
            martinize_flags = martinize_flags+' '+key+' '+options[key]

    clean_flags = ''

    return martinize_flags, clean_flags

