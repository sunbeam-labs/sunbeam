"""
varsub - variable substitution of parsed yaml objects

from varsub import varsub
import yaml

config = yaml.load("config.yaml")

or in snakefile

configfile: "config.yaml"

varsub(config)

or to suppress error messages

varsub(config, False)

"""
import sys
import re

verbal = True
version = ''
warns = {}

def varsub(ys, v = True):
    global verbal
    global version
    verbal = v
    version = sys.version[0]

    # definition must be the top level 1

    hit = 1
 
    while(hit):
        hit = 0
        objs = [ys]

        while(len(objs)):
            obj = objs.pop()
            t = str(type(obj))
            if isinstance(obj, dict):
                keys = obj.keys()

                # print(keys)

                for key in keys:
                    # print(key)
                    # print(obj[key])
                    #print(key + ': ' + obj[key])
                    t = str(type(obj[key]))
                    # print(t)
                    if isinstance(obj[key], str):
                        #print("str")
                        if __scan4vars(ys, obj, key):
                            hit = 1
                    elif isinstance(obj[key], list):
                        # print("list")
                        objs.append(obj[key])
                    elif isinstance(obj[key], dict):
                        # print("dict")
                        objs.append(obj[key])
                    #else:  # may be int, float, etc.
                    #    print(t + " unknown")
            else:  # must be list
                for idx in range(len(obj)):
                    el = obj[idx]
                    if isinstance(el, str):
                        if __scan4vars(ys, obj, idx):
                            hit = 1
                    elif isinstance(el, list):
                        objs.append(el)
                    elif isinstance(el, dict):
                        objs.append(el)
                    #else:  # may be int, float, etc.
                    #    print(t + " unknown", file=sys.stderr)
    if verbal and warns:
       for k in warns: 
           print(warns[k], file=sys.stderr)

def __scan4vars(ys, obj, key):
    pat = re.compile('\${?([\w\d]+)}?')
    full = re.compile('\${?([\w\d]+)}?$')
    ps = {}
    hit = 0
    for m in re.findall(pat, obj[key]):
        ps[m] = 1  # handle dulicated matches
    if ps:
        if re.match(full, obj[key]): # not embreded substitution
            k = list(ps.keys())[0]
            if k in ys:
               obj[key] = ys[k]
               hit = 1
            elif verbal:  # warning, quiting?
                warns[k] = k + ' in "%s" not defined' % (obj[key])
        else:
            for k in ps.keys():
                if k in ys:
                    #if not isinstance(ys[k], (list, dict, tuple, complex)):
                    p = '\${?%s}?' % k
                    if isinstance(ys[k], str):
                       obj[key] = re.sub(p, ys[k], obj[key])
                       hit = 1
                    elif isinstance(ys[k], (int, float)):
                       obj[key] = re.sub(p, str(ys[k]), obj[key])
                       hit = 1
                    elif version == 2 and isinstance(ys[k], long):
                       obj[key] = re.sub(p, str(ys[k]), obj[key])
                       hit = 1
                    elif verbal:
                       warns[k] = k + ' is not a valid type for substituion in "%s"' % (obj[key])
                elif verbal:  # warning, quiting?
                    warns[k] = k + ' in "%s" not defined' % (obj[key])
        #print(obj[key])
    return(hit)
