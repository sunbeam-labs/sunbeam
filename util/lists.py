
class onechromosome:
    def __init__(self, dlocs):
        self.dlocs = dlocs
    def cmp(a,b):
        if self.dlocs[a][0] == self.dlocs[b][0]:
           if self.dlocs[a][1] > self.dlocs[b][1]:
               return 1
           else:
               return -1
        else:
           if self.dlocs[a][0] > self.dlocs[b][0]:
               return 1
           else:
               return -1

def combine(inputf,outputf):
    "combine indels in the individual lists into a single list"
    indels = dict()
    dlocs = dict()
    for name in inputf:
        fin = open(name, "r")
        for line in fin.readlines():
            line = re.sub('\n', '', line)
            chr, loc = line.split(':')
            if chr in indels:
                if loc not in indels[chr]:
                    locs = loc.split('-')
                    if len(locs) == 1:
                        locs.append(locs[0])
                    indels[chr][loc] = [int(locs[0]), int(locs[1])]
            else:
                indels[chr] = dict()
                locs = loc.split('-')
                if len(locs) == 1:
                    locs.append(locs[0])
                indels[chr][loc] = [int(locs[0]), int(locs[1])]
                    
    # how to sniff version?
    chrs = list(indels.keys())  # works for both version 2 and 3
    chrs.sort()
    with open(outputf, "w") as out:
        for chr in chrs:
            dlocs = indels[chr]
            onechr = onechromosome(dlocs)
            ks = list(dlocs.keys())
            # ks.sort(key=lamda x:dlocs[x[0])
            ks.sort(key=cmp_to_key(onechr.cmp))
            for loc in ks:
                out.write(chr + ':' + loc + "\n")
