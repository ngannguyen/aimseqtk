
import os
import sys

def read_clonesize(file):
    s2clones = {}
    f = open(file, 'r')
    f.readline()
    for line in f:
        items = line.strip().split('\t')
        sample = items[0]
        clones = int(float(items[1]))
        s2clones[sample] = clones
    f.close()
    return s2clones

def get_np_pc(p_s2clones, np_s2clones, outfile):
    f = open(outfile, 'w')
    f.write("#Sample\t%%productive\t%%non_productive\n")
    for s, p in p_s2clones.iteritems():
        np = 0
        if s in np_s2clones:
            np = np_s2clones[s]
        total = p + np
        if total > 0:
            f.write("%s\t%f\t%f\n" % (s, 100.0*p/total, 100.0*np/total))
    f.close()

def main():
    pfile = sys.argv[1]
    npfile = sys.argv[2]
    outfile = sys.argv[3]
    p_s2clones = read_clonesize(pfile)
    np_s2clones = read_clonesize(npfile)
    get_np_pc(p_s2clones, np_s2clones, outfile)

if __name__ == '__main__':
    main()
