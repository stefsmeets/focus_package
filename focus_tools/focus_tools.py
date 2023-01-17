import os, sys

from operator import itemgetter
from itertools import groupby


drc = os.path.abspath(os.path.dirname( __file__ )) # get path of program


def move(fname, target):
    if os.path.exists(target):
        os.remove(target)
    os.rename(fname, target)


def parse_coseq(f, neighbours=100):
    if isinstance(f, (str, str)):
        f = open(f)
    
    prev = ""
    db  = {}
    rdb = {}

    for line in f:  
        inp = line.split()
        if not inp:
            continue

        name  = inp[0]
        coseq = tuple(map(int,inp[1:neighbours+1]))

        if coseq in rdb:
            if name not in rdb[coseq]:
                rdb[coseq].append(name)
            else:
                continue
        else:
            rdb[coseq] = [name]


        if name == prev:
            db[name].append(coseq)
        else:
            db[name] = [coseq]
            prev = name

    return db,rdb

coseqdb, coseqrdb = parse_coseq(os.path.join(drc, "..", "resources", "coseq"))   # coseq length = 10


def uniquify_list(lst):
    return list(set(lst))

def uniquify_invert_dict(d):
    uniq_map = {}
    counts = {}
    keys = list(d.keys())
    keys.sort()
    for key in keys:
        values = tuple(d[key])
        if values in uniq_map:
            pass
        else:
            uniq_map[values] = key
    
    return uniq_map

def invert_dict(d):
    inv_map = {}
    for key, values in d.items():
        for val in values:
            inv_map[val] = inv_map.get(val, [])
            inv_map[val].append(key)
    return inv_map

def section(f, intro="F", key=None):
    if isinstance(f, (str, str)):
        f = open(f)
    
    start = f">Begin {key}"
    end   = f">End {key}"
    
    out = False
    num = intro+"0000"

    for line in f:
        if line.startswith(start):
            out = True
            num = intro+line.split(":")[-1].strip()
            continue
    
        if line.startswith(end):
            out = False
            continue       
            
        if out:
            yield num+" "+line.strip()


def cut(line, col, label=None):
    inp = line.split()
    inp.pop(col)
    if label:
        return inp[0] + f"-{label} " + " ".join(inp[1:])
    else:
        return " ".join(inp)

def coseq_reduce(lines):
    last_key = ""
    seq = []
    d = {}
    key = None
    for line in lines:
        inp = line.split()
        key = inp.pop(0)
        if key != last_key:
            seq = uniquify_list(seq)
            seq.sort()
            d[last_key] = seq
            seq = []
            last_key = key
        try:
            inp = tuple([int(item) for item in inp])
        except ValueError:
            continue
        
        assert len(inp) == 10
        seq.append(inp)
    # remove empty key
    if "" in d:
        del d[""]
    # ensure last FW is added too
    if key:
        seq = uniquify_list(seq)
        seq.sort()
        d[key] = seq
    return d

def get_coseq(fn):
    if isinstance(fn, list) and len(fn) == 1:
        fn = fn[0]

    if isinstance(fn, list):
        lines = []
        for i, f in enumerate(fn):
            s = section(f, key="coseq")
            l = [cut(line, 1, label=str(i)) for line in s]
            lines.extend(l)
    else:
        lines = section(fn, key="coseq")
        lines = [cut(line, 1) for line in lines]
    
    coseq = coseq_reduce(lines)
    return coseq

def coseq_cmp(d, coseqdb):
    di = uniquify_invert_dict(d)
    rdb = uniquify_invert_dict(coseqdb)

    dcounts = {}
    dmap = {}
    for k,v in list(d.items()):
        v = tuple(v)
        try:
            fw = di[v]
        except KeyError:
            print("KeyError:", v)
            fw = None
        else:
            try:
                dcounts[fw] += 1
            except KeyError:
                dcounts[fw] = 1
                
    for key in list(dcounts.keys()):
        v = tuple(d[key])
        try:
            ftc = rdb[v]
        except KeyError:
            ftc = None
        else:
            dmap[key] = ftc
    
    return dcounts, dmap

def fo2hist(f):
    d = get_coseq(f)

    dcounts, dmap = coseq_cmp(d, coseqdb)
    items = sorted(list(dcounts.items()), key=itemgetter(1), reverse=True)
    total = sum(dcounts.values())
    
    if isinstance(f, list):
        print(f"Summary of {len(f)} files")
    else:
        print(f)
    print("fw              #      %  Natoms  ftc")
    for key, count in items:
        print("{:6} {:10d} {:6.3f} {:7d}  {}".format(key, count, float(count)/total, len(d[key]), dmap.get(key, "-")))
    print(f"Total: {total}")
    print()

def fo2strudat(f, fout=None, unique_only=True):
    lines = section(f, key="Framework")
    if unique_only:
        d = get_coseq(f)
        dcounts, dmap = coseq_cmp(d, coseqdb)
        keys = list(dcounts.keys())

        for line in lines:
            try:
                key, title = line.split(None, 1)
            except ValueError:    # if no Title is given in foc.inp
                key = line.strip()
                title = ""
            if key in keys:
                print(title, file=fout)  
    else:
        for line in lines:
            print(line.split(None, 1)[0], file=fout)

def get_R_from_dls76_out():
    f = open("dls76.out")

    for line in f:
        if line.startswith("0R=SQRT( SUM(W*(DO-D))**2 / SUM(W*DO)**2 )="):
            rval = line[44:57]
    return float(rval)


def nfilea2cif(nfilea="nfilea.inp", cif="structure.cif", out=None):
    if isinstance(nfilea, (str, str)):
        nfilea = open("nfilea.inp")
    if isinstance(cif, (str, str)):
        cif = open(cif)
    if out:
        if isinstance(out, (str, str)):
            out = open(out, "w")

    for line in cif:
        if not (line.startswith("O") or line.startswith("T")):
            print(line, file=out)

    for line in nfilea:
        if line.startswith("NATOM"):
            print(f"{line[7:13]}  {line[13:21]}  {line[21:29]}  {line[29:37]}", file=out)


def enable_cdls():
    fin = open("dlsinp")
    fout = open("dlsinp_new", "w")

    for line in fin:
        line = line.strip()
        if line == "DLS-76  -5    30     0         2                1":
            line = "DLS-76  -5    30     1         2                1"

        print(line, file=fout)
    fin.close()
    fout.close()

    move(fout.name, fin.name)


def dlsall(threshold=0.01, refine_cell=False):
    from . import pykriber
    from . import pydls

    keys = pykriber.extract_all_keys_from_strudat()

    print(f"Threshold = {threshold}")
    print()
    print("framework    Rval")
    for key in keys:
        pykriber.strudat2dls(args=["addo"], keys=key)

        if refine_cell:
            enable_cdls()

        pydls.dls76(args=["dlsinp"])
        rval = get_R_from_dls76_out()
        marker  = "**" if rval < threshold else "  \n"
        print(f"{key:10s} {rval:.4f} {marker}", end=' ')

        if rval < threshold:
            pykriber.strudat2cif(args=["addo"], keys=key, verbose=False)

            nfilea2cif(cif=key+".cif", out=key+"_dls.cif")
            print(" >> Wrote file {}".format(key+"_dls.cif"))


def cdlsall(threshold=0.01):
    dlsall(threshold=threshold, refine_cell=True)


def dlsall_entry():
    args = sys.argv[1:]
    if len(args) == 1:
        threshold = float(args[0])
        dlsall(threshold=threshold)
    else:
        dlsall()


def cdlsall_entry():
    args = sys.argv[1:]
    if len(args) == 1:
        threshold = float(args[0])
        cdlsall(threshold=threshold)
    else:
        cdlsall() 


def fo2cif_entry():
    fns = sys.argv[1:]
    if not fns:
        print("No files given. \n \n >> Usage: fo2cif foc.out [...]")
        sys.exit()

    strudat = open("strudat", "w")
    for fn in fns:
        fo2strudat(fn, fout=strudat)
    strudat.close()

    from . import pykriber
    pykriber.strudat2cif()


def fo2hist_entry():
    fns = sys.argv[1:]
    if not fns:
        print("No files given. \n \n >> Usage: fo2hist foc.out [...]")
        sys.exit()

    for fn in fns:
        fo2hist(fn)

    print("-------------------------------------")
    fo2hist(fns)


def fo2strudat_entry():
    fns = sys.argv[1:]
    if not fns:
        print("No files given. \n \n >> Usage: fo2strudat foc.out [...]")
        sys.exit()

    strudat = open("strudat", "w")
    for fn in fns:
        fo2strudat(fn, fout=strudat)


if __name__ == '__main__':
    fo2hist_entry()
    # fo2strudat_entry()
    # dlsall()
    # cdlsall()