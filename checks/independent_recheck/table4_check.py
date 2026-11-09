import sys
for fn in sys.argv[1:]:
    rows=[[float(t) for t in l.split()] for l in open(fn,encoding='utf-8') if l.strip() and not l.startswith('%')]
    r=rows[0]; print(fn.split('/single_gap/')[-1], 'rows=',len(rows), 'cols=',len(r))
    print('  last col = %.16e ; recomputed sout/in = %.16e ; recomputed imbalance = %.6e (file %.6e)'%(r[4], r[2]/r[0], (r[0]-r[1]-r[2])/r[0], r[3]))
