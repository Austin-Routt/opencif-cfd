import sys, numpy as np
def load_points(fn):
    rows=[]
    with open(fn, encoding='utf-8', errors='replace') as f:
        for line in f:
            if line.startswith('%') or not line.strip(): continue
            rows.append([float(t) for t in line.split()])
    return np.array(rows)
for fn in sys.argv[1:]:
    a=load_points(fn)
    print('==',fn.split('/comsol/')[-1], a.shape)
    for j,nm in enumerate('xyz'):
        u=np.unique(a[:,j]); print(' ',nm,'n=',len(u),'min',u.min(),'max',u.max(), 'steps', np.unique(np.round(np.diff(u),6))[:6])
    nanrows=np.isnan(a[:,3:]).any(axis=1)
    print('  NaN rows',nanrows.sum(), ' all-NaN rows', np.isnan(a[:,3:]).all(axis=1).sum())
    # check grid order
    print('  first 3 rows', a[:3,:3].tolist(), ' row 4000', a[4000,:3].tolist())
