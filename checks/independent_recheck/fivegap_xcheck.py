"""Cross-check table pressures against the mid-depth plane; locate max |w|. Usage: fivegap_xcheck.py <five_gap_dir>"""
import sys, os, glob, numpy as np
ROOT=sys.argv[1]
for m in ["V1_coarse","V1_extra_coarse","V3_coarse_resolved","V3_coarse_stored_May7","V3_extra_coarse"]:
    d=os.path.join(ROOT,m)
    b=np.loadtxt(glob.glob(os.path.join(d,"*middepth_z70*.txt"))[0],comments="%",encoding="utf-8")
    def tab(n):
        p=[q for q in glob.glob(os.path.join(d,"*")) if os.path.basename(q).lower().startswith(f"table{n}")][0]
        L=[l for l in open(p,encoding='utf-8') if not l.startswith('%') and l.strip()]
        hdr=[l for l in open(p,encoding='utf-8') if 'Point:' in l][0]
        import re
        pts=[tuple(map(float,t.split(','))) for t in re.findall(r'Point: \(([^)]*)\)',hdr)]
        return pts,[float(t) for t in L[0].split()]
    print('==',m)
    for n in (4,5):
        pts,vals=tab(n)
        for (x,y,z),val in zip(pts,vals):
            i=np.where((abs(b[:,0]-x)<1e-6)&(abs(b[:,1]-y)<1e-6))[0]
            pp=b[i[0],6] if len(i) else np.nan
            print(f'  T{n} ({x:g},{y:g},{z:g}) table={val:.9f} plane={pp:.9f} diff={val-pp:.3e}')
    ok=~np.isnan(b[:,3]); B=b[ok]
    k=np.argmax(abs(B[:,5])); print('  max|w| at x,y=',B[k,:2],'w=',B[k,5],' speed max at', B[np.argmax(np.sqrt((B[:,3:6]**2).sum(1))),:2])
    print('  y extent of fluid', B[:,1].min(), B[:,1].max(), ' u at inlet x=0 max', B[B[:,0]==0,3].max())
