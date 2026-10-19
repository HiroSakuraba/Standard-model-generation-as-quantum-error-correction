"""Reproducible exhaustive 2^5 mask placement audit for SU3 selective Z3 echoes."""
import json
from itertools import combinations
from scipy.sparse.linalg import expm_multiply
from su3_native_control_compilation import construct,source_data,syndrome_decoder,evolution,toggled_matrices,entanglement_fidelity,T


def main():
    states,K,d,W,H=construct()
    V,_,gram=syndrome_decoder(*source_data())
    assert gram.rank()==5
    rows=[]
    for r in range(6):
        for subset in combinations(range(5),r):
            out=(evolution(K,toggled_matrices(d,W,subset),T,8) if subset
                 else expm_multiply(-1j*T*H,K,traceA=0))
            row={'subset':list(subset),'cost_links':r,
                 'F_raw':entanglement_fidelity([out],K),
                 'F_recovered':entanglement_fidelity([out],K,V)}
            rows.append(row)
            print('MASK',json.dumps(row,sort_keys=True),flush=True)
    assert len(rows)==32
    best={r:max((x for x in rows if x['cost_links']==r),key=lambda x:x['F_recovered']) for r in range(6)}
    assert [r['subset'] for r in best.values()]==[[],[4],[0,4],[0,1,4],[0,2,3,4],[0,1,2,3,4]]
    assert best[3]['F_recovered']<.99 and best[4]['F_recovered']>.99999
    assert abs(best[5]['F_recovered']-.999997303395923)<1e-9
    results={'n_plaquettes':5,'dimension':243,'cycles':8,'g2':1,'time':T,'rows':rows,'best_per_cost':{str(k):v for k,v in best.items()}}
    with open('su3_control_mask_scan_results.json','w') as f:json.dump(results,f,indent=2,sort_keys=True)
    print('ALL 32 MASKS AUDITED: best four-site mask',best[4]['subset'],'F_e',best[4]['F_recovered'],flush=True)

if __name__=='__main__':main()
