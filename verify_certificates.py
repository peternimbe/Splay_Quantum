from splay_core import *
import json
from pathlib import Path
import sys
def verify_payload(payload):
    n=payload['n']; g=build_graph(n)
    assert payload['states']==[list(s) for s in g['states']]
    assert payload['edges']==g['edges'].tolist()
    assert payload['diagonal_state_ids']==g['diagonal']
    p,q=payload['K']; assert isinstance(p,int) and isinstance(q,int) and q>0
    H=payload['H']; assert len(H)==len(g['states']) and all(type(h)==int for h in H)
    for u,v,a,b,x,d in payload['edges']:
        assert H[v]-H[u] >= q*b-p*a
    cyc=payload['cycle_edge_ids']; assert cyc
    first=payload['edges'][cyc[0]][0]; state=first; A=B=0
    for eid in cyc:
        u,v,a,b,x,d=payload['edges'][eid]; assert u==state
        A+=a; B+=b; state=v
    assert state==first and q*B==p*A and A>0
    prefix,diag=prefix_to(g,first); assert diag in g['diagonal']
    return {'n':n,'K':str(Fraction(p,q)),'states':len(g['states']),'edges_verified':len(g['edges']),'status':'PASS'}

if __name__=='__main__':
    directory=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).parent/'splay_outputs'
    paths=sorted(directory.glob('certificate_n*.json'))
    if not paths: raise SystemExit('No certificate files found')
    for path in paths: print(verify_payload(json.loads(path.read_text())))
