"""Structural splay potentials with rational, all-edge verification."""
from splay_core import *
import math
import pandas as pd

FEATURE_NAMES=[
 'depth_sum_X','depth_sum_Y','depth_square_sum_X','depth_square_sum_Y',
 'subtree_square_sum_X','subtree_square_sum_Y','integer_rank_sum_X','integer_rank_sum_Y',
 'leaf_count_X','leaf_count_Y','imbalance_sum_X','imbalance_sum_Y',
 'depth_absolute_difference','subtree_absolute_difference','ancestor_disagreement',
 'parent_disagreement','root_gap','same_root']

@lru_cache(None)
def tree_statistics(t):
    n=len(inorder(t)); depths=np.zeros(n,dtype=np.int64); sizes=np.zeros(n,dtype=np.int64)
    parents=np.full(n,-1,dtype=np.int64); ancestors=np.zeros((n,n),dtype=bool)
    leaves=0; imbalance=0
    def walk(z,depth=0,parent=-1,path=()):
        nonlocal leaves,imbalance
        if z is None:return 0
        k,l,r=z; depths[k]=depth; parents[k]=parent
        for a in path:ancestors[a,k]=True
        ls=walk(l,depth+1,k,path+(k,)); rs=walk(r,depth+1,k,path+(k,))
        sizes[k]=1+ls+rs; leaves+=int(l is None and r is None); imbalance+=abs(ls-rs)
        return int(sizes[k])
    walk(t)
    stats=(int(depths.sum()),int(depths@depths),int(sizes@sizes),
           sum(int(z).bit_length()-1 for z in sizes),leaves,imbalance)
    return depths,sizes,parents,ancestors,stats

def structural_features(tx,ty):
    dx,sx,px,ax,fx=tree_statistics(tx); dy,sy,py,ay,fy=tree_statistics(ty)
    paired=[v for pair in zip(fx,fy) for v in pair]
    paired += [int(np.abs(dx-dy).sum()),int(np.abs(sx-sy).sum()),int(np.count_nonzero(ax!=ay)),
               int(np.count_nonzero(px!=py)),abs(tx[0]-ty[0]),int(tx[0]==ty[0])]
    return np.array(paired,dtype=np.int64)

def feature_table(g):
    return np.array([structural_features(g['trees'][i],g['trees'][j]) for i,j in g['states']],dtype=np.int64)

def feature_problem(graphs,ns):
    differences=[]; aa=[]; bb=[]; labels=[]
    for n in ns:
        g=graphs[n]; F=feature_table(g); E=g['edges']
        differences.append(F[E[:,1]]-F[E[:,0]])
        aa.append(E[:,2]); bb.append(E[:,3])
        labels.extend((n,eid) for eid in range(len(E)))
    return np.vstack(differences),np.concatenate(aa),np.concatenate(bb),labels

def propose_structural(D,a,b,chosen):
    ix=np.array(sorted(chosen),dtype=int); d=D.shape[1]
    A=np.column_stack([-a[ix],D[ix]])
    obj=np.zeros(d+1);obj[0]=1
    fit=linprog(obj,A_ub=A,b_ub=-b[ix],bounds=[(0,None)]+[(None,None)]*d,method='highs')
    if not fit.success:raise RuntimeError(fit.message)
    return float(fit.x[0]),fit.x[1:]

def sparse_formula(D,a,b,Kfloat):
    d=D.shape[1]
    fit=linprog(np.ones(2*d),A_ub=np.column_stack([D,-D]),b_ub=Kfloat*a-b+1e-8,
                bounds=[(0,None)]*(2*d),method='highs')
    if not fit.success:raise RuntimeError(fit.message)
    raw=fit.x[:d]-fit.x[d:]
    coeff=[Fraction(0) if abs(c)<1e-7 else Fraction(float(c)).limit_denominator(1_000) for c in raw]
    L=math.lcm(*(c.denominator for c in coeff))
    numer=np.array([c.numerator*(L//c.denominator) for c in coeff],dtype=object)
    vals=np.asarray(D,dtype=object)@numer
    # Recompute the required coefficient exactly after rational reconstruction.
    K=max(Fraction(L*int(bi)+int(vi),L*int(ai)) for ai,bi,vi in zip(a,b,vals))
    return K,coeff

def check_structural(g,K,coeff,F=None):
    F=feature_table(g) if F is None else F
    E=g['edges']; L=math.lcm(K.denominator,*(c.denominator for c in coeff))
    cn=np.array([c.numerator*(L//c.denominator) for c in coeff],dtype=object)
    P=np.asarray(F,dtype=object)@cn
    kn=K.numerator*(L//K.denominator)
    residual=L*np.asarray(E[:,3],dtype=object)-kn*np.asarray(E[:,2],dtype=object)+P[E[:,1]]-P[E[:,0]]
    bad=np.flatnonzero(np.array([int(z)>0 for z in residual]))
    return {'passed':len(bad)==0,'bad_edges':bad,'residual':residual,'potential_numerators':P,
            'scale':L,'max_residual':Fraction(int(max(residual)),L)}

def learn_structural(graphs,ns,seed):
    D,a,b,labels=feature_problem(graphs,ns)
    local_rng=np.random.default_rng(seed)
    chosen=set(local_rng.choice(len(a),size=min(32,len(a)),replace=False).tolist())
    records=[]
    for iteration in range(len(a)+1):
        K,c=propose_structural(D,a,b,chosen)
        residual=b-K*a+D@c; bad=np.flatnonzero(residual>1e-7)
        records.append({'training_sizes':str(tuple(ns)),'iteration':iteration,'active_constraints':len(chosen),
                        'proposed_K':K,'violations':len(bad),'largest_residual':float(residual.max())})
        if len(bad)==0:break
        worst=bad[np.argsort(residual[bad])[-min(32,len(bad)):]]
        before=len(chosen);chosen.update(worst.tolist());assert len(chosen)>before
    else:raise RuntimeError('Structural refinement did not terminate')
    Kexact,coeff=sparse_formula(D,a,b,K)
    for n in ns:assert check_structural(graphs[n],Kexact,coeff)['passed']
    return Kexact,coeff,pd.DataFrame(records)

def structural_payload(g,K,coeff,stage):
    check=check_structural(g,K,coeff); E=g['edges']; P=check['potential_numerators']; L=check['scale']
    B=max(Fraction(int(P[s]-min(P)),L) for s in g['diagonal'])
    worst=int(np.argmax(check['residual']))
    u,v,a,b,x,d=map(int,E[worst]); i,j=g['states'][u]
    prefix,diag=prefix_to(g,u)
    initial_index=g['states'][diag][0]
    witness={'prefix_actions':E[prefix,4:6].tolist(),'initial_tree':g['trees'][initial_index],
             'edge_id':worst,'full_tree':g['trees'][i],'subsequence_tree':g['trees'][j],
             'key':x,'keep':d,'a':a,'b':b,'residual':str(check['max_residual'])}
    return {'stage':stage,'n':g['n'],'K':[K.numerator,K.denominator],
            'feature_names':FEATURE_NAMES,'coefficients':[[c.numerator,c.denominator] for c in coeff],
            'passed':bool(check['passed']),'violation_count':len(check['bad_edges']),
            'max_residual':str(check['max_residual']),
            'B_max_if_passed':str(B) if check['passed'] else None,'worst_edge':witness,
            'cost_model':'bottom-up splay; pre-access depth plus one'}

def verify_structural_payload(payload,g=None):
    assert payload['feature_names']==FEATURE_NAMES
    g=build_graph(payload['n']) if g is None else g
    K=Fraction(*payload['K']); coeff=[Fraction(*z) for z in payload['coefficients']]
    assert len(coeff)==len(FEATURE_NAMES)
    fresh=structural_payload(g,K,coeff,payload['stage'])
    for key in ('passed','violation_count','max_residual','B_max_if_passed'):
        assert fresh[key]==payload[key],key
    # JSON normalizes nested immutable tree tuples to lists.
    import json
    assert json.loads(json.dumps(fresh['worst_edge']))==payload['worst_edge']
    return {'stage':payload['stage'],'n':payload['n'],'classification':'PASS' if fresh['passed'] else 'REJECTED',
            'edges_checked':len(g['edges']),'claim_reproduced':True}
