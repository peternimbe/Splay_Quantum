from functools import lru_cache
from fractions import Fraction
import numpy as np
import networkx as nx
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
@lru_cache(None)
def trees(lo,hi):
    if lo>=hi: return (None,)
    return tuple((k,l,r) for k in range(lo,hi)
                 for l in trees(lo,k) for r in trees(k+1,hi))

def inorder(t):
    return inorder(t[1])+(t[0],)+inorder(t[2]) if t else ()

def balanced(lo,hi):
    if lo>=hi: return None
    k=(lo+hi)//2
    return (k,balanced(lo,k),balanced(k+1,hi))

@lru_cache(None)
def access(t,x):
    keys=inorder(t); n=len(keys)
    if keys!=tuple(range(n)) or x not in keys: raise ValueError('keys must be 0..n-1')
    left=[-1]*n; right=[-1]*n; parent=[-1]*n
    def load(z,p=-1):
        if z is None: return
        k,l,r=z; parent[k]=p
        if l: left[k]=l[0]; load(l,k)
        if r: right[k]=r[0]; load(r,k)
    load(t); root=t[0]; v=root; c=1
    while v!=x:
        v=left[v] if x<v else right[v]; c+=1
    rotations=0
    def rotate(z):
        nonlocal root,rotations
        p=parent[z]; g=parent[p]
        if left[p]==z: b=right[z]; left[p]=b; right[z]=p
        else: b=left[z]; right[p]=b; left[z]=p
        if b!=-1: parent[b]=p
        parent[p]=z; parent[z]=g
        if g==-1: root=z
        elif left[g]==p: left[g]=z
        else: right[g]=z
        rotations+=1
    while parent[x]!=-1:
        p=parent[x]; g=parent[p]
        if g==-1: rotate(x)
        elif (left[g]==p)==(left[p]==x): rotate(p); rotate(x)
        else: rotate(x); rotate(x)
    def freeze(k):
        return None if k==-1 else (k,freeze(left[k]),freeze(right[k]))
    result=freeze(root)
    assert rotations==c-1 and root==x and inorder(result)==keys
    return result,c

def ref_rotate(t,direction):
    k,l,r=t
    if direction=='R':
        a,b,c=l; return (a,b,(k,c,r))
    a,b,c=r; return (a,(k,l,b),c)

def ref_at(t,path,direction):
    if not path: return ref_rotate(t,direction)
    k,l,r=t
    return (k,ref_at(l,path[1:],direction),r) if path[0]==0 else (k,l,ref_at(r,path[1:],direction))

def reference_access(t,x):
    initial=t; c=0
    while True:
        path=[]; z=t
        while z[0]!=x:
            j=0 if x<z[0] else 1; path.append(j); z=z[j+1]
        if c==0: c=len(path)+1
        if not path: return t,c
        if len(path)==1: t=ref_at(t,(), 'R' if path[-1]==0 else 'L')
        elif path[-2]==path[-1]:
            d='R' if path[-1]==0 else 'L'
            t=ref_at(t,tuple(path[:-2]),d); t=ref_at(t,tuple(path[:-2]),d)
        else:
            t=ref_at(t,tuple(path[:-1]),'R' if path[-1]==0 else 'L')
            t=ref_at(t,tuple(path[:-2]),'R' if path[-2]==0 else 'L')

def execute(t,seq):
    trace=[]; states=[t]
    for x in seq:
        t,c=access(t,int(x)); trace.append(c); states.append(t)
    return sum(trace),trace,states

def mirror(t,n):
    return (n-1-t[0],mirror(t[2],n),mirror(t[1],n)) if t else None

def build_graph(n):
    ts=trees(0,n); ti={t:i for i,t in enumerate(ts)}; N=len(ts)
    nxt=np.empty((N,n),dtype=int); cost=np.empty((N,n),dtype=int)
    for i,t in enumerate(ts):
        for x in range(n):
            z,c=access(t,x); nxt[i,x]=ti[z]; cost[i,x]=c
    states=[(i,i) for i in range(N)]; ids={s:i for i,s in enumerate(states)}
    edges=[]; parent=[None]*N; pos=0
    while pos<len(states):
        i,j=states[pos]
        for x in range(n):
            for keep in (0,1):
                z=(int(nxt[i,x]),int(nxt[j,x]) if keep else j)
                if z not in ids:
                    ids[z]=len(states); states.append(z); parent.append(len(edges))
                edges.append((pos,ids[z],int(cost[i,x]),int(cost[j,x]) if keep else 0,x,keep))
        pos+=1
    return {'n':n,'trees':ts,'states':states,'edges':np.array(edges,dtype=np.int64),
            'parent':parent,'diagonal':list(range(N)),'nxt':nxt,'cost':cost}

def solve_potential(g,edge_ids=None):
    E=g['edges']; V=len(g['states'])
    if edge_ids is not None: E=E[np.array(edge_ids,dtype=int)]
    m=len(E); rr=np.repeat(np.arange(m),3)
    cc=np.column_stack([np.zeros(m,dtype=int),1+E[:,0],1+E[:,1]]).ravel()
    vv=np.column_stack([-E[:,2],-np.ones(m),np.ones(m)]).ravel()
    A=coo_matrix((vv,(rr,cc)),shape=(m,V+1)).tocsr()
    objective=np.zeros(V+1); objective[0]=1
    fit=linprog(objective,A_ub=A,b_ub=-E[:,3],bounds=[(0,g['n'])]+[(None,None)]*V,method='highs')
    if not fit.success: raise RuntimeError(fit.message)
    return fit.x[0],fit.x[1:]

def exact_certificate(g,K):
    # Integer longest-path inequalities. If no positive-weight cycle exists,
    # a maximizing walk can have all repeated-vertex cycles removed.
    E=g['edges']; u,v,a,b=E[:,:4].T; V=len(g['states'])
    p,q=K.numerator,K.denominator; w=q*b-p*a
    H=np.zeros(V,dtype=np.int64)
    for iteration in range(V):
        new=H.copy(); np.maximum.at(new,v,H[u]+w)
        if np.array_equal(new,H):
            slack=H[v]-H[u]-w
            assert int(slack.min())>=0
            return H,slack,iteration
        H=new
    raise ValueError('Candidate K lacks an exact certificate: positive cycle or insufficient rational precision')

def tight_cycle(g,K,H):
    E=g['edges']; p,q=K.numerator,K.denominator
    graph=nx.DiGraph()
    for eid,(u,v,a,b,x,d) in enumerate(E):
        if H[v]-H[u]==q*b-p*a: graph.add_edge(int(u),int(v),eid=eid)
    try: cyc=nx.find_cycle(graph)
    except nx.NetworkXNoCycle: raise ValueError('No tight cycle: candidate is not proven optimal')
    ids=[graph[u][v]['eid'] for u,v in cyc]
    e=E[ids]; assert Fraction(int(e[:,3].sum()),int(e[:,2].sum()))==K
    return ids

def prefix_to(g,state):
    ids=[]
    while g['parent'][state] is not None:
        eid=g['parent'][state]; ids.append(eid); state=int(g['edges'][eid,0])
    return list(reversed(ids)),state
