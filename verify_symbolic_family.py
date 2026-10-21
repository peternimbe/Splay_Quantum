from splay_core import *
from itertools import product
def switches(seq):return sum(x!=y for x,y in zip(seq,seq[1:]))
def adjacent_identity(t,seq):
    if not seq:return 0,0,0
    total,trace,_=execute(t,seq)
    first_switch=next((i for i in range(1,len(seq)) if seq[i]!=seq[i-1]),None)
    alpha=trace[0]-1+(trace[first_switch]-2 if first_switch is not None else 0)
    assert total==len(seq)+switches(seq)+alpha
    return total,alpha,switches(seq)
adjacency_checks=identity_checks=deletion_checks=0
for nn in range(2,7):
    for t in trees(0,nn):
        x=t[0]
        for y in (x-1,x+1):
            if 0<=y<nn:
                result,c=access(t,y);child=result[1] if x<y else result[2]
                assert result[0]==y and child and child[0]==x
                adjacency_checks+=1
for nn in range(2,5):
    for t in trees(0,nn):
        for aa in range(nn-1):
            for length in range(1,9):
                for bits in product((0,1),repeat=length):
                    seq=tuple(aa+b for b in bits);c,alpha,sw=adjacent_identity(t,seq)
                    assert 0<=alpha<=2*nn-3;identity_checks+=1
                    if length<=5:
                        for mask in range(1<<length):
                            sub=tuple(x for j,x in enumerate(seq) if mask>>j&1)
                            d,_,_=adjacent_identity(t,sub)
                            assert switches(sub)<=sw and d<=c+2*nn-3;deletion_checks+=1
family_check_counts={'adjacency_lemma_checks':adjacency_checks,'exact_cost_identities':identity_checks,'deletion_checks':deletion_checks}
print(family_check_counts)
def adjacent_initial_allowance(t,a):
    # Four structural initialization cases determine every nonempty execution.
    b=a+1
    alpha_cases=[execute(t,(a,))[0]-1,execute(t,(b,))[0]-1,
                 execute(t,(a,b))[0]-3,execute(t,(b,a))[0]-3]
    return max(alpha_cases)-min(alpha_cases)

print("Symbolic-family implementation checks passed.")
