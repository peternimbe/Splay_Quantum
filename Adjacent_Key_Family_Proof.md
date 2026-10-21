## 14. Symbolic theorem for an arbitrary-size access family

We now give a separate, modest symbolic result for an explicitly defined infinite family. This elementary proof demonstrates the proof-producing objective; it is not claimed as a new theorem in the literature or as a consequence automatically learned by the feature-fitting code.

**Theorem (two adjacent accessed keys).** For any initial BST \(T\) on \(n\ge2\) keys, adjacent keys \(a,b\) in sorted order, any sequence \(X\) over \(\{a,b\}\), and subsequence \(Y\preceq X\),
\[
\boxed{C_T(Y)\le C_T(X)+2n-3.}
\]
The asymptotic coefficient 1 is sharp because \(Y=X\) is permitted. The additive allowance is conservative.

**Adjacency lemma.** If the root is \(a\) and the next access is its successor \(b\), the search path goes right and then only left, because \(b\) is the minimum of the right subtree. Bottom-up zig-zig rotations below the root shorten that left chain; the final step is a zig or zig-zag at the original root. The resulting root is \(b\), with \(a\) its immediate left child. For a predecessor the argument is reflected. After both keys have been accessed, repeated accesses cost 1 and switches cost 2. That root/child configuration is preserved by each switch.

Let \(m=|X|\), and let \(s(X)\) count switches between consecutive keys. Let \(d_1\) be the first access cost, and \(d_*\) the cost at the first switch, when one exists. For nonempty \(X\),
\[
 C_T(X)=m+s(X)+\alpha_T(X),\qquad
 \alpha_T(X)=d_1-1+\begin{cases}d_*-2,&s(X)>0,\\0,&s(X)=0.\end{cases}
\]
Indeed, after the first access a repetition leaves the tree unchanged; after the first switch the adjacency lemma gives the stated unit/repetition and two-unit/switch costs. Since \(d_1\le n\) and \(2\le d_*\le n\),
\[
 0\le\alpha_T(X)\le(n-1)+(n-2)=2n-3.
\]
For a sequence with no switch, \(\alpha_T(X)\le n-1\le2n-3\).

Deleting positions cannot increase the number of switches: each surviving switch spans a disjoint interval containing an original switch. Thus \(|Y|\le|X|\) and \(s(Y)\le s(X)\). The cost identity then gives
\[
 C_T(Y)-C_T(X)
 =(|Y|-|X|)+(s(Y)-s(X))+(\alpha_T(Y)-\alpha_T(X))\le2n-3.
\]
An empty subsequence is immediate, completing the proof for all sizes and lengths in this restricted family.

This structural formula has a bulk term \(m+s(X)\) and an explicitly bounded initialization remainder. It is independent of the learned finite-state potentials. It does not address access sequences over unrestricted key sets. Exhaustive and sampled checks below validate the implementation; the written proof, rather than those checks, supports the arbitrary-size statement.

**Tree-specific refinement.** Every nonempty execution has one of four initialization remainders: those of the single-key starts \((a)\), \((b)\), and the first-switch orders \((a,b)\), \((b,a)\). Set
\[
\mathcal A_T=\{C_T(a)-1,C_T(b)-1,C_T(a,b)-3,C_T(b,a)-3\},
\quad B_T=\max\mathcal A_T-\min\mathcal A_T.
\]
Then the same argument proves \(C_T(Y)\le C_T(X)+B_T\), and \(0\le B_T\le2n-3\). This is an explicit compact structural allowance obtained from four short executions, valid for arbitrary sequence lengths in the family. It is conservative, not claimed to be the smallest possible allowance.
