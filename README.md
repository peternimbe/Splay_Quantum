# Exact Finite-State Certification of Splay Subsequence Bounds with Quantum-Assisted Counterexample Search

Open `Quantum_Splay_Certification.ipynb` and choose **Restart & Run All**. The notebook is self-contained and includes executed outputs. No API keys or QPU are required.

## Setup

Use Python 3.10 or later:

```bash
python -m pip install -r requirements.txt
python -m pip install jupyterlab
jupyter lab Quantum_Splay_Certification.ipynb
```

Once dependencies are installed, computations run offline. If a separate kernel cannot be launched, `python run_inprocess.py` executes the cells through IPython and saves outputs into the notebook.

## Extended research workflow

1. Enumerate paired-tree transitions and compute exact state-table bounds for n=1..5.
2. Match rational upper certificates with repeatable cycle witnesses.
3. Learn an 18-feature structural potential on n=1..4; minimize coefficient L1 norm as a sparsity heuristic; reconstruct rational coefficients and verify every training edge with Python integers.
4. Freeze the formula and test size 5. Preserve failed transitions as executable prefix/action witnesses.
5. Refine through size 5, freeze again, and verify all 209,088 reachable transitions at the new size 6.
6. Record failure honestly, then calibrate only K while retaining the refined structural coefficients. Verify that post-test certificate through n=6.
7. Compare feature ablations, run ideal quantum-search simulations, and account explicitly for oracle preprocessing.
8. Prove an arbitrary-size theorem for access sequences supported on two adjacent keys; validate the implementation exhaustively and with longer sampled sequences.

The state-table coefficients are 1,1,1,3/2,8/5 at sizes 1..5, with finite additive allowances. Their optimality is certified by matching cycle witnesses under the specified node-visit cost model.

The initial structural formula has 8 nonzero terms and coefficient 26/17, but fails on 352 size-5 transitions. The refined formula has 14 nonzero terms and coefficient 40/21, but fails on 1,538 size-6 transitions. A post-test calibration to 64/21, with those same 14 coefficients, gives a verified upper certificate through size 6. This calibrated coefficient is not claimed optimal, and its verification is not a held-out prediction.

## Files and checks

- `splay_core.py`: base transition and graph implementation.
- `structural_model.py`: features, compact-template fitting, rational verification, witness export.
- `verify_certificates.py`: rebuild and verify state-table certificates and equal-ratio cycles.
- `verify_structural_certificates.py`: rebuild and verify all exported structural PASS/REJECTED classifications, without LP fitting.
- `verify_symbolic_family.py`: exhaustive implementation checks for the adjacent-key lemma, cost identity, and deletion bound.
- `Adjacent_Key_Family_Proof.md`: complete elementary arbitrary-size proof.
- `splay_outputs/`: 16 figures in PNG and vector PDF, result CSVs, state-table and structural JSON certificates/witnesses, environment metadata and SHA-256 manifest.

Run:

```bash
python verify_certificates.py splay_outputs
python verify_structural_certificates.py splay_outputs
python verify_symbolic_family.py
```

The symbolic theorem states C_T(Y) <= C_T(X)+2n-3 for arbitrary initial n-key BSTs when the accesses use only two adjacent keys. A tighter tree-specific allowance is computed from four short initialization executions. This is a restricted, elementary family result; no first-in-literature claim is made.

## Scope and reproducibility

Costs are pre-access depth plus one for standard bottom-up splaying. Every finite-size upper certificate allows an explicitly bounded additive initialization term. Default exact state-table enumeration stops at five keys; structural evaluation builds size six. State-space size grows exponentially, so increase limits deliberately.

Quantum experiments use ideal NumPy statevectors and classically precomputed phase-oracle tables. They do not implement a reversible procedural splay oracle or establish an end-to-end hardware speedup. The global-depolarization sensitivity model is illustrative, not calibrated hardware noise.

Passing a finite-size certificate does not prove a uniform constant for unrestricted splay executions. Rejecting a potential inequality is not automatically a counterexample to the underlying all-path cost bound. The notebook distinguishes these claims, preserves failed transfer outcomes, and makes no claim to solve the dynamic optimality conjecture. Publication-level originality remains to be established by a broader prior-work comparison.
