# Progress log

Status of the PQEC-demo project and what has been built so far.

## 1. Core PQEC primitive (baseline)

The purification-based QEC primitive from Raghoonanan & Byrnes,
*Quantum Error Correction by Purification* (arXiv:2603.11568), implemented as a
genuine PennyLane circuit and verified numerically.

- **`pqec.py`** — SWAP gadget as an explicit `default.mixed` circuit
  (H–CSWAP–H), the purification map `P(ρ) = ρ²/Tr[ρ²]`, multi-round
  purification, and the depolarizing / dephasing noise channels.
- **`verify_pqec.py`** — checks Eqs. (5)–(10), (7), (34) to ~1e-16 over 2000
  random states; reproduces fidelity-vs-rounds curves and the error thresholds
  (depolarizing p=3/4, dephasing p=1/2). Saves `pqec_purification_map.png` and
  `pqec_thresholds.png`.
- **`draw_circuits.py`** — draws the elementary gadget (M=1, M=2) and the
  depth-efficient binary tree (ℓ=2, N=4 copies).

**Verified equations**

| Quantity | Paper Eq. | Max error (2000 states) |
|----------|-----------|-------------------------|
| `P± = (1 ± Tr ρ²)/2` | (5) | 5.6e-16 |
| `ρ± = (ρ ± ρ²)/2P±` | (6) | 5.6e-16 |
| `P₊ρ₊ + P₋ρ₋ = ρ` | (8) | 5.6e-16 |
| `P₊ρ₊ − P₋ρ₋ = ρ²` | (9) | 5.6e-16 |
| circuit output `= ρ²/Tr ρ²` | (7) | 3.3e-16 |
| Bloch rescaling `r → 2r/(1+|r|²)` | (34) | 5.0e-16 |

**Measured thresholds** (ℓ=1 vs ℓ=3 crossing): depolarizing **0.750** (3/4),
dephasing **0.500** (1/2).

## 2. PQEC applied to an algorithm — the Deutsch algorithm (new)

Plugged the purification gadget into a real computation to show it corrects the
output of an actual quantum algorithm, not just abstract density matrices.

- **`deutsch_pqec.py`** — runs the 2-qubit Deutsch algorithm on the mixed-state
  simulator with a depolarizing channel on the output query qubit, then purifies
  that qubit with the genuine SWAP-gadget circuit from `pqec.py`.
- **`draw_deutsch_pqec.py`** — draws the full circuit for one purification
  round: two noisy Deutsch runs feeding a SWAP gadget between the two query
  qubits (`deutsch_pqec_circuit.png`).

**Key idea.** Ideally the Deutsch query qubit ends pure — `|0⟩` (constant) or
`|1⟩` (balanced) — read with certainty. Noise makes it mixed and the answer
uncertain. The correct answer is the *dominant eigenvector* of the noisy state,
so purification (`ρ → ρ²/Tr[ρ²]`) concentrates weight back onto it, with no
knowledge of which answer is correct.

**Results**

- All four oracles (constant f=0, f=1; balanced f=x, f=1−x) classify correctly
  with no noise.
- At depolarizing p=0.30 the success probability recovers over purification
  rounds: `0.80 → 0.94 → 0.996 → 1.000` (ℓ = 0,1,2,3).
- The same **p=3/4** threshold reappears: below it purification drives
  P(correct) → 1, above it it amplifies the wrong answer → 0; all round-count
  curves cross exactly at p=3/4 (measured 0.750).

Figures: `deutsch_pqec.png` (recovery + threshold), `deutsch_pqec_circuit.png`
(full circuit).

## How to run

```bash
pip install -r requirements.txt
python verify_pqec.py         # core equation checks + figures
python draw_circuits.py       # gadget / binary-tree circuits
python deutsch_pqec.py        # PQEC applied to the Deutsch algorithm
python draw_deutsch_pqec.py   # full Deutsch + PQEC circuit diagram
```

## Possible next steps

- Dephasing version of the Deutsch demo (threshold p=1/2).
- ℓ=2 binary-tree variant that actually consumes N=2^ℓ fresh copies per output.
- Extend to a slightly larger algorithm (e.g. 2-bit Deutsch–Jozsa) with a
  multi-qubit SWAP gadget.
