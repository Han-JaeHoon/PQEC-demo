# Quantum Error Correction by Purification (PQEC) — PennyLane verification

PennyLane-based numerical verification of the purification-based quantum error
correction (**PQEC**) primitive from

> J. Raghoonanan & T. Byrnes, *Quantum Error Correction by Purification*,
> arXiv:2603.11568 (2026).

The core primitive is the **SWAP gadget**: a SWAP test applied to two identical
noisy copies `ρ ⊗ ρ` with an ancilla. Reading out the ancilla and post-processing
the outcomes with a parity sign extracts the *purified* component

```
P(ρ) = ρ² / Tr[ρ²] = Σ λ_i² |i⟩⟨i| / Σ λ_i²
```

which concentrates weight on the dominant eigenvector, raising the state's purity
and fidelity. Repeating over ℓ rounds consumes `N = 2^ℓ` copies and produces
`ρ^N / Tr[ρ^N]`, driving any non-maximally-mixed state toward a pure state.

This repo implements the gadget as a **genuine PennyLane circuit** on the
mixed-state simulator and checks every key equation of the paper to machine
precision, then reproduces the error-correction and threshold behaviour.

## Files

| File | Description |
|------|-------------|
| [`pqec.py`](pqec.py) | Core library: SWAP gadget as an explicit `default.mixed` circuit (H–CSWAP–H), purification map, noise channels, fidelity helpers |
| [`verify_pqec.py`](verify_pqec.py) | Verifies Eqs. (5)–(10), (7), (34); reproduces fidelity-vs-rounds and error thresholds; saves figures |
| [`draw_circuits.py`](draw_circuits.py) | Draws the full quantum circuits (Figs. 1 & 2) |

## Setup & run

```bash
python -m venv pqec_env
source pqec_env/bin/activate          # Windows: pqec_env\Scripts\activate
pip install -r requirements.txt

python verify_pqec.py     # numerical checks + result figures
python draw_circuits.py   # circuit diagrams
```

## What is verified

The SWAP-test circuit reproduces the paper's equations to ~1e-16:

| Quantity | Paper Eq. | Max error (2000 random states) |
|----------|-----------|-------------------------------|
| `P± = (1 ± Tr ρ²)/2` | (5) | 5.6e-16 |
| `ρ± = (ρ ± ρ²)/2P±` | (6) | 5.6e-16 |
| `P₊ρ₊ + P₋ρ₋ = ρ` (average → input) | (8) | 5.6e-16 |
| `P₊ρ₊ − P₋ρ₋ = ρ²` (subtract → purified) | (9) | 5.6e-16 |
| circuit output `= ρ²/Tr ρ²` | (7) | 3.3e-16 |
| Bloch rescaling `r → 2r/(1+|r|²)` | (34) | 5.0e-16 |

Measured error thresholds (crossing point of ℓ=1 vs ℓ=3 curves):

| Channel | Measured | Paper |
|---------|----------|-------|
| Local depolarizing | **0.750** | 3/4 |
| Local dephasing | **0.500** | 1/2 |

## Circuits

**Elementary SWAP gadget** (Fig. 1b): `H` on ancilla, controlled-SWAP (Fredkin),
`H`, measure. For an M-qubit register the CSWAP becomes M parallel Fredkin gates.

![SWAP gadget](circuit_gadget_M1.png)

**Depth-efficient binary tree**, ℓ=2, N=4 copies (Fig. 1a):

![Binary tree](circuit_tree_ell2.png)

## Results

**Purification as a radial Bloch rescaling** — fixed points at |r|=0 (mixed) and
|r|=1 (pure); repeated rounds drive purity → 1:

![Purification map](pqec_purification_map.png)

**Fidelity vs cycles and error thresholds** — curves for different round counts ℓ
cross exactly at the thresholds p=3/4 (depolarizing) and p=1/2 (dephasing):

![Thresholds](pqec_thresholds.png)

## Example: restoring a Bell state from noisy copies (M=2 PQEC)

`bell_pqec.py` generalizes PQEC to a **2-qubit register**: a Bell factory
(`H`–`CNOT` → `|Φ⁺⟩`) emits noisy mixed copies, and the SWAP gadget becomes a
SWAP *test* between two 2-qubit registers (**M=2 = two parallel Fredkin gates**
+ ancilla). Reading the ancilla still extracts `ρ²/Tr[ρ²]` (dimension-
independent), restoring both fidelity and entanglement (Wootters concurrence).

Local depolarizing p=0.30 — one copy is barely entangled, purification restores
it over `N=2^ℓ` copies:

| rounds ℓ | copies N | fidelity | concurrence |
|:--------:|:--------:|:--------:|:-----------:|
| 0 | 1 | 0.520 | 0.040 |
| 1 | 2 | 0.779 | 0.558 |
| 2 | 4 | 0.974 | 0.948 |
| 3 | 8 | 1.000 | 0.999 |

Threshold structure by channel:

- **Local depolarizing (both qubits):** `|Φ⁺⟩` stays the dominant Bell component
  for every `p`, so it is restored everywhere **except the fully-mixed point
  `p=3/4`** (`ρ=I/4`, a fixed point).
- **Bit-flip (one qubit):** clean threshold **`p=1/2`** — above it purification
  amplifies the *wrong* Bell state `|Ψ⁺⟩`.

```bash
python bell_pqec.py        # fidelity + concurrence recovery, saves bell_pqec.png
python draw_bell_pqec.py   # full noisy-Bell x2 + M=2 gadget circuit
```

![Bell restoration](bell_pqec.png)

![Bell + PQEC circuit](bell_pqec_circuit.png)
