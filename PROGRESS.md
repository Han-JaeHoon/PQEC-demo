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
  simulator with a noise channel on the output query qubit, then purifies that
  qubit with the genuine SWAP-gadget circuit from `pqec.py`. Parametrized by the
  noise channel via `run_demo(...)`; default is depolarizing (`p_th=3/4`).
- **`deutsch_pqec_bitflip.py`** — reuses `run_demo(...)` with a **bit-flip**
  channel on the output qubit (`p_th=1/2`), saving `deutsch_pqec_bitflip.png`.
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
- The **threshold depends on the output-qubit noise channel**:

  | Channel      | Recovery (p=0.30)        | Measured `p_th` | Theory |
  |--------------|--------------------------|-----------------|--------|
  | Depolarizing | `0.80 → 1.000`           | 0.750           | 3/4    |
  | Bit-flip     | `0.70 → 1.000`           | 0.500           | 1/2    |

  Below threshold purification drives P(correct) → 1; above it, it amplifies the
  wrong answer → 0. All round-count curves cross exactly at `p_th`.

Figures: `deutsch_pqec.png`, `deutsch_pqec_bitflip.png` (recovery + threshold),
`deutsch_pqec_circuit.png` (full circuit).

## 3. PQEC on a 2-qubit register — restoring a Bell state (new)

Generalized the primitive from one qubit to a **2-qubit register** and used it
to recover an entangled Bell state from many noisy copies.

- **`bell_pqec.py`** — a noisy Bell factory (`H`–`CNOT` → `|Φ⁺⟩` + noise) emits
  mixed 4×4 copies; the SWAP gadget becomes a genuine 5-wire SWAP *test* between
  two 2-qubit registers (M=2 = two parallel Fredkin gates). Tracks fidelity and
  Wootters concurrence. Saves `bell_pqec.png`.
- **`draw_bell_pqec.py`** — draws the full circuit: two noisy Bell factories
  feeding the M=2 SWAP gadget (`bell_pqec_circuit.png`).

**Key idea.** The SWAP test extracts `ρ²/Tr[ρ²]` in any dimension (verified to
2e-16 on 500 random 2-qubit states), concentrating weight on the dominant
eigenvector. Below threshold that is `|Φ⁺⟩`, so purification restores both the
fidelity and the entanglement of the Bell state.

**Results**

- Local depolarizing p=0.30: one copy is barely entangled (`F=0.52`,
  concurrence `0.04`); purification restores both to 1 over `N=2^ℓ` copies
  (`F: 0.52 → 0.78 → 0.97 → 1.00`; concurrence `0.04 → 0.56 → 0.95 → 1.00`).
- Threshold structure by channel:

  | Channel (on the Bell state)   | Behaviour |
  |-------------------------------|-----------|
  | Local depolarizing, both qubits | `|Φ⁺⟩` restored for all `p` except the fully-mixed point `p=3/4` (`ρ=I/4`) |
  | Bit-flip, one qubit           | clean threshold `p=1/2`; above it purification amplifies the wrong Bell state `|Ψ⁺⟩` |

Figures: `bell_pqec.png` (fidelity + concurrence recovery, both thresholds),
`bell_pqec_circuit.png` (full circuit).

## 4. The paper's actual protocol — purified observables (new)

The earlier demos read the purified state off the density matrix
(`ρ² = blk₊ − blk₋`), which is an algebraic shortcut. The paper's real protocol
is an **observable-estimation** scheme: keep the SWAP-gadget ancilla, measure
**Z on the ancilla** together with the target observable **O**, and take the
ratio of correlators.

- **`pqec_observable.py`** — after the gadget the joint (ancilla, register)
  state is `½(I⊗ρ + Z⊗ρ²)`, so `⟨Z⊗O⟩ = Tr(Oρ²)`, `⟨Z⊗I⟩ = Tr(ρ²)`, and
  `⟨O⟩_purified = ⟨Z⊗O⟩/⟨Z⊗I⟩ = Tr(Oρ²)/Tr(ρ²)`. For `ℓ` rounds the sign is the
  total parity `Ω = Πᵢ Z_i` of all ancillas: `⟨O⟩_ℓ = ⟨Ω⊗O⟩/⟨Ω⟩ =
  Tr(Oρ^N)/Tr(ρ^N)`, `N=2^ℓ`. Implemented as genuine circuits (ℓ=1 on 5 wires,
  ℓ=2 binary tree on 11 wires) and matched to the analytic value to ~1e-16.

**Scenario (depolarized Bell, `ρ_ε=(1-ε)|Φ⁺⟩⟨Φ⁺|+ε I/D`, `O=|Φ⁺⟩⟨Φ⁺|`, ε=0.40):**

| | `⟨O⟩` | effective `ε′` |
|--|:-----:|:--------------:|
| ideal | 1.000 | 0 |
| no QEC | 0.700 | 0.400 |
| PQEC ℓ=1 | 0.942 | 0.077 |
| PQEC ℓ=2 | 0.999 | 0.002 |
| PQEC ℓ=3 | 1.000 | 0.000 |

The measured observable improves toward its ideal value; PQEC returns the
expectation value of a much less noisy effective state (`ε′ ≪ ε`). This is the
key correction over §2–3: the purified value is **measured** (ancilla-parity
correlator), not obtained by subtracting density matrices.

Figure: `pqec_observable.png`.

## 5. Faulty PQEC — noise on the gadget operations (new)

Earlier demos assumed a perfect gadget. Here the **gadget gates themselves are
noisy** (the realistic fault-tolerance question).

- **`pqec_gate_noise.py`** — inserts a depolarizing channel after every gadget
  gate (H, and each Fredkin, strength `g`) plus an optional ancilla readout
  bit-flip `r`, then measures `⟨O⟩ = ⟨Z⊗O⟩/⟨Z⊗I⟩` on the genuine circuit.

**Findings**

- **Gate-error threshold `g*`:** gate noise degrades `⟨O⟩_PQEC`; beyond `g*` the
  purifier is worse than no QEC. At `ε=0.40`, `g*≈0.145`.
- `g*` grows with input noise `ε` (`0.2→0.075`, `0.4→0.145`, `0.6→0.205`): a
  nearly-clean input is easily spoiled by a faulty purifier; a noisy input
  tolerates a sloppier gadget.
- **Readout error self-mitigates:** a symmetric ancilla flip scales `⟨Z⊗O⟩` and
  `⟨Z⊗I⟩` by the same `1−2r`, cancelling in the ratio.

Figure: `pqec_gate_noise.png`.

## How to run

```bash
pip install -r requirements.txt
python verify_pqec.py         # core equation checks + figures
python draw_circuits.py       # gadget / binary-tree circuits
python deutsch_pqec.py         # PQEC on Deutsch, depolarizing noise (p_th=3/4)
python deutsch_pqec_bitflip.py # PQEC on Deutsch, bit-flip noise    (p_th=1/2)
python draw_deutsch_pqec.py    # full Deutsch + PQEC circuit diagram
python bell_pqec.py            # restore a Bell state from noisy copies
python draw_bell_pqec.py       # full noisy-Bell x2 + M=2 gadget circuit
python pqec_observable.py      # purified observable via ancilla-parity correlator
python pqec_gate_noise.py      # faulty gadget: gate-error threshold g*
```

## Possible next steps

- ℓ=2 binary-tree variant that actually consumes N=2^ℓ fresh copies per output
  (currently modelled by recursively re-purifying the same ρ).
- GHZ / larger entangled states with the M-qubit gadget.
- Extend to a slightly larger algorithm (e.g. 2-bit Deutsch–Jozsa).
