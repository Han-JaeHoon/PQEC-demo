"""
Applying PQEC purification to a tiny quantum algorithm: the Deutsch algorithm.
==============================================================================

The 2-qubit Deutsch algorithm decides whether a one-bit function
f: {0,1} -> {0,1} is *constant* (f(0)=f(1)) or *balanced* (f(0)!=f(1)) with a
single oracle query.  Circuit (query wire q, ancilla wire a):

        q: |0> --H--------------------H--  measure  -> 0: constant, 1: balanced
                      |  U_f  |
        a: |1> --H----|       |-----------

    U_f : |x>|y> -> |x>|y (+) f(x)>.

Ideally the query qubit ends in the *pure* state |0> (constant) or |1>
(balanced), so a single measurement gives the answer with certainty.

Under noise the query qubit becomes a *mixed* state and the answer is no longer
certain.  This is exactly what the PQEC purification primitive fixes: the
purification map

    P(rho) = rho^2 / Tr[rho^2]         (see pqec.py, Eq. 7)

concentrates weight on the dominant eigenvector of rho.  As long as the noise
is below threshold the *correct* answer is still the dominant eigenvector, so
purifying the output drives the success probability back toward 1 -- without
the algorithm ever knowing which answer is correct.

The threshold depends on the channel:

  * local depolarizing on the output qubit ->  p_th = 3/4
  * bit-flip on the output qubit           ->  p_th = 1/2   (deutsch_pqec_bitflip.py)

The whole pipeline (noisy algorithm + purification) reuses the genuine
PennyLane SWAP-gadget circuit from pqec.py, so it is a real mixed-state
simulation, not an algebraic shortcut.  The demo below is parametrized by the
noise channel so it can be reused for both cases.

Run:  python deutsch_pqec.py            # depolarizing, threshold 3/4
      python deutsch_pqec_bitflip.py    # bit-flip,     threshold 1/2
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

from pqec import swap_gadget, purity


# ---------------------------------------------------------------------------
# The four possible oracles U_f  (query = wire 0, ancilla = wire 1)
# ---------------------------------------------------------------------------
def oracle_const0(q, a):
    """f(x) = 0  -> constant.  Identity (do nothing)."""
    qml.Identity(wires=q)


def oracle_const1(q, a):
    """f(x) = 1  -> constant.  Flip the ancilla unconditionally."""
    qml.PauliX(wires=a)


def oracle_balanced_id(q, a):
    """f(x) = x  -> balanced.  CNOT from query to ancilla."""
    qml.CNOT(wires=[q, a])


def oracle_balanced_not(q, a):
    """f(x) = 1 - x  -> balanced.  CNOT then flip the ancilla."""
    qml.CNOT(wires=[q, a])
    qml.PauliX(wires=a)


ORACLES = {
    "constant  f=0":   (oracle_const0,       "constant"),
    "constant  f=1":   (oracle_const1,       "constant"),
    "balanced  f=x":   (oracle_balanced_id,  "balanced"),
    "balanced  f=1-x": (oracle_balanced_not, "balanced"),
}

# The query qubit's ideal (noiseless) answer state and the measurement outcome.
#   constant -> |0>  (bit 0),   balanced -> |1>  (bit 1)
ANSWER_STATE = {"constant": np.array([1, 0], dtype=complex),
                "balanced": np.array([0, 1], dtype=complex)}
ANSWER_BIT = {"constant": 0, "balanced": 1}


# ---------------------------------------------------------------------------
# Noisy Deutsch algorithm on the mixed-state simulator
# ---------------------------------------------------------------------------
_dev = qml.device("default.mixed", wires=2)


@qml.qnode(_dev)
def deutsch_query_rho(oracle, p_noise, channel=qml.DepolarizingChannel):
    """Run the Deutsch algorithm with a single-qubit noise `channel` of strength
    p_noise on the output query qubit, and return that qubit's 2x2 density
    matrix (the reduced state we read the answer from).

    `channel` is any PennyLane single-qubit channel with signature
    (p, wires=...):  qml.DepolarizingChannel (default) or qml.BitFlip.
    """
    qml.PauliX(wires=1)          # ancilla -> |1>
    qml.Hadamard(wires=0)        # query   -> |+>
    qml.Hadamard(wires=1)        # ancilla -> |->
    oracle(0, 1)                 # apply U_f
    qml.Hadamard(wires=0)        # final interference on the query qubit
    channel(p_noise, wires=0)    # gate/readout noise on the output qubit
    return qml.density_matrix(wires=0)


def success_prob(rho, kind):
    """P(read the correct answer) = <answer| rho |answer>."""
    psi = ANSWER_STATE[kind]
    return float(np.real(np.vdot(psi, rho @ psi)))


def purify_once(rho):
    """One purification round via the genuine SWAP-gadget circuit (pqec.py).
    An already-pure state is a fixed point, so we skip the round to avoid the
    harmless 0/0 in the (discarded) minus block when P_- -> 0."""
    if purity(rho) > 1 - 1e-12:
        return rho
    return swap_gadget(rho)["purified"]


def purify_rounds_circuit(rho, ell):
    """Apply `ell` purification rounds using the *genuine* SWAP-gadget circuit
    from pqec.py.  Each round consumes two copies of the current state and
    returns rho^2/Tr[rho^2], extracted from an actual H-CSWAP-H simulation."""
    for _ in range(ell):
        rho = purify_once(rho)
    return rho


# ===========================================================================
# Parametrized demo -- reused for depolarizing (3/4) and bit-flip (1/2).
# ===========================================================================
def run_demo(channel=qml.DepolarizingChannel, channel_name="depolarizing",
             p_th=0.75, p_th_label="3/4", p_example=0.30,
             outfile="deutsch_pqec.png"):
    print("=" * 74)
    print(f" PART 1 -- Noiseless Deutsch algorithm classifies all four oracles")
    print("=" * 74)
    for name, (oracle, kind) in ORACLES.items():
        rho = deutsch_query_rho(oracle, 0.0, channel)
        bit = int(np.argmax(np.real(np.diag(rho))))
        ok = "OK" if bit == ANSWER_BIT[kind] else "FAIL"
        print(f"  {name:16s}  ->  measured bit {bit}  ({kind:8s})  [{ok}]")
    print("  -> query bit 0 = constant, 1 = balanced, each read with certainty.")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(" PART 2 -- Noise degrades the answer; PQEC purification restores it")
    print("=" * 74)

    oracle, kind = ORACLES["balanced  f=x"]   # ideal output |1>
    rho_noisy = deutsch_query_rho(oracle, p_example, channel)

    print(f"\n  Oracle: balanced f=x   {channel_name} p = {p_example}")
    print(f"  noisy output query qubit rho =\n{rho_noisy}")
    print(f"    purity(rho)          = {purity(rho_noisy):.4f}")
    print(f"    P(correct = '1')     = {success_prob(rho_noisy, kind):.4f}   (no QEC)")

    print("\n  Applying PQEC purification rounds (genuine SWAP-gadget circuit):")
    rho = rho_noisy
    for ell in range(1, 5):
        rho = purify_once(rho)
        print(f"    after ell = {ell}:  purity = {purity(rho):.4f}   "
              f"P(correct) = {success_prob(rho, kind):.4f}")
    print("  -> below threshold, purification drives the success probability -> 1.")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(f" PART 3 -- Threshold: purification helps iff p < {p_th_label} ({channel_name})")
    print("=" * 74)

    ps = np.linspace(0.0, 0.99, 100)
    curves = {}
    for ell in [0, 1, 2, 3, 5]:
        succ = []
        for p in ps:
            r = deutsch_query_rho(oracle, p, channel)
            r = purify_rounds_circuit(r, ell)
            succ.append(success_prob(r, kind))
        curves[ell] = np.array(succ)

    # Crossing of the ell=1 and ell=3 curves = empirical threshold (skip the
    # p~0 region where all curves coincide at 1).
    mask = ps > 0.4
    i = int(np.argmin(np.abs(curves[1] - curves[3])[mask]) + np.sum(~mask))
    print(f"\n  measured threshold (ell=1 vs ell=3 crossing): p ~= {ps[i]:.3f}"
          f"   (theory: {p_th})")

    # ---- Figure ----------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    # (a) success probability vs cumulative rounds at fixed noise
    ax = axes[0]
    for p in [0.2, 0.4, 0.6, 0.8]:
        r = deutsch_query_rho(oracle, p, channel)
        succ = [success_prob(r, kind)]
        for _ in range(6):
            r = purify_once(r)
            succ.append(success_prob(r, kind))
        ax.plot(range(7), succ, "-o", ms=4, label=f"p = {p}")
    ax.axhline(0.5, color="k", ls=":", lw=1)
    ax.set_xlabel("purification rounds  $\\ell$")
    ax.set_ylabel("P(correct answer)")
    ax.set_title(f"(a) Deutsch output recovers under purification\n({channel_name} noise)")
    ax.set_ylim(0.0, 1.02)
    ax.legend(frameon=False)

    # (b) success probability vs noise p, several depths -> cross at p_th
    ax = axes[1]
    for ell in [0, 1, 2, 3, 5]:
        ax.plot(ps, curves[ell], ("r:" if ell == 0 else "-"),
                label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
    ax.axvline(p_th, color="k", ls="--", lw=1.2)
    ax.text(p_th, 0.06, f" $p_{{th}}={p_th_label}$", fontsize=10)
    ax.axhline(0.5, color="k", ls=":", lw=1)
    ax.set_xlabel("physical error rate  p")
    ax.set_ylabel("P(correct answer)")
    ax.set_title(f"(b) Threshold: purification helps iff  p < {p_th_label}")
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig(outfile, dpi=140)
    print(f"\n  saved  {outfile}")
    print("\nDone.")
    return curves


if __name__ == "__main__":
    # Default: local depolarizing channel on the output qubit, threshold 3/4.
    run_demo()
