"""
Recovering a Bell state from many noisy copies with PQEC.
=========================================================

A Bell-state factory runs the standard circuit

        q0: |0> --H--*--            |Phi+> = (|00> + |11>) / sqrt(2)
                     |
        q1: |0> -----X--

but every copy comes out *noisy* -- a mixed two-qubit state rho.  Given many
such copies we purify them with the PQEC primitive, generalised from one qubit
to a two-qubit register: the SWAP gadget becomes a SWAP *test* between two
2-qubit registers (two parallel Fredkin gates + an ancilla), and reading the
ancilla extracts

    P(rho) = rho^2 / Tr[rho^2]        (pqec.py, Eq. 7, dimension-independent).

This concentrates weight on the dominant eigenvector of rho.  While the noise is
below threshold that eigenvector is still |Phi+>, so purification drives both
the fidelity F = <Phi+|rho|Phi+> and the entanglement (concurrence) back toward
1 -- the Bell state is restored from copies that were barely entangled.

Two noise models on the freshly made Bell state:

  * local depolarizing on BOTH qubits -> |Phi+> stays the dominant Bell
    component for every p, so purification restores it everywhere except the
    single maximally-mixed point p = 3/4 (there rho = I/4 is a fixed point).

  * bit-flip on ONE qubit -> mixes |Phi+> with |Psi+> = (|01>+|10>)/sqrt(2).
    |Phi+> is dominant only for p < 1/2; above 1/2 purification locks onto the
    WRONG Bell state |Psi+>.  Clean threshold p_th = 1/2.

Everything uses genuine PennyLane mixed-state circuits (Bell prep + noise, and
the 5-wire SWAP-test gadget), not algebraic shortcuts.

Run:  python bell_pqec.py       # console checks + bell_pqec.png
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# Target Bell states and figures of merit
# ---------------------------------------------------------------------------
PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)   # (|00>+|11>)/v2
PSI_PLUS = np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2)   # (|01>+|10>)/v2
_sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
_YY = np.kron(_sy, _sy)


def fidelity(rho, psi):
    """F = <psi| rho |psi> for a pure target |psi>."""
    return float(np.real(np.vdot(psi, rho @ psi)))


def purity(rho):
    return float(np.real(np.trace(rho @ rho)))


def concurrence(rho):
    """Wootters concurrence of a 2-qubit state (0 = separable, 1 = max. entangled)."""
    R = rho @ _YY @ rho.conj() @ _YY
    ev = np.sort(np.sqrt(np.clip(np.real(np.linalg.eigvals(R)), 0, None)))[::-1]
    return max(0.0, ev[0] - ev[1] - ev[2] - ev[3])


# ---------------------------------------------------------------------------
# Genuine noisy Bell-state factory on the mixed-state simulator
# ---------------------------------------------------------------------------
_bell_dev = qml.device("default.mixed", wires=2)


@qml.qnode(_bell_dev)
def bell_rho(p, noise="depol"):
    """Prepare |Phi+> and apply noise; return the 4x4 density matrix.

    noise = "depol"   : local depolarizing of strength p on BOTH qubits
    noise = "bitflip" : bit-flip of strength p on qubit 0 only
    """
    qml.Hadamard(0)
    qml.CNOT(wires=[0, 1])
    if noise == "depol":
        qml.DepolarizingChannel(p, wires=0)
        qml.DepolarizingChannel(p, wires=1)
    elif noise == "bitflip":
        qml.BitFlip(p, wires=0)
    else:
        raise ValueError(noise)
    return qml.density_matrix(wires=[0, 1])


# ---------------------------------------------------------------------------
# PQEC purification for a 2-qubit register: genuine SWAP-test gadget
#   wire 0 = ancilla,  wires [1,2] = register A (kept),  [3,4] = register B
#   controlled-SWAP of the two registers = two parallel Fredkin gates
# ---------------------------------------------------------------------------
_gadget_dev = qml.device("default.mixed", wires=5)


@qml.qnode(_gadget_dev)
def _swap_gadget_2q(rho_AB):
    """SWAP test on two 2-qubit copies; return joint (ancilla, register A)."""
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2, 3, 4])   # two noisy copies
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])          # Fredkin on qubit pair 0
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])          # Fredkin on qubit pair 1
    qml.Hadamard(0)
    return qml.density_matrix(wires=[0, 1, 2])           # ancilla + register A


def purify_once(rho):
    """One PQEC round on a 2-qubit state via the genuine SWAP-test circuit.
    Extracts rho^2 = P_+ rho_+ - P_- rho_- (Eq. 9) from the ancilla blocks and
    normalises.  Already-pure states are fixed points (skip to avoid 0/0)."""
    if purity(rho) > 1 - 1e-12:
        return rho
    d = rho.shape[0]
    sigma = _swap_gadget_2q(np.kron(rho, rho))    # joint (anc, A), 8x8
    rho_sq = sigma[:d, :d] - sigma[d:, d:]        # ancilla |0> block minus |1> block
    return rho_sq / np.trace(rho_sq)


def purify_rounds(rho, ell):
    for _ in range(ell):
        rho = purify_once(rho)
    return rho


# ===========================================================================
def main():
    print("=" * 74)
    print(" PART 1 -- Sanity: the SWAP-test gadget purifies a 2-qubit register")
    print("=" * 74)
    rng = np.random.default_rng(2)
    err = 0.0
    for _ in range(500):
        A = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
        rho = A @ A.conj().T
        rho = rho / np.trace(rho)
        circ = purify_once(rho)
        alg = (rho @ rho) / np.trace(rho @ rho)
        err = max(err, np.abs(circ - alg).max())
    print(f"  circuit output == rho^2/Tr[rho^2]   max error (500 states): {err:.2e}")

    rho0 = bell_rho(0.0)
    print(f"  noiseless Bell state:  F = {fidelity(rho0, PHI_PLUS):.4f}   "
          f"concurrence = {concurrence(rho0):.4f}")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(" PART 2 -- Restoring |Phi+> from noisy copies (depolarizing p=0.30)")
    print("=" * 74)
    p = 0.30
    rho = bell_rho(p, "depol")
    print(f"\n  one noisy copy (local depol p={p}) rho =\n{np.real(rho)}")
    print(f"    F(|Phi+>)   = {fidelity(rho, PHI_PLUS):.4f}   "
          f"concurrence = {concurrence(rho):.4f}   purity = {purity(rho):.4f}")
    print("\n  PQEC purification rounds (each consumes 2 copies of the current state):")
    r = rho
    for ell in range(1, 6):
        r = purify_once(r)
        print(f"    ell={ell}  (N=2^{ell}={2**ell} copies):  "
              f"F = {fidelity(r, PHI_PLUS):.4f}   "
              f"concurrence = {concurrence(r):.4f}   purity = {purity(r):.4f}")
    print("  -> fidelity and entanglement both driven back to 1: Bell state restored.")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(" PART 3 -- Bit-flip on one qubit: threshold p=1/2, wrong state above")
    print("=" * 74)
    print("\n   p     F(Phi+)->   F(Psi+)->   (after 8 rounds)")
    for p in [0.2, 0.4, 0.5, 0.6, 0.8]:
        r = purify_rounds(bell_rho(p, "bitflip"), 8)
        print(f"  {p:.2f}    {fidelity(r, PHI_PLUS):.4f}      "
              f"{fidelity(r, PSI_PLUS):.4f}")
    print("  -> below 1/2 purification restores |Phi+>; above 1/2 it locks onto |Psi+>.")

    # ---- Figure ----------------------------------------------------------
    ps = np.linspace(0.0, 0.99, 100)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))

    # (a) depolarizing p=0.3: fidelity & concurrence vs rounds
    ax = axes[0, 0]
    r = bell_rho(0.30, "depol")
    F = [fidelity(r, PHI_PLUS)]
    C = [concurrence(r)]
    rr = r
    for _ in range(6):
        rr = purify_once(rr)
        F.append(fidelity(rr, PHI_PLUS))
        C.append(concurrence(rr))
    ax.plot(range(7), F, "-o", ms=4, label="fidelity  F")
    ax.plot(range(7), C, "-s", ms=4, label="concurrence")
    ax.axhline(1.0, color="k", ls=":", lw=1)
    ax.set_xlabel("purification rounds  $\\ell$")
    ax.set_ylabel("value")
    ax.set_title("(a) Depolarizing p=0.30: Bell state restored")
    ax.set_ylim(0.0, 1.05)
    ax.legend(frameon=False)

    # (b) depolarizing: fidelity vs p for several depths
    ax = axes[0, 1]
    for ell in [0, 1, 2, 3, 5]:
        Fp = [fidelity(purify_rounds(bell_rho(p, "depol"), ell), PHI_PLUS) for p in ps]
        ax.plot(ps, Fp, ("r:" if ell == 0 else "-"),
                label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
    ax.axvline(0.75, color="k", ls="--", lw=1.2)
    ax.text(0.75, 0.06, " p=3/4 (I/4)", fontsize=9)
    ax.set_xlabel("physical error rate  p")
    ax.set_ylabel("fidelity to |Phi+>")
    ax.set_title("(b) Depolarizing: restored except at the fully-mixed point")
    ax.legend(frameon=False)

    # (c) bit-flip: fidelity to |Phi+> vs p -> clean threshold 1/2
    ax = axes[1, 0]
    for ell in [0, 1, 2, 3, 5]:
        Fp = [fidelity(purify_rounds(bell_rho(p, "bitflip"), ell), PHI_PLUS) for p in ps]
        ax.plot(ps, Fp, ("r:" if ell == 0 else "-"),
                label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
    ax.axvline(0.5, color="k", ls="--", lw=1.2)
    ax.text(0.5, 0.06, " $p_{th}=1/2$", fontsize=10)
    ax.set_xlabel("physical error rate  p")
    ax.set_ylabel("fidelity to |Phi+>")
    ax.set_title("(c) Bit-flip: threshold p=1/2")
    ax.legend(frameon=False)

    # (d) bit-flip at ell=5: which Bell state wins
    ax = axes[1, 1]
    Fphi = [fidelity(purify_rounds(bell_rho(p, "bitflip"), 5), PHI_PLUS) for p in ps]
    Fpsi = [fidelity(purify_rounds(bell_rho(p, "bitflip"), 5), PSI_PLUS) for p in ps]
    ax.plot(ps, Fphi, "-", label="F(|Phi+>) (correct)")
    ax.plot(ps, Fpsi, "-", label="F(|Psi+>) (wrong)")
    ax.axvline(0.5, color="k", ls="--", lw=1.2)
    ax.set_xlabel("physical error rate  p")
    ax.set_ylabel("fidelity after $\\ell=5$")
    ax.set_title("(d) Bit-flip: purification amplifies the dominant Bell state")
    ax.legend(frameon=False)

    fig.tight_layout()
    fig.savefig("bell_pqec.png", dpi=140)
    print("\n  saved  bell_pqec.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
