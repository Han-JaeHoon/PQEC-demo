"""
PQEC the way the paper actually does it: purified OBSERVABLE expectation values.
================================================================================

The SWAP gadget does **not** hand you a purified density matrix, and it is
**not** postselection.  Following the paper (and the whiteboard derivation),
after the gadget the joint (ancilla, kept register A) state is

    P_+ |0><0| (x) rho_+  +  P_- |1><1| (x) rho_-  =  (1/2)( I (x) rho + Z (x) rho^2 ).

So for ANY observable O on the register, measuring Z on the ancilla together
with O on the register gives, on average,

    <Z_anc (x) O>  =  Tr(O rho^2),        <Z_anc (x) I>  =  Tr(rho^2),

and therefore the expectation value with respect to the *purified* state
P(rho) = rho^2/Tr(rho^2) is obtained as a ratio of two physically measured
correlators:

    <O>_purified  =  <Z_anc (x) O> / <Z_anc (x) I>  =  Tr(O rho^2) / Tr(rho^2).

For ell rounds (binary tree, N = 2^ell copies) the sign is the TOTAL PARITY
Omega = prod_i Z_(anc,i) of all 2^ell - 1 ancillas:

    <O>_ell  =  <Omega (x) O> / <Omega>  =  Tr(O rho^N) / Tr(rho^N).

This is a genuine measurement protocol (Pauli-Z on the ancillas correlated with
O on the target), NOT an algebraic subtraction of density matrices.  Below we
verify it on a real PennyLane circuit and then show that PQEC improves the
observable of a noisy Bell state toward its ideal value.

Scenario (whiteboard example (2)):
    ideal      : |Phi+> = (|00>+|11>)/sqrt2,  O = |Phi+><Phi+|,  <O> = 1
    depolarized: rho_eps = (1-eps)|Phi+><Phi+| + eps I/D,  <O> = 1 - 3eps/4
    after PQEC : effective rho_eps' with eps' < eps,  <O> -> 1

Run:  python pqec_observable.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

np.set_printoptions(precision=4, suppress=True)

# ---------------------------------------------------------------------------
# Single- and two-qubit operators
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
_PAULIS = [I2, X, Y, Z]
_TWOQ_PAULIS = [np.kron(a, b) for a in _PAULIS for b in _PAULIS]   # 16, [0]=I⊗I

PHI = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)           # |Phi+>
O_FID = np.outer(PHI, PHI.conj())                                 # |Phi+><Phi+|
O_ZZ = np.kron(Z, Z)                                              # a generic Pauli obs
D = 4                                                             # 2-qubit dimension


def global_depolarizing_kraus(eps):
    """Kraus operators of the 2-qubit global depolarizing channel
        E(rho) = (1-eps) rho + eps I/D,   D = 4.
    Built from the 16 two-qubit Paulis (a valid CPTP set)."""
    K = [np.sqrt(1 - eps + eps / 16) * _TWOQ_PAULIS[0]]
    K += [np.sqrt(eps / 16) * P for P in _TWOQ_PAULIS[1:]]
    return K


def noisy_bell_rho(eps):
    """The depolarized Bell state rho_eps = (1-eps)|Phi+><Phi+| + eps I/4."""
    return (1 - eps) * O_FID + eps * np.eye(4) / D


def obs_purified(O, rho, ell):
    """Analytic purified observable Tr(O rho^N)/Tr(rho^N), N = 2^ell.
    This is exactly what the parity-correlator measurement below estimates."""
    rN = np.linalg.matrix_power(rho, 2 ** ell)
    return float(np.real(np.trace(O @ rN)) / np.real(np.trace(rN)))


def effective_eps(O_fid_value):
    """Invert <O_fid> = 1 - 3 eps/4 to report the effective error rate."""
    return 4 * (1 - O_fid_value) / 3


# ---------------------------------------------------------------------------
# Genuine circuits: measure <Omega (x) O> and <Omega> and take the ratio
# ---------------------------------------------------------------------------
def _prep_noisy_bell(q0, q1, eps):
    qml.Hadamard(q0)
    qml.CNOT(wires=[q0, q1])
    qml.QubitChannel(global_depolarizing_kraus(eps), wires=[q0, q1])


_dev1 = qml.device("default.mixed", wires=5)     # ell=1: anc + 2 copies


@qml.qnode(_dev1)
def _gadget_ell1(eps, O):
    """One purification round. wires: 0=anc, [1,2]=copy A (kept), [3,4]=copy B."""
    _prep_noisy_bell(1, 2, eps)
    _prep_noisy_bell(3, 4, eps)
    qml.Hadamard(0)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    qml.Hadamard(0)
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),   # <Z (x) O>
            qml.expval(qml.PauliZ(0)))                                    # <Z (x) I>


_dev2 = qml.device("default.mixed", wires=11)    # ell=2: 3 anc + 4 copies


@qml.qnode(_dev2)
def _gadget_ell2(eps, O):
    """Two rounds (binary tree, N=4).  ancillas 0,1,2 ; copies [3,4],[5,6],[7,8],[9,10]."""
    c0, c1, c2, c3 = [3, 4], [5, 6], [7, 8], [9, 10]
    for c in (c0, c1, c2, c3):
        _prep_noisy_bell(c[0], c[1], eps)
    # round 1: two gadgets in parallel, keep c0 and c2
    for anc, A, B in [(0, c0, c1), (1, c2, c3)]:
        qml.Hadamard(anc)
        qml.ctrl(qml.SWAP, control=anc)(wires=[A[0], B[0]])
        qml.ctrl(qml.SWAP, control=anc)(wires=[A[1], B[1]])
        qml.Hadamard(anc)
    # round 2: merge c0, c2, keep c0
    qml.Hadamard(2)
    qml.ctrl(qml.SWAP, control=2)(wires=[c0[0], c2[0]])
    qml.ctrl(qml.SWAP, control=2)(wires=[c0[1], c2[1]])
    qml.Hadamard(2)
    parity = qml.PauliZ(0) @ qml.PauliZ(1) @ qml.PauliZ(2)   # total parity Omega
    return (qml.expval(parity @ qml.Hermitian(O, wires=c0)),
            qml.expval(parity))


def obs_purified_circuit(eps, O, ell):
    """Purified observable from the GENUINE circuit: <Omega (x) O> / <Omega>."""
    zO, zI = (_gadget_ell1 if ell == 1 else _gadget_ell2)(eps, O)
    return float(zO) / float(zI)


# ===========================================================================
def main():
    eps = 0.40

    print("=" * 74)
    print(" PART 1 -- The observable IS measured: <Omega (x) O> / <Omega>")
    print("=" * 74)
    rho = noisy_bell_rho(eps)
    print(f"\n  Noisy Bell state, depolarizing eps = {eps}   (rho_eps = (1-eps)|Phi+><Phi+| + eps I/4)")
    print(f"  Observable O = |Phi+><Phi+|   (ideal <O> = 1)\n")

    zO, zI = _gadget_ell1(eps, O_FID)
    print(f"  ell=1 circuit:  <Z_anc (x) O> = {float(zO):.4f}  = Tr(O rho^2) = {np.real(np.trace(O_FID@rho@rho)):.4f}")
    print(f"                  <Z_anc (x) I> = {float(zI):.4f}  = Tr(rho^2)   = {np.real(np.trace(rho@rho)):.4f}")
    o1_circ = float(zO) / float(zI)
    o1_ana = obs_purified(O_FID, rho, 1)
    print(f"     -> <O>_PQEC(ell=1) = {o1_circ:.4f}   (analytic Tr(O rho^2)/Tr(rho^2) = {o1_ana:.4f},  match {abs(o1_circ-o1_ana):.1e})")

    o2_circ = obs_purified_circuit(eps, O_FID, 2)   # genuine 11-wire tree (~20s)
    o2_ana = obs_purified(O_FID, rho, 2)
    print(f"  ell=2 tree circuit (parity Z1 Z2 Z3):  <O>_PQEC(ell=2) = {o2_circ:.4f}   (analytic {o2_ana:.4f},  match {abs(o2_circ-o2_ana):.1e})")
    print("\n  -> the purified expectation value is obtained from a REAL measurement")
    print("     (ancilla Z correlated with O), not from subtracting density matrices.")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 74)
    print(" PART 2 -- Does PQEC improve the observable?  (whiteboard scenario)")
    print("=" * 74)
    for name, O, ideal in [("O = |Phi+><Phi+|  (fidelity)", O_FID, 1.0),
                           ("O = Z(x)Z          (parity)  ", O_ZZ, 1.0)]:
        print(f"\n  {name}   ideal <O> = {ideal:.0f}")
        noisy = float(np.real(np.trace(O @ rho)))
        print(f"    no QEC              <O> = {noisy:.4f}")
        for ell in [1, 2, 3]:
            val = obs_purified(O, rho, ell)
            print(f"    PQEC ell={ell} (N={2**ell:>2}) <O> = {val:.4f}")

    print("\n  Effective error rate (from <O_fid> = 1 - 3 eps'/4):")
    print(f"    no QEC : eps' = {eps:.3f}")
    for ell in [1, 2, 3]:
        epsp = effective_eps(obs_purified(O_FID, rho, ell))
        print(f"    ell={ell}  : eps' = {epsp:.3f}   (< eps = {eps})")
    print("  -> PQEC yields the expectation value of a much less noisy Bell state.")

    # ---- Figure ----------------------------------------------------------
    epss = np.linspace(0.0, 1.0, 101)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    ax = axes[0]
    ax.axhline(1.0, color="k", ls=":", lw=1, label="ideal  <O>=1")
    rhos = [noisy_bell_rho(e) for e in epss]
    ax.plot(epss, [float(np.real(np.trace(O_FID @ r))) for r in rhos], "r--",
            label="no QEC")
    for ell in [1, 2, 3]:
        ax.plot(epss, [obs_purified(O_FID, r, ell) for r in rhos], "-",
                label=f"PQEC $\\ell={ell}$")
    ax.set_xlabel("depolarizing strength  $\\varepsilon$")
    ax.set_ylabel(r"$\langle O_{\Phi}\rangle = $ fidelity to $|\Phi^+\rangle$")
    ax.set_title("(a) PQEC improves the measured observable toward 1")
    ax.set_ylim(0.2, 1.02)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    ax.plot(epss, epss, "r--", label="no QEC  ($\\varepsilon'=\\varepsilon$)")
    for ell in [1, 2, 3]:
        epsp = [effective_eps(obs_purified(O_FID, r, ell)) for r in rhos]
        ax.plot(epss, epsp, "-", label=f"PQEC $\\ell={ell}$")
    ax.plot([0, 1], [0, 1], "k:", lw=0.8)
    ax.set_xlabel("physical error  $\\varepsilon$")
    ax.set_ylabel("effective error  $\\varepsilon'$  after PQEC")
    ax.set_title(r"(b) Effective noise $\varepsilon' < \varepsilon$")
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("pqec_observable.png", dpi=140)
    print("\n  saved  pqec_observable.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
