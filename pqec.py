"""
Purification Quantum Error Correction (PQEC)
=============================================

PennyLane-based verification of

    Raghoonanan & Byrnes, "Quantum Error Correction by Purification",
    arXiv:2603.11568v1 (2026).

The paper's central primitive is the *SWAP gadget*: apply a SWAP test to two
identical noisy copies rho (x) rho with an ancilla, and read out the ancilla.
Conditioning on the ancilla outcome projects the first copy onto

    rho_+/-  =  (rho +/- rho^2) / (2 P_+/-)                              (Eq. 6)

with outcome probabilities  P_+/-  =  (1 +/- Tr[rho^2]) / 2             (Eq. 5).

Averaging the two conditional states just returns the input,

    P_+ rho_+ + P_- rho_-  =  rho                                       (Eq. 8)

but subtracting them (using the ancilla/parity sign) extracts the *purified*
component

    P_+ rho_+ - P_- rho_-  =  rho^2                                     (Eq. 9)

whose normalised form is the purification map

    P(rho)  =  rho^2 / Tr[rho^2]  =  sum_i lambda_i^2 |i><i| / sum_i lambda_i^2.  (Eq. 7)

This concentrates weight on the dominant eigenvector -> higher purity/fidelity.
For a single qubit rho = (I + r.sigma)/2 this is a radial Bloch rescaling

    r  ->  2 r / (1 + |r|^2)                                            (Eq. 34)

with fixed points |r| = 0 (maximally mixed) and |r| = 1 (pure), so purification
drives any non-maximally-mixed qubit toward a pure state.

This module implements the SWAP gadget as a genuine PennyLane circuit on the
mixed-state simulator and checks every one of the equations above numerically,
then studies the multi-round error-correction behaviour and thresholds.
"""

import numpy as np
import pennylane as qml

# ---------------------------------------------------------------------------
# Single-qubit helpers
# ---------------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = (X, Y, Z)


def bloch_to_rho(r):
    """Bloch vector r=(x,y,z) -> density matrix rho = (I + r.sigma)/2."""
    x, y, z = r
    return 0.5 * (I2 + x * X + y * Y + z * Z)


def rho_to_bloch(rho):
    """density matrix -> Bloch vector."""
    return np.real(np.array([np.trace(rho @ P) for P in PAULI]))


def purity(rho):
    return np.real(np.trace(rho @ rho))


def fidelity_pure(rho, psi):
    """F = <psi| rho |psi> for a pure target |psi>."""
    return float(np.real(np.vdot(psi, rho @ psi)))


# ---------------------------------------------------------------------------
# The SWAP gadget as an explicit PennyLane circuit  (Fig. 1b)
# ---------------------------------------------------------------------------
# wire 0 : ancilla       wire 1 : copy A (kept)     wire 2 : copy B (discarded)
_dev_gadget = qml.device("default.mixed", wires=3)


@qml.qnode(_dev_gadget)
def _swap_gadget_circuit(rho_AB):
    """Run the SWAP test on two copies and return the joint (ancilla, A) state.

    rho_AB is the 4x4 density matrix of the two copies (usually kron(rho, rho)).
    """
    qml.QubitDensityMatrix(rho_AB, wires=[1, 2])   # prepare the two noisy copies
    qml.Hadamard(0)                                 # ancilla into |+>
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 2])     # controlled-SWAP (Fredkin)
    qml.Hadamard(0)                                 # interfere
    return qml.density_matrix(wires=[0, 1])         # keep ancilla + copy A


def swap_gadget(rho):
    """Apply one SWAP-gadget purification round to a single-qubit state rho.

    Returns a dict with the conditional states rho_+/-, their probabilities
    P_+/-, and the purified state P(rho) = rho^2/Tr[rho^2], all obtained from
    the explicit circuit above.
    """
    sigma = _swap_gadget_circuit(np.kron(rho, rho))   # joint (ancilla, A), 4x4

    # Ancilla |0> outcome -> top-left 2x2 block;  |1> -> bottom-right block.
    blk_plus = sigma[0:2, 0:2]
    blk_minus = sigma[2:4, 2:4]
    P_plus = float(np.real(np.trace(blk_plus)))
    P_minus = float(np.real(np.trace(blk_minus)))
    rho_plus = blk_plus / P_plus
    rho_minus = blk_minus / P_minus

    # Eq. 9: subtract with the parity sign to isolate rho^2, then normalise.
    rho_sq = P_plus * rho_plus - P_minus * rho_minus   # == rho @ rho
    purified = rho_sq / np.trace(rho_sq)

    return dict(
        rho_plus=rho_plus, rho_minus=rho_minus,
        P_plus=P_plus, P_minus=P_minus,
        rho_sq=rho_sq, purified=purified,
    )


def purify_map(rho):
    """Direct algebraic form of the purification map, Eq. 7: rho^2 / Tr[rho^2]."""
    rho_sq = rho @ rho
    return rho_sq / np.trace(rho_sq)


def purify_rounds(rho, ell):
    """Apply ell rounds of purification.  ell rounds consume N = 2^ell copies
    and produce rho^N / Tr[rho^N]  (Eq. 22)."""
    for _ in range(ell):
        rho = purify_map(rho)
    return rho


# ---------------------------------------------------------------------------
# Single-qubit noise channels
# ---------------------------------------------------------------------------
def depolarizing(rho, p):
    """Local depolarizing channel, Eq. 52:
        (1-p) rho + (p/3)(X rho X + Y rho Y + Z rho Z).
    Bloch vector contracts as r -> (1 - 4p/3) r.  Threshold p = 3/4."""
    return (1 - p) * rho + (p / 3) * (X @ rho @ X + Y @ rho @ Y + Z @ rho @ Z)


def dephasing(rho, p):
    """Local dephasing (phase-flip) channel:
        (1-p) rho + p Z rho Z.
    Bloch vector: (x,y) -> (1-2p)(x,y),  z -> z.  Threshold p = 1/2."""
    return (1 - p) * rho + p * (Z @ rho @ Z)
