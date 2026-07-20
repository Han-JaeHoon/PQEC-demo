"""
Draw the faulty PQEC circuit: the gadget with noise after every gate.

Same topology as the clean gadget (two noisy Bell copies -> SWAP test), but a
depolarizing channel (strength g) is inserted after each gate that the purifier
actually runs:
  * after the ancilla H gates
  * after each Fredkin (controlled-SWAP), on all three qubits it touches
and an optional bit-flip (r) models ancilla readout error before the Z measurement.

The input noise on the two Bell copies (strength eps) is shown as the
"Depol(eps)" boxes right after each H-CNOT preparation.

Produces:
  pqec_gate_noise_circuit.png
and prints the text drawing to the console.
"""

import matplotlib.pyplot as plt
import pennylane as qml

# wires: 0 = ancilla,  [1,2] = Bell copy A (kept),  [3,4] = Bell copy B (discarded)
P_ANC = 0
REG_A = [1, 2]
REG_B = [3, 4]
N = 5
EPS = 0.30     # input noise on the Bell copies
G = 0.05       # gate noise after every gadget gate
R = 0.05       # ancilla readout error


def noisy_bell(q0, q1):
    """Noisy |Phi+> factory: H, CNOT, then input depolarizing noise (eps)."""
    qml.Hadamard(q0)
    qml.CNOT(wires=[q0, q1])
    qml.DepolarizingChannel(EPS, wires=q0)
    qml.DepolarizingChannel(EPS, wires=q1)


def gate_noise(*wires):
    for w in wires:
        qml.DepolarizingChannel(G, wires=w)


def faulty_gadget():
    noisy_bell(*REG_A)                         # copy A  (input noise eps)
    noisy_bell(*REG_B)                         # copy B
    qml.Barrier(wires=range(N))
    qml.Hadamard(P_ANC);                 gate_noise(P_ANC)          # H + gate noise
    qml.ctrl(qml.SWAP, control=P_ANC)(wires=[REG_A[0], REG_B[0]])   # Fredkin 1
    gate_noise(P_ANC, REG_A[0], REG_B[0])
    qml.ctrl(qml.SWAP, control=P_ANC)(wires=[REG_A[1], REG_B[1]])   # Fredkin 2
    gate_noise(P_ANC, REG_A[1], REG_B[1])
    qml.Hadamard(P_ANC);                 gate_noise(P_ANC)          # H + gate noise
    qml.BitFlip(R, wires=P_ANC)                                     # readout error
    qml.measure(P_ANC)                                             # measure Z (x) O


dev = qml.device("default.mixed", wires=N)
qnode = qml.QNode(faulty_gadget, dev)

title = "Faulty PQEC gadget: Depol(g) after every gate + readout BitFlip(r)"
print("\n" + "=" * 70)
print(" " + title)
print("=" * 70)
print(qml.draw(qnode, wire_order=list(range(N)), show_all_wires=True)())

fig, ax = qml.draw_mpl(
    qnode, wire_order=list(range(N)), show_all_wires=True,
    wire_options={"color": "#333"}, style="pennylane",
)()
ax.set_title(title, fontsize=11)
fig.savefig("pqec_gate_noise_circuit.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("\n  saved  pqec_gate_noise_circuit.png")
print("\nDone.")
