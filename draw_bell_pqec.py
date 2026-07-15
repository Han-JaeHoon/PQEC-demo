"""
Draw the full Bell-state + PQEC quantum circuit (one purification round).

One round purifies two noisy Bell copies, so the complete circuit is:

    two noisy Bell factories  (each: H, CNOT, depolarizing on both qubits)
    -> a 2-qubit SWAP gadget   (H, two parallel Fredkin gates, H, measure)

The controlled-SWAP of the two 2-qubit registers is M=2 parallel Fredkin gates.
Reading the ancilla and post-processing (bell_pqec.py) yields rho^2/Tr[rho^2].

Produces:
  bell_pqec_circuit.png
and prints the text drawing to the console.
"""

import matplotlib.pyplot as plt
import pennylane as qml

# Wires:
#   0 = purification ancilla
#   1,2 = Bell copy A (kept, purified output)     3,4 = Bell copy B (discarded)
P_ANC = 0
REG_A = [1, 2]
REG_B = [3, 4]
N_WIRES = 5
P_NOISE = 0.30


def noisy_bell(q0, q1):
    """Noisy |Phi+> factory on qubits (q0, q1): H, CNOT, local depolarizing."""
    qml.Hadamard(q0)
    qml.CNOT(wires=[q0, q1])
    qml.DepolarizingChannel(P_NOISE, wires=q0)
    qml.DepolarizingChannel(P_NOISE, wires=q1)


def bell_pqec():
    noisy_bell(*REG_A)          # copy A
    noisy_bell(*REG_B)          # copy B
    qml.Barrier(wires=range(N_WIRES))
    # 2-qubit SWAP gadget = two parallel Fredkin gates + ancilla interference.
    qml.Hadamard(P_ANC)
    qml.ctrl(qml.SWAP, control=P_ANC)(wires=[REG_A[0], REG_B[0]])
    qml.ctrl(qml.SWAP, control=P_ANC)(wires=[REG_A[1], REG_B[1]])
    qml.Hadamard(P_ANC)
    qml.measure(P_ANC)


dev = qml.device("default.mixed", wires=N_WIRES)
qnode = qml.QNode(bell_pqec, dev)

title = "Noisy Bell factory x2 + one PQEC purification round (M=2 SWAP gadget)"
print("\n" + "=" * 70)
print(" " + title)
print("=" * 70)
print(qml.draw(qnode, wire_order=list(range(N_WIRES)), show_all_wires=True)())

fig, ax = qml.draw_mpl(
    qnode, wire_order=list(range(N_WIRES)), show_all_wires=True,
    wire_options={"color": "#333"}, style="pennylane",
)()
ax.set_title(title, fontsize=12)
fig.savefig("bell_pqec_circuit.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("\n  saved  bell_pqec_circuit.png")
print("\nDone.")
