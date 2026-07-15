"""
Draw the full Deutsch + PQEC quantum circuit.

One purification round needs two identical noisy copies of the Deutsch output,
so the complete circuit is:

    two independent runs of the Deutsch algorithm  (each: X, H, H, U_f, H, noise)
    -> a SWAP gadget between the two query qubits   (H, CSWAP, H, measure)

Reading the purification-ancilla outcome and post-processing (pqec.py) yields
the purified query state rho^2/Tr[rho^2].  Here we draw the M=1 case with the
balanced oracle f(x)=x (CNOT from query to ancilla).

Produces:
  deutsch_pqec_circuit.png   - full Deutsch + one-round PQEC circuit
and prints the text drawing to the console.
"""

import matplotlib.pyplot as plt
import pennylane as qml

# Wires:
#   0 = purification ancilla
#   1 = query A  (kept, purified output)      2 = Deutsch ancilla A
#   3 = query B  (second copy, discarded)     4 = Deutsch ancilla B
P_ANC, QA, DA, QB, DB = 0, 1, 2, 3, 4
N_WIRES = 5
P_NOISE = 0.30


def deutsch(query, danc, oracle):
    """One run of the Deutsch algorithm with depolarizing noise on the output.
    `oracle` acts on (query, danc)."""
    qml.PauliX(wires=danc)          # ancilla -> |1>
    qml.Hadamard(wires=query)       # query   -> |+>
    qml.Hadamard(wires=danc)        # ancilla -> |->
    oracle(query, danc)             # apply U_f
    qml.Hadamard(wires=query)       # final interference
    qml.DepolarizingChannel(P_NOISE, wires=query)   # output noise


def oracle_balanced_id(query, danc):
    """f(x) = x  ->  balanced.  CNOT from query to ancilla."""
    qml.CNOT(wires=[query, danc])


def deutsch_pqec():
    # Two identical noisy copies of the Deutsch output.
    deutsch(QA, DA, oracle_balanced_id)
    deutsch(QB, DB, oracle_balanced_id)
    qml.Barrier(wires=range(N_WIRES))
    # SWAP-gadget purification round on the two query qubits.
    qml.Hadamard(P_ANC)
    qml.ctrl(qml.SWAP, control=P_ANC)(wires=[QA, QB])   # Fredkin
    qml.Hadamard(P_ANC)
    qml.measure(P_ANC)


dev = qml.device("default.mixed", wires=N_WIRES)
qnode = qml.QNode(deutsch_pqec, dev)

title = "Deutsch algorithm + one PQEC purification round (balanced f=x)"
print("\n" + "=" * 70)
print(" " + title)
print("=" * 70)
print(qml.draw(qnode, wire_order=list(range(N_WIRES)), show_all_wires=True)())

fig, ax = qml.draw_mpl(
    qnode, wire_order=list(range(N_WIRES)), show_all_wires=True,
    wire_options={"color": "#333"}, style="pennylane",
)()
ax.set_title(title, fontsize=12)
fig.savefig("deutsch_pqec_circuit.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print("\n  saved  deutsch_pqec_circuit.png")
print("\nDone.")
