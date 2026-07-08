"""
Draw the full PQEC quantum circuits (arXiv:2603.11568, Figs. 1 & 2).

Produces:
  circuit_gadget_M1.png   - elementary SWAP gadget, 1-qubit registers  (Fig. 1b)
  circuit_gadget_M2.png   - SWAP gadget, 2-qubit registers = 2 parallel Fredkins
  circuit_tree_ell2.png   - depth-efficient binary tree, ell=2 rounds   (Fig. 1a)
and prints the text drawings to the console.
"""

import matplotlib.pyplot as plt
import pennylane as qml


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------
def prep_noisy_copy(reg, theta=0.9):
    """Stand-in 'prepare a noisy copy |psi>' box on a register (list of wires).
    (In the real protocol this is the output of the algorithm + noise; here it
    is a placeholder rotation just so the circuit topology is visible.)"""
    for w in reg:
        qml.RY(theta, wires=w)


def swap_gadget(anc, regA, regB):
    """One SWAP-test purification gadget (Fig. 1b).

    anc  : ancilla wire
    regA : register kept       (list of M wires)
    regB : register discarded  (list of M wires)
    For M-qubit registers the controlled-SWAP is M parallel Fredkin gates.
    """
    qml.Hadamard(anc)
    for a, b in zip(regA, regB):
        qml.ctrl(qml.SWAP, control=anc)(wires=[a, b])   # Fredkin on qubit pair
    qml.Hadamard(anc)
    qml.measure(anc)                                     # read out the ancilla


# ---------------------------------------------------------------------------
# 1) Elementary gadget, M = 1  (wires: anc, A, B)
# ---------------------------------------------------------------------------
def gadget_M1():
    prep_noisy_copy([1]); prep_noisy_copy([2])
    swap_gadget(anc=0, regA=[1], regB=[2])


# ---------------------------------------------------------------------------
# 2) Gadget with M = 2 registers -> two parallel Fredkin gates
#    wires: 0 = ancilla, [1,2] = register A, [3,4] = register B
# ---------------------------------------------------------------------------
def gadget_M2():
    prep_noisy_copy([1, 2]); prep_noisy_copy([3, 4])
    swap_gadget(anc=0, regA=[1, 2], regB=[3, 4])


# ---------------------------------------------------------------------------
# 3) Depth-efficient binary tree, ell = 2 rounds, M = 1  (Fig. 1a)
#    4 noisy copies consumed, N = 2^ell = 4, one purified output.
#    data wires : d0=1, d1=2, d2=3, d3=4      ancillas : a=0, b=5, c=6
#      round 1 : gadget(a; d0,d1)  and  gadget(b; d2,d3)   -> keep d0, d2
#      round 2 : gadget(c; d0,d2)                          -> output on d0
# ---------------------------------------------------------------------------
def tree_ell2():
    for d in (1, 2, 3, 4):
        prep_noisy_copy([d])
    qml.Barrier(wires=range(7))
    swap_gadget(anc=0, regA=[1], regB=[2])      # merge copies 0,1 -> keep d0
    swap_gadget(anc=5, regA=[3], regB=[4])      # merge copies 2,3 -> keep d2
    qml.Barrier(wires=range(7))
    swap_gadget(anc=6, regA=[1], regB=[3])      # merge the two -> output on d0


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
configs = [
    ("gadget_M1", gadget_M1, 3,
     {0: "anc", 1: "A |psi>", 2: "B |psi>  (discard)"},
     "SWAP gadget, M=1 register  (Fig. 1b)"),
    ("gadget_M2", gadget_M2, 5,
     {0: "anc", 1: "A q0", 2: "A q1", 3: "B q0 (disc)", 4: "B q1 (disc)"},
     "SWAP gadget, M=2 registers = 2 parallel Fredkin gates"),
    ("tree_ell2", tree_ell2, 7,
     {0: "anc r1a", 1: "copy0 -> out", 2: "copy1 (disc)", 3: "copy2",
      4: "copy3 (disc)", 5: "anc r1b", 6: "anc r2"},
     "Depth-efficient binary tree, ell=2 (N=4 copies)  (Fig. 1a)"),
]

for name, fn, nwires, labels, title in configs:
    dev = qml.device("default.qubit", wires=nwires)
    qnode = qml.QNode(fn, dev)

    print("\n" + "=" * 70)
    print(" " + title)
    print("=" * 70)
    print(qml.draw(qnode, wire_order=list(range(nwires)), show_all_wires=True)())

    fig, ax = qml.draw_mpl(
        qnode, wire_order=list(range(nwires)), show_all_wires=True,
        wire_options={"color": "#333"}, style="pennylane",
    )()
    # relabel wires
    ax.set_title(title, fontsize=12)
    for txt in ax.texts:
        pass
    fig.savefig(f"circuit_{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved  circuit_{name}.png")

print("\nDone.")
