"""
PQEC on the Deutsch algorithm under a BIT-FLIP noise channel.
=============================================================

Companion to deutsch_pqec.py (which uses a depolarizing channel).  Here the
output query qubit is hit by a bit-flip channel

    rho -> (1 - p) rho + p X rho X

which flips the answer bit 0 <-> 1 with probability p.  Because the Deutsch
answer is a *computational-basis* state (|0> = constant, |1> = balanced), the
noisy state is diagonal in that basis,

    rho = (1 - p) |ans><ans| + p |flip><flip|,

so the correct answer stays the *dominant* eigenvector exactly while p < 1/2.
The PQEC purification map rho -> rho^2/Tr[rho^2] then concentrates weight back
onto it -- like a coherent majority vote -- driving P(correct) -> 1.

Hence the threshold is p_th = 1/2 (not 3/4 as for depolarizing): below 1/2
purification recovers the answer, above 1/2 it amplifies the flipped (wrong)
answer.  Everything reuses the genuine SWAP-gadget circuit from pqec.py via
run_demo() in deutsch_pqec.py.

Run:  python deutsch_pqec_bitflip.py
Produces deutsch_pqec_bitflip.png.
"""

import pennylane as qml

from deutsch_pqec import run_demo

if __name__ == "__main__":
    run_demo(
        channel=qml.BitFlip,
        channel_name="bit-flip",
        p_th=0.5,
        p_th_label="1/2",
        p_example=0.30,
        outfile="deutsch_pqec_bitflip.png",
    )
