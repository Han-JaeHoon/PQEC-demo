"""
Faulty PQEC: noise on the gadget operations themselves (H, Fredkin, readout).
=============================================================================

So far noise lived only on the INPUT copies; the SWAP gadget (H, controlled-SWAP,
H, ancilla measurement) was assumed perfect.  In reality those gates are noisy
too -- especially the 3-qubit controlled-SWAP (Fredkin) -- so the purifier can
inject as much error as it removes.  This asks the fault-tolerance question:

    does PQEC still help when the correction hardware is itself faulty,
    and is there a threshold on the GATE error rate?

How to simulate it: insert noise channels **after every gate** in the gadget.
Here (a standard local model):
  * after each single-qubit gate (H on the ancilla):        DepolarizingChannel(g)
  * after each Fredkin (3-qubit gate on anc + 2 data):      DepolarizingChannel(g)
                                                             on each of the 3 wires
  * ancilla readout error before the Z measurement:         BitFlip(r)

We measure the purified observable exactly as in pqec_observable.py,
<O>_PQEC = <Z_anc (x) O> / <Z_anc (x) I> = Tr(O rho^2)/Tr(rho^2) when g=r=0, and
watch it degrade as the gadget becomes noisy.

Findings (see console + figure):
  * gate noise g degrades <O>_PQEC; beyond a gate-error threshold g* the purifier
    does WORSE than no QEC at all (it injects more error than it removes).
  * the smaller the input noise eps, the smaller g* -- a nearly-clean input has
    little to gain and is easily spoiled by a faulty purifier.
  * symmetric ancilla READOUT error cancels in the ratio (both <Z(x)O> and
    <Z(x)I> pick up the same factor 1-2r), so it does not bias <O>_PQEC.

Run:  python pqec_gate_noise.py
"""

import numpy as np
import matplotlib.pyplot as plt
import pennylane as qml

np.set_printoptions(precision=4, suppress=True)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
_TWOQ = [np.kron(a, b) for a in (I2, X, Y, Z) for b in (I2, X, Y, Z)]
PHI = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
O_FID = np.outer(PHI, PHI.conj())


def _global_depol_kraus(eps):
    K = [np.sqrt(1 - eps + eps / 16) * _TWOQ[0]]
    K += [np.sqrt(eps / 16) * P for P in _TWOQ[1:]]
    return K


def _prep_noisy_bell(q0, q1, eps):
    qml.Hadamard(q0)
    qml.CNOT(wires=[q0, q1])
    qml.QubitChannel(_global_depol_kraus(eps), wires=[q0, q1])


def _dep(wire, g):
    if g > 0:
        qml.DepolarizingChannel(g, wires=wire)


_dev = qml.device("default.mixed", wires=5)   # 0=anc, [1,2]=copy A, [3,4]=copy B


@qml.qnode(_dev)
def _gadget(eps, g, r, O):
    """One purification round with gate noise g (after every gate) and ancilla
    readout error r.  Returns (<Z (x) O>, <Z (x) I>)."""
    _prep_noisy_bell(1, 2, eps)
    _prep_noisy_bell(3, 4, eps)
    qml.Hadamard(0);                     _dep(0, g)
    qml.ctrl(qml.SWAP, control=0)(wires=[1, 3])
    for w in (0, 1, 3):                  _dep(w, g)          # Fredkin gate noise
    qml.ctrl(qml.SWAP, control=0)(wires=[2, 4])
    for w in (0, 2, 4):                  _dep(w, g)
    qml.Hadamard(0);                     _dep(0, g)
    if r > 0:
        qml.BitFlip(r, wires=0)                              # readout error
    return (qml.expval(qml.PauliZ(0) @ qml.Hermitian(O, wires=[1, 2])),
            qml.expval(qml.PauliZ(0)))


def obs_pqec(eps, g, r=0.0, O=O_FID):
    zO, zI = _gadget(eps, g, r, O)
    return float(zO) / float(zI)


def no_qec(eps):
    """<O_fid> with no purification = Tr(O rho_eps) = 1 - 3 eps/4."""
    return 1 - 3 * eps / 4


# ===========================================================================
def main():
    print("=" * 74)
    print(" PART 1 -- Gate noise g degrades the purifier; a gate-error threshold g*")
    print("=" * 74)
    eps = 0.40
    print(f"\n  input noise eps = {eps},  no-QEC <O> = {no_qec(eps):.4f}")
    print(f"  observable O = |Phi+><Phi+|,  one purification round (ell=1)\n")
    print(f"  {'g (gate)':>9} {'<O>_PQEC':>9} {'vs no-QEC':>10}   verdict")
    g_prev, v_prev = 0.0, obs_pqec(eps, 0.0)
    gstar = None
    for g in [0.0, 0.02, 0.05, 0.10, 0.12, 0.14, 0.16, 0.20, 0.30]:
        v = obs_pqec(eps, g)
        verdict = "PQEC helps" if v > no_qec(eps) else "PQEC HURTS"
        print(f"  {g:>9.2f} {v:>9.4f} {v-no_qec(eps):>+10.4f}   {verdict}")
        if gstar is None and v < no_qec(eps):
            gstar = g_prev + (g - g_prev) * (v_prev - no_qec(eps)) / (v_prev - v)
        g_prev, v_prev = g, v
    print(f"\n  -> gate-error threshold  g* ~= {gstar:.3f}"
          f"   (above it, purifying is worse than doing nothing)")

    print("\n" + "=" * 74)
    print(" PART 2 -- Ancilla readout error cancels in the ratio <Z(x)O>/<Z(x)I>")
    print("=" * 74)
    print(f"\n  {'readout r':>10} {'g=0.00':>9} {'g=0.05':>9}")
    for r in [0.0, 0.05, 0.10, 0.20]:
        print(f"  {r:>10.2f} {obs_pqec(eps,0.0,r):>9.4f} {obs_pqec(eps,0.05,r):>9.4f}")
    print("  -> <O>_PQEC is independent of r: both correlators scale by (1-2r),")
    print("     which cancels in the ratio (self-mitigating readout error).")

    # ---- Figure ----------------------------------------------------------
    gs = np.linspace(0.0, 0.30, 46)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    ax = axes[0]
    for e, c in [(0.2, "C0"), (0.4, "C1"), (0.6, "C2")]:
        vals = [obs_pqec(e, g) for g in gs]
        ax.plot(gs, vals, "-", color=c, label=f"PQEC, $\\varepsilon$={e}")
        ax.axhline(no_qec(e), color=c, ls=":", lw=1)
    ax.text(0.305, no_qec(0.2), " no-QEC", va="center", fontsize=8)
    ax.set_xlabel("gate error rate  $g$")
    ax.set_ylabel(r"$\langle O\rangle$")
    ax.set_title("(a) Faulty gadget: PQEC (solid) vs no-QEC (dotted)")
    ax.set_xlim(0, 0.33)
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    for e, c in [(0.2, "C0"), (0.4, "C1"), (0.6, "C2")]:
        adv = [obs_pqec(e, g) - no_qec(e) for g in gs]
        ax.plot(gs, adv, "-", color=c, label=f"$\\varepsilon$={e}")
        # mark the gate-error threshold (advantage crosses 0)
        adv = np.array(adv)
        s = np.where(np.diff(np.sign(adv)))[0]
        if len(s):
            i = s[0]
            gstar = gs[i] + (gs[i+1]-gs[i]) * adv[i]/(adv[i]-adv[i+1])
            ax.plot(gstar, 0, "o", color=c, ms=6)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("gate error rate  $g$")
    ax.set_ylabel(r"PQEC advantage  $\langle O\rangle_{PQEC}-\langle O\rangle_{noQEC}$")
    ax.set_title("(b) Gate-error threshold $g^*$ (advantage = 0)")
    ax.legend(frameon=False, fontsize=9)

    fig.tight_layout()
    fig.savefig("pqec_gate_noise.png", dpi=140)
    print("\n  saved  pqec_gate_noise.png")
    print("\nDone.")


if __name__ == "__main__":
    main()
