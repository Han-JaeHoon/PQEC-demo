"""
Numerical verification of the PQEC paper (arXiv:2603.11568).

Run:  python verify_pqec.py
Produces console checks of Eqs. (5)-(10), (34), (7) and saves two figures:
  pqec_purification_map.png   - Bloch-sphere radial rescaling & purity growth
  pqec_thresholds.png         - fidelity vs rounds and error thresholds
"""

import numpy as np
import matplotlib.pyplot as plt

from pqec import (
    bloch_to_rho, rho_to_bloch, purity, fidelity_pure,
    swap_gadget, purify_map, purify_rounds,
    depolarizing, dephasing,
)

np.set_printoptions(precision=4, suppress=True)
rng = np.random.default_rng(0)


def random_qubit(rmax=0.99):
    """A random mixed single-qubit state with |r| <= rmax."""
    v = rng.normal(size=3)
    v = v / np.linalg.norm(v) * rng.uniform(0.0, rmax)
    return bloch_to_rho(v)


# ===========================================================================
print("=" * 74)
print(" PART 1 -- The SWAP gadget circuit reproduces Eqs. (5), (6), (8), (9)")
print("=" * 74)

max_err = dict(P=0.0, avg=0.0, sub=0.0, rhoplus=0.0)
for _ in range(2000):
    rho = random_qubit()
    g = swap_gadget(rho)

    # Eq. 5:  P_+/- = (1 +/- Tr rho^2)/2
    tr2 = purity(rho)
    max_err["P"] = max(max_err["P"], abs(g["P_plus"] - 0.5 * (1 + tr2)))

    # Eq. 6:  rho_+ = (rho + rho^2)/(2 P_+)
    rho_plus_pred = (rho + rho @ rho) / (2 * g["P_plus"])
    max_err["rhoplus"] = max(max_err["rhoplus"],
                             np.abs(g["rho_plus"] - rho_plus_pred).max())

    # Eq. 8:  P_+ rho_+ + P_- rho_- = rho     (averaging removes purification)
    avg = g["P_plus"] * g["rho_plus"] + g["P_minus"] * g["rho_minus"]
    max_err["avg"] = max(max_err["avg"], np.abs(avg - rho).max())

    # Eq. 9:  P_+ rho_+ - P_- rho_- = rho^2   (subtracting extracts purified part)
    max_err["sub"] = max(max_err["sub"], np.abs(g["rho_sq"] - rho @ rho).max())

print(f"  Eq.5  P_+ = (1+Tr rho^2)/2      max error over 2000 states: {max_err['P']:.2e}")
print(f"  Eq.6  rho_+ = (rho+rho^2)/2P_+  max error: {max_err['rhoplus']:.2e}")
print(f"  Eq.8  P_+rho_+ + P_-rho_- = rho max error: {max_err['avg']:.2e}")
print(f"  Eq.9  P_+rho_+ - P_-rho_- = rho^2 max err: {max_err['sub']:.2e}")
print("  -> the SWAP-test circuit implements the purification primitive exactly.")

# One concrete example
rho = bloch_to_rho([0.3, 0.0, 0.5])
g = swap_gadget(rho)
print("\n  Example  r = (0.3, 0.0, 0.5),  |r| = %.4f" % np.linalg.norm([0.3, 0, 0.5]))
print("    P_+ = %.4f   P_- = %.4f" % (g["P_plus"], g["P_minus"]))
print("    purity(rho)         = %.4f" % purity(rho))
print("    purity(P(rho))      = %.4f   (higher -> purified)" % purity(g["purified"]))


# ===========================================================================
print("\n" + "=" * 74)
print(" PART 2 -- Purification map, Eq. (7) & Bloch rescaling, Eq. (34)")
print("=" * 74)

# circuit-extracted purified state  ==  rho^2/Tr rho^2  (Eq. 7)
err_map = 0.0
err_bloch = 0.0
for _ in range(2000):
    rho = random_qubit()
    g = swap_gadget(rho)
    err_map = max(err_map, np.abs(g["purified"] - purify_map(rho)).max())

    # Eq. 34:  r -> 2 r /(1 + |r|^2)
    r = rho_to_bloch(rho)
    r_pred = 2 * r / (1 + np.dot(r, r))
    r_out = rho_to_bloch(g["purified"])
    err_bloch = max(err_bloch, np.abs(r_out - r_pred).max())

print(f"  Eq.7   circuit output == rho^2/Tr rho^2   max error: {err_map:.2e}")
print(f"  Eq.34  Bloch: r -> 2r/(1+|r|^2)           max error: {err_bloch:.2e}")
print("  Fixed points of the radial map |r|=0 (mixed) and |r|=1 (pure):")
for rr in [0.0, 0.3, 0.6, 0.9, 1.0]:
    print(f"    |r| = {rr:.2f}  ->  |r'| = {2*rr/(1+rr**2):.4f}")


# ===========================================================================
print("\n" + "=" * 74)
print(" PART 3 -- Error correction: fidelity vs rounds, and thresholds")
print("=" * 74)

# One PQEC "cycle" = apply the error channel once, then ell purification rounds.
# We track F = <psi0| rho |psi0> over many cycles.  The paper's Definition 1
# calls p_th the largest physical error for which the logical error -> 0 as the
# number of rounds grows.  For local depolarizing p_th = 3/4, dephasing = 1/2.

def run_cycles(channel, p, ell, psi0, n_cycles=40):
    """Return the fidelity trajectory over n_cycles of [channel -> ell rounds]."""
    rho = np.outer(psi0, psi0.conj())
    fids = [fidelity_pure(rho, psi0)]
    for _ in range(n_cycles):
        rho = channel(rho, p)
        rho = purify_rounds(rho, ell)
        fids.append(fidelity_pure(rho, psi0))
    return np.array(fids)


# |psi0> = |+>  (eigenstate stressed by both depolarizing and dephasing)
plus = np.array([1, 1], dtype=complex) / np.sqrt(2)

print("\n  Depolarizing p=0.5, |+>: steady-state fidelity vs rounds ell")
for ell in [0, 1, 2, 3]:
    f = run_cycles(depolarizing, 0.5, ell, plus)[-1]
    print(f"    ell = {ell}:  F_inf = {f:.4f}")

print("\n  Dephasing p=0.3, |+>: steady-state fidelity vs rounds ell")
for ell in [0, 1, 2, 3]:
    f = run_cycles(dephasing, 0.3, ell, plus)[-1]
    print(f"    ell = {ell}:  F_inf = {f:.4f}")


# ---- Figure 1: purification map (Bloch rescaling + purity growth) ----------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

# (a) radial map r -> r'
rr = np.linspace(0, 1, 400)
axes[0].plot(rr, 2 * rr / (1 + rr ** 2), lw=2.2, label=r"$r'=2r/(1+r^2)$ (Eq. 34)")
axes[0].plot(rr, rr, "k--", lw=1, label="identity")
axes[0].scatter([0, 1], [0, 1], c="crimson", zorder=5, s=45, label="fixed points")
axes[0].set_xlabel("Bloch radius  |r|  (before)")
axes[0].set_ylabel("Bloch radius  |r'|  (after one round)")
axes[0].set_title("(a) Purification = radial rescaling toward the pure shell")
axes[0].legend(frameon=False)
axes[0].set_aspect("equal")

# (b) purity vs rounds, starting from several mixed states
for r0 in [0.2, 0.4, 0.6, 0.8]:
    rho = bloch_to_rho([r0, 0, 0])
    pur = [purity(rho)]
    for _ in range(8):
        rho = purify_map(rho)
        pur.append(purity(rho))
    axes[1].plot(range(9), pur, "-o", ms=4, label=f"|r0| = {r0}")
axes[1].axhline(1.0, color="k", ls=":", lw=1)
axes[1].set_xlabel("purification rounds  $\\ell$")
axes[1].set_ylabel(r"purity  Tr$[\rho^2]$")
axes[1].set_title("(b) Repeated rounds drive purity -> 1")
axes[1].legend(frameon=False)

fig.tight_layout()
fig.savefig("pqec_purification_map.png", dpi=140)
print("\n  saved  pqec_purification_map.png")


# ---- Figure 2: fidelity vs rounds and error thresholds ---------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))

# (a) depolarizing: fidelity vs cycles, p=0.5, various ell   (cf. Fig. 5a)
ax = axes[0, 0]
for ell in [0, 1, 2, 3]:
    f = run_cycles(depolarizing, 0.5, ell, plus)
    style = "r:o" if ell == 0 else "-o"
    ax.plot(range(len(f)), f, style, ms=3,
            label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
ax.set_xlabel("PQEC cycles  t"); ax.set_ylabel("fidelity  F")
ax.set_title("(a) Depolarizing p=0.5, |+>"); ax.legend(frameon=False); ax.set_ylim(0.4, 1.02)

# (b) dephasing: fidelity vs cycles, p=0.3, various ell
ax = axes[0, 1]
for ell in [0, 1, 2, 3]:
    f = run_cycles(dephasing, 0.3, ell, plus)
    style = "r:o" if ell == 0 else "-o"
    ax.plot(range(len(f)), f, style, ms=3,
            label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
ax.set_xlabel("PQEC cycles  t"); ax.set_ylabel("fidelity  F")
ax.set_title("(b) Dephasing p=0.3, |+>"); ax.legend(frameon=False); ax.set_ylim(0.4, 1.02)

# (c) DEPOLARIZING threshold: F after one cycle vs p, various ell -> cross at 3/4
ax = axes[1, 0]
ps = np.linspace(0.01, 0.99, 99)
for ell in [0, 1, 2, 3, 5]:
    F = [run_cycles(depolarizing, p, ell, plus, n_cycles=1)[-1] for p in ps]
    ax.plot(ps, F, ("r:" if ell == 0 else "-"),
            label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
ax.axvline(0.75, color="k", ls="--", lw=1.2)
ax.text(0.75, 0.05, " $p_{th}=3/4$", fontsize=10)
ax.set_xlabel("physical error rate  p"); ax.set_ylabel("fidelity after 1 cycle")
ax.set_title("(c) Depolarizing threshold  p=3/4"); ax.legend(frameon=False)

# (d) DEPHASING threshold: F after one cycle vs p, various ell -> cross at 1/2
ax = axes[1, 1]
for ell in [0, 1, 2, 3, 5]:
    F = [run_cycles(dephasing, p, ell, plus, n_cycles=1)[-1] for p in ps]
    ax.plot(ps, F, ("r:" if ell == 0 else "-"),
            label=("no QEC" if ell == 0 else f"$\\ell={ell}$"))
ax.axvline(0.5, color="k", ls="--", lw=1.2)
ax.text(0.5, 0.05, " $p_{th}=1/2$", fontsize=10)
ax.set_xlabel("physical error rate  p"); ax.set_ylabel("fidelity after 1 cycle")
ax.set_title("(d) Dephasing threshold  p=1/2"); ax.legend(frameon=False)

fig.tight_layout()
fig.savefig("pqec_thresholds.png", dpi=140)
print("  saved  pqec_thresholds.png")

# Numerically locate the crossing points as a check of the thresholds.
def crossing(channel):
    ps_fine = np.linspace(0.3, 0.95, 651)
    f1 = np.array([run_cycles(channel, p, 1, plus, 1)[-1] for p in ps_fine])
    f3 = np.array([run_cycles(channel, p, 3, plus, 1)[-1] for p in ps_fine])
    i = np.argmin(np.abs(f1 - f3))
    return ps_fine[i]

print("\n  measured crossing (ell=1 vs ell=3 curves):")
print(f"    depolarizing p_th ~= {crossing(depolarizing):.3f}   (paper: 0.75)")
print(f"    dephasing    p_th ~= {crossing(dephasing):.3f}   (paper: 0.50)")
print("\nDone.")
