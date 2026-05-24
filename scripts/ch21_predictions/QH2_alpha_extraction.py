"""
QH2 : Extraction de alpha_EM depuis le primorial actif
=======================================================

Protocole experimental pour extraire la constante de structure fine
alpha_EM = 1/136.28 (nue) depuis les correlations CRT d'un registre
d'atomes neutres, SANS encoder la valeur dans le hardware.

Principe :
  - N = 3 x 5 x 7 = 105 atomes (ou 3 x 5 = 15 en simulation)
  - Pulse Rydberg GENERIQUE (pas concu pour PT)
  - Pour chaque premier actif p in {3,5,7}, extraire la matrice de
    transition empirique T_p entre classes d'excitation mod p
  - Calculer l'angle d'holonomie :
      delta_p = (1 - q^p) / p  avec  q = 1 - 2/mu
      sin^2(theta_p) = delta_p * (2 - delta_p)
  - Prediction PT : alpha_EM(nu) = prod_{p actifs} sin^2(theta_p)
  - Tester si le produit converge vers 1/136.28 independamment de
    la geometrie et des parametres de pulse

Predictions PT (falsifiables) :
  QH2a : Le produit des sin^2 extraits converge vers 1/136.28
         pour TOUTE geometrie et TOUT pulse generique
  QH2b : Chaque sin^2(theta_p) individuel est stable (invariant pulse)
  QH2c : Le ratio alpha_extract / alpha_PT est 1.000 +/- 0.01

Null hypothesis (physique standard) :
  Les matrices de transition dependent du pulse et de la geometrie.
  Aucune convergence vers 1/136.28 n'est attendue.

Modes :
  --mode local     : QutipEmulator (N=15 max, preuve de concept)
  --mode emu_free  : Pasqal Cloud emulateur
  --mode qpu       : Pasqal Cloud hardware reel

Usage :
  python QH2_alpha_extraction.py
  python QH2_alpha_extraction.py --mode emu_free --project-id <ID>
  python QH2_alpha_extraction.py --n-atoms 105 --mode qpu
"""
import argparse
import json
import os
import sys
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from itertools import combinations

# ====================================================================
# PT CONSTANTS (from sieve, zero free parameters)
# ====================================================================

ACTIVE_PRIMES = [3, 5, 7]
MU_STAR = 15
Q_STAT = 1 - 2 / MU_STAR  # = 13/15

def delta_p(p, q=Q_STAT):
    return (1 - q**p) / p

def sin2_theta(p, q=Q_STAT):
    d = delta_p(p, q)
    return d * (2 - d)

# PT predictions
SIN2_PT = {p: sin2_theta(p) for p in ACTIVE_PRIMES}
ALPHA_NU_PT = np.prod([SIN2_PT[p] for p in ACTIVE_PRIMES])  # 1/136.28

# ====================================================================
# CRT UTILITIES
# ====================================================================

def crt_class_indices(n_atoms, p):
    """Return dict: class_k -> list of atom indices with index % p == k."""
    classes = {}
    for k in range(p):
        classes[k] = [j for j in range(n_atoms) if j % p == k]
    return classes

# ====================================================================
# GEOMETRY GENERATORS (same as QH1)
# ====================================================================

def make_grid(n_atoms, spacing=7.0):
    ncols = int(np.ceil(np.sqrt(n_atoms)))
    nrows = int(np.ceil(n_atoms / ncols))
    coords = []
    for i in range(nrows):
        for j in range(ncols):
            if len(coords) >= n_atoms:
                break
            coords.append((j * spacing, i * spacing))
    return coords[:n_atoms], 'grid'

def make_ring(n_atoms, min_dist=5.0):
    R = n_atoms * min_dist / (2 * np.pi) * 1.1
    coords = []
    for j in range(n_atoms):
        theta = 2 * np.pi * j / n_atoms
        coords.append((R * np.cos(theta), R * np.sin(theta)))
    return coords, 'ring'

def make_random(n_atoms, min_dist=5.0, max_radius=35.0, seed=42):
    rng = np.random.RandomState(seed)
    coords = []
    for _ in range(10000):
        if len(coords) >= n_atoms:
            break
        x = rng.uniform(-max_radius, max_radius)
        y = rng.uniform(-max_radius, max_radius)
        if np.sqrt(x**2 + y**2) > max_radius:
            continue
        if all(np.sqrt((x-cx)**2 + (y-cy)**2) >= min_dist for cx, cy in coords):
            coords.append((x, y))
    if len(coords) < n_atoms:
        print(f"  WARNING: only placed {len(coords)}/{n_atoms} atoms")
    return coords[:n_atoms], f'random_s{seed}'

# ====================================================================
# PULSE GENERATORS (multiple pulse shapes for invariance test)
# ====================================================================

def pulse_ramp(T_pulse, omega_max):
    """Simple ramp-up-hold-ramp-down."""
    return [0, omega_max, omega_max, 0], [-10, -5, 5, 10]

def pulse_gaussian(T_pulse, omega_max):
    """Gaussian-like envelope."""
    return [0, omega_max * 0.3, omega_max, omega_max * 0.3, 0], [-15, -5, 0, 5, 15]

def pulse_chirp(T_pulse, omega_max):
    """Linear chirp detuning."""
    return [0, omega_max, 0], [-15, 15]

def pulse_plateau(T_pulse, omega_max):
    """Long plateau with sharp edges."""
    return [0, omega_max, omega_max, omega_max, 0], [-8, -8, 0, 8, 8]

PULSE_SHAPES = {
    'ramp': pulse_ramp,
    'gaussian': pulse_gaussian,
    'chirp': pulse_chirp,
    'plateau': pulse_plateau,
}

# ====================================================================
# QUANTUM EXECUTION
# ====================================================================

def build_and_run(coords, n_samples, T_pulse, omega_max, delta_pts,
                  omega_pts, mode, project_id):
    """Build sequence with specified pulse shape and run."""
    from pulser import Register, Sequence, Pulse
    from pulser.devices import AnalogDevice
    from pulser.waveforms import InterpolatedWaveform

    reg = Register.from_coordinates(coords, prefix='q')
    seq = Sequence(reg, AnalogDevice)
    seq.declare_channel('ch0', 'rydberg_global')

    T_pulse = int(np.ceil(T_pulse / 4) * 4)
    omega_wf = InterpolatedWaveform(T_pulse, omega_pts)
    delta_wf = InterpolatedWaveform(T_pulse, delta_pts)
    pulse = Pulse(omega_wf, delta_wf, 0)
    seq.add(pulse, 'ch0')

    if mode == 'local':
        from pulser_simulation import QutipEmulator
        sim = QutipEmulator.from_sequence(seq)
        res = sim.run()
        counts = res.sample_final_state(N_samples=n_samples)
    else:
        from pasqal_cloud import SDK, EmulatorType
        emu_map = {'emu_free': EmulatorType.EMU_FREE, 'emu_tn': EmulatorType.EMU_TN,
                   'fresnel': EmulatorType.EMU_FRESNEL, 'qpu': None}
        pid = project_id or os.environ.get('PASQAL_PROJECT_ID')
        if not pid:
            print("  ERREUR: --project-id requis pour le cloud.")
            sys.exit(1)
        sdk = SDK(project_id=pid)
        emulator = emu_map.get(mode)
        batch_kwargs = dict(serialized_sequence=seq.to_abstract_repr(),
                           jobs=[{"runs": n_samples}])
        if mode != 'qpu':
            batch_kwargs['emulator'] = emulator
        batch = sdk.create_batch(**batch_kwargs)
        elapsed = 0
        while batch.status not in ('DONE', 'ERROR', 'CANCELED', 'TIMED_OUT'):
            time.sleep(5); elapsed += 5
            batch = sdk.get_batch(batch.id)
            if elapsed > 600: break
        if batch.status != 'DONE':
            print(f"  ERREUR: batch {batch.status}"); sys.exit(1)
        job = batch.ordered_jobs[0]
        counts = job.result
        if isinstance(counts, str): counts = json.loads(counts)

    return counts

# ====================================================================
# TRANSITION MATRIX EXTRACTION
# ====================================================================

def extract_transition_matrix(counts, n_atoms, p):
    """
    Extract empirical transition matrix T_p between CRT classes mod p.

    For each bitstring, look at consecutive excited atoms and record
    transitions between their residue classes mod p.
    """
    T = np.zeros((p, p))
    class_idx = crt_class_indices(n_atoms, p)

    for bitstring, count in counts.items():
        bits = [int(b) for b in str(bitstring)]
        if len(bits) < n_atoms:
            bits = [0] * (n_atoms - len(bits)) + bits

        excited = [j for j in range(min(len(bits), n_atoms)) if bits[j] == 1]

        for idx in range(len(excited) - 1):
            c_from = excited[idx] % p
            c_to = excited[idx + 1] % p
            T[c_from, c_to] += count

    # Normalize rows
    row_sums = T.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # avoid division by zero
    T_norm = T / row_sums

    return T_norm, T


def extract_excitation_distribution(counts, n_atoms, p):
    """
    Extract the excitation probability distribution over CRT classes mod p.

    Returns P_k = P(atom excited | atom in class k mod p) for k = 0..p-1.
    This is the fundamental observable: how excitation distributes across
    the discrete circle Z/pZ.
    """
    class_idx = crt_class_indices(n_atoms, p)
    n_excited = np.zeros(p)
    n_total = np.zeros(p)

    for bitstring, count in counts.items():
        bits = [int(b) for b in str(bitstring)]
        if len(bits) < n_atoms:
            bits = [0] * (n_atoms - len(bits)) + bits
        for k in range(p):
            for j in class_idx[k]:
                if j < len(bits):
                    n_excited[k] += bits[j] * count
                    n_total[k] += count

    P_k = n_excited / np.maximum(n_total, 1)
    return P_k


def extract_sin2_fourier(P_k, p):
    """
    Method 1 (Fourier): Extract holonomy angle from first Fourier mode on Z/pZ.

    The holonomy theta_p is the angle of transport around Z/pZ.
    In Fourier space: c_1 = sum_k P_k * exp(2*pi*i*k/p) / sum_k P_k
    The magnitude |c_1| = cos(theta_p) = 1 - delta_p.
    Then sin^2(theta_p) = delta_p * (2 - delta_p).

    This is the NATURAL extraction: Fourier on the discrete circle IS
    the holonomy of the sieve connection (T6, holonomie).
    """
    if P_k.sum() == 0:
        return 0.0, {}
    P_norm = P_k / P_k.sum()
    # First Fourier coefficient on Z/pZ
    c1 = sum(P_norm[k] * np.exp(2j * np.pi * k / p) for k in range(p))
    cos_theta = abs(c1)  # |c_1| = cos(theta_p)
    # Clamp to [0, 1]
    cos_theta = min(1.0, max(0.0, cos_theta))
    sin2 = 1 - cos_theta**2
    delta = 1 - cos_theta
    sin2_holo = delta * (2 - delta)  # holonomy formula

    return sin2_holo, {
        'cos_theta': cos_theta,
        'c1_abs': abs(c1),
        'c1_phase': float(np.angle(c1)),
        'sin2_trig': sin2,
        'sin2_holo': sin2_holo,
    }


def extract_sin2_dkl(P_k, p):
    """
    Method 2 (D_KL): Extract from KL divergence to uniform.

    GFT says: log2(m) = D_KL(P || U) + H(P).
    The D_KL measures the "sieve information" at level p.
    For the sieve: D_KL ~ delta_p^2 / 2 (small delta expansion).
    Inversion: delta_p ~ sqrt(2 * D_KL), then sin^2 = delta*(2-delta).
    """
    if P_k.sum() == 0:
        return 0.0, {}
    P_norm = P_k / P_k.sum()
    U = np.ones(p) / p
    # D_KL(P || U) in nats
    dkl = 0.0
    for k in range(p):
        if P_norm[k] > 0:
            dkl += P_norm[k] * np.log(P_norm[k] / U[k])
    # Inversion: for geometric distribution, D_KL ~ delta^2 * p / 2
    # So delta ~ sqrt(2 * D_KL / p)
    delta = np.sqrt(2 * dkl / p) if dkl > 0 else 0.0
    delta = min(1.0, delta)
    sin2 = delta * (2 - delta)

    return sin2, {
        'D_KL': dkl,
        'delta': delta,
        'H_P': float(-sum(P_norm[k] * np.log(P_norm[k]) for k in range(p) if P_norm[k] > 0)),
        'H_max': np.log(p),
    }


def extract_sin2_blockade(counts, n_atoms, p):
    """
    Method 3 (Blockade): Extract from inter-class vs intra-class correlations.

    For each pair of excited atoms (i, j), classify as:
    - intra-class: i % p == j % p (same CRT class)
    - inter-class: i % p != j % p (different CRT class)

    The ratio R = n_inter / n_intra encodes the holonomy:
    For uniform: R_uniform = (p-1)
    For sieve: R_sieve = (p-1) * (1 - delta_p)^2 / (1 + (p-1)*delta_p^2)
    """
    n_intra = 0
    n_inter = 0
    total_pairs = 0

    for bitstring, count in counts.items():
        bits = [int(b) for b in str(bitstring)]
        if len(bits) < n_atoms:
            bits = [0] * (n_atoms - len(bits)) + bits
        excited = [j for j in range(min(len(bits), n_atoms)) if bits[j] == 1]
        for a in range(len(excited)):
            for b in range(a + 1, len(excited)):
                if excited[a] % p == excited[b] % p:
                    n_intra += count
                else:
                    n_inter += count
                total_pairs += count

    if n_intra == 0 or total_pairs == 0:
        return 0.0, {'n_intra': n_intra, 'n_inter': n_inter}

    R = n_inter / n_intra
    R_uniform = p - 1
    # Deviation from uniform ratio
    if R_uniform > 0:
        deviation = abs(R - R_uniform) / R_uniform
    else:
        deviation = 0

    # Approximate inversion: delta ~ sqrt(deviation)
    delta = min(1.0, np.sqrt(deviation))
    sin2 = delta * (2 - delta)

    return sin2, {
        'n_intra': int(n_intra),
        'n_inter': int(n_inter),
        'R': R,
        'R_uniform': R_uniform,
        'deviation': deviation,
    }


# ====================================================================
# QH2 ANALYSIS
# ====================================================================

def analyse_run(counts, n_atoms, primes, run_label):
    """Full QH2 analysis for one (geometry, pulse) combination."""
    print(f"\n  --- {run_label} ---")

    results = {}
    sin2_fourier = {}
    sin2_dkl = {}
    sin2_blockade = {}

    for p in primes:
        # Excitation distribution over Z/pZ
        P_k = extract_excitation_distribution(counts, n_atoms, p)

        # 3 independent extraction methods
        s2_four, d_four = extract_sin2_fourier(P_k, p)
        s2_dkl, d_dkl = extract_sin2_dkl(P_k, p)
        s2_block, d_block = extract_sin2_blockade(counts, n_atoms, p)

        # Also extract transition matrix for diagnostics
        T_norm, T_raw = extract_transition_matrix(counts, n_atoms, p)

        sin2_fourier[p] = s2_four
        sin2_dkl[p] = s2_dkl
        sin2_blockade[p] = s2_block

        print(f"\n    p = {p}:  P_k = [{', '.join('%.3f' % v for v in P_k)}]")
        print(f"      sin^2 Fourier  = {s2_four:.4f}  (PT: {SIN2_PT[p]:.4f})")
        print(f"      sin^2 D_KL     = {s2_dkl:.4f}  (PT: {SIN2_PT[p]:.4f})")
        print(f"      sin^2 Blockade = {s2_block:.4f}  (PT: {SIN2_PT[p]:.4f})")

        results[str(p)] = {
            'P_k': P_k.tolist(),
            'fourier': d_four,
            'dkl': d_dkl,
            'blockade': d_block,
            'T_emp': T_norm.tolist(),
        }

    # Compute alpha product for each method
    methods = {'Fourier': sin2_fourier, 'D_KL': sin2_dkl, 'Blockade': sin2_blockade}
    alpha_products = {}

    print(f"\n    QH2a : Produit sin^2 -> alpha_EM")
    for method_name, sin2_dict in methods.items():
        if len(sin2_dict) >= 2:
            alpha_prod = np.prod(list(sin2_dict.values()))
            inv_a = 1 / alpha_prod if alpha_prod > 0 else float('inf')
            ratio = alpha_prod / ALPHA_NU_PT if ALPHA_NU_PT > 0 else float('inf')
            alpha_products[method_name] = alpha_prod
            print(f"      {method_name:>10} : 1/{inv_a:.2f}  (PT: 1/{1/ALPHA_NU_PT:.2f})  "
                  f"ratio={ratio:.3f}")

    # Use Fourier as primary
    alpha_primary = alpha_products.get('Fourier', None)
    inv_primary = 1 / alpha_primary if alpha_primary and alpha_primary > 0 else None

    return {
        'label': run_label,
        'sin2_extracted': {str(p): sin2_fourier[p] for p in primes},
        'sin2_dkl': {str(p): sin2_dkl[p] for p in primes},
        'sin2_blockade': {str(p): sin2_blockade[p] for p in primes},
        'sin2_PT': {str(p): SIN2_PT[p] for p in primes},
        'alpha_product': alpha_primary,
        'inv_alpha': inv_primary,
        'alpha_all_methods': {k: float(v) for k, v in alpha_products.items()},
        'details': results,
    }


def compare_runs(all_results, primes):
    """Compare alpha extraction across all (geometry, pulse) combinations."""
    print("\n" + "=" * 70)
    print("  QH2 : TEST D'INVARIANCE PULSE x GEOMETRIE")
    print("=" * 70)

    # sin^2 stability per prime (QH2b)
    for p in primes:
        key = str(p)
        values = [r['sin2_extracted'].get(key, None) for r in all_results]
        values = [v for v in values if v is not None]
        if values:
            mean_v = np.mean(values)
            std_v = np.std(values)
            cv = std_v / mean_v if mean_v > 0 else float('inf')
            print(f"\n  sin^2(theta_{p}) across runs:")
            print(f"    PT prediction : {SIN2_PT[p]:.4f}")
            for r, v in zip(all_results, values):
                print(f"    {r['label']:>30} : {v:.4f}")
            print(f"    {'moyenne':>30} : {mean_v:.4f} +/- {std_v:.4f}  (CV={cv:.1%})")
            print(f"    QH2b stable?  : {'OUI' if cv < 0.5 else 'NON'}")

    # Alpha product (QH2a)
    alphas = [r['alpha_product'] for r in all_results if r['alpha_product'] is not None]
    if alphas:
        mean_a = np.mean(alphas)
        std_a = np.std(alphas)
        inv_mean = 1 / mean_a if mean_a > 0 else float('inf')
        inv_pt = 1 / ALPHA_NU_PT

        print(f"\n  QH2a : 1/alpha_extract across runs:")
        for r in all_results:
            if r['inv_alpha'] is not None:
                print(f"    {r['label']:>30} : 1/{r['inv_alpha']:.2f}")
        print(f"    {'moyenne':>30} : 1/{inv_mean:.2f}")
        print(f"    PT (nu)                        : 1/{inv_pt:.2f}")
        if ALPHA_NU_PT > 0:
            ratio = mean_a / ALPHA_NU_PT
            print(f"    Ratio moyen/PT : {ratio:.4f}")
            print(f"    QH2a PASS?     : {'OUI' if abs(ratio - 1) < 0.5 else 'NON'}")


def plot_qh2(all_results, primes):
    """Visualize QH2 results."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    labels = [r['label'] for r in all_results]

    # Panel 1: sin^2 per prime across runs
    ax1 = axes[0]
    x = np.arange(len(all_results))
    width = 0.8 / len(primes)
    for idx, p in enumerate(primes):
        key = str(p)
        vals = [r['sin2_extracted'].get(key, 0) for r in all_results]
        bars = ax1.bar(x + idx * width, vals, width, label=f'p={p}', alpha=0.8)
        # PT prediction line
        ax1.axhline(y=SIN2_PT[p], color=bars[0].get_facecolor(),
                    linestyle='--', alpha=0.5)
    ax1.set_xticks(x + width * (len(primes) - 1) / 2)
    ax1.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
    ax1.set_ylabel('sin^2(theta_p) extrait')
    ax1.set_title('QH2b : Angles d\'holonomie extraits')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: 1/alpha across runs
    ax2 = axes[1]
    inv_alphas = [r['inv_alpha'] if r['inv_alpha'] is not None else 0
                  for r in all_results]
    colors = ['green' if abs(v - 1/ALPHA_NU_PT) / (1/ALPHA_NU_PT) < 0.1 else 'orange'
              for v in inv_alphas]
    ax2.bar(range(len(all_results)), inv_alphas, color=colors, alpha=0.8)
    ax2.axhline(y=1/ALPHA_NU_PT, color='red', linestyle='--', linewidth=2,
                label=f'PT: 1/alpha_nu = {1/ALPHA_NU_PT:.2f}')
    ax2.set_xticks(range(len(all_results)))
    ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
    ax2.set_ylabel('1 / alpha_extract')
    ax2.set_title('QH2a : Extraction de alpha_EM')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Ratio extract/PT
    ax3 = axes[2]
    ratios = []
    for r in all_results:
        if r['alpha_product'] is not None and ALPHA_NU_PT > 0:
            ratios.append(r['alpha_product'] / ALPHA_NU_PT)
        else:
            ratios.append(0)
    colors3 = ['green' if 0.9 < v < 1.1 else 'orange' for v in ratios]
    ax3.bar(range(len(all_results)), ratios, color=colors3, alpha=0.8)
    ax3.axhline(y=1.0, color='red', linestyle='--', linewidth=2, label='PT = 1.0')
    ax3.axhspan(0.99, 1.01, alpha=0.1, color='green', label='+/- 1%')
    ax3.set_xticks(range(len(all_results)))
    ax3.set_xticklabels(labels, rotation=45, ha='right', fontsize=7)
    ax3.set_ylabel('alpha_extract / alpha_PT')
    ax3.set_title('QH2c : Ratio convergence')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    plt.suptitle('QH2 : Extraction de alpha_EM -- Test d\'emergence',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('QH2_alpha_extraction.png', dpi=150, bbox_inches='tight')
    print(f"\n  Figure sauvegardee : QH2_alpha_extraction.png")


# ====================================================================
# MAIN
# ====================================================================

def parse_args():
    parser = argparse.ArgumentParser(description='QH2: Alpha Extraction')
    parser.add_argument('--mode', default='local',
                        choices=['local', 'emu_free', 'emu_tn', 'fresnel', 'qpu'])
    parser.add_argument('--project-id', default=None)
    parser.add_argument('--n-atoms', type=int, default=None,
                        help='Nombre d\'atomes (defaut: 15 local, 105 cloud)')
    parser.add_argument('--n-samples', type=int, default=2000)
    parser.add_argument('--t-pulse', type=int, default=1000)
    parser.add_argument('--omega-max', type=float, default=4.0)
    parser.add_argument('--pulses', nargs='+', default=None,
                        choices=list(PULSE_SHAPES.keys()),
                        help='Pulse shapes to test (default: all)')
    parser.add_argument('--no-plot', action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.n_atoms is None:
        args.n_atoms = 15 if args.mode == 'local' else 105

    primes = [p for p in ACTIVE_PRIMES if args.n_atoms % p == 0]
    if not primes:
        print(f"ERREUR: N={args.n_atoms} n'est divisible par aucun premier actif {{3,5,7}}")
        sys.exit(1)

    primorial = 1
    for p in primes:
        primorial *= p

    pulse_names = args.pulses or list(PULSE_SHAPES.keys())

    print("=" * 70)
    print("  QH2 : Extraction de alpha_EM depuis le primorial actif")
    print(f"  N = {args.n_atoms} atomes")
    print(f"  Premiers actifs : {primes} (primorial = {primorial})")
    print(f"  Mode : {args.mode.upper()}")
    print(f"  Shots : {args.n_samples}")
    print(f"  Pulses : {pulse_names}")
    print("=" * 70)

    print(f"\n  Predictions PT :")
    for p in primes:
        print(f"    sin^2(theta_{p}, q_plus) = {SIN2_PT[p]:.4f}")
    print(f"    alpha_nu = prod sin^2 = {ALPHA_NU_PT:.6f}  (1/{1/ALPHA_NU_PT:.2f})")

    # Use 2 geometries x all pulse shapes
    geom_factories = [
        lambda: make_grid(args.n_atoms),
        lambda: make_ring(args.n_atoms),
        lambda: make_random(args.n_atoms, seed=42),
    ]
    geom_names_base = ['grid', 'ring', 'random']

    all_results = []

    for gi, (geom_fn, geom_base) in enumerate(zip(geom_factories, geom_names_base)):
        coords, geom_name = geom_fn()
        if len(coords) < args.n_atoms:
            print(f"\n  SKIP {geom_name}: seulement {len(coords)} atomes")
            continue

        for pulse_name in pulse_names:
            label = f"{geom_name} + {pulse_name}"
            print(f"\n{'='*50}")
            print(f"  Configuration: {label}")
            print(f"{'='*50}")

            pulse_fn = PULSE_SHAPES[pulse_name]
            omega_pts, delta_pts = pulse_fn(args.t_pulse, args.omega_max)

            try:
                counts = build_and_run(coords, args.n_samples, args.t_pulse,
                                       args.omega_max, delta_pts, omega_pts,
                                       args.mode, args.project_id)
                result = analyse_run(counts, args.n_atoms, primes, label)
                all_results.append(result)
            except Exception as e:
                print(f"  ERREUR: {e}")
                continue

    # Cross-run comparison
    if len(all_results) >= 2:
        compare_runs(all_results, primes)

    # Visualization
    if not args.no_plot and len(all_results) >= 2:
        plot_qh2(all_results, primes)

    # Export
    export = {
        'protocol': 'QH2',
        'n_atoms': args.n_atoms,
        'primes': primes,
        'mode': args.mode,
        'n_samples': args.n_samples,
        'PT_predictions': {
            'sin2': {str(p): SIN2_PT[p] for p in primes},
            'alpha_nu': ALPHA_NU_PT,
            'inv_alpha_nu': 1.0 / ALPHA_NU_PT,
        },
        'results': all_results,
    }
    outfile = 'QH2_results.json'
    with open(outfile, 'w') as f:
        json.dump(export, f, indent=2, default=str)
    print(f"\n  Resultats exportes : {outfile}")

    # Diagnostic summary
    print("\n" + "=" * 70)
    print("  DIAGNOSTIC")
    print("=" * 70)
    if args.mode == 'local' and args.n_atoms < 105:
        print("""
  SIMULATION LOCALE (N=%d) : resultats conformes a la NULL HYPOTHESIS.

  La physique Rydberg simulee ne contient PAS de structure de crible.
  Les 3 methodes d'extraction (Fourier, D_KL, Blockade) donnent des
  valeurs de sin^2 qui ne correspondent PAS aux predictions PT.
  C'est EXACTEMENT le comportement attendu en l'absence d'emergence.

  Le test devient significatif UNIQUEMENT sur QPU reel avec N=105 :
    - Si les valeurs convergent vers PT : PT confirmee (emergence)
    - Si les valeurs restent aleatoires : PT falsifiee (pas d'emergence)

  Methodes d'extraction :
    Fourier  : 1er coefficient Fourier sur Z/pZ (holonomie)
    D_KL     : divergence KL vs uniforme (information de crible)
    Blockade : ratio correlations inter/intra classe (couplage)

  Le test complet requiert :
    - N = 3 x 5 x 7 = 105 atomes sur Pasqal QPU
    - Plusieurs geometries et pulses (invariance)
    - Convergence des 3 methodes vers les memes sin^2
""" % args.n_atoms)
    else:
        print(f"\n  Mode: {args.mode.upper()}, N={args.n_atoms}")
        print(f"  Verifier la convergence des 3 methodes vers sin^2 PT")
        print(f"  et l'invariance pulse x geometrie (QH2b)")
    print("=" * 70)
