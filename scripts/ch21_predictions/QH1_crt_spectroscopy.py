"""
QH1 : Spectroscopie CRT sur le primorial actif
================================================

Protocole experimental pour tester si la structure CRT du crible
EMERGE d'un systeme quantique a atomes neutres, sans etre encodee.

Principe :
  - N = 3 x 5 x 7 = 105 atomes (ou 3 x 5 = 15 en simulation)
  - Chaque atome j a une decomposition CRT unique :
    j <-> (j mod 3, j mod 5, j mod 7)
  - Pulse Rydberg GENERIQUE (pas concu pour PT)
  - Mesure des correlations entre facteurs CRT
  - Repetition sur PLUSIEURS geometries

Predictions PT (falsifiables) :
  QH1a : L'information mutuelle entre facteurs CRT est INVARIANTE
         par changement de geometrie
  QH1b : Le ratio de couplage R_35/R_37 = sin^2(theta_5)/sin^2(theta_7)
         = 0.1940/0.1726 = 1.124
  QH1c : La transition T[1][1] = 0 (mod 3) est supprimee (signature T0)

Null hypothesis (physique standard) :
  Les correlations dependent de la geometrie et du potentiel Rydberg.
  Aucun ratio universel lie aux angles du crible n'est attendu.

Modes :
  --mode local     : QutipEmulator (N=15 max, preuve de concept)
  --mode emu_free  : Pasqal Cloud emulateur
  --mode qpu       : Pasqal Cloud hardware reel

Usage :
  python QH1_crt_spectroscopy.py
  python QH1_crt_spectroscopy.py --mode emu_free --project-id <ID>
  python QH1_crt_spectroscopy.py --n-atoms 105 --mode qpu
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
SIN2 = {p: sin2_theta(p) for p in ACTIVE_PRIMES}
RATIO_35_37 = SIN2[5] / SIN2[7]  # QH1b prediction: 1.124

# ====================================================================
# CRT UTILITIES
# ====================================================================

def crt_decompose(j, primes):
    """Decompose j into CRT residues mod each prime."""
    return tuple(j % p for p in primes)

def crt_class_indices(n_atoms, p):
    """Return dict: class_k -> list of atom indices with index % p == k."""
    classes = {}
    for k in range(p):
        classes[k] = [j for j in range(n_atoms) if j % p == k]
    return classes

# ====================================================================
# GEOMETRY GENERATORS
# ====================================================================

def make_grid(n_atoms, spacing=7.0):
    """Rectangular grid geometry."""
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
    """Circular ring geometry."""
    R = n_atoms * min_dist / (2 * np.pi) * 1.1
    coords = []
    for j in range(n_atoms):
        theta = 2 * np.pi * j / n_atoms
        coords.append((R * np.cos(theta), R * np.sin(theta)))
    return coords, 'ring'

def make_random(n_atoms, min_dist=5.0, max_radius=35.0, seed=42):
    """Random geometry with minimum distance constraint."""
    rng = np.random.RandomState(seed)
    coords = []
    max_attempts = 10000
    for _ in range(max_attempts):
        if len(coords) >= n_atoms:
            break
        x = rng.uniform(-max_radius, max_radius)
        y = rng.uniform(-max_radius, max_radius)
        if np.sqrt(x**2 + y**2) > max_radius:
            continue
        ok = True
        for cx, cy in coords:
            if np.sqrt((x-cx)**2 + (y-cy)**2) < min_dist:
                ok = False
                break
        if ok:
            coords.append((x, y))
    if len(coords) < n_atoms:
        print(f"  WARNING: only placed {len(coords)}/{n_atoms} atoms in random geometry")
    return coords[:n_atoms], f'random_s{seed}'

def make_zigzag(n_atoms, spacing=7.0, amplitude=5.0, max_radius=35.0):
    """Zigzag line geometry, centered to fit within max_radius."""
    coords = []
    for j in range(n_atoms):
        x = j * spacing
        y = amplitude * (j % 2)
        coords.append((x, y))
    # Center the array so it fits within max_radius
    cx = np.mean([c[0] for c in coords])
    cy = np.mean([c[1] for c in coords])
    coords = [(x - cx, y - cy) for x, y in coords]
    # Check if it fits; if not, reduce spacing
    max_r = max(np.sqrt(x**2 + y**2) for x, y in coords)
    if max_r > max_radius:
        scale = max_radius / max_r * 0.95
        coords = [(x * scale, y * scale) for x, y in coords]
        # Ensure min distance respected
        min_d = min(np.sqrt((coords[i][0]-coords[j][0])**2 + (coords[i][1]-coords[j][1])**2)
                    for i in range(len(coords)) for j in range(i+1, len(coords)))
        if min_d < 5.0:
            print(f"  WARNING: zigzag scaled to fit, min_dist={min_d:.1f} um < 5 um")
    return coords, 'zigzag'

# ====================================================================
# QUANTUM EXECUTION
# ====================================================================

def build_and_run(coords, n_samples, T_pulse, omega_max, mode, project_id):
    """Build sequence and run on specified backend."""
    from pulser import Register, Sequence, Pulse
    from pulser.devices import AnalogDevice
    from pulser.waveforms import InterpolatedWaveform

    reg = Register.from_coordinates(coords, prefix='q')
    seq = Sequence(reg, AnalogDevice)
    seq.declare_channel('ch0', 'rydberg_global')

    T_pulse = int(np.ceil(T_pulse / 4) * 4)
    omega_wf = InterpolatedWaveform(T_pulse, [0, omega_max, 0])
    delta_wf = InterpolatedWaveform(T_pulse, [-10, 10])
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
# QH1 ANALYSIS
# ====================================================================

def compute_crt_fractions(counts, n_atoms, primes, n_samples):
    """Compute CRT-projected excitation fractions per shot."""
    # For each bitstring, compute n_k^(p) for each prime p and class k
    fractions = {p: {k: [] for k in range(p)} for p in primes}
    class_indices = {p: crt_class_indices(n_atoms, p) for p in primes}

    for bitstring, count in counts.items():
        bits = [int(b) for b in str(bitstring)]
        if len(bits) < n_atoms:
            bits = [0] * (n_atoms - len(bits)) + bits

        for p in primes:
            for k in range(p):
                indices = class_indices[p][k]
                n_in_class = len(indices)
                if n_in_class > 0:
                    frac = sum(bits[i] for i in indices if i < len(bits)) / n_in_class
                else:
                    frac = 0.0
                for _ in range(count):
                    fractions[p][k].append(frac)

    return fractions


def compute_mutual_information(fractions, p1, p2):
    """Compute mutual information I(n^(p1); n^(p2)) via binning."""
    # Use class-averaged excitation as observable
    x = np.array([np.mean([fractions[p1][k][i] for k in range(p1)])
                   for i in range(len(fractions[p1][0]))])
    y = np.array([np.mean([fractions[p2][k][i] for k in range(p2)])
                   for i in range(len(fractions[p2][0]))])

    # Binned MI estimation
    n_bins = min(20, int(np.sqrt(len(x))))
    if n_bins < 2:
        return 0.0

    hist_xy, _, _ = np.histogram2d(x, y, bins=n_bins)
    hist_xy = hist_xy / hist_xy.sum()

    hist_x = hist_xy.sum(axis=1)
    hist_y = hist_xy.sum(axis=0)

    mi = 0.0
    for i in range(n_bins):
        for j in range(n_bins):
            if hist_xy[i, j] > 0 and hist_x[i] > 0 and hist_y[j] > 0:
                mi += hist_xy[i, j] * np.log(hist_xy[i, j] / (hist_x[i] * hist_y[j]))
    return mi


def compute_coupling_ratios(fractions, primes):
    """Compute variance-based coupling ratios R_pq."""
    variances = {}
    for p in primes:
        all_fracs = []
        for k in range(p):
            all_fracs.extend(fractions[p][k])
        variances[p] = np.var(all_fracs) if len(all_fracs) > 0 else 0.0

    ratios = {}
    for p1, p2 in combinations(primes, 2):
        key = f"R_{p1}{p2}"
        ratios[key] = variances[p1] / variances[p2] if variances[p2] > 0 else float('inf')

    return variances, ratios


def check_T0_suppression(counts, n_atoms, n_samples):
    """Check if T[1][1] = 0 (mod 3) is suppressed in excitation patterns."""
    # For consecutive atom pairs (j, j+1), count transitions between classes mod 3
    transitions = np.zeros((3, 3))
    total = 0

    for bitstring, count in counts.items():
        bits = [int(b) for b in str(bitstring)]
        if len(bits) < n_atoms:
            bits = [0] * (n_atoms - len(bits)) + bits

        excited = [j for j in range(min(len(bits), n_atoms)) if bits[j] == 1]

        for idx in range(len(excited) - 1):
            c_from = excited[idx] % 3
            c_to = excited[idx + 1] % 3
            transitions[c_from, c_to] += count
            total += count

    if total > 0:
        T_emp = transitions / transitions.sum(axis=1, keepdims=True)
        T_emp = np.nan_to_num(T_emp)
    else:
        T_emp = np.zeros((3, 3))

    return T_emp, transitions


# ====================================================================
# MAIN ANALYSIS
# ====================================================================

def analyse_geometry(counts, n_atoms, primes, n_samples, geom_name):
    """Full QH1 analysis for one geometry."""
    print(f"\n  --- Geometrie: {geom_name} ({n_atoms} atomes) ---")

    # CRT fractions
    fractions = compute_crt_fractions(counts, n_atoms, primes, n_samples)

    # Mutual information (QH1a)
    print(f"\n  QH1a : Information mutuelle entre facteurs CRT")
    mi_values = {}
    for p1, p2 in combinations(primes, 2):
        mi = compute_mutual_information(fractions, p1, p2)
        mi_values[f"I({p1},{p2})"] = mi
        print(f"    I(mod {p1} ; mod {p2}) = {mi:.4f} nats")

    # Coupling ratios (QH1b)
    print(f"\n  QH1b : Ratios de couplage")
    variances, ratios = compute_coupling_ratios(fractions, primes)
    for p in primes:
        print(f"    Var(mod {p}) = {variances[p]:.6f}")
    for key, val in ratios.items():
        print(f"    {key} = {val:.4f}")

    if 5 in primes and 7 in primes and 3 in primes:
        r35 = ratios.get('R_35', ratios.get('R_53', 0))
        r37 = ratios.get('R_37', ratios.get('R_73', 0))
        if r37 > 0:
            measured_ratio = r35 / r37
            print(f"    R_35/R_37 = {measured_ratio:.4f}  (PT predit: {RATIO_35_37:.4f})")
            print(f"    Ecart: {abs(measured_ratio - RATIO_35_37) / RATIO_35_37 * 100:.1f}%")

    # T0 suppression (QH1c)
    print(f"\n  QH1c : Suppression T0 (transitions mod 3)")
    T_emp, T_counts = check_T0_suppression(counts, n_atoms, n_samples)
    print(f"    Matrice de transition empirique (mod 3):")
    for i in range(3):
        row = '    '.join(f'{T_emp[i,j]:.3f}' for j in range(3))
        print(f"      [{row}]")
    if T_counts.sum() > 0:
        t11_frac = T_emp[1, 1]
        t11_count = int(T_counts[1, 1])
        t11_total = int(T_counts[1].sum())
        print(f"    T[1][1] = {t11_frac:.4f} ({t11_count}/{t11_total})")
        print(f"    T0 suppression: {'OUI' if t11_frac < 0.2 else 'NON'}")

    return {
        'geometry': geom_name,
        'n_atoms': n_atoms,
        'mi': mi_values,
        'variances': {str(k): v for k, v in variances.items()},
        'ratios': ratios,
        'T_emp_mod3': T_emp.tolist(),
    }


def compare_geometries(all_results, primes):
    """Compare CRT correlations across geometries (QH1a test)."""
    print("\n" + "=" * 70)
    print("  QH1a : TEST D'INVARIANCE GEOMETRIQUE")
    print("=" * 70)

    # Collect MI values across geometries
    mi_keys = list(all_results[0]['mi'].keys())
    for key in mi_keys:
        values = [r['mi'][key] for r in all_results]
        mean_val = np.mean(values)
        std_val = np.std(values)
        cv = std_val / mean_val if mean_val > 0 else float('inf')
        geom_names = [r['geometry'] for r in all_results]

        print(f"\n  {key}:")
        for name, val in zip(geom_names, values):
            print(f"    {name:>12} : {val:.4f}")
        print(f"    {'moyenne':>12} : {mean_val:.4f} +/- {std_val:.4f}")
        print(f"    {'CV':>12} : {cv:.2%}")
        invariant = cv < 0.5  # threshold: 50% relative variation
        print(f"    Invariant?   : {'OUI' if invariant else 'NON'} (seuil CV < 50%)")

    # QH1b ratio across geometries
    print(f"\n  QH1b : RATIO sin^2(theta_5)/sin^2(theta_7) across geometries")
    print(f"  PT predit: {RATIO_35_37:.4f}")

    ratio_key = 'R_35'
    ratio_key2 = 'R_37'
    if all(ratio_key in r['ratios'] and ratio_key2 in r['ratios'] for r in all_results):
        measured = []
        for r in all_results:
            r35 = r['ratios'][ratio_key]
            r37 = r['ratios'][ratio_key2]
            ratio = r35 / r37 if r37 > 0 else float('inf')
            measured.append(ratio)
            print(f"    {r['geometry']:>12} : {ratio:.4f}")
        mean_r = np.mean(measured)
        std_r = np.std(measured)
        print(f"    {'moyenne':>12} : {mean_r:.4f} +/- {std_r:.4f}")
        print(f"    Ecart vs PT  : {abs(mean_r - RATIO_35_37) / RATIO_35_37 * 100:.1f}%")


def plot_qh1(all_results, primes):
    """Visualize QH1 results."""
    n_geom = len(all_results)
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    geom_names = [r['geometry'] for r in all_results]

    # Panel 1: MI across geometries
    ax1 = axes[0]
    mi_keys = list(all_results[0]['mi'].keys())
    x = np.arange(n_geom)
    width = 0.8 / len(mi_keys)
    for idx, key in enumerate(mi_keys):
        values = [r['mi'][key] for r in all_results]
        ax1.bar(x + idx * width, values, width, label=key, alpha=0.8)
    ax1.set_xticks(x + width * (len(mi_keys) - 1) / 2)
    ax1.set_xticklabels(geom_names, rotation=30, ha='right')
    ax1.set_ylabel('Information mutuelle (nats)')
    ax1.set_title('QH1a : MI entre facteurs CRT')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    # Panel 2: T0 suppression
    ax2 = axes[1]
    t11_values = [r['T_emp_mod3'][1][1] for r in all_results]
    colors = ['green' if v < 0.2 else 'red' for v in t11_values]
    ax2.bar(geom_names, t11_values, color=colors, alpha=0.8)
    ax2.axhline(y=1/3, color='gray', linestyle='--', label='Uniforme (1/3)')
    ax2.axhline(y=0, color='red', linestyle='-', alpha=0.3, label='PT predit (0)')
    ax2.set_ylabel('T[1][1] empirique')
    ax2.set_title('QH1c : Suppression T0 (mod 3)')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    plt.setp(ax2.get_xticklabels(), rotation=30, ha='right')

    # Panel 3: Coupling ratio R35/R37
    ax3 = axes[2]
    ratio_vals = []
    for r in all_results:
        r35 = r['ratios'].get('R_35', 0)
        r37 = r['ratios'].get('R_37', 1)
        ratio_vals.append(r35 / r37 if r37 > 0 else 0)
    ax3.bar(geom_names, ratio_vals, color='steelblue', alpha=0.8)
    ax3.axhline(y=RATIO_35_37, color='red', linestyle='--', linewidth=2,
                label=f'PT: sin2(th5)/sin2(th7) = {RATIO_35_37:.3f}')
    ax3.set_ylabel('R_35 / R_37')
    ax3.set_title('QH1b : Ratio de couplage CRT')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)
    plt.setp(ax3.get_xticklabels(), rotation=30, ha='right')

    plt.suptitle('QH1 : Spectroscopie CRT -- Test d\'emergence', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig('QH1_crt_spectroscopy.png', dpi=150, bbox_inches='tight')
    print(f"\n  Figure sauvegardee : QH1_crt_spectroscopy.png")


# ====================================================================
# MAIN
# ====================================================================

def parse_args():
    parser = argparse.ArgumentParser(description='QH1: CRT Spectroscopy')
    parser.add_argument('--mode', default='local',
                        choices=['local', 'emu_free', 'emu_tn', 'fresnel', 'qpu'])
    parser.add_argument('--project-id', default=None)
    parser.add_argument('--n-atoms', type=int, default=None,
                        help='Nombre d\'atomes (defaut: 15 local, 105 cloud)')
    parser.add_argument('--n-samples', type=int, default=2000)
    parser.add_argument('--t-pulse', type=int, default=1000)
    parser.add_argument('--omega-max', type=float, default=4.0)
    parser.add_argument('--no-plot', action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    if args.n_atoms is None:
        args.n_atoms = 15 if args.mode == 'local' else 105

    # Determine active primes for this register size
    primes = [p for p in ACTIVE_PRIMES if args.n_atoms % p == 0]
    if not primes:
        print(f"ERREUR: N={args.n_atoms} n'est divisible par aucun premier actif {{3,5,7}}")
        sys.exit(1)

    primorial = 1
    for p in primes:
        primorial *= p

    print("=" * 70)
    print("  QH1 : Spectroscopie CRT sur le primorial actif")
    print(f"  N = {args.n_atoms} atomes")
    print(f"  Premiers actifs testes : {primes} (primorial = {primorial})")
    print(f"  Mode : {args.mode.upper()}")
    print(f"  Shots : {args.n_samples}")
    print("=" * 70)

    print(f"\n  Predictions PT :")
    for p in primes:
        print(f"    sin^2(theta_{p}, q_plus) = {SIN2[p]:.4f}")
    if 5 in primes and 7 in primes:
        print(f"    Ratio QH1b = sin^2(th5)/sin^2(th7) = {RATIO_35_37:.4f}")
    print(f"    T[1][1] mod 3 = 0 (QH1c)")

    # Generate multiple geometries
    if args.n_atoms <= 20:
        spacing = 7.0
    else:
        spacing = 5.5

    geometries = [
        make_grid(args.n_atoms, spacing=spacing),
        make_ring(args.n_atoms),
        make_random(args.n_atoms, seed=42),
        make_zigzag(args.n_atoms, spacing=spacing),
    ]

    all_results = []

    for coords, geom_name in geometries:
        if len(coords) < args.n_atoms:
            print(f"\n  SKIP {geom_name}: seulement {len(coords)} atomes places")
            continue

        print(f"\n{'='*50}")
        print(f"  Geometrie: {geom_name}")
        print(f"{'='*50}")

        counts = build_and_run(coords, args.n_samples, args.t_pulse,
                               args.omega_max, args.mode, args.project_id)

        result = analyse_geometry(counts, args.n_atoms, primes,
                                  args.n_samples, geom_name)
        all_results.append(result)

    # Cross-geometry comparison
    if len(all_results) >= 2:
        compare_geometries(all_results, primes)

    # Visualization
    if not args.no_plot and len(all_results) >= 2:
        plot_qh1(all_results, primes)

    # Export
    export = {
        'protocol': 'QH1',
        'n_atoms': args.n_atoms,
        'primes': primes,
        'mode': args.mode,
        'n_samples': args.n_samples,
        'PT_predictions': {
            'sin2': {str(p): SIN2[p] for p in primes},
            'ratio_QH1b': RATIO_35_37 if (5 in primes and 7 in primes) else None,
            'T11_mod3': 0.0,
        },
        'results': all_results,
    }
    outfile = 'QH1_results.json'
    with open(outfile, 'w') as f:
        json.dump(export, f, indent=2)
    print(f"\n  Resultats exportes : {outfile}")

    print("\n" + "=" * 70)
    if args.mode == 'local' and args.n_atoms < 105:
        print(f"  NOTE: Simulation locale avec N={args.n_atoms} (preuve de concept)")
        print(f"  Le test complet QH1 requiert N=105 atomes sur Pasqal QPU")
    print("=" * 70)
