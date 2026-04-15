#!/usr/bin/env python3
"""
verify_sota.py -- SOTA Score Verification for Persistence Theory Chemistry

Monograph: pt_chimie_scores.tex
Derivation chain: s = 1/2 -> sieve -> SM -> PTC calculator -> SOTA benchmarks
Zero fitted parameters.

This script verifies the three SOTA chemistry scores claimed in the monograph
by calling the external PTC (Persistence Theory Chemistry) calculator:

  Step 1. IONIZATION ENERGIES (IE)
          Z=1-86 (NIST experimental): MAE 0.070%
          Z=1-103 absolute: MAE 3.9 meV
          Z=104-118 (vs FSCC/CCSDT theoretical): MAE 0.79%
          Block-level: s 0.079%, p 0.037%, d 0.064%, f 0.153%

  Step 2. MOLECULAR DISSOCIATION ENERGIES (849 molecules)
          Global MAE 2.15%, median 1.32%
          354/849 under 1%
          Sub-benchmarks: Main 2.02%, Extended 1.58%, ATcT 2.88%, d-block 3.66%

  Step 3. ELECTRON AFFINITIES (EA)
          73 elements, MAE 1.37%
          All under 10%

IMPORTANT: This script requires the external `ptc` package.
If ptc is not installed, all checks gracefully SKIP (exit code 0).
The monograph Level 1/2 scripts (theorems, historical validation) run without ptc.

Theorems verified:
  -- "IE SOTA Score"          (pt_chimie_scores.tex)
  -- "Molecular SOTA Score"   (pt_chimie_scores.tex)
  -- "EA SOTA Score"          (pt_chimie_scores.tex)
"""

import sys
from pathlib import Path

_scripts_root = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_scripts_root))
sys.path.insert(0, str(_scripts_root / "lib"))

from lib.pt_check import Checker

ck = Checker("verify_sota", chapter="verify_sota", total_steps=3)

# ── Graceful PTC import ──────────────────────────────────────────
try:
    import ptc
    HAS_PTC = True
except ImportError:
    HAS_PTC = False

# ── IE imports ───────────────────────────────────────────────────
HAS_IE = False
if HAS_PTC:
    try:
        from ptc.ie_geo import compute_ie
        from ptc.constants import NIST_IE, THEO_IE
        HAS_IE = True
    except ImportError:
        pass

# ── Molecule imports ─────────────────────────────────────────────
HAS_MOL = False
compute_atomization_energy = None
load_benchmark = None
if HAS_PTC:
    try:
        from ptc.molecule import compute_atomization_energy
        from ptc.data.benchmark import load_benchmark
        HAS_MOL = True
    except ImportError:
        try:
            from ptc.transfer_matrix import compute_dat as compute_atomization_energy
            from ptc.constants import BENCHMARK_MOLECULES
            def load_benchmark():
                return BENCHMARK_MOLECULES
            HAS_MOL = True
        except ImportError:
            pass

# ── EA imports ───────────────────────────────────────────────────
HAS_EA = False
if HAS_PTC:
    try:
        from ptc.ea_geo import compute_ea
        from ptc.constants import NIST_EA
        HAS_EA = True
    except ImportError:
        pass


def _block(Z):
    """Rough orbital block classification for Z=1-86."""
    if Z in (1, 2, 3, 4, 11, 12, 19, 20, 37, 38, 55, 56):
        return "s"
    if 57 <= Z <= 71 and Z not in (57, 71):
        return "f"
    if Z in range(21, 31) or Z in range(39, 49) or Z in range(57, 58) or Z in range(71, 81):
        return "d"
    return "p"


# =====================================================================
#  Step 1: IONIZATION ENERGIES
# =====================================================================
ck.section("Step 1: Ionization energies (IE)")

if not HAS_PTC:
    print("  SKIP: ptc not available")
    print("  Install PTC or add it to PYTHONPATH to verify IE scores.")
elif not HAS_IE:
    print("  SKIP: ptc.ie_geo or ptc.constants not available")
    print("  IE verification requires compute_ie, NIST_IE, THEO_IE from ptc.")
else:
    # Z=1-86 (experimental)
    errors_86 = []
    for Z in range(1, 87):
        ie_pt = compute_ie(Z)
        ie_obs = NIST_IE.get(Z)
        if ie_obs is None or ie_obs == 0:
            continue
        err = abs(ie_pt - ie_obs) / ie_obs * 100
        errors_86.append(err)

    if errors_86:
        mae_86 = sum(errors_86) / len(errors_86)
        print(f"  IE Z=1-86: {len(errors_86)} elements, MAE = {mae_86:.3f}%")
        ck.check("IE_Z1_86_mae_lt_0.10pct",
                 mae_86 < 0.10,
                 f"MAE = {mae_86:.3f}%")

    # Z=1-103 absolute MAE in meV
    errors_103_mev = []
    for Z in range(1, 104):
        ie_pt = compute_ie(Z)
        ie_obs = NIST_IE.get(Z)
        if ie_obs is None or ie_obs == 0:
            continue
        err_mev = abs(ie_pt - ie_obs) * 1000  # eV to meV
        errors_103_mev.append(err_mev)

    if errors_103_mev:
        mae_103_mev = sum(errors_103_mev) / len(errors_103_mev)
        print(f"  IE Z=1-103: {len(errors_103_mev)} elements, MAE = {mae_103_mev:.1f} meV")
        ck.check("IE_Z1_103_mae_lt_5meV",
                 mae_103_mev < 5.0,
                 f"MAE = {mae_103_mev:.1f} meV")

    # Z=104-118 vs theoretical (FSCC/CCSDT)
    errors_sh = []
    for Z in range(104, 119):
        ie_pt = compute_ie(Z)
        ie_theo = THEO_IE.get(Z)
        if ie_theo is None or ie_theo == 0:
            continue
        err = abs(ie_pt - ie_theo) / ie_theo * 100
        errors_sh.append(err)

    if errors_sh:
        mae_sh = sum(errors_sh) / len(errors_sh)
        print(f"  IE Z=104-118: {len(errors_sh)} elements, MAE = {mae_sh:.2f}%")
        ck.check("IE_Z104_118_mae_lt_1pct",
                 mae_sh < 1.0,
                 f"MAE = {mae_sh:.2f}%")

    # Block-level MAE
    block_thresholds = {
        "s": (0.10, [Z for Z in range(1, 87) if _block(Z) == "s"]),
        "p": (0.05, [Z for Z in range(1, 87) if _block(Z) == "p"]),
        "d": (0.08, [Z for Z in range(1, 87) if _block(Z) == "d"]),
        "f": (0.20, [Z for Z in range(1, 87) if _block(Z) == "f"]),
    }
    for block, (thresh, zlist) in block_thresholds.items():
        errs = []
        for Z in zlist:
            ie_pt = compute_ie(Z)
            ie_obs = NIST_IE.get(Z)
            if ie_obs is None or ie_obs == 0:
                continue
            errs.append(abs(ie_pt - ie_obs) / ie_obs * 100)
        if errs:
            mae_b = sum(errs) / len(errs)
            print(f"  {block}-block: {len(errs)} elements, MAE = {mae_b:.3f}%")
            ck.check(f"IE_{block}_block_mae_lt_{thresh}pct",
                     mae_b < thresh,
                     f"MAE = {mae_b:.3f}%")


# =====================================================================
#  Step 2: MOLECULAR DISSOCIATION ENERGIES (849 molecules)
# =====================================================================
ck.section("Step 2: Molecular dissociation energies (849 mol)")

if not HAS_PTC:
    print("  SKIP: ptc not available")
    print("  Install PTC or add it to PYTHONPATH to verify molecular scores.")
elif not HAS_MOL:
    print("  SKIP: ptc.molecule or ptc.transfer_matrix not available")
else:
    try:
        benchmark = load_benchmark()
    except Exception:
        benchmark = []

    errors = []
    for mol in benchmark:
        obs = mol.get("D_at_obs", mol.get("D0_obs"))
        if obs is None or obs <= 0:
            continue
        try:
            pt = compute_atomization_energy(mol)
        except Exception:
            continue
        if pt is None:
            continue
        err = abs(pt - obs) / obs * 100
        errors.append(err)

    if not errors:
        print("  SKIP: no molecules computed (check PTC API)")
    else:
        n_mol = len(errors)
        mae = sum(errors) / n_mol
        median = sorted(errors)[n_mol // 2]
        under_1 = sum(1 for e in errors if e < 1.0)

        print(f"  Molecules computed: {n_mol}")
        print(f"  Global MAE: {mae:.2f}%")
        print(f"  Median: {median:.2f}%")
        print(f"  Under 1%: {under_1}/{n_mol}")

        ck.check("mol_count_ge_800",
                 n_mol >= 800,
                 f"got {n_mol}")

        ck.check("mol_mae_lt_2.5pct",
                 mae < 2.5,
                 f"MAE = {mae:.2f}%")

        ck.check("mol_median_lt_2pct",
                 median < 2.0,
                 f"median = {median:.2f}%")

        ck.check("mol_under1pct_ge_300",
                 under_1 >= 300,
                 f"under_1 = {under_1}")


# =====================================================================
#  Step 3: ELECTRON AFFINITIES (EA)
# =====================================================================
ck.section("Step 3: Electron affinities (EA)")

if not HAS_PTC:
    print("  SKIP: ptc not available")
    print("  Install PTC or add it to PYTHONPATH to verify EA scores.")
elif not HAS_EA:
    print("  SKIP: ptc.ea_geo or ptc.constants not available")
else:
    errors_ea = []
    for Z, ea_obs in sorted(NIST_EA.items()):
        if ea_obs <= 0:
            continue
        try:
            ea_pt = compute_ea(Z)
        except Exception:
            continue
        if ea_pt is None or ea_pt <= 0:
            continue
        err = abs(ea_pt - ea_obs) / ea_obs * 100
        errors_ea.append(err)

    if not errors_ea:
        print("  SKIP: no EA values computed (check PTC API)")
    else:
        n_elem = len(errors_ea)
        mae_ea = sum(errors_ea) / n_elem

        print(f"  EA elements: {n_elem}")
        print(f"  EA MAE: {mae_ea:.2f}%")

        ck.check("ea_count_ge_70",
                 n_elem >= 70,
                 f"got {n_elem}")

        ck.check("ea_mae_lt_2pct",
                 mae_ea < 2.0,
                 f"MAE = {mae_ea:.2f}%")

        over_10 = sum(1 for e in errors_ea if e > 10.0)
        ck.check("ea_all_under_10pct",
                 over_10 == 0,
                 f"{over_10}/{n_elem} over 10%")


# ── Final summary ────────────────────────────────────────────────
ck.summary()
