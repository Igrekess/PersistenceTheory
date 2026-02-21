"""
run_demos.py - Script maitre pour les demonstrations de la Theorie de la Persistance

Usage:
    python scripts/demonstrations/run_demos.py          # menu interactif
    python scripts/demonstrations/run_demos.py --all    # batch + JSON auto
    python scripts/demonstrations/run_demos.py --all --json results/run.json

Auteur : Yan Senez (2026)
"""

import argparse
import datetime
import json
import os
import subprocess
import sys
import time

# ---------------------------------------------------------------------------
# Registre des 24 demonstrations (hardcode, ordre canonique)
# ---------------------------------------------------------------------------
DEMOS = [
    {
        "id": "L0",
        "title_fr": "Unicite de q_stat",
        "title_en": "Uniqueness of q_stat",
        "script": "L0_uniqueness_q_stat.py",
    },
    {
        "id": "D00",
        "title_fr": "Transitions interdites mod 3",
        "title_en": "Forbidden transitions mod 3",
        "script": "D00_forbidden_transitions.py",
    },
    {
        "id": "D01",
        "title_fr": "Theoreme de conservation",
        "title_en": "Conservation theorem",
        "script": "D01_conservation_theorem.py",
    },
    {
        "id": "D02",
        "title_fr": "GFT et operateur de Ruelle",
        "title_en": "GFT and Ruelle operator",
        "script": "D02_gft_ruelle.py",
    },
    {
        "id": "D03",
        "title_fr": "Sieve variationnel",
        "title_en": "Variational sieve",
        "script": "D03_variational_sieve.py",
    },
    {
        "id": "D04",
        "title_fr": "Formule maitre et point fixe",
        "title_en": "Master formula and fixed point",
        "script": "D04_master_formula_fp.py",
    },
    {
        "id": "D05",
        "title_fr": "Loi de Mertens",
        "title_en": "Mertens law",
        "script": "D05_mertens_law.py",
    },
    {
        "id": "D06",
        "title_fr": "Preuve q positif",
        "title_en": "Proof q positive",
        "script": "D06_proof_q_positive.py",
    },
    {
        "id": "D07",
        "title_fr": "Identite sin^2",
        "title_en": "sin^2 identity",
        "script": "D07_sin2_identity.py",
    },
    {
        "id": "D08",
        "title_fr": "Point fixe mu* = 15",
        "title_en": "Fixed point mu* = 15",
        "script": "D08_fixed_point_mu15.py",
    },
    {
        "id": "D09",
        "title_fr": "Constante de structure fine alpha_EM",
        "title_en": "Fine structure constant alpha_EM",
        "script": "D09_alpha_em.py",
    },
    {
        "id": "D10",
        "title_fr": "Metrique de Bianchi",
        "title_en": "Bianchi metric",
        "script": "D10_bianchi_metric.py",
    },
    {
        "id": "D11",
        "title_fr": "Potentiel de persistance",
        "title_en": "Persistence potential",
        "script": "D11_persistence_potential.py",
    },
    {
        "id": "D12",
        "title_fr": "Equations d'Einstein",
        "title_en": "Einstein equations",
        "script": "D12_einstein_equations.py",
    },
    {
        "id": "D13",
        "title_fr": "Spin foam U(1)",
        "title_en": "Spin foam U(1)",
        "script": "D13_spin_foam_u1.py",
    },
    {
        "id": "D14",
        "title_fr": "Unification GFT",
        "title_en": "GFT unification",
        "script": "D14_unification_gft.py",
    },
    {
        "id": "D15",
        "title_fr": "Charge topologique",
        "title_en": "Topological charge",
        "script": "D15_charge_topological.py",
    },
    {
        "id": "D16",
        "title_fr": "Chaine causale",
        "title_en": "Causal chain",
        "script": "D16_causal_chain.py",
    },
    {
        "id": "D17",
        "title_fr": "Profondeur exactement 2",
        "title_en": "Depth exactly 2",
        "script": "D17_depth_exactly_2.py",
    },
    {
        "id": "D17b",
        "title_fr": "Modulation de Catalan",
        "title_en": "Catalan modulation",
        "script": "D17b_catalan_modulation.py",
    },
    {
        "id": "D18",
        "title_fr": "Hardy-Littlewood et persistance",
        "title_en": "Hardy-Littlewood and persistence",
        "script": "D18_hardy_littlewood.py",
    },
    {
        "id": "D19",
        "title_fr": "Correction a une boucle",
        "title_en": "One-loop correction",
        "script": "D19_one_loop_correction.py",
    },
    {
        "id": "D27",
        "title_fr": "Longueur de coherence",
        "title_en": "Coherence length",
        "script": "D27_coherence_length.py",
    },
    {
        "id": "D28",
        "title_fr": "Bruit de fond spatio-temporel PT",
        "title_en": "Spacetime background noise PT",
        "script": "D28_spacetime_noise.py",
    },
]

TIMEOUT_S = 60  # timeout par script en secondes

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_scripts_dir():
    """Repertoire contenant ce fichier."""
    return os.path.dirname(os.path.abspath(__file__))


def run_script(script_path, stream_output=False):
    """
    Execute un script Python.
    Retourne (return_code, stdout, stderr, duration_s).
    """
    cmd = [sys.executable, script_path]
    t0 = time.time()
    try:
        if stream_output:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
                errors="replace",
            )
            stdout_lines = []
            stderr_lines = []
            # Lire stdout en temps reel
            for line in proc.stdout:
                print(line, end="", flush=True)
                stdout_lines.append(line)
            proc.stdout.close()
            proc.wait(timeout=TIMEOUT_S)
            stderr_data = proc.stderr.read()
            proc.stderr.close()
            rc = proc.returncode
            return rc, "".join(stdout_lines), stderr_data, time.time() - t0
        else:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=TIMEOUT_S,
            )
            return result.returncode, result.stdout, result.stderr, time.time() - t0
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT (>{} s)".format(TIMEOUT_S), time.time() - t0
    except Exception as exc:
        return -2, "", "ERROR: {}".format(exc), time.time() - t0


def determine_status(return_code, stderr):
    """Determine PASS / FAIL / ERROR / TIMEOUT a partir du code de retour."""
    if return_code == -1:
        return "TIMEOUT"
    if return_code == -2:
        return "ERROR"
    if return_code == 0:
        return "PASS"
    return "FAIL"


# ---------------------------------------------------------------------------
# Mode interactif
# ---------------------------------------------------------------------------

def print_menu():
    print("")
    print("=== Demonstrations de la Theorie de la Persistance ===")
    for idx, demo in enumerate(DEMOS):
        print(" {:2d}. {:5s} - {}".format(idx, demo["id"], demo["title_fr"]))
    print("")
    print("  q. Quitter")
    print("")


def run_interactive():
    while True:
        print_menu()
        try:
            raw = input("Choisir (0-{}) ou 'q' : ".format(len(DEMOS) - 1)).strip()
        except (EOFError, KeyboardInterrupt):
            print("\nAu revoir.")
            break

        if raw.lower() == "q":
            print("Au revoir.")
            break

        try:
            idx = int(raw)
        except ValueError:
            print("Entree invalide.")
            continue

        if idx < 0 or idx >= len(DEMOS):
            print("Numero hors plage (0-{}).".format(len(DEMOS) - 1))
            continue

        demo = DEMOS[idx]
        script_path = os.path.join(get_scripts_dir(), demo["script"])

        if not os.path.isfile(script_path):
            print("Script introuvable : {}".format(script_path))
            continue

        print("")
        print("--- {} : {} ---".format(demo["id"], demo["title_fr"]))
        print("")

        rc, _stdout, stderr, duration = run_script(script_path, stream_output=True)
        status = determine_status(rc, stderr)

        if stderr.strip():
            print("")
            print("[stderr]")
            print(stderr.rstrip())

        print("")
        print(
            "--- Fin {} | status={} | duree={:.2f}s ---".format(
                demo["id"], status, duration
            )
        )
        print("")

        if status != "PASS":
            print("Code de retour : {}".format(rc))


# ---------------------------------------------------------------------------
# Mode batch
# ---------------------------------------------------------------------------

def run_batch(json_path=None):
    scripts_dir = get_scripts_dir()
    results = []

    total = len(DEMOS)
    print("")
    print("=== Batch : {} demonstrations ===".format(total))
    print("")

    for idx, demo in enumerate(DEMOS):
        script_path = os.path.join(scripts_dir, demo["script"])
        prefix = "[{:2d}/{:2d}] {:5s}".format(idx + 1, total, demo["id"])

        if not os.path.isfile(script_path):
            print("{} SKIP (script introuvable)".format(prefix))
            results.append(
                {
                    "id": demo["id"],
                    "title_fr": demo["title_fr"],
                    "title_en": demo["title_en"],
                    "script": demo["script"],
                    "status": "SKIP",
                    "return_code": -3,
                    "duration_s": 0.0,
                    "output": "",
                    "error": "Script not found",
                }
            )
            continue

        rc, stdout, stderr, duration = run_script(script_path, stream_output=False)
        status = determine_status(rc, stderr)

        print(
            "{} {} ({:.2f}s)".format(prefix, status, duration)
        )

        results.append(
            {
                "id": demo["id"],
                "title_fr": demo["title_fr"],
                "title_en": demo["title_en"],
                "script": demo["script"],
                "status": status,
                "return_code": rc,
                "duration_s": round(duration, 3),
                "output": stdout,
                "error": stderr,
            }
        )

    # --- Recapitulatif ---
    counts = {"PASS": 0, "FAIL": 0, "ERROR": 0, "TIMEOUT": 0, "SKIP": 0}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1

    passed = counts["PASS"]
    failed = total - passed

    print("")
    print("=== Recapitulatif ===")
    print("  Total   : {}".format(total))
    print("  PASS    : {}".format(counts["PASS"]))
    print("  FAIL    : {}".format(counts["FAIL"]))
    print("  ERROR   : {}".format(counts["ERROR"]))
    print("  TIMEOUT : {}".format(counts["TIMEOUT"]))
    print("  SKIP    : {}".format(counts["SKIP"]))
    print("")

    # Scripts echoues
    failed_list = [r for r in results if r["status"] != "PASS"]
    if failed_list:
        print("Scripts non PASS :")
        for r in failed_list:
            print("  {:5s} {} - {}".format(r["id"], r["status"], r["title_fr"]))
        print("")

    # --- Ecriture JSON ---
    run_date = datetime.datetime.now()
    timestamp = run_date.strftime("%Y%m%d_%H%M%S")

    if json_path is None:
        results_dir = os.path.join(
            os.path.dirname(os.path.dirname(scripts_dir)), "results"
        )
        os.makedirs(results_dir, exist_ok=True)
        json_path = os.path.join(
            results_dir, "demo_results_{}.json".format(timestamp)
        )
    else:
        os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)

    payload = {
        "run_date": run_date.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": "{}.{}.{}".format(
            sys.version_info.major,
            sys.version_info.minor,
            sys.version_info.micro,
        ),
        "total": total,
        "passed": passed,
        "failed": failed,
        "results": results,
    }

    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False, indent=2)

    print("Resultats JSON ecrits dans : {}".format(json_path))
    return passed == total


# ---------------------------------------------------------------------------
# Point d'entree
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Script maitre des demonstrations de la Theorie de la Persistance"
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Mode batch : execute tous les scripts et produit un JSON",
    )
    parser.add_argument(
        "--json",
        metavar="PATH",
        default=None,
        help="Chemin du fichier JSON de sortie (mode --all uniquement)",
    )
    args = parser.parse_args()

    if args.all:
        success = run_batch(json_path=args.json)
        sys.exit(0 if success else 1)
    else:
        run_interactive()


if __name__ == "__main__":
    main()
