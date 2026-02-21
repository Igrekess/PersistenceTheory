#!/usr/bin/env python3
"""
Run ALL test scripts and collect results in a JSON file.
Extracts scores, PASS/FAIL, and key numerical values.
"""
import subprocess
import json
import os
import sys
import time
import re

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR = os.path.dirname(SCRIPTS_DIR)

# Collect all .py scripts (excluding this runner)
def find_all_scripts():
    # Skip entire directories
    SKIP_DIRS = {'universality', '_supporting', '__pycache__'}
    # Skip specific file prefixes
    SKIP_FILES = [
        'generate_gaps_streaming', 'generate_gaps_streaming_cuda',
        'test_TP_complete_donnees_reelles',
    ]
    scripts = []
    for root, dirs, files in os.walk(SCRIPTS_DIR):
        # Prune skipped directories
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for f in files:
            if f.endswith('.py') and f != 'run_all_tests.py':
                if any(f.startswith(s) for s in SKIP_FILES):
                    continue
                full = os.path.join(root, f)
                rel = os.path.relpath(full, REPO_DIR)
                scripts.append((rel, full))
    scripts.sort()
    return scripts

def extract_score(output):
    """Extract score patterns like 46/46, 10/10, 7/7, 25/25, etc."""
    # Look for "Score: X/Y" or "X/Y PASS" or "SCORE: X/Y"
    patterns = [
        r'[Ss]core\s*[:=]\s*(\d+)/(\d+)',
        r'(\d+)/(\d+)\s*(?:PASS|tests?\s+pass)',
        r'TOTAL\s*:\s*(\d+)/(\d+)',
        r'(\d+)/(\d+)\s*(?:etapes?|steps?)',
    ]
    scores = []
    for pat in patterns:
        for m in re.finditer(pat, output):
            scores.append(f"{m.group(1)}/{m.group(2)}")
    return scores[-1] if scores else None

def extract_verdict(output):
    """Extract overall verdict."""
    last_lines = output[-2000:]
    if 'FAIL' in last_lines and 'PASS' not in last_lines[-500:]:
        return 'FAIL'
    if re.search(r'(?:PROUVE|PROVEN|PASS|SUCCESS)', last_lines[-500:], re.IGNORECASE):
        return 'PASS'
    if re.search(r'Error|Traceback|Exception', last_lines[-1000:]):
        return 'ERROR'
    return 'COMPLETED'

def extract_key_values(output):
    """Extract key numerical results."""
    values = {}
    # alpha_EM
    m = re.search(r'1/alpha\s*[=:]\s*([\d.]+)', output)
    if m:
        values['1/alpha_EM'] = float(m.group(1))
    # Error percentages
    for m in re.finditer(r'(\w[\w\s]*?)\s*[:=]\s*([\d.]+)\s*%\s*(?:err|error|ecart)', output, re.IGNORECASE):
        key = m.group(1).strip()[:40]
        values[f'err_{key}'] = float(m.group(2))
    return values

def run_script(rel_path, full_path, timeout=120):
    """Run a single script and return results."""
    result = {
        'script': rel_path.replace('\\', '/'),
        'status': None,
        'score': None,
        'verdict': None,
        'duration_s': None,
        'error': None,
    }

    start = time.time()
    try:
        proc = subprocess.run(
            [sys.executable, full_path],
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=REPO_DIR,
            encoding='utf-8',
            errors='replace'
        )
        duration = time.time() - start
        result['duration_s'] = round(duration, 1)

        output = proc.stdout + proc.stderr

        if proc.returncode == 0:
            result['status'] = 'OK'
            result['score'] = extract_score(output)
            result['verdict'] = extract_verdict(output)
        else:
            result['status'] = f'EXIT_{proc.returncode}'
            result['verdict'] = 'ERROR'
            # Get last few lines of error
            err_lines = (proc.stderr or proc.stdout or '').strip().split('\n')
            result['error'] = '\n'.join(err_lines[-3:])[:300]

    except subprocess.TimeoutExpired:
        result['status'] = 'TIMEOUT'
        result['duration_s'] = timeout
        result['verdict'] = 'TIMEOUT'
    except Exception as e:
        result['status'] = 'EXCEPTION'
        result['error'] = str(e)[:200]

    return result

def main():
    scripts = find_all_scripts()
    print(f"Found {len(scripts)} scripts to run")
    print("=" * 80)

    results = []
    passed = 0
    failed = 0
    errors = 0
    timeouts = 0

    for i, (rel, full) in enumerate(scripts, 1):
        short = rel.replace('scripts\\', '').replace('scripts/', '')
        print(f"[{i:3d}/{len(scripts)}] {short:60s} ", end='', flush=True)

        r = run_script(rel, full, timeout=120)
        results.append(r)

        status_str = r['status']
        if r['score']:
            status_str += f" ({r['score']})"

        if r['status'] == 'OK':
            passed += 1
            score_str = r.get('score') or ''
            dur = r.get('duration_s') or 0
            print(f"OK  {score_str:<10s} {dur:.1f}s")
        elif r['status'] == 'TIMEOUT':
            timeouts += 1
            print(f"TIMEOUT")
        else:
            failed += 1
            err_short = (r.get('error', '') or '')[:60].replace('\n', ' ')
            err_short = err_short.encode('ascii', 'replace').decode('ascii')
            print(f"FAIL  {err_short}")

    print("=" * 80)
    print(f"TOTAL: {len(scripts)} scripts")
    print(f"  OK:      {passed}")
    print(f"  FAIL:    {failed}")
    print(f"  TIMEOUT: {timeouts}")
    print()

    # Summary
    summary = {
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        'total_scripts': len(scripts),
        'passed': passed,
        'failed': failed,
        'timeouts': timeouts,
        'results': results
    }

    out_path = os.path.join(REPO_DIR, 'gaps_data', 'test_results.json')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {out_path}")

    # Print failures
    if failed > 0:
        print("\n--- FAILURES ---")
        for r in results:
            if r['status'] not in ('OK', 'TIMEOUT'):
                print(f"  {r['script']}: {r.get('error', '')[:100]}")

    if timeouts > 0:
        print("\n--- TIMEOUTS ---")
        for r in results:
            if r['status'] == 'TIMEOUT':
                print(f"  {r['script']}")

if __name__ == '__main__':
    main()
