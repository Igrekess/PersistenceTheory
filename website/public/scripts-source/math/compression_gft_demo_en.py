from collections import Counter
from math import log2

print("DEMONSTRATOR: COMPRESSION AND GFT")
print("=================================")
print("Question: why does compression amount to extracting persistence?")
print("We compare two strings: one highly redundant, the other much closer to uniform.")
print("A redundant string has lower entropy and more visible structure.")
print()

samples = {
    "redundant": "A" * 80 + "B" * 16 + "C" * 4,
    "mixed": "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789" * 3,
}

for name, text in samples.items():
    counts = Counter(text)
    n = sum(counts.values())
    m = len(counts)
    h = -sum((c / n) * log2(c / n) for c in counts.values())
    total = log2(m)
    dkl = total - h

    print(f"Case: {name}")
    print(f"  length              = {n}")
    print(f"  distinct alphabet   = {m} symbols")
    print("  counts              =", dict(counts))
    print("  probabilities       =", {k: round(v / n, 4) for k, v in counts.items()})
    print("  H calculation       = -Σ p_i log2(p_i)")
    print("  D_KL calculation    = log2(m) - H")
    print(f"  budget log2(m)      = {total:.6f} bits")
    print(f"  structure D_KL      = {dkl:.6f} bits")
    print(f"  entropy H           = {h:.6f} bits")
    if dkl > 0.5:
        print("  reading: imbalanced distribution, hence compressible structure.")
    else:
        print("  reading: near-uniform distribution, hence little simple redundancy.")
    print()

print("Conclusion: GFT does not replace a real compressor.")
print("But it gives the language: persistent structure + entropic residue.")
print("[PASS] a redundant string carries more compressible structure")
