from math import prod

print("CRT DEMONSTRATOR: WHY PRIME CHANNELS?")
print("=====================================")
print("Question: why are primes the right base channels?")
print("Idea: modulo 30, a position can be read separately modulo 2, 3, and 5.")
print("This is the intuition behind the Chinese Remainder Theorem.")
print()

channels = [2, 3, 5]
modulus = prod(channels)

def signature(n):
    return tuple(n % p for p in channels)

residues = [n for n in range(modulus) if all(n % p != 0 for p in channels)]
signatures = {n: signature(n) for n in residues}

print("Independent prime channels:", channels)
print("Recombined modulus        :", modulus)
print("Calculation: M = 2 x 3 x 5 =", modulus)
print()
print("Each survivor carries a separated signature:")
for n, sig in signatures.items():
    checks = [f"{n} mod {p} = {n % p}" for p in channels]
    print(f"  n={n:2d}  ->  {', '.join(checks)}  -> signature {sig}")

print()
print("Signature uniqueness check:")
print("  number of survivors           =", len(signatures))
print("  number of distinct signatures =", len(set(signatures.values())))
print()
print("PT reading:")
print("  - a composite channel is not fundamental: it factorizes;")
print("  - a prime channel adds an irreducible constraint;")
print("  - CRT explains why prime channels can be read separately and recomposed.")
print()
assert len(set(signatures.values())) == len(signatures)
print("[PASS] every survivor carries a CRT signature")
print("[PASS] prime channels are read separately and recomposed")
