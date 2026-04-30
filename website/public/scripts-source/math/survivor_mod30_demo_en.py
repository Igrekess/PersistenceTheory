from math import gcd

print("SURVIVORS MODULO 30 DEMONSTRATOR")
print("================================")
print("Question: what remains after imposing the constraints 2, 3, and 5?")
print("In PT, positions eliminated by no channel are persistence points of the sieve.")
print("The gaps are not added afterwards: they are distances between those survivors.")
print()

modulus = 2 * 3 * 5
survivors = [n for n in range(1, modulus + 1) if gcd(n, modulus) == 1]
gaps = []
for i, n in enumerate(survivors):
    nxt = survivors[(i + 1) % len(survivors)]
    if i == len(survivors) - 1:
        nxt += modulus
    gaps.append(nxt - n)

print(f"Full modulus M = 2 x 3 x 5 = {modulus}")
print("Tested positions: 1, 2, ..., 30")
print("Survival rule   : gcd(n, 30) = 1")
print()
print("Local calculation on the 30 positions:")
for n in range(1, modulus + 1):
    status = "SURVIVES" if gcd(n, modulus) == 1 else "removed"
    print(f"  n={n:2d}  gcd({n:2d},30)={gcd(n, modulus):2d}  -> {status}")
print()
print("Survivors found :", survivors)
print("Circular gaps   :", gaps)
print("Exact density   :", f"{len(survivors)}/{modulus}", "=", f"{len(survivors)/modulus:.6f}")
print("Gap calculation :")
for n, gap in zip(survivors, gaps):
    nxt = n + gap
    if nxt > modulus:
        nxt -= modulus
    print(f"  from {n:2d} to {nxt:2d}: gap = {gap}")
print()
print("PT reading:")
print("  - survivors are points that persist under constraints;")
print("  - gaps are produced by their circular geometry;")
print("  - the pattern 6,4,2,4,2,4,6,2 appears with no free parameter.")
print()
assert survivors == [1, 7, 11, 13, 17, 19, 23, 29]
assert gaps == [6, 4, 2, 4, 2, 4, 6, 2]
print("[PASS] coprime residues modulo 30")
print("[PASS] circular gaps produced by survivors")
