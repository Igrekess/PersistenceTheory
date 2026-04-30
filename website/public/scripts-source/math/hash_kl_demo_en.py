from hashlib import sha256
from math import log2

print("DEMONSTRATOR: HASH AND DIVERGENCE FROM UNIFORMITY")
print("=================================================")
print("Question: how can a hash function be read with GFT?")
print("We take 512 toy inputs, then look only at the first SHA-256 half-byte.")
print("If outputs are close to uniform, D_KL should be small.")
print()

def kl_to_uniform(counts):
    total_count = sum(counts)
    m = len(counts)
    dkl = 0.0
    h = 0.0
    for c in counts:
        if c:
            p = c / total_count
            dkl += p * log2(m * p)
            h -= p * log2(p)
    return dkl, h, log2(m)

counts = [0] * 16
for nonce in range(512):
    digest = sha256(f"PT:{nonce}".encode()).digest()
    counts[digest[0] >> 4] += 1

print("Examples of class calculation:")
for nonce in range(5):
    digest = sha256(f"PT:{nonce}".encode()).hexdigest()
    print(f"  input PT:{nonce:<2d} -> SHA-256 {digest[:12]}... -> first half-byte = {digest[0]}")
print()

dkl, h, total = kl_to_uniform(counts)
residual = abs(total - dkl - h)

print("Observed hexadecimal-prefix classes:")
for i, c in enumerate(counts):
    print(f"  class {i:02x} : {c:3d}")
print()
print(f"Total budget log2(16) = {total:.6f} bits")
print("D_KL calculation      = Σ p_i log2(16 p_i)")
print("H calculation         = -Σ p_i log2(p_i)")
print(f"D_KL to uniform        = {dkl:.6f} bits")
print(f"Entropy H             = {h:.6f} bits")
print(f"GFT residual          = {residual:.2e}")
print()
print("PT reading:")
print("  - here D_KL is small: classes are close to uniform;")
print("  - information about the input is strongly dispersed;")
print("  - this toy illustrates the PT language, not a cryptographic security proof.")
assert residual < 1e-12
print("[PASS] GFT partition on toy SHA-256 prefixes")
