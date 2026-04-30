from hashlib import sha256
from math import log2

print("DEMONSTRATEUR : HASH ET DIVERGENCE À L'UNIFORME")
print("================================================")
print("Question : comment lire une fonction de hachage avec GFT ?")
print("On prend 512 entrées jouets, puis on regarde seulement le premier demi-octet")
print("du SHA-256. Si les sorties sont proches de l'uniforme, D_KL doit être faible.")
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

print("Exemples de calcul de classe :")
for nonce in range(5):
    digest = sha256(f"PT:{nonce}".encode()).hexdigest()
    print(f"  entree PT:{nonce:<2d} -> SHA-256 {digest[:12]}... -> premier demi-octet = {digest[0]}")
print()

dkl, h, total = kl_to_uniform(counts)
residual = abs(total - dkl - h)

print("Classes observées du préfixe hexadécimal :")
for i, c in enumerate(counts):
    print(f"  classe {i:02x} : {c:3d}")
print()
print(f"Budget total log2(16) = {total:.6f} bits")
print("Calcul D_KL           = Σ p_i log2(16 p_i)")
print("Calcul H              = -Σ p_i log2(p_i)")
print(f"D_KL à l'uniforme     = {dkl:.6f} bits")
print(f"Entropie H            = {h:.6f} bits")
print(f"Résidu GFT            = {residual:.2e}")
print()
print("Lecture PT :")
print("  - ici D_KL est faible : les classes sont proches de l'uniforme ;")
print("  - l'information sur l'entrée est fortement dispersée ;")
print("  - ce jouet illustre la langue PT, pas une preuve de sécurité cryptographique.")
assert residual < 1e-12
print("[PASS] partition GFT sur préfixes SHA-256 jouets")
