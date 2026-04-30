from math import gcd

print("DEMONSTRATEUR DES SURVIVANTS MODULO 30")
print("======================================")
print("Question : que reste-t-il si l'on impose les contraintes 2, 3 et 5 ?")
print("En PT, les positions qui ne sont éliminées par aucun canal sont les")
print("points de persistance du crible. Les gaps ne sont pas ajoutés ensuite :")
print("ils sont les distances entre ces points survivants.")
print()

modulus = 2 * 3 * 5
survivors = [n for n in range(1, modulus + 1) if gcd(n, modulus) == 1]
gaps = []
for i, n in enumerate(survivors):
    nxt = survivors[(i + 1) % len(survivors)]
    if i == len(survivors) - 1:
        nxt += modulus
    gaps.append(nxt - n)

print(f"Module complet M = 2 x 3 x 5 = {modulus}")
print("Positions testées : 1, 2, ..., 30")
print("Règle de survie   : gcd(n, 30) = 1")
print()
print("Calcul local sur les 30 positions :")
for n in range(1, modulus + 1):
    status = "SURVIT" if gcd(n, modulus) == 1 else "elimine"
    print(f"  n={n:2d}  gcd({n:2d},30)={gcd(n, modulus):2d}  -> {status}")
print()
print("Survivants trouvés :", survivors)
print("Gaps circulaires   :", gaps)
print("Densité exacte     :", f"{len(survivors)}/{modulus}", "=", f"{len(survivors)/modulus:.6f}")
print("Calcul des gaps    :")
for n, gap in zip(survivors, gaps):
    nxt = n + gap
    if nxt > modulus:
        nxt -= modulus
    print(f"  depuis {n:2d} vers {nxt:2d} : gap = {gap}")
print()
print("Lecture PT :")
print("  - les survivants sont les points qui persistent sous les contraintes ;")
print("  - les gaps sont produits par leur géométrie circulaire ;")
print("  - on obtient déjà le motif 6,4,2,4,2,4,6,2 sans paramètre libre.")
print()
assert survivors == [1, 7, 11, 13, 17, 19, 23, 29]
assert gaps == [6, 4, 2, 4, 2, 4, 6, 2]
print("[PASS] résidus copremiers modulo 30")
print("[PASS] gaps circulaires produits par les survivants")
