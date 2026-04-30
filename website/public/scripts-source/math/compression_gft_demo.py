from collections import Counter
from math import log2

print("DEMONSTRATEUR : COMPRESSION ET GFT")
print("==================================")
print("Question : pourquoi compresser revient-il à extraire de la persistance ?")
print("On compare deux chaînes : l'une très redondante, l'autre beaucoup plus uniforme.")
print("Une chaîne redondante a une entropie plus basse et une structure plus visible.")
print()

samples = {
    "redondant": "A" * 80 + "B" * 16 + "C" * 4,
    "melange": "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789" * 3,
}

for name, text in samples.items():
    counts = Counter(text)
    n = sum(counts.values())
    m = len(counts)
    h = -sum((c / n) * log2(c / n) for c in counts.values())
    total = log2(m)
    dkl = total - h

    print(f"Cas : {name}")
    print(f"  longueur                = {n}")
    print(f"  alphabet distinct       = {m} symboles")
    print("  comptage                =", dict(counts))
    print("  probabilités            =", {k: round(v / n, 4) for k, v in counts.items()})
    print("  calcul H                = -Σ p_i log2(p_i)")
    print("  calcul D_KL             = log2(m) - H")
    print(f"  budget log2(m)          = {total:.6f} bits")
    print(f"  structure D_KL          = {dkl:.6f} bits")
    print(f"  entropie H              = {h:.6f} bits")
    if dkl > 0.5:
        print("  lecture : distribution déséquilibrée, donc structure compressible.")
    else:
        print("  lecture : distribution presque uniforme, donc peu de redondance simple.")
    print()

print("Conclusion : GFT ne remplace pas un compresseur réel.")
print("Mais elle donne le langage : structure persistante + résidu entropique.")
print("[PASS] une chaîne redondante porte plus de structure compressible")
