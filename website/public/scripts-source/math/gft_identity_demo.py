from math import log2

print("DEMONSTRATEUR GFT")
print("=================")
print("Question : si un système a m possibilités, où va son budget d'information ?")
print("Réponse PT : il se partage exactement en deux parts :")
print("  - D_KL : structure persistante, écart à l'uniforme")
print("  - H    : entropie, incertitude restante")
print("La revendication testée ici est l'identité exacte : log2(m) = D_KL + H.")
print()

distributions = [
    ("uniforme", [0.25, 0.25, 0.25, 0.25]),
    ("legerement structuree", [0.40, 0.30, 0.20, 0.10]),
    ("fortement structuree", [0.70, 0.10, 0.10, 0.10]),
    ("presque un seul etat", [0.97, 0.01, 0.01, 0.01]),
]

for label, p in distributions:
    m = len(p)
    h = -sum(x * log2(x) for x in p if x > 0)
    dkl = sum(x * log2(m * x) for x in p if x > 0)
    total = log2(m)
    residual = abs(total - (dkl + h))

    print(f"Cas : {label}")
    print("  distribution P =", [round(x, 3) for x in p])
    print("  calcul H        = -Σ p_i log2(p_i)")
    for i, x in enumerate(p, start=1):
        print(f"    i={i} : -{x:.3f} log2({x:.3f}) = {-x * log2(x):.6f}")
    print("  calcul D_KL     = Σ p_i log2(m p_i)")
    for i, x in enumerate(p, start=1):
        print(f"    i={i} : {x:.3f} log2({m} x {x:.3f}) = {x * log2(m * x):+.6f}")
    print(f"  budget total log2(m)        = {total:.6f} bits")
    print(f"  structure persistante D_KL  = {dkl:.6f} bits")
    print(f"  entropie restante H         = {h:.6f} bits")
    print(f"  verification D_KL + H       = {dkl + h:.6f} bits")
    print(f"  residu numerique            = {residual:.2e}")
    if dkl < 1e-9:
        print("  lecture : tout est entropie ; aucune structure ne se distingue.")
    elif h < 0.5:
        print("  lecture : la masse est concentree ; la persistance domine.")
    else:
        print("  lecture : la structure augmente quand P s'eloigne de l'uniforme.")
    print()
    assert residual < 1e-12

print("Conclusion : ce script ne fait aucun ajustement.")
print("Il montre que le budget total reste fixe et que seule sa partition change.")
print("[PASS] identite GFT verifiee sur toutes les distributions")
