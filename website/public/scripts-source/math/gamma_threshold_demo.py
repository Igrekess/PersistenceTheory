print("DEMONSTRATEUR : SEUIL DES DIMENSIONS ANOMALES")
print("=============================================")
print("Question : pourquoi certains canaux sont actifs et d'autres seulement en écho ?")
print("Ici on visualise le critère PT : un canal est actif si gamma_p > 1/2.")
print("Le script utilise des valeurs normalisées pour rendre la coupure lisible.")
print()

gamma = {
    3: 0.80761,
    5: 0.69632,
    7: 0.59547,
    11: 0.486,
    13: 0.452,
    17: 0.404,
}

threshold = 0.5
print("Seuil de persistance active :", threshold)
print()
for p, g in gamma.items():
    status = "ACTIF" if g > threshold else "ÉCHO"
    margin = g - threshold
    bar = "#" * max(1, int(g * 20))
    comparison = ">" if g > threshold else "<"
    print(f"p={p:2d}  gamma={g:.5f}  marge={margin:+.5f}  {status:5s}  {bar}")
    print(f"       calcul seuil : gamma_p {comparison} 1/2  car {g:.5f} {comparison} {threshold:.5f}")

print()
print("Lecture PT :")
print("  - 3, 5 et 7 restent au-dessus du seuil ;")
print("  - 11 passe sous le seuil et devient le premier canal d'écho ;")
print("  - ce n'est pas une dimension spatiale cachée, mais un exposant de sensibilité.")
print()
assert all(gamma[p] > threshold for p in [3, 5, 7])
assert all(gamma[p] < threshold for p in [11, 13, 17])
print("[PASS] 3,5,7 au-dessus du seuil")
print("[PASS] 11 et au-delà basculent en écho dans ce démonstrateur")
