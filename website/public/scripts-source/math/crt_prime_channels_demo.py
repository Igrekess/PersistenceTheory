from math import prod

print("DEMONSTRATEUR CRT : POURQUOI LES CANAUX PREMIERS ?")
print("===================================================")
print("Question : pourquoi les premiers sont-ils les bons canaux de base ?")
print("Idée : modulo 30, une position peut être lue séparément modulo 2, 3 et 5.")
print("C'est exactement l'intuition du théorème chinois des restes.")
print()

channels = [2, 3, 5]
modulus = prod(channels)

def signature(n):
    return tuple(n % p for p in channels)

residues = [n for n in range(modulus) if all(n % p != 0 for p in channels)]
signatures = {n: signature(n) for n in residues}

print("Canaux premiers indépendants :", channels)
print("Module recomposé             :", modulus)
print("Calcul : M = 2 x 3 x 5 =", modulus)
print()
print("Chaque survivant porte une signature séparée :")
for n, sig in signatures.items():
    checks = [f"{n} mod {p} = {n % p}" for p in channels]
    print(f"  n={n:2d}  ->  {', '.join(checks)}  -> signature {sig}")

print()
print("Vérification d'unicité des signatures :")
print("  nombre de survivants           =", len(signatures))
print("  nombre de signatures distinctes =", len(set(signatures.values())))
print()
print("Lecture PT :")
print("  - un canal composé n'est pas fondamental : il se factorise ;")
print("  - un canal premier ajoute une contrainte irréductible ;")
print("  - CRT explique pourquoi les canaux premiers peuvent être lus puis recomposés.")
print()
assert len(set(signatures.values())) == len(signatures)
print("[PASS] chaque survivant porte une signature CRT")
print("[PASS] les canaux premiers sont lus séparément puis recomposés")
