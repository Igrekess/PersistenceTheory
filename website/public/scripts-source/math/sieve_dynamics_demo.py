from math import prod

print("DEMONSTRATEUR : LE CRIBLE COMME DYNAMIQUE")
print("=========================================")
print("Question : que fait une contrainte nouvelle au système des survivants ?")
print("On ajoute les canaux 2, puis 3, puis 5, puis 7.")
print("À chaque étape, la densité restante doit être le produit exact des pertes locales.")
print()

primes = [2, 3, 5, 7]
active = []

for p in primes:
    active.append(p)
    modulus = prod(active)
    survivors = [n for n in range(1, modulus + 1) if all(n % q != 0 for q in active)]
    density_exact = len(survivors) / modulus
    density_product = prod(1 - 1 / q for q in active)
    removed = 1 - density_exact
    product_terms = " x ".join(f"(1 - 1/{q})" for q in active)

    print(f"Étape : ajout du canal p={p}")
    print(f"  canaux actifs       = {active}")
    print(f"  période complète M  = {modulus}")
    print(f"  survivants          = {len(survivors)}")
    print(f"  calcul densité      = {len(survivors)}/{modulus}")
    print(f"  produit local       = {product_terms}")
    print(f"  densité observée    = {density_exact:.6f}")
    print(f"  produit théorique   = {density_product:.6f}")
    print(f"  part dissipée       = {removed:.6f}")
    print(f"  premiers survivants = {survivors[:12]}")
    print()
    assert abs(density_exact - density_product) < 1e-12

print("Lecture PT : le crible n'est pas seulement un tri.")
print("C'est une dynamique de contraintes : chaque canal retire une direction,")
print("et les survivants restants portent la structure persistante.")
print("[PASS] la densité suit le produit exact des contraintes")
print("[PASS] ajouter un canal réorganise les survivants")
