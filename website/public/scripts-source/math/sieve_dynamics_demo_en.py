from math import prod

print("DEMONSTRATOR: THE SIEVE AS A DYNAMICS")
print("=====================================")
print("Question: what does a new constraint do to the survivor system?")
print("We add channels 2, then 3, then 5, then 7.")
print("At each step, the remaining density must equal the exact product of local losses.")
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

    print(f"Step: add channel p={p}")
    print(f"  active channels    = {active}")
    print(f"  full period M      = {modulus}")
    print(f"  survivors          = {len(survivors)}")
    print(f"  density calculation= {len(survivors)}/{modulus}")
    print(f"  local product      = {product_terms}")
    print(f"  observed density   = {density_exact:.6f}")
    print(f"  theoretical product= {density_product:.6f}")
    print(f"  dissipated share   = {removed:.6f}")
    print(f"  first survivors    = {survivors[:12]}")
    print()
    assert abs(density_exact - density_product) < 1e-12

print("PT reading: the sieve is not just a sorting procedure.")
print("It is a constraint dynamics: each channel removes one direction,")
print("and the remaining survivors carry persistent structure.")
print("[PASS] density follows the exact product of constraints")
print("[PASS] adding a channel reorganizes the survivors")
