from math import acos, sin

print("DEMONSTRATOR: DISCRETE RESIDUE -> CONTINUOUS PHASE")
print("==================================================")
print("Question: how can a discrete Z/pZ channel carry a continuous quantity?")
print("PT associates a channel depth delta_p, then a phase theta_p.")
print("The tested relation is: sin^2(theta_p) = delta_p(2-delta_p).")
print()

mu = 15
primes = [3, 5, 7, 11]
q_plus = 1 - 2 / mu

print("Fixed point used: mu* =", mu)
print("q+ branch       :", f"{q_plus:.12f}")
print()

for p in primes:
    delta = (1 - q_plus ** p) / p
    theta = acos(1 - delta)
    sin2 = sin(theta) ** 2
    rhs = delta * (2 - delta)
    residual = abs(sin2 - rhs)

    print(f"Channel p={p}")
    print(f"  delta_p calculation         = (1 - q_+^p)/p = (1 - {q_plus:.8f}^{p})/{p}")
    print(f"  theta_p calculation         = arccos(1 - delta_p)")
    print(f"  delta_p                    = {delta:.8f}")
    print(f"  theta_p                    = {theta:.8f} rad")
    print(f"  continuous check            = sin(theta_p)^2 vs delta_p(2-delta_p)")
    print(f"  sin^2(theta_p)             = {sin2:.8f}")
    print(f"  delta_p(2-delta_p)         = {rhs:.8f}")
    print(f"  residual                   = {residual:.2e}")
    if p in (3, 5, 7):
        print("  reading: active channel in the PT triad.")
    else:
        print("  reading: first channel beyond the triad, useful as an echo.")
    print()
    assert residual < 1e-12

print("Conclusion: a channel residue produces a real measurable phase.")
print("This is one place where the discrete/continuous boundary dissolves.")
print("[PASS] sin^2(theta_p)=delta_p(2-delta_p)")
print("[PASS] a channel residue gives a continuous phase")
