from math import cos, sin, sqrt

print("DEMONSTRATOR: PRIME SPIRAL")
print("==========================")
print("Question: why are prime spirals pedagogically useful?")
print("They do not prove a law, but they make a structure visible:")
print("sieve survivors do not look like uniform noise.")
print()

def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True

points = []
for n in range(2, 80):
    if is_prime(n):
        r = sqrt(n)
        theta = n * 0.65
        points.append((n, r, theta, r * cos(theta), r * sin(theta)))

print("We place each prime p on a simplified Archimedean spiral:")
print("  radius r = sqrt(p)")
print("  angle theta = 0.65 p")
print()
for n, r, theta, x, y in points[:18]:
    print(f"p={n:2d}  r=sqrt({n})={r:7.4f}  theta=0.65*{n}={theta:7.4f}  x={x:8.4f}  y={y:8.4f}")

print()
print("PT reading:")
print("  - the spiral is a visualization, not a proof;")
print("  - it helps show that survivors have geometry;")
print("  - the demonstrative step then comes from modular classes and the sieve.")
print("[PASS] coordinates generated for prime points")
