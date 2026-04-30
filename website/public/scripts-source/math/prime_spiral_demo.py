from math import cos, sin, sqrt

print("DEMONSTRATEUR : SPIRALE DES PREMIERS")
print("====================================")
print("Question : pourquoi les spirales de premiers sont-elles pédagogiques ?")
print("Elles ne prouvent pas une loi, mais elles rendent visible une structure :")
print("les survivants du crible ne ressemblent pas à un bruit uniforme.")
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

print("On place chaque premier p sur une spirale d'Archimède simplifiée :")
print("  rayon r = sqrt(p)")
print("  angle theta = 0.65 p")
print()
for n, r, theta, x, y in points[:18]:
    print(f"p={n:2d}  r=sqrt({n})={r:7.4f}  theta=0.65*{n}={theta:7.4f}  x={x:8.4f}  y={y:8.4f}")

print()
print("Lecture PT :")
print("  - la spirale est une visualisation, pas une démonstration ;")
print("  - elle aide à voir que les survivants ont une géométrie ;")
print("  - l'étape démonstrative vient ensuite par classes modulaires et crible.")
print("[PASS] coordonnées générées pour les points premiers")
