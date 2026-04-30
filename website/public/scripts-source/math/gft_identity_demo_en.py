from math import log2

print("GFT DEMONSTRATOR")
print("================")
print("Question: if a system has m possible states, where does its information budget go?")
print("PT answer: it splits exactly into two parts:")
print("  - D_KL: persistent structure, distance from uniformity")
print("  - H   : entropy, remaining uncertainty")
print("The claim tested here is the exact identity: log2(m) = D_KL + H.")
print()

distributions = [
    ("uniform", [0.25, 0.25, 0.25, 0.25]),
    ("slightly structured", [0.40, 0.30, 0.20, 0.10]),
    ("strongly structured", [0.70, 0.10, 0.10, 0.10]),
    ("almost one state", [0.97, 0.01, 0.01, 0.01]),
]

for label, p in distributions:
    m = len(p)
    h = -sum(x * log2(x) for x in p if x > 0)
    dkl = sum(x * log2(m * x) for x in p if x > 0)
    total = log2(m)
    residual = abs(total - (dkl + h))

    print(f"Case: {label}")
    print("  distribution P =", [round(x, 3) for x in p])
    print("  H calculation   = -Σ p_i log2(p_i)")
    for i, x in enumerate(p, start=1):
        print(f"    i={i}: -{x:.3f} log2({x:.3f}) = {-x * log2(x):.6f}")
    print("  D_KL calculation= Σ p_i log2(m p_i)")
    for i, x in enumerate(p, start=1):
        print(f"    i={i}: {x:.3f} log2({m} x {x:.3f}) = {x * log2(m * x):+.6f}")
    print(f"  total budget log2(m)     = {total:.6f} bits")
    print(f"  persistent structure D_KL = {dkl:.6f} bits")
    print(f"  remaining entropy H       = {h:.6f} bits")
    print(f"  check D_KL + H            = {dkl + h:.6f} bits")
    print(f"  numerical residual        = {residual:.2e}")
    if dkl < 1e-9:
        print("  reading: everything is entropy; no structure stands out.")
    elif h < 0.5:
        print("  reading: probability mass is concentrated; persistence dominates.")
    else:
        print("  reading: structure grows as P moves away from uniformity.")
    print()
    assert residual < 1e-12

print("Conclusion: this script uses no fitting.")
print("It shows that the total budget is fixed and only its partition changes.")
print("[PASS] GFT identity verified on all distributions")
