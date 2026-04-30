print("DEMONSTRATOR: ANOMALOUS-DIMENSION THRESHOLD")
print("===========================================")
print("Question: why are some channels active while others are only echoes?")
print("Here we visualize the PT criterion: a channel is active if gamma_p > 1/2.")
print("The script uses normalized values to make the cutoff readable.")
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
print("Active-persistence threshold:", threshold)
print()
for p, g in gamma.items():
    status = "ACTIVE" if g > threshold else "ECHO"
    margin = g - threshold
    bar = "#" * max(1, int(g * 20))
    comparison = ">" if g > threshold else "<"
    print(f"p={p:2d}  gamma={g:.5f}  margin={margin:+.5f}  {status:6s}  {bar}")
    print(f"       threshold calculation: gamma_p {comparison} 1/2 because {g:.5f} {comparison} {threshold:.5f}")

print()
print("PT reading:")
print("  - 3, 5, and 7 remain above the threshold;")
print("  - 11 falls below it and becomes the first echo channel;")
print("  - this is not a hidden spatial dimension, but a sensitivity exponent.")
print()
assert all(gamma[p] > threshold for p in [3, 5, 7])
assert all(gamma[p] < threshold for p in [11, 13, 17])
print("[PASS] 3,5,7 above the threshold")
print("[PASS] 11 and beyond switch to echo in this demonstrator")
