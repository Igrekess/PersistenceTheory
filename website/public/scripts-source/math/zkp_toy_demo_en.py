from hashlib import sha256

print("DEMONSTRATOR: TOY PROOF WITHOUT REVEALING")
print("========================================")
print("Question: how can one prove a property without revealing the whole witness?")
print("This script is not a real cryptographic ZKP protocol.")
print("It only illustrates the PT separation: persistent property vs hidden content.")
print()

secret = "hidden-witness"
public_statement = sha256(secret.encode()).hexdigest()

def prove(witness):
    return sha256(witness.encode()).hexdigest()

def verify(statement, proof):
    return statement == proof

proof = prove(secret)
bad_proof = prove("wrong-witness")

print("Private witness      : not shown to the verifier")
print("Public calculation   : proof = SHA-256(witness)")
print("Public statement     :", public_statement[:24] + "...")
print("Proof sent           :", proof[:24] + "...")
print("Bad proof tested     :", bad_proof[:24] + "...")
print()
print("Verification calculation: compare the public statement and the received proof.")
print("Good proof verifies  :", verify(public_statement, proof))
print("Bad proof verifies   :", verify(public_statement, bad_proof))
print()
print("PT reading:")
print("  - the verifiable property persists in the proof;")
print("  - the complete witness is not transmitted;")
print("  - real ZKP cryptography adds randomness, interaction, and formal guarantees.")

assert verify(public_statement, proof)
assert not verify(public_statement, bad_proof)
print("[PASS] property verified without displaying the witness")
