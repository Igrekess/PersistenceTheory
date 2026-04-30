from hashlib import sha256

print("DEMONSTRATEUR : PREUVE JOUET SANS RÉVÉLATION")
print("============================================")
print("Question : comment prouver une propriété sans révéler tout le témoin ?")
print("Ce script n'est pas un vrai protocole ZKP cryptographique.")
print("Il illustre seulement la séparation PT : propriété persistante vs contenu caché.")
print()

secret = "temoin-cache"
public_statement = sha256(secret.encode()).hexdigest()

def prove(witness):
    return sha256(witness.encode()).hexdigest()

def verify(statement, proof):
    return statement == proof

proof = prove(secret)
bad_proof = prove("mauvais-temoin")

print("Témoin privé         : non affiché au vérificateur")
print("Calcul public        : preuve = SHA-256(témoin)")
print("Énoncé public        :", public_statement[:24] + "...")
print("Preuve envoyée       :", proof[:24] + "...")
print("Mauvaise preuve test :", bad_proof[:24] + "...")
print()
print("Calcul de vérification : comparer l'énoncé public et la preuve reçue.")
print("Vérification bonne preuve    :", verify(public_statement, proof))
print("Vérification mauvaise preuve :", verify(public_statement, bad_proof))
print()
print("Lecture PT :")
print("  - la propriété vérifiable persiste dans la preuve ;")
print("  - le témoin complet n'est pas transmis ;")
print("  - la vraie cryptographie ZKP ajoute l'aléa, l'interaction et les garanties formelles.")

assert verify(public_statement, proof)
assert not verify(public_statement, bad_proof)
print("[PASS] propriété vérifiée sans afficher le témoin")
