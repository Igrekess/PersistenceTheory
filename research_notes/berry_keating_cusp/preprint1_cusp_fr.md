# L'opérateur de Berry–Keating sur un cusp hyperbolique

## Spectre, extensions auto-adjointes et formule de trace de type Selberg

**Yan Senez** · 2026-05-17 · *draft 1*

---

## Résumé

L'opérateur de Berry–Keating `H = -i(x ∂_x + 1/2)`, candidat usuel à une réalisation spectrale de la conjecture d'Hilbert–Polya, est ici construit canoniquement sur le cusp parabolique d'une courbe spectrale hyperbolique `Σ` de rayon `R = √3` et de courbure scalaire `K = -1/3`. La présence d'un poids cuspidal `e^{-n η / R}` (où `n` est la dimension du sous-groupe parabolique) impose un terme géométrique additif `-n η / (2 R)` sans lequel l'opérateur n'est pas symétrique. On démontre que cet opérateur cuspidal `H_n^{\text{cusp}}` est unitairement équivalent à `-i ∂_ξ` sur un intervalle de longueur `r = \log(η_{\max} / η_{\min})`, avec indices de défaut `(1, 1)` et famille `U(1)` d'extensions auto-adjointes. Une condition de symétrie spectrale anti-unitaire restreint cette famille à deux choix, et l'exclusion du mode marginal sélectionne l'extension antipériodique. Sous une troncature dynamique `r(γ) = \log(γ / (2 π))`, le nombre de valeurs propres positives satisfait `N_+(γ) = (γ / 2π) \log(γ / (2π)) + O(1)`, soit la densité asymptotique de Riemann–von Mangoldt. On dérive enfin une formule de trace de type Selberg en trois termes, dans laquelle la matrice de scattering du cusp prend une forme produit `φ(s) = ζ_+(s) ζ_-(s)` portée par deux branches dynamiques de la courbe spectrale.

---

## 1. Introduction

L'opérateur formel `H = -i(x ∂_x + 1/2)` sur `L²(ℝ_+, dx)` est le candidat le plus simple pour une réalisation spectrale de la conjecture d'Hilbert–Polya, telle qu'examinée par Berry et Keating [1]. C'est l'opérateur de dilatation symétrisé, l'analogue mathématique d'un produit `x p` quantifié de Weyl. Trois obstacles s'opposent à son utilisation directe :

1. **Domaine** : `H` n'est pas essentiellement auto-adjoint sur `C_c^∞(ℝ_+)`. Une famille d'extensions à un paramètre existe, mais aucun argument de premier principe ne distingue l'une d'entre elles.

2. **Spectre continu** : sur la demi-droite non tronquée, `H` a un spectre continu. Pour produire un spectre discret comparable aux zéros de `ζ`, une troncature est nécessaire ; le choix de cette troncature n'est pas dicté par l'opérateur lui-même.

3. **Formule de trace** : l'identification du spectre discret avec les ordonnées des zéros non triviaux de `ζ` exige une formule de trace explicite reliant les valeurs propres aux orbites primitives de la dynamique de dilatation, longueurs `log p` pour `p` parcourant les premiers.

Cette note construit un cadre géométrique pour `H` qui répond aux trois points. La construction part d'une courbe spectrale hyperbolique `Σ` de rayon `R = √3`, dont le bord asymptotique est un cusp parabolique. L'origine de cette courbe est externe : elle provient d'un programme de construction d'une variété de moduli spectrale plus large [2], dont les détails ne sont pas requis ici. Ce qui importe pour la présente note est l'**énoncé géométrique minimal** : une métrique hyperbolique près du cusp, une dimension de sous-groupe parabolique `n`, et un choix de coordonnée canonique sur la coordonnée radiale du cusp.

Sur cette donnée, on construit en §3 un opérateur cuspidal `H_n^{\text{cusp}}` qui généralise `H` en absorbant la courbure du cusp via un terme de transport parallèle. §4 calcule les indices de défaut et exhibe la famille `U(1)` des extensions auto-adjointes. §5 dégage le rôle de la troncature et donne la densité spectrale asymptotique. §6 dérive une formule de trace de type Selberg. §7 discute les limites et perspectives.

**Théorèmes principaux.**

- *Théorème 3.1 (symétrie cuspidale).* `H_n^{\text{cusp}}` sur le domaine `C_c^∞(η_{\min}, η_{\max})` est un opérateur symétrique densément défini.
- *Théorème 4.1 (indices de défaut).* Les indices de défaut sont `(n_+, n_-) = (1, 1)`. La famille des extensions auto-adjointes est `\{\tilde H_θ\}_{θ ∈ [0, 2π)}`, paramétrée par une condition de bord twisted.
- *Théorème 4.2 (extension distinguée).* Sous (i) symétrie spectrale anti-unitaire et (ii) exclusion du mode marginal `λ = 0`, l'unique extension est l'antipériodique `θ = π`.
- *Théorème 5.1 (densité de Riemann–von Mangoldt).* Sous la troncature dynamique `r(γ) = \log(γ/(2π))`, le compteur des valeurs propres positives jusqu'à `γ` satisfait `N_+(γ) = (γ / 2π) \log(γ / (2π)) + O(1)`.
- *Théorème 6.1 (formule de trace).* La distribution `\sum_n h(γ_n)` (où `γ_n` parcourt le demi-spectre régularisé) se décompose en quatre termes (Weyl, géométrique, parabolique, archimédien), dans laquelle la matrice de scattering du cusp admet la forme produit `φ(s) = ζ_+(s) ζ_-(s)`.

---

## 2. Géométrie du cusp

Soit `Σ` une surface hyperbolique avec un cusp parabolique au bord asymptotique. Près de ce cusp, la métrique prend la **forme demi-espace de Poincaré généralisée**
```
ds² = dη² + e^{-2 η / R} Σ_{a = 1}^{n} dX_a²,                       (2.1)
```
où :
- `η > 0` est la coordonnée radiale, avec `η → ∞` au cusp ;
- `X_a ∈ ℝ / L_a ℤ` (`a = 1, …, n`) sont les coordonnées horocycliques, compactifiées par le sous-groupe parabolique du cusp ;
- `R > 0` est le rayon hyperbolique du cusp ;
- `n` est la dimension du sous-groupe parabolique (et du tore horocyclique).

La courbure sectionnelle de (2.1) est constante `K = -1 / R²`. Pour `R = √3`, `K = -1/3`. Pour `R = 1`, on retrouve le demi-espace standard `ℍ^{n+1}` de courbure `-1`. Dans la suite on travaillera avec `R = √3` et `n = 3`, valeurs spécifiques au cadre [2], mais les énoncés des sections 3-5 valent pour tout choix `(R, n)` avec `R > 0` et `n ∈ ℕ`.

L'**élément de volume cuspidal** est
```
d \text{Vol} = e^{-n η / R} dη \prod_{a = 1}^{n} dX_a.
```
L'espace de Hilbert `L²(Σ, d\text{Vol})` se décompose en série de Fourier sur les directions horocycliques :
```
f(η, X) = \sum_{k ∈ ℤ^n} f_k(η) e^{2 π i \, k · X / L},
```
avec norme séparable
```
\|f\|² = \sum_k \|f_k\|²_{ℋ_n},\quad
ℋ_n := L²(ℝ_η, e^{-n η / R} dη).                                  (2.2)
```
Le **sous-espace cuspidal radial** `ℋ_n` (mode `k = 0`) est celui sur lequel l'opérateur de dilatation agit naturellement.

---

## 3. L'opérateur cuspidal de Berry–Keating

### 3.1 Définition

On considère sur `ℋ_n` l'opérateur
```
H_n^{\text{cusp}} := -i \bigl(η \partial_η + 1/2 - n η / (2 R)\bigr).  (3.1)
```
Le premier terme `-i(η ∂_η + 1/2)` est la version symétrisée standard du générateur de dilatation. Le terme additif `-n η / (2 R)` est une **correction géométrique** qui restore la symétrie en présence du poids cuspidal `e^{-n η / R}` de `ℋ_n`. Pour `n = 0`, on récupère exactement l'opérateur de Berry–Keating de [1].

Le domaine test canonique est
```
\mathcal{D}_0 := C_c^∞\bigl((η_{\min}, η_{\max})\bigr),
```
où `η_{\min}, η_{\max}` sont des bornes à fixer en §5. `\mathcal{D}_0` est dense dans `ℋ_n` (théorème d'approximation standard).

### 3.2 Symétrie

**Théorème 3.1.** *Pour tout `n ≥ 0`, l'opérateur `H_n^{\text{cusp}}` sur `\mathcal{D}_0` est symétrique : pour tous `φ, ψ ∈ \mathcal{D}_0`,*
```
\langle φ, H_n^{\text{cusp}} ψ \rangle_{ℋ_n}
= \langle H_n^{\text{cusp}} φ, ψ \rangle_{ℋ_n}.
```

*Démonstration.* Notons `w(η) := e^{-n η / R}` le poids cuspidal. On a
```
\langle φ, H_n^{\text{cusp}} ψ \rangle
= -i \int \bar φ \bigl(η ψ' + ψ/2 - (n η / (2 R)) ψ\bigr) w(η) dη.
```
Une intégration par parties sur le terme `η ψ'`, en utilisant que `(η w)' = w - (n η / R) w`, donne après simplification
```
\langle φ, H_n^{\text{cusp}} ψ \rangle
- \langle H_n^{\text{cusp}} φ, ψ \rangle
= 0,
```
le terme additif `-i (n η / (2 R))` étant exactement le facteur qui compense le défaut généré par la dérivation de `w`. Les termes de bord s'annulent car `φ, ψ` sont à support compact dans l'ouvert `(η_{\min}, η_{\max})`. ∎

### 3.3 Fonctions propres formelles

L'équation aux valeurs propres `H_n^{\text{cusp}} f = λ f` donne (pour `λ ∈ ℝ`)
```
η f'(η) + f(η)/2 - (n η / (2 R)) f(η) = i λ f(η),
```
dont la solution est
```
f_λ^{(n)}(η) = η^{i λ - 1/2} e^{n η / (2 R)},                       (3.2)
```
avec norme cuspidale
```
\|f_λ^{(n)}\|²_{ℋ_n}
= \int_{η_{\min}}^{η_{\max}} η^{-1} e^{n η / R} e^{-n η / R} dη
= \log(η_{\max} / η_{\min}).                                         (3.3)
```
La norme est **indépendante de `λ` et de `n`**. Le facteur de courbure `e^{n η / (2 R)}` compense exactement le poids cuspidal `e^{-n η / R}` via un demi-facteur, rétablissant la **mesure de Haar multiplicative** `dη / η`. C'est la raison structurelle du shift `+1/2` : il est dicté par cette invariance.

### 3.4 Équivalence unitaire avec la translation

Soit le changement de variable Mellin
```
ξ = \log(η / η_{\min}) \in [0, r],\quad r = \log(η_{\max} / η_{\min}),
```
et la jauge `U_n : ℋ_n → L²([0, r], dξ)` définie par
```
(U_n f)(ξ) = η(ξ)^{1/2} e^{-n η(ξ) / (2 R)} f(η(ξ)),\quad η(ξ) = η_{\min} e^ξ.
```
**Proposition 3.2.** *`U_n` est unitaire, et*
```
U_n H_n^{\text{cusp}} U_n^{-1} = -i \partial_ξ
\quad \text{sur } L²([0, r], dξ).                                   (3.4)
```
*Démonstration.* Le changement de variable `η = η_{\min} e^ξ` donne `dη = η dξ`, donc
```
\int |f|² e^{-n η / R} dη = \int |f|² η e^{-n η / R} dξ
= \int η^{1/2} e^{-n η / (2 R)} |f|² η^{1/2} e^{-n η / (2 R)} dξ
= \int |U_n f|² dξ,
```
ce qui prouve l'unitarité. L'action de `H_n^{\text{cusp}}` se réécrit alors par calcul direct comme `-i ∂_ξ` sur les images ; le terme `+1/2` provient du Jacobien de la transformation, et le terme `-n η / (2 R)` est absorbé exactement par la jauge `e^{-n η / (2 R)}`. ∎

**Corollaire 3.3.** *Le spectre de `H_n^{\text{cusp}}` est invariant en `n` et égal au spectre de l'opérateur de translation `-i ∂_ξ` sur `L²([0, r], dξ)`. La correction de courbure `-n η / (2 R)` ne fait que **passer entre jauges** ; elle n'agit pas sur la classe spectrale.*

---

## 4. Extensions auto-adjointes

Par la Proposition 3.2, il suffit d'étudier `-i ∂_ξ` sur `L²([0, r], dξ)` avec domaine `\mathcal{D}_0 = C_c^∞((0, r))`.

### 4.1 Indices de défaut

L'adjoint `(-i ∂_ξ)^*` agit également comme `-i ∂_ξ` mais sur le domaine maximal `H^1([0, r])` sans condition de bord. Les espaces de défaut sont
```
K_± = \ker((-i ∂_ξ)^* \mp i\, \text{Id})
= \{g ∈ L²([0,r]) : g'(ξ) = ∓ g(ξ)\}.
```
Les solutions sont `g(ξ) = c \, e^{∓ ξ}`, qui appartiennent à `L²([0, r])` car l'intervalle est borné. Donc `\dim K_+ = \dim K_- = 1`.

**Théorème 4.1.** *L'opérateur `-i ∂_ξ` sur `\mathcal{D}_0` a indices de défaut `(n_+, n_-) = (1, 1)`. Il n'est donc pas essentiellement auto-adjoint, et la famille des extensions auto-adjointes est paramétrée par `U(1)`.*

### 4.2 Extensions paramétrées par une phase

Par la formule de von Neumann (Reed–Simon vol. II, Thm. X.2 [3]), l'extension `\tilde H_θ` correspondant à la phase `e^{i θ}` a pour domaine
```
\mathcal{D}(\tilde H_θ) = \{g ∈ H^1([0, r]) : g(r) = e^{i θ} g(0)\}. (4.1)
```
Pour `θ = 0` la condition est périodique ; pour `θ = π`, antipériodique.

### 4.3 Spectre des extensions

Les fonctions propres `\tilde H_θ g = λ g` sont `g(ξ) = e^{i λ ξ}` ; la condition (4.1) donne `e^{i λ r} = e^{i θ}`, soit
```
σ(\tilde H_θ) = \bigl\{ (θ + 2 π k) / r : k ∈ ℤ \bigr\}.            (4.2)
```

### 4.4 Sélection par symétrie spectrale anti-unitaire

Sur `L²([0, r], dξ)` on dispose de deux involutions naturelles :
- la **conjugaison complexe** `T : g(ξ) ↦ \overline{g(ξ)}`, antilinéaire, `T² = \text{Id}` ;
- la **réflexion spatiale** `R : g(ξ) ↦ g(r - ξ)`, linéaire, `R² = \text{Id}`.

Sur l'opérateur `-i ∂_ξ` :
- `T` flippe `-i ↔ +i` (conjugaison) sans agir sur `∂_ξ`, donc `T \, H = -H \, T` (**anti-commutation**) ;
- `R` flippe `∂_ξ ↔ -∂_ξ` (chain rule) sans agir sur `-i`, donc `R \, H = -H \, R` (**anti-commutation**) ;
- leur produit `J := T ∘ R` cumule les deux signes, donc `J \, H = +H \, J` (**commutation** : J préserve le spectre sans le flipper, et **n'impose aucune contrainte** sur `θ`).

C'est `T` (ou de façon équivalente `R`) qui réalise la symétrie spectrale `λ → -λ`. Sur les fonctions propres `g_λ(ξ) = e^{i λ ξ}` :
```
T g_λ = e^{-i λ ξ} = g_{-λ}.                                          (4.3)
```

Pour qu'une extension `\tilde H_θ` soit `T`-invariante (au sens où `T \, D(\tilde H_θ) = D(\tilde H_θ)`), la condition de bord doit être stable sous `T`. Or
```
T g(r) = \overline{g(r)} = \overline{e^{i θ} g(0)} = e^{-i θ} \, T g(0),
```
donc `T g ∈ D(\tilde H_{-θ})`. L'invariance impose `θ = -θ \bmod 2π`, soit
```
θ ∈ \{0, π\}.                                                         (4.4)
```

### 4.5 Élimination du mode marginal

L'extension `θ = 0` admet `λ_0 = 0` comme valeur propre, avec fonction propre constante. Or :
1. Sous la jauge inverse `U_n^{-1}`, le mode `g_0 = \text{cste}` correspond à `f_0(η) = c \, η^{-1/2} e^{n η / (2 R)}` qui croît exponentiellement pour `n > 0` ; sur intervalle tronqué il reste admissible, mais croît en module.
2. Dans le cadre d'application visé (zéros non triviaux de `ζ`, cf. §5 et [2]), aucun `γ_n` n'est nul : par le théorème classique `ζ(1/2) ≠ 0` (Titchmarsh §10.2, valeur explicite `ζ(1/2) ≈ -1{,}4604`), il n'existe pas de zéro de `ζ` à `s = 1/2`. Sous l'identification spectrale (établie hors de cette note), le mode `λ = 0` n'a pas de contrepartie, et l'extension périodique introduirait un mode parasite.

**Théorème 4.2.** *Sous (i) symétrie spectrale `λ → -λ` (involution `T`, eq. (4.3-4.4)) et (ii) exclusion de `λ = 0` (équivalente à `ζ(1/2) ≠ 0`), l'unique extension auto-adjointe distinguée est l'antipériodique*
```
\tilde H := \tilde H_{θ = π},\quad
σ(\tilde H) = \bigl\{ (2k + 1) π / r : k ∈ ℤ \bigr\}.                 (4.5)
```

---

## 5. Spectre discret et densité de Riemann–von Mangoldt

### 5.1 Densité brute

Le nombre de valeurs propres de `\tilde H` dans `[-Λ, Λ]` est
```
N(Λ) = \#\{k ∈ ℤ : |(2k+1)π/r| ≤ Λ\} \sim r Λ / π,
```
soit densité constante `ρ = r / π` indépendante de `Λ`. Le **demi-spectre** positif `λ > 0` a densité `r / (2π)`.

### 5.2 Troncature dynamique et densité asymptotique

Pour comparer le demi-spectre aux ordonnées `\{γ_n\}` des zéros non triviaux de `ζ` (qui ont densité asymptotique `\log γ / (2π)`), on impose que `r` croisse logarithmiquement avec une échelle d'observation `γ`. Une motivation est la suivante : sur le cusp de `Σ`, une **cellule de phase canonique** (analogue cuspidal de la cellule de Planck) impose une borne supérieure adaptative
```
η_{\max}(γ) = √R \, γ / √{2π}                                       (5.1)
```
en fonction de l'échelle `γ`. Pour `η_{\min}` constant (par exemple `√{2 π R}`), on obtient
```
r(γ) = \log(η_{\max}(γ) / η_{\min}) = \log(γ / (2π)) + O(1).        (5.2)
```
La motivation physique de (5.1) est externe à cette note (cf. [2]) ; on peut aussi l'adopter comme **postulat de troncature** dont les conséquences spectrales sont l'objet de l'étude.

**Théorème 5.1.** *Sous la troncature (5.2), le nombre de valeurs propres positives du spectre antipériodique (4.5) jusqu'à `γ` satisfait*
```
N_+(γ) = (γ / (2π)) \log(γ / (2π)) + O(1).                          (5.3)
```
C'est la densité asymptotique de Riemann–von Mangoldt pour les zéros non triviaux de `ζ` (Titchmarsh [4], chap. 9).

*Démonstration.* `N_+(γ) = r(γ) γ / (2π) = \log(γ/(2π)) γ / (2π)`, modulo un terme `O(1)` provenant de la borne `η_{\min}` constante et de la parité du spectre. ∎

### 5.3 Le terme constant et la phase de Maslov

La densité classique de Riemann–von Mangoldt a une correction sous-dominante :
```
N_{\text{RvM}}(γ) = (γ / (2π)) \log(γ / (2π e)) + O(1)
= (γ / (2π)) \log(γ / (2π)) - γ / (2π) + O(1).
```
Le terme `-γ / (2π)` est la **phase de Maslov** d'une trajectoire classique fermée. Dans le contexte présent, elle est reliée à l'antipériodicité (`θ = π`) : un demi-tour de phase par point de retournement, deux points de retournement (au cusp et à la coupure inférieure), donne `θ_{\text{Maslov}} = -π = θ - 0`. Cette identification est cohérente mais ne donne pas le coefficient exact `-γ/(2π)` ; on renvoie à [2] pour une analyse Bohr–Sommerfeld détaillée.

---

## 6. Une formule de trace de type Selberg

### 6.1 Forme générale (Iwaniec [5], chap. 10)

Pour un opérateur scalaire sur une variété hyperbolique à cusp de rang 1, la formule de trace de Selberg s'écrit, pour une fonction test admissible `h(λ)` (paire, analytique dans `|\Im λ| < 1/2 + ε`, décroissant plus vite que `1 / λ^2`) :
```
\sum_n h(λ_n) = T_W(h) + T_G(h) + T_P(h) + T_A(h),                  (6.1)
```
avec :
- `T_W(h) = \text{Vol} · \int h(λ) ρ_W(λ) dλ` : terme de Weyl, densité d'état d'aire ;
- `T_G(h) = \sum_{\text{orbites primitives } γ} (\ell(γ) / (1 - e^{-\ell(γ)})) \hat h(\ell(γ))` : terme géométrique, somme sur orbites périodiques de longueur `\ell` ;
- `T_P(h) = (1 / 4π) \int \hat h(τ) \, (φ'/φ)(1/2 + i τ) dτ` : terme parabolique, intégrale logarithmique de la **matrice de scattering du cusp** `φ(s)` ;
- `T_A(h)` : terme archimédien, Γ-corrections.

### 6.2 Adaptation à `\tilde H`

Le flot Hamiltonien associé à `\tilde H` (en représentation Mellin-conjuguée à `H_n^{\text{cusp}}`) est le flot de dilatation `(x, t) ↦ (e^t x, t)`. Ses orbites primitives correspondent aux multiples logarithmiques d'un sous-ensemble distingué de réels positifs. Dans le cadre [2], ce sous-ensemble est l'ensemble `\{\log p : p \text{ premier}\}`, soit les longueurs `\ell_p = \log p`. Le terme géométrique prend alors la forme
```
T_G(h) = \sum_p (\log p / \sqrt p) \, \hat h(\log p),                (6.2)
```
qui est, modulo normalisation, la **somme de von Mangoldt** transformée de Fourier de `h`. Sa version analytique compacte est
```
T_G(s) = -\frac{d \log ζ(s)}{ds} = \sum_p \frac{\log p}{p^s - 1}.   (6.3)
```

### 6.3 Matrice de scattering du cusp

Le cadre [2] munit la courbe `Σ` de **deux branches dynamiques distinguées**, indexées par deux paramètres `q_+, q_-` issus d'une bifurcation interne à la dynamique cuspidale. À chaque branche est associée une **fonction de coefficients** sur l'ensemble des premiers :
```
\delta_\pm(p) = (1 - q_\pm^p) / p,\quad
a_p^\pm = \delta_\pm(p) \, (2 - \delta_\pm(p)),                       (6.4)
```
et une **fonction zêta associée** définie par produit eulérien
```
ζ_\pm(s) := \exp\!\Bigl(\, \sum_p a_p^\pm \, p^{-s}\,\Bigr).        (6.5)
```
Ces fonctions `ζ_\pm` sont les analogues PT-internes de la fonction zêta de Riemann ; elles sont régulières dans `\Re s > 1` et admettent une continuation analytique aux pôles près. Leur produit
```
φ(s) := ζ_+(s) \, ζ_-(s)                                             (6.6)
```
joue le rôle de la **matrice de scattering** du cusp de `Σ`.

### 6.4 Décomposition Selberg-PT

**Théorème 6.1.** *Pour `h` test admissible, le spectre régularisé de `\tilde H` satisfait*
```
\sum_n h(γ_n) = T_W(h) + T_G(h) + T_P(h) + T_A(h),                  (6.7)
```
*avec :*

- *terme de Weyl `T_W(h) = (r / 2π) \int h(λ) dλ`,*
- *terme géométrique `T_G` donné par (6.2)/(6.3),*
- *terme parabolique*
```
T_P(h) = \frac{1}{4π} \int_{-∞}^{∞} \hat h(τ) \, \frac{φ'(1/2 + i τ)}{φ(1/2 + i τ)} \, dτ
       = \sum_p A_p \, \frac{\log p}{\sqrt p} \, \hat h(\log p),
```
*où `A_p = a_p^+ + a_p^-`,*
- *terme archimédien `T_A` provenant des Γ-corrections de `ζ_\pm` et du pôle de `ζ` en `s = 1`.*

**Lemme 6.2 (identité analytique).** *Pour `\Re s > 1`,*
```
\sum_p \log p \, \Bigl[ \frac{1}{p^s - 1} - A_p \, p^{-s} \Bigr]
= -\frac{d \log ζ(s)}{ds} + \frac{d \log φ(s)}{ds},
```
*qui est l'expression compacte de `T_G + T_P` au sens distributionnel.*

*Démonstration de 6.1.* Suit de la dérivation logarithmique de `R(s) := ζ(s) / φ(s)`, des conventions de Fourier pour `h, \hat h`, et de la régularisation Hadamard du contour de Hankel autour des zéros et pôles de `R`. Les détails techniques sont reportés à [2]. ∎

### 6.5 Identification

La structure (6.6) avec deux branches `q_\pm` est spécifique au cadre [2] ; mais l'**identification structurelle** scattering matrix `=` produit de fonctions zêta apparentées est, à notre connaissance, la **première réalisation explicite** d'une matrice de scattering cuspidale en termes de produits eulériens dérivés. Elle suggère une route concrète vers la formule explicite de Riemann–Weil : intégrer (6.7) contre des fonctions test bien choisies devrait reproduire les sommes classiques sur premiers.

### 6.6 Décomposition par tours de spirale log-polaire

Le terme géométrique `T_G(h)` de (6.2) admet une **décomposition canonique par tours de la spirale d'Archimède log-polaire** qui fait apparaître une cascade de seuils universels en `√(π(k+1))`. Cette structure géométrique, qui n'utilise que la densité PNT et la géométrie de la spirale, complète la formule de trace par une signature de la **période fondamentale `2π`** des coordonnées Mellin.

#### 6.6.1 Coordonnées log-polaires et tours

En coordonnées log-polaires `(r, θ) ↦ (\ln r, θ)`, la suite des entiers `n` est placée à angle `θ = \ln n`, et la spirale dégénère en une représentation linéaire des **décades multiplicatives** sous forme additive. La période angulaire `2π` engendre une partition canonique des premiers en **tours** :
```
\text{turn}_k := \{p \text{ premier} : \ln p \in [2π k, 2π (k+1))\},   \quad k \in \mathbb N.   (6.8)
```

Les premiers tours sont :
- `turn_0` : `p \in [1, e^{2π}) = [1, 535)`, soit 99 premiers ;
- `turn_1` : `p \in [535, e^{4π}) = [535, 286\,751)`, soit 24 877 premiers ;
- `turn_2` : `p \in [286\,751, e^{6π}) \approx [2.9 \times 10^5, 1.5 \times 10^8)`, soit `\approx 8.6 \times 10^6` premiers ;
- `turn_k` (k grand) : contient `\sim (e^{2π} - 1) e^{2πk}/(2π(k+1)) \approx 85 \cdot e^{2πk}/(k+1)` premiers par PNT.

**Observation visuelle** : le seuil `\ln 23 \approx 3.135 \approx π` (à 0.2 %, via `e^π \approx 23.14`) place les premiers `{3, 5, 7, 11, 13, 17, 19, 23}` exactement dans le premier **demi-tour** `θ \in [0, π]`. Pour le cadre PT, ceci recouvre la cascade complète {actif ∪ écho ∪ super-écho}.

#### 6.6.2 Décomposition par tours du terme géométrique

Pour une fonction test admissible `h`, le terme géométrique se décompose en somme orthogonale par tours :
```
T_G(h) = 2 \sum_p \frac{\log p}{\sqrt p} \, \hat h(\log p)
       = \sum_{k \geq 0} T_G^{turn_k}(h),
\quad T_G^{turn_k}(h) := 2 \sum_{p \in turn_k} \frac{\log p}{\sqrt p} \, \hat h(\log p).
\tag{6.9}
```

(Pour la formule explicite de Riemann–Weil complète, on étend la somme aux puissances de premiers `p^j`, `j \geq 1` ; pour `h` à support de transformée bornée et primes `p \in turn_k` avec `k \geq 1`, les contributions `j \geq 2` sont négligeables. Nous écrivons (6.9) pour `j = 1` ; l'extension est immédiate.)

#### 6.6.3 Théorème géométrique 6.3

**Théorème 6.3 (Identité de spirale).** *Soit `h(γ) = e^{-σ^2 γ^2}` (`σ > 0`) la famille des fonctions test gaussiennes admissibles, dont la transformée de Fourier est `\hat h(t) = (\sqrt π / σ) e^{-t^2/(4σ^2)}`. Dans l'approximation de densité PNT (`\# \{p \leq x\} \sim x / \log x`), on a*
```
T_G^{turn_k}(σ) = T_G^{turn_0}(σ) \quad \Longleftrightarrow \quad σ^2 = π (k + 1).      (6.10)
```

*En particulier, le seuil critique où la `k`-ième tour de la spirale égale en contribution la tour zéro est*
```
\boxed{\quad σ_{\mathrm{crit}}^{(k)} = \sqrt{π (k+1)} \quad}                          (6.11)
```

*— une cascade arithmétique en `π` propre à la période fondamentale de la spirale log-polaire.*

**Démonstration.** Dans la limite continue PNT, en substituant `u = \log p`, `dp/\log p \to e^u du / u`, le facteur intégrand `(\log p / \sqrt p) \cdot dp/\log p` devient `e^{u/2} du`. Pour la transformée gaussienne `\hat h(u) = (\sqrt π / σ) e^{-u^2/(4σ^2)}`, on a donc :
```
T_G^{turn_k}(σ) \approx \frac{2}{σ \sqrt π} \int_{2π k}^{2π (k+1)} e^{u/2 - u^2/(4σ^2)} du.
```
Complétant le carré dans l'exposant, `u/2 - u^2/(4σ^2) = σ^2/4 - (u - σ^2)^2/(4σ^2)`, d'où :
```
T_G^{turn_k}(σ) = \frac{e^{σ^2/4}}{σ \sqrt π} \int_{2π k}^{2π (k+1)} e^{-(u - σ^2)^2/(4σ^2)} du.   (6.12)
```
Posons `I_k(σ) := \int_{2π k}^{2π (k+1)} e^{-(u - σ^2)^2/(4σ^2)} du`. Le ratio `T_G^{turn_k} / T_G^{turn_0} = I_k / I_0`. Sur substitution `v = u - σ^2`, `I_0(σ) = \int_{-σ^2}^{2π - σ^2} e^{-v^2/(4σ^2)} dv` et `I_k(σ) = \int_{2π k - σ^2}^{2π(k+1) - σ^2} e^{-v^2/(4σ^2)} dv`. Par parité de la gaussienne `e^{-v^2/(4σ^2)}`, on a `I_0(σ) = I_k(σ)` ssi les intervalles d'intégration sont reliés par la symétrie `v \leftrightarrow -v`. Cette condition se réduit à
```
-σ^2 = -(2π (k+1) - σ^2) \quad \text{et} \quad 2π - σ^2 = -(2π k - σ^2),
```
les deux donnant `σ^2 = π (k+1)`. ∎

#### 6.6.4 Triple validation numérique

L'identité (6.10)-(6.11) a été validée numériquement sur trois ordres distincts `k \in \{1, 2, 3\}` avec couverture totale des tours respectives (turn 0 toujours inclus, turn `k` enumérée intégralement). Les écarts relatifs entre prédiction et observation décroissent avec `k`, cohérents avec une convergence vers l'identité continue dans la limite PNT.

| `k` | `σ_{\mathrm{crit}}^{(k)} = \sqrt{π(k+1)}` | `σ` observé | écart | `\#` premiers turn `k` |
|---|---|---|---|---|
| `1` | `\sqrt{2π} \approx 2.5066` | `2.4909` | `-0.626 \%` | `24\,877` |
| `2` | `\sqrt{3π} \approx 3.0700` | `3.0630` | `-0.228 \%` | `\approx 8.6 \times 10^6` |
| `3` | `\sqrt{4π} \approx 3.5449` | `3.5419` | `-0.085 \%` | `\approx 3.4 \times 10^9` |

(Scripts `test_w7bis_extended_sieve.py` pour `k = 1, 2` ; `test_w7ter_primesieve_turn3.py` pour `k = 3` via primesieve CLI streamé, 12.4 min, 83 chunks de range `10^9`.)

La **décroissance des écarts** (`6.3 \times 10^{-3} \to 2.3 \times 10^{-3} \to 8.5 \times 10^{-4}`) suit approximativement `O(1/\log(\#\text{primes}))`, cohérent avec la fluctuation Mertens-Riemann au-delà de PNT pur.

#### 6.6.5 Lien avec le cadre cuspidal

Le scale `σ_{\mathrm{crit}}^{(1)} = \sqrt{2π}` n'utilise que la géométrie log-polaire et PNT — il ne dépend ni du paramètre Berry-Keating `r`, ni de la cascade PT `{q_+, q_-}`. Il s'agit donc d'un **invariant universel de la formule explicite de Weil** sous la décomposition par tours de spirale. La structure (6.11) imprime sur la somme géométrique `T_G` une **cascade arithmétique en π** parallèle à la cascade arithmétique en entiers de la formule explicite (orbites primitives `\ell_p = \log p`, longueurs additives `k \log p` pour les puissances).

**Conséquence pour la formule de trace (6.7)** : la somme `T_G(h)` est, pour `h` gaussienne de paramètre `σ`, un **paquet gaussien en index de tour `k`** centré sur
```
k_{\mathrm{peak}}(σ) = \frac{σ^2}{2π} - \frac{1}{2}.   (6.13)
```
Pour `σ \to ∞`, `k_{\mathrm{peak}} \to ∞`, donc la masse de `T_G` se déplace vers les hautes tours de la spirale. Cette diffusion non-bornée reproduit, sous prisme géométrique, la **non-localité fondamentale** de la condition de positivité de Weil [9] : aucune somme finie de tours ne borne `T_G(σ)` pour `σ` arbitraire, et aucune information locale (cascade `{3, 5, 7}` du cœur PT, ou tours initiales de la spirale) ne contrôle l'asymptotique.

**Statut épistémique** : le Théorème 6.3 est une identité **inconditionnelle** dans la limite continue PNT, prouvée par parité gaussienne et triple-validée numériquement. Il ne franchit pas RH (la non-localité subsiste), mais fournit une **signature géométrique propre** de la période `2π` dans la formule explicite. Il est par construction **PT-compatible** mais **non PT-dépendant** : la décomposition par tours est définie par la seule géométrie log-polaire des premiers.

---

## 7. Discussion et perspectives

### 7.1 Le gap résiduel

Cette note construit un cadre spectral propre pour l'opérateur de Berry–Keating sur un cusp hyperbolique. **Elle ne démontre pas la conjecture d'Hilbert–Polya**. Le spectre brut de `\tilde H` est arithmétiquement régulier `\{(2k+1)π/r\}`, alors que les ordonnées `\{γ_n\}` des zéros de `ζ` sont irrégulières (statistiques GUE conjecturées). Le pont entre les deux exige une **régularisation** au sens distributionnel de la formule de trace (6.7), suivie d'une identification stricte des pôles régularisés avec `\{γ_n\}`. Cette étape reste le gap classique du programme d'Hilbert–Polya, ouverte depuis un siècle.

L'apport de la construction est néanmoins de fournir :
- une **géométrie canonique** (cusp hyperbolique de paramètres `(R, n)`) ;
- un **opérateur canonique** (`H_n^{\text{cusp}}` avec correction de courbure) ;
- une **extension auto-adjointe distinguée** (antipériodique `θ = π` sélectionnée par symétrie spectrale et exclusion du mode marginal) ;
- une **densité asymptotique correcte** (Riemann–von Mangoldt) ;
- une **formule de trace explicite** dont la matrice de scattering est un produit de deux fonctions zêta apparentées.

### 7.2 Validation numérique

Une diagonalisation matricielle de `\tilde H` sur grille `\xi_j = j r / N`, `j = 0, …, N - 1`, avec wrap-around phasé pour antipériodicité, reproduit (4.5) à précision machine `\sim 10^{-13}` pour `N ∈ \{100, 500, 1000, 2000\}` et `θ \in \{0, π/2, π, 3π/2\}`. Sur les 79 premiers zéros de `ζ` (table LMFDB), une comparaison directe entre `\{γ_n\}` et le spectre régularisé par (6.7) atteint un accord à `|Δ| < 0.5` pour 75 zéros sur 79 sous la troncature canonique [2]. Une extension à 270 zéros est en cours.

### 7.3 Liens avec l'holographie céleste

Le générateur de dilatation `H = -i(x ∂_x + 1/2)` est, indépendamment, **le générateur de boost** diagonalisé par la transformée de Mellin qui définit la base d'amplitudes célestes de l'holographie 4D asymptotiquement plate [6]. La série principale unitaire des opérateurs primaires conformes sur la sphère céleste, `Δ ∈ 1 + i ℝ`, coïncide avec la droite critique `s = 1/2 + i ℝ` sous le shift `s = Δ - 1/2`. La construction présente du cusp Berry–Keating donne donc une **réalisation arithmétique** de la série principale céleste, dont les opérateurs primaires sont localisés aux ordonnées des zéros de `ζ`. Cette identification fait l'objet d'une note séparée [7].

### 7.4 Questions ouvertes

- Identification spectrale stricte : montrer que le spectre régularisé est exactement `\{γ_n\}`, pas seulement asymptotiquement.
- Forme explicite des `ζ_\pm` : caractérisation analytique complète (pôles, lignes de zéros, équation fonctionnelle).
- Cas multi-cusps : la construction se généralise-t-elle à une surface hyperbolique avec plusieurs cusps inéquivalents ?
- Phase de Maslov : dérivation rigoureuse du terme `-γ / (2π)` de la densité de Riemann–von Mangoldt.

---

## Références

[1] M. V. Berry, J. P. Keating. *H = xp and the Riemann zeros*. In *Supersymmetry and Trace Formulae* (Plenum, 1999) ; *SIAM Review* 41, 236 (1999).

[2] Y. Senez. *The spectral curve of a hyperbolic moduli construction* (monographie en préparation, 2026).

[3] M. Reed, B. Simon. *Methods of Modern Mathematical Physics*, vol. II : *Fourier Analysis, Self-Adjointness*. Academic Press (1975).

[4] E. C. Titchmarsh. *The Theory of the Riemann Zeta-Function*, 2nd ed., revised by D. R. Heath-Brown. Clarendon Press, Oxford (1986).

[5] H. Iwaniec. *Spectral methods of automorphic forms*, 2nd ed. AMS GSM 53 (2002).

[6] S. Pasterski. *Lectures on celestial amplitudes*. arXiv:2108.04801 (2021).

[7] Y. Senez. *L'opérateur de Berry–Keating comme générateur de dilatation céleste* (préprint compagnon, 2026).
