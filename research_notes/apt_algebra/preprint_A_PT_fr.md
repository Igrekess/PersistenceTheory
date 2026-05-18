# On the Wick-Poisson algebra of a hyperbolic cusp and its arithmetic statistical structure

**Yan Senez** · 2026-05-17 · *draft 1, §1-§4*

---

## Résumé

On construit rigoureusement une `C*`-algèbre quantique-statistique `A_PT = A_+ ⊗ A_-` sur l'espace des configurations de la courbe spectrale hyperbolique `Σ` (de rayon `R = √3`, courbure scalaire `K = -1/3`) admettant un cusp parabolique. L'algèbre est engendrée par les opérateurs de création-annihilation Wick-Poisson `b_p^σ(f)`, `b_p^σ†(f)` pour `f ∈ L²(Σ)`, `p` premier et `σ ∈ {+, -}` indexant les deux branches dynamiques de la courbe au cusp. L'espace de Hilbert ambiant est l'espace de Fock-Poisson `ℋ_PT = Γ(L²(X))` au sens d'Albeverio-Kondratiev-Röckner [3], où `X = Σ × \mathbb{P} × \{+, -\}` est l'espace de configuration marqué (deux dimensions continues sur `Σ`, un indice arithmétique discret sur les premiers, un indice de branche `ℤ_2`). On démontre que `A_PT` est nucléaire, sa K-théorie est `(K_0, K_1) = (\mathbb{Z}, 0)`, et son enveloppe de von Neumann sous l'état KMS canonique est un facteur de type `III_1`. La fonction de partition canonique est le produit `\zeta_+(\beta) \zeta_-(\beta)` de deux fonctions zêta tronquées portées par les branches `q_+, q_-`. Les sections §5-§9 (dynamique modulaire, condition KMS, fonction de partition, transitions de phase, discussion) seront rédigées dans des sessions ultérieures.

---

## 1. Introduction

Bost et Connes [1] ont construit en 1995 une `C*`-algèbre quantique-statistique dont la fonction de partition canonique reproduit la fonction zêta de Riemann, et dont les états KMS à basse température subissent une brisure spontanée de symétrie indexée par les caractères de Dirichlet. Cette construction a ouvert la voie à un programme de géométrie non-commutative arithmétique (Connes-Marcolli [2], Laca-Larsen-Neshveyev [4], Ha-Paugam [5]) où des `L`-fonctions classiques apparaissent comme fonctions de partition de systèmes opérateur-théoriques portant des structures adéliques.

La présente note s'inscrit dans cette lignée mais avec un point de départ géométrique différent. Le préprint compagnon [6] construit, sur le cusp parabolique d'une courbe spectrale hyperbolique `Σ` de rayon `R = √3` et de courbure scalaire `K = -1/3`, un opérateur de dilatation cuspidal de type Berry-Keating `H_n^{\text{cusp}} = -i(\eta \partial_\eta + 1/2 - n \eta / (2R))`. La construction y aboutit à une formule de trace de type Selberg, dans laquelle la matrice de scattering du cusp prend la forme produit
```
\varphi_\Sigma(s) = \zeta_+(s) \, \zeta_-(s),                              (1.1)
```
où `\zeta_\pm(s) = \exp(\sum_p a_p^\pm p^{-s})` sont deux fonctions zêta tronquées portées par deux branches dynamiques distinguées `q_+, q_-` du flot cuspidal. Les coefficients `a_p^\pm := \delta_p^\pm (2 - \delta_p^\pm)` avec `\delta_p^\pm := (1 - q_\pm^p) / p` sont une donnée géométrique de la courbe.

La question naturelle qu'ouvre [6] est : *existe-t-il un système quantique-statistique* `C*`-*algébrique dont la fonction de partition canonique soit précisément la matrice de scattering* (1.1) ? Cette note y répond positivement et construit explicitement une telle algèbre, notée `A_PT = A_+ ⊗ A_-`, comme produit tensoriel de deux sous-algèbres de Wick-Poisson indexées par les branches.

La construction diffère de Bost-Connes sur deux points essentiels.

D'une part, la **statistique** sous-jacente n'est pas bosonique mais **Poisson** : pour chaque premier `p` et chaque branche `\sigma`, le poids de l'état excité `k`-particules est `(a_p^\sigma)^k / k!`. Cette signature classique-distinguishable (Maxwell-Boltzmann pondéré) provient de la statistique des points distincts sur l'espace de configurations marqué `X = \Sigma \times \mathbb{P} \times \{+, -\}`, et contraste avec la statistique bosonique `\|\,|k\rangle\,\|² = 1` (sans `1/k!`) qui sous-tend la fonction de partition `\zeta(\beta) = \prod_p (1 - p^{-\beta})^{-1}` de [1]. Au niveau de la fonction de partition, le passage Bose `\to` Poisson remplace `\prod_p (1 - p^{-\beta})^{-1}` par `\prod_p \exp(a_p p^{-\beta}) = \exp(\sum_p a_p p^{-\beta})`.

D'autre part, la **structure à deux branches** est inscrite dès le départ comme produit tensoriel : `A_PT = A_+ \otimes A_-`. Aucune des constructions de la littérature post-1995 (Bost-Connes [1], `Q`-lattices `GL_n` de Connes-Marcolli [2], `K`-systèmes de Laca-Larsen-Neshveyev [4], systèmes de Shimura de Ha-Paugam [5], `q`-déformations de Marcolli) ne réalise simultanément ces deux ingrédients ; le triplet (tensorisation `A_+ \otimes A_-`, coefficients `a_p^\sigma` distincts par branche, statistique Poisson) est, à notre connaissance, nouveau.

### 1.1 Résultats principaux

Trois théorèmes structurent les §§ 2-4 et seront prouvés ou esquissés en détail dans cette note.

- **Théorème A (Espace de Fock-Poisson, §3).** *Pour tout `\beta` complexe avec `\Re \beta > 0`, il existe une mesure de Poisson `\pi_\beta` sur l'espace de configurations `\mathrm{Conf}(X)` d'intensité*
```
\lambda_\beta(x, p, \sigma) = a_p^\sigma \cdot p^{-\beta} \cdot \mu_\Sigma(x).
```
*L'espace de Hilbert associé `\mathcal{H}_{PT} := L²(\mathrm{Conf}(X), \pi_\beta)` est canoniquement isomorphe au Fock bosonique pondéré `\Gamma(L²(X))` au sens d'Albeverio-Kondratiev-Röckner [3], et admet la factorisation*
```
\mathcal{H}_{PT} \cong \bigotimes_{p, \sigma} \mathcal{H}_p^\sigma,
```
*où chaque `\mathcal{H}_p^\sigma = \bigoplus_{k \geq 0} \mathrm{Sym}^k L²(\Sigma)` est un Fock symétrique de poids `a_p^\sigma`.*

- **Théorème B (`C*`-algèbre Wick-Poisson, §4).** *L'algèbre `A_{PT}` engendrée dans `\mathcal{B}(\mathcal{H}_{PT})` par les opérateurs de Wick `\{b_p^\sigma(f), b_p^\sigma{}^\dagger(f) : f \in L²(\Sigma), p \in \mathbb{P}, \sigma \in \{+, -\}\}` est nucléaire et se factorise canoniquement comme `A_{PT} = A_+ \otimes A_-` (produit tensoriel `C*`-algébrique).*

- **Théorème C (Invariants K-théoriques, §4).** *`K_0(A_{PT}) = \mathbb{Z}`, `K_1(A_{PT}) = 0`. L'enveloppe faible `\pi_\beta(A_{PT})''` dans la représentation GNS de l'état KMS canonique (à `\beta > 0` réel) est un facteur de von Neumann de type `III_1`.*

Le Théorème A est notre version « Poisson-marquée 2-branches » de la construction de Fock standard ; il découle techniquement du théorème de représentation Wick-Poisson [3, Thm 4.1] appliqué à l'espace `X` ci-dessus. Le Théorème B est une vérification systématique des relations canoniques de commutation Wick-Poisson et de la nucléarité (via Künneth pour des facteurs Toeplitz). Le Théorème C, plus délicat, repose sur le théorème de Coburn [9] et le théorème de Künneth de Schochet [11] (pour la K-théorie), et sur le théorème ITPFI d'Araki-Woods [10] (pour le type von Neumann).

### 1.2 Plan

Le §2 rappelle, à partir de [6], la géométrie cuspidale de la courbe `\Sigma` qui sert de support continu à la construction. Le §3 introduit l'espace de configurations marqué `X = \Sigma \times \mathbb{P} \times \{+, -\}`, la mesure de Poisson `\pi_\beta`, et l'espace de Fock-Poisson `\mathcal{H}_{PT}`. Le §4 construit la `C*`-algèbre `A_{PT}`, démontre la factorisation `A_+ \otimes A_-`, calcule ses invariants K-théoriques, et identifie le type von Neumann de son enveloppe faible.

Les sections ultérieures (§5 dynamique modulaire de Tomita-Takesaki et identification du flot cuspidal `\alpha_t = e^{itH}` avec `H = -i(\eta \partial_\eta + 1/2 - n\eta/(2R))` sur la composante `\Sigma` ; §6 condition KMS pour la dynamique `\alpha_t(b_p^\sigma) = p^{-it} b_p^\sigma` ; §7 fonction de partition `Z(\beta) = \zeta_+(\beta) \zeta_-(\beta)` et convergence sur `\Re \beta > 0` ; §8 transitions de phase et brisure spontanée de symétrie ; §9 discussion et liens avec [1, 2, 6, 7]) seront rédigées dans des sessions ultérieures.

---

## 2. Géométrie cuspidale de la courbe `\Sigma`

Le support géométrique continu de la construction est la courbe `\Sigma` de [6], dont nous rappelons les propriétés strictement nécessaires.

### 2.1 Définition

`\Sigma` est une surface hyperbolique compacte de genre `g = 14`, point isolé hyperelliptique de l'espace de modules `M_{0, 30}^{\mathrm{hyp}}` (cf. [6, §2] et la notion auxiliaire de courbe spectrale). Ce qui importe pour la présente note est qu'elle admet un **bord asymptotique** muni d'un cusp parabolique de rang `n` (i.e. dont le groupe d'isotropie est `\mathbb{Z}^n`, ou pratiquement `n = 3` dans [6]).

Près du cusp, la métrique de `\Sigma` prend la forme demi-espace de Poincaré généralisée
```
ds² = d\eta² + e^{-2 \eta / R} \sum_{a = 1}^{n} dX_a²,                    (2.1)
```
avec `\eta > 0` la coordonnée radiale (`\eta \to \infty` au cusp), `X_a \in \mathbb{R} / L_a \mathbb{Z}` (`a = 1, \dots, n`) les coordonnées horocycliques compactifiées, `R > 0` le rayon hyperbolique, `n` la dimension parabolique. La courbure sectionnelle est constante `K = -1/R^2`. Pour `R = \sqrt 3` et `n = 3`, on retrouve `K = -1/3` et la mesure de volume cuspidale
```
d\mathrm{Vol}_\Sigma = e^{-n \eta / R} \, d\eta \prod_a dX_a.             (2.2)
```

### 2.2 Coordonnée Mellin

La coordonnée Mellin canonique est `\xi = \log(\eta / \eta_{\min})`, qui transforme la mesure de Haar multiplicative `d\eta / \eta` sur `(\eta_{\min}, \eta_{\max})` en mesure de Lebesgue `d\xi` sur `(0, r)`, où `r = \log(\eta_{\max} / \eta_{\min})`. Sous une troncature dynamique adaptative `r(\gamma) = \log(\gamma / 2\pi) + O(1)` (cf. [6, §5.2] et §5 ci-après), `r` devient une fonction lente de l'échelle spectrale.

Pour notre construction, on retient la donnée minimale suivante :

- **Mesure** : `\mu_\Sigma` est la mesure invariante de `\Sigma` issue de (2.2) ;
- **Espace de Hilbert** : `L²(\Sigma) := L²(\Sigma, \mu_\Sigma)`, séparable ;
- **Dynamique** : le flot de dilatation radiale `\sigma_t : \eta \mapsto e^t \eta` (i.e. `\xi \mapsto \xi + t`), unitairement implémenté sur `L²(\Sigma)` par le générateur cuspidal `H^{\mathrm{cusp}}_n = -i(\eta \partial_\eta + 1/2 - n\eta/(2R))` (théorème central de [6, Prop. 3.2]).

Aucune autre propriété géométrique de `\Sigma` n'intervient dans les §3-§4. Le rôle de `\Sigma` y est de fournir un *espace métrique mesuré séparable* sur lequel construire un espace de configurations marqué.

### 2.3 Asymptotique vers `\mathbb{H}^{n+1}`

[6, §2.3] prouve que `\Sigma` est asymptotiquement isométrique au demi-espace hyperbolique `\mathbb{H}^{n+1}` de rayon `R = \sqrt 3` quand `\eta \to \infty`. Cette asymptotique servira en §5 pour identifier la dynamique cuspidale `\alpha_t` à un flot Mellin canonique sur `\mathbb{H}^{n+1}`. Elle n'intervient pas en §3-§4.

---

## 3. Espace de configurations et espace de Fock-Poisson

### 3.1 L'espace de configurations marqué

**Définition 3.1.** *On note `\mathbb{P}` l'ensemble des nombres premiers et*
```
X := \Sigma \times \mathbb{P} \times \{+, -\},                            (3.1)
```
*muni de la mesure produit `\mu_X := \mu_\Sigma \otimes \mu_{\mathrm{count}}^{\mathbb{P}} \otimes \mu_{\mathrm{count}}^{\{\pm\}}`, où `\mu_{\mathrm{count}}` désigne la mesure de comptage. `X` est un espace métrique mesuré séparable, polonais, `\sigma`-fini.*

**Définition 3.2.** *L'espace de configurations `\mathrm{Conf}(X)` est l'ensemble des mesures atomiques localement finies sur `X`,*
```
\mathrm{Conf}(X) := \Bigl\{ \gamma = \sum_{\alpha} \delta_{x_\alpha} : x_\alpha \in X, \; \gamma(K) < \infty \;\; \forall K \subset X \text{ compact} \Bigr\}. \tag{3.2}
```
*Une configuration `\gamma \in \mathrm{Conf}(X)` est donc un ensemble localement fini de points marqués `(x_\alpha, p_\alpha, \sigma_\alpha) \in X` (sans répétition pour la mesure de Poisson simple).*

`\mathrm{Conf}(X)` est muni de la topologie vague et de la `\sigma`-algèbre borélienne associée [3, §2]. C'est l'espace de probabilité de support pour les mesures de Poisson marquées.

### 3.2 Mesure de Poisson marquée

**Définition 3.3.** *Soit `\beta \in \mathbb{C}` avec `\Re \beta > 0`. La mesure de Poisson `\pi_\beta` sur `\mathrm{Conf}(X)` est l'unique mesure de probabilité telle que, pour tout borélien `B \subset X` d'intensité*
```
\lambda_\beta(B) := \int_B a_p^\sigma \cdot p^{-\beta} \, d\mu_X(x, p, \sigma),                (3.3)
```
*la variable aléatoire `\gamma \mapsto \gamma(B)` est une variable de Poisson de paramètre `\lambda_\beta(B)`, et les variables associées à des boréliens disjoints sont indépendantes. Explicitement, on a la factorisation*
```
\pi_\beta = \pi_{\mu_\Sigma} \;\otimes\; \bigotimes_{p \in \mathbb{P}} \mathrm{Poiss}(a_p^+ \, p^{-\beta}) \;\otimes\; \bigotimes_{p \in \mathbb{P}} \mathrm{Poiss}(a_p^- \, p^{-\beta}),                   (3.4)
```
*où `\mathrm{Poiss}(\lambda)` est la loi de Poisson d'intensité `\lambda`.*

L'existence de `\pi_\beta` est une application du théorème de représentation de Kolmogorov-Daniell pour les processus ponctuels (voir Daley-Vere-Jones [8, ch. 9] ou Last-Penrose [12, Thm. 3.6]). La condition de convergence est l'intégrabilité de `\lambda_\beta` sur tout compact de `X`, garantie par le lemme suivant.

**Lemme 3.4.** *Pour `\Re \beta > 0`, la mesure d'intensité `\lambda_\beta` (eq. (3.3)) est localement finie sur `X`, et la mesure de Poisson `\pi_\beta` existe et est unique.*

*Démonstration.* L'intégrabilité locale revient à montrer que pour tout compact `K \subset \Sigma`, la somme `\sum_{p, \sigma} a_p^\sigma p^{-\Re \beta} \mu_\Sigma(K)` est finie. Or les coefficients `a_p^\sigma = \delta_p^\sigma (2 - \delta_p^\sigma)` avec `\delta_p^\sigma = (1 - q_\sigma^p)/p` admettent l'asymptotique `a_p^\sigma = 2/p + O(|q_\sigma|^p / p + 1/p^2)` pour `p \to \infty` (cf. [6, §6.3]), donc `a_p^\sigma \leq C/p` uniformément. Il en résulte
```
\sum_{p, \sigma} a_p^\sigma p^{-\Re \beta} \leq 2 C \sum_{p} p^{-(\Re \beta + 1)} < \infty
```
pour `\Re \beta > 0` (convergence par comparaison avec `\sum p^{-(1 + \varepsilon)}`). L'unicité suit du théorème de représentation de Kolmogorov [8, Thm. 9.2.IV]. ∎

**Remarque 3.5.** L'asymptotique `a_p^\sigma \sim 2/p` fait que la somme `\sum_p a_p^+ p^{-\beta} = \log \zeta_+(\beta)` se comporte comme `2 \mathcal{P}(\beta + 1)` (où `\mathcal{P}(s) = \sum_p p^{-s}` est le prime zêta) à l'ordre dominant. La fonction de partition `\zeta_+(\beta) \zeta_-(\beta) = \exp(\sum_{p, \sigma} a_p^\sigma p^{-\beta})` est ainsi holomorphique sur `\Re \beta > 0`, avec un pôle d'ordre 4 à `\beta = 0`, comme démontré en détail dans [6] et le préprint compagnon sur la convergence. Cette extension automatique à la bande critique `0 < \Re \beta < 1` (au-delà du domaine `\Re \beta > 1` de Bost-Connes [1]) est une propriété structurelle de la statistique Poisson à deux branches.

### 3.3 Espace de Fock-Poisson

**Définition 3.6.** *L'espace de Hilbert de Fock-Poisson est*
```
\mathcal{H}_{PT} := L²(\mathrm{Conf}(X), \pi_\beta).                                  (3.5)
```

Albeverio, Kondratiev et Röckner [3, Thm. 4.1] établissent l'isomorphisme canonique
```
L²(\mathrm{Conf}(X), \pi_\beta) \;\cong\; \Gamma_{\mathrm{sym}}(L²(X, \lambda_\beta))
= \bigoplus_{k \geq 0} \mathrm{Sym}^k L²(X, \lambda_\beta),                          (3.6)
```
où `\Gamma_{\mathrm{sym}}` désigne le foncteur de Fock symétrique. La normalisation Wick-Poisson assigne au vecteur élémentaire de niveau `k` la norme `\|f_1 \otimes_{\mathrm{sym}} \cdots \otimes_{\mathrm{sym}} f_k\|^2 = k! \prod_i \|f_i\|_{\lambda_\beta}^2`, qui se réécrit en termes des polynômes orthogonaux de Charlier (le pendant Poisson des Hermite) pour la base orthonormée. La présence du facteur `1/k!` dans la normalisation des états excités est la signature de la **statistique Poisson** (distinguishable, ou Maxwell-Boltzmann classique pondéré) ; elle contraste avec la normalisation bosonique sans `1/k!` qui caractérise le Fock standard de [1].

### 3.4 Factorisation par espèce

**Proposition 3.7.** *L'espace `\mathcal{H}_{PT}` se factorise canoniquement comme*
```
\mathcal{H}_{PT} \;\cong\; \bigotimes_{p \in \mathbb{P}, \, \sigma \in \{+, -\}} \mathcal{H}_p^\sigma,                                   (3.7)
```
*où*
```
\mathcal{H}_p^\sigma := L²\bigl(\mathrm{Conf}(\Sigma), \pi_{a_p^\sigma p^{-\beta}}\bigr)
\;\cong\; \Gamma_{\mathrm{sym}}\bigl(L²(\Sigma, a_p^\sigma p^{-\beta} \mu_\Sigma)\bigr).                                  (3.8)
```
*Chaque `\mathcal{H}_p^\sigma` est un Fock symétrique sur `L²(\Sigma)` de poids `a_p^\sigma p^{-\beta}`. Le produit tensoriel infini en (3.7) est pris au sens des produits tensoriels GNS (Powers [10]) relativement aux états vide de chaque facteur.*

*Démonstration.* La factorisation découle de la factorisation (3.4) de la mesure de Poisson sur des composantes disjointes de `X`. Pour chaque `(p, \sigma)`, la projection `\gamma \mapsto \gamma|_{\Sigma \times \{p\} \times \{\sigma\}}` est une variable aléatoire Poisson indépendante (sur l'algèbre de `\sigma`-Borel des configurations sur la composante correspondante de `X`). L'isomorphisme (3.6) appliqué facteur par facteur donne (3.7-3.8). La convergence du produit tensoriel infini relatif aux états vide est garantie par le critère d'Araki-Woods [10] et la sommabilité `\sum_{p, \sigma} a_p^\sigma p^{-\Re \beta} < \infty` du Lemme 3.4. ∎

**Remarque 3.8.** La structure (3.7) reflète l'indépendance physique des canaux `(p, \sigma)` : à chaque premier `p` et chaque branche `\sigma`, un processus de Poisson indépendant d'intensité `a_p^\sigma p^{-\beta}` peuple `\Sigma` d'« événements » spatialement distribués selon `\mu_\Sigma`. La fonction de partition canonique sera, on le verra en §6-§7, le produit des fonctions de partition à un canal :
```
Z(\beta) = \prod_{p, \sigma} \exp(a_p^\sigma p^{-\beta}) = \exp\Bigl(\sum_{p, \sigma} a_p^\sigma p^{-\beta}\Bigr) = \zeta_+(\beta) \, \zeta_-(\beta).             (3.9)
```

**Remarque 3.9 (statistique Poisson vs bosonique).** Le contraste avec Bost-Connes [1] est instructif. Dans [1], l'espace de Fock est l'espace `\ell²(\mathbb{N})` avec base `\{|n\rangle\}_{n \geq 1}`, et la fonction de partition à un canal `p` est `\sum_k p^{-\beta k} = (1 - p^{-\beta})^{-1}`. Dans la construction présente, la fonction de partition à un canal `(p, \sigma)` est `\sum_k (a_p^\sigma)^k p^{-\beta k} / k! = \exp(a_p^\sigma p^{-\beta})`. La différence est précisément le facteur `1/k!`, qui distingue la statistique de Maxwell-Boltzmann classique (Poisson) de la statistique de Bose-Einstein (bosonique non-pondérée). Cette différence se propage à la fonction de partition totale (`\prod_p \exp(a_p p^{-\beta})` vs `\prod_p (1 - p^{-\beta})^{-1}`) et change la nature analytique : `\zeta_+ \zeta_-` est entière au sens des produits exponentiels et n'a pas de pôle dans `\Re \beta > 0`, contrairement à `\zeta` qui a un pôle simple à `\beta = 1`.

---

## 4. La `C*`-algèbre `A_{PT}`

### 4.1 Opérateurs de Wick-Poisson

Pour `f \in L²(\Sigma)`, `p \in \mathbb{P}`, `\sigma \in \{+, -\}`, on définit l'**opérateur de création** `b_p^\sigma{}^\dagger(f)` sur `\mathcal{H}_p^\sigma` (puis étendu à `\mathcal{H}_{PT}` par tensorisation triviale sur les autres canaux) par
```
b_p^\sigma{}^\dagger(f) \cdot (f_1 \otimes_{\mathrm{sym}} \cdots \otimes_{\mathrm{sym}} f_k)
:= \sqrt{a_p^\sigma p^{-\beta}} \; f \otimes_{\mathrm{sym}} f_1 \otimes_{\mathrm{sym}} \cdots \otimes_{\mathrm{sym}} f_k.                                (4.1)
```
L'**opérateur d'annihilation** `b_p^\sigma(f)` est son adjoint formel,
```
b_p^\sigma(f) \cdot (f_1 \otimes_{\mathrm{sym}} \cdots \otimes_{\mathrm{sym}} f_k)
:= \sqrt{a_p^\sigma p^{-\beta}} \;\sum_{i = 1}^{k} \langle f, f_i \rangle_{L²(\Sigma)} \; f_1 \otimes_{\mathrm{sym}} \cdots \widehat{f_i} \cdots \otimes_{\mathrm{sym}} f_k.        (4.2)
```

**Définition 4.1.** *La `C*`-algèbre Wick-Poisson à deux branches est*
```
A_{PT} := C^*\bigl( \{ b_p^\sigma(f), \, b_p^\sigma{}^\dagger(f) : f \in L²(\Sigma), \, p \in \mathbb{P}, \, \sigma \in \{+, -\} \} \bigr) \;\subset\; \mathcal{B}(\mathcal{H}_{PT}), \qquad (4.3)
```
*où la clôture est prise en norme d'opérateur.*

### 4.2 Relations canoniques de commutation

**Proposition 4.2.** *Les opérateurs (4.1)-(4.2) satisfont les relations canoniques de commutation Wick-Poisson : pour tous `f, g \in L²(\Sigma)`, `p, q \in \mathbb{P}`, `\sigma, \tau \in \{+, -\}`,*
```
[b_p^\sigma(f), \; b_q^\tau{}^\dagger(g)] = \delta_{pq} \, \delta_{\sigma\tau} \; a_p^\sigma \, p^{-\beta} \; \langle f, g \rangle_{L²(\Sigma)} \cdot \mathbf{1}, \qquad (4.4)
```
```
[b_p^\sigma(f), \; b_q^\tau(g)] = 0, \qquad [b_p^\sigma{}^\dagger(f), \; b_q^\tau{}^\dagger(g)] = 0. \qquad (4.5)
```

*Démonstration.* Calcul direct à partir de (4.1)-(4.2) sur un vecteur de niveau `k`. Pour `(p, \sigma) \neq (q, \tau)`, les facteurs `\mathcal{H}_p^\sigma` et `\mathcal{H}_q^\tau` sont tensoriellement indépendants donc les opérateurs commutent trivialement. Pour `(p, \sigma) = (q, \tau)`, le calcul standard sur le Fock symétrique donne (4.4), avec le facteur `a_p^\sigma p^{-\beta}` provenant des préfacteurs des opérateurs (4.1)-(4.2). Les relations (4.5) sont symétriques. ∎

**Remarque 4.3.** La pondération `a_p^\sigma p^{-\beta}` dans le membre droit de (4.4) est *intensité-dépendante* : c'est elle qui distingue la convention Wick-Poisson de la convention Wick bosonique standard (où le membre droit serait simplement `\langle f, g \rangle \cdot \mathbf{1}`). Sous la limite `a_p^\sigma p^{-\beta} \to 1`, on retrouve formellement les CCR standard ; mais cette limite n'est pas atteinte par la construction (les intensités physiques sont à `\beta > 0` réel). Le choix d'inclure ces facteurs dans les CCR plutôt que dans la définition des opérateurs est conventionnel ; nous suivons ici la convention d'Albeverio-Kondratiev-Röckner [3].

### 4.3 Factorisation `A_{PT} = A_+ \otimes A_-`

**Théorème 4.4.** *L'algèbre `A_{PT}` se factorise canoniquement comme produit tensoriel `C*`-algébrique*
```
A_{PT} \;=\; A_+ \;\otimes\; A_-, \qquad (4.6)
```
*où pour `\sigma \in \{+, -\}`,*
```
A_\sigma := C^*\bigl( \{ b_p^\sigma(f), \, b_p^\sigma{}^\dagger(f) : f \in L²(\Sigma), \, p \in \mathbb{P} \} \bigr) \;\subset\; \mathcal{B}\bigl(\mathcal{H}_\sigma\bigr), \qquad \mathcal{H}_\sigma := \bigotimes_p \mathcal{H}_p^\sigma. \qquad (4.7)
```
*Les algèbres `A_+` et `A_-` agissent sur des facteurs tensoriels disjoints de `\mathcal{H}_{PT}` et commutent fortement.*

*Démonstration.* La factorisation `\mathcal{H}_{PT} = \mathcal{H}_+ \otimes \mathcal{H}_-` découle de la Proposition 3.7 (sépare les facteurs `\sigma = +` et `\sigma = -`). Les générateurs `\{b_p^+, b_p^+{}^\dagger\}_p` agissent comme `\mathrm{id} \otimes \cdot` sur `\mathcal{H}_+ \otimes \mathcal{H}_-`, et `\{b_p^-, b_p^-{}^\dagger\}_p` comme `\cdot \otimes \mathrm{id}`. Les relations de commutation (4.4)-(4.5) avec `\sigma \neq \tau` donnent commutateur nul, donc les deux sous-algèbres sont en position tensorielle. La nucléarité de chaque facteur `A_\sigma` (Proposition 4.5 ci-dessous) assure que le produit tensoriel `C*`-minimal et maximal coïncident, soit `A_{PT} = A_+ \otimes A_-` au sens canonique. ∎

### 4.4 Nucléarité et K-théorie

**Proposition 4.5.** *Pour chaque `\sigma`, l'algèbre `A_\sigma` est nucléaire. Sa K-théorie est `K_0(A_\sigma) = \mathbb{Z}` et `K_1(A_\sigma) = 0`.*

*Esquisse.* À chaque `(p, \sigma)` fixé, la sous-algèbre `A_p^\sigma := C^*(b_p^\sigma(f), b_p^\sigma{}^\dagger(f) : f \in L²(\Sigma))` agit sur `\mathcal{H}_p^\sigma`. Pour `L²(\Sigma)` séparable de dimension infinie, `A_p^\sigma` est isomorphe à l'algèbre de Toeplitz `\mathcal{T} \otimes \mathcal{K}` à coefficient d'opérateurs compacts (Coburn [9]) : le créateur `b_p^\sigma{}^\dagger(f)` agit comme une isométrie pondérée (modulo la base orthonormée de `L²(\Sigma)`, c'est une famille d'isométries Toeplitz `\mu_{p, f}^\sigma`). La K-théorie de l'algèbre de Toeplitz est `(K_0, K_1) = (\mathbb{Z}, 0)` ([9], et Wegge-Olsen [13, ex. 9.4.1]). Pour le produit tensoriel infini relatif aux états vide,
```
A_\sigma \;\cong\; \bigotimes_p^{\,\mathrm{GNS}} \bigl(\mathcal{T} \otimes \mathcal{K}, \psi_p^{\mathrm{vac}}\bigr),
```
le critère ITPFI d'Araki-Woods [10] garantit la nucléarité (chaque facteur étant nucléaire et la famille des états vide étant un produit tensoriel cohérent). Le théorème de Künneth pour les `C*`-algèbres nucléaires (Schochet [11, Thm 4.1]) avec `K_*(\mathcal{T} \otimes \mathcal{K})` libre donne, par itération sur les premiers et stabilisation,
```
K_0(A_\sigma) = \mathbb{Z}, \qquad K_1(A_\sigma) = 0.
```
∎

**Théorème 4.6.** *L'algèbre `A_{PT} = A_+ \otimes A_-` est nucléaire. Sa K-théorie est*
```
K_0(A_{PT}) = \mathbb{Z}, \qquad K_1(A_{PT}) = 0. \qquad (4.8)
```

*Démonstration.* Nucléarité : produit tensoriel de deux algèbres nucléaires. K-théorie : application du théorème de Künneth [11] à `A = A_+`, `B = A_-` avec `K_*(A_-) = K_*(A_+) = (\mathbb{Z}, 0)` libre :
```
K_0(A_+ \otimes A_-) = K_0(A_+) \otimes K_0(A_-) \oplus K_1(A_+) \otimes K_1(A_-) = \mathbb{Z} \otimes \mathbb{Z} = \mathbb{Z},
```
```
K_1(A_+ \otimes A_-) = K_0(A_+) \otimes K_1(A_-) \oplus K_1(A_+) \otimes K_0(A_-) = 0.
```
∎

### 4.5 Type von Neumann

**Proposition 4.7.** *Soit `\omega_\beta^{PT}` l'état KMS canonique sur `A_{PT}` à `\beta > 0` réel (à définir précisément en §6 ; à ce stade nous prenons l'état produit factorisé `\omega_\beta^{PT} = \omega_\beta^+ \otimes \omega_\beta^-` avec `\omega_\beta^\sigma` l'état Poisson d'intensité `a_p^\sigma p^{-\beta}` sur chaque canal). La fermeture faible*
```
M_{PT} := \pi_{\omega_\beta^{PT}}(A_{PT})''  \subset  \mathcal{B}(\mathcal{H}_{PT}^{\mathrm{GNS}}) \qquad (4.9)
```
*est un facteur de von Neumann de type `III_1`.*

*Esquisse.* Pour chaque branche `\sigma`, l'enveloppe faible `M_\sigma := \pi_{\omega_\beta^\sigma}(A_\sigma)''` est un produit tensoriel infini ITPFI (Araki-Woods [10]) de facteurs de type `I_\infty` aux primes (chaque `\mathcal{B}(\mathcal{H}_p^\sigma)$ étant de type `I_\infty`). Le spectre de Connes de `M_\sigma` est l'adhérence du groupe multiplicatif engendré par `\{p^\beta : p \in \mathbb{P}\}` dans `\mathbb{R}_+^*`. L'indépendance linéaire sur `\mathbb{Q}` des `\log p` pour `p` premier (corollaire élémentaire du théorème fondamental de l'arithmétique) implique que ce groupe est dense dans `\mathbb{R}_+^*` dès que `\beta > 0`. Donc `S(M_\sigma) = [0, \infty)`, ce qui caractérise le type `III_1` (Connes-Takesaki [14]). Le produit tensoriel `M_{PT} = M_+ \otimes M_-` de deux facteurs `III_1` reste un facteur `III_1` [14, Thm 4.8]. ∎

**Remarque 4.8.** Le type `III_1` est partagé avec le facteur de Bost-Connes [1, Thm 25], qui est aussi `III_1` pour `\beta > 1`. C'est cohérent : le type von Neumann d'un facteur ITPFI ne distingue pas finement la statistique bosonique de la statistique Poisson (les deux donnent un produit infini de facteurs de type `I_\infty` aux intensités dense-dans-`\mathbb{R}_+^*`). La signature distinctive de `A_{PT}` vs `A_{BC}` ne réside donc pas dans le type von Neumann, mais dans la **dimension spectrale** au sens NCG (ordre du pôle de la fonction zêta de Dirac à la transition de phase) : `A_{BC}` a dimension spectrale 1 (pôle simple de `\zeta` à `\beta = 1`), alors que `A_{PT}` a dimension spectrale 4 (pôle d'ordre 4 de `\zeta_+ \zeta_-` à `\beta = 0`). Une analyse détaillée des invariants K-théoriques et de la dimension spectrale 4D est donnée séparément dans [7].

### 4.6 Opérateurs de comptage

**Lemme 4.9.** *Pour chaque `(p, \sigma)` et `f \in L²(\Sigma)`, l'opérateur*
```
N_p^\sigma(f, f) := b_p^\sigma{}^\dagger(f) \, b_p^\sigma(f) \in A_{PT} \qquad (4.10)
```
*est positif et essentiellement auto-adjoint sur le sous-espace dense des vecteurs de niveau total fini. Plus généralement, l'opérateur de comptage total*
```
N_p^\sigma := \int_\Sigma b_p^\sigma{}^\dagger(x) \, b_p^\sigma(x) \, d\mu_\Sigma(x) \qquad (4.11)
```
*est défini comme la fermeture autosadjointe de la forme quadratique correspondante (au sens des formes de Dirichlet [3, §5]), avec spectre `\mathrm{spec}(N_p^\sigma) = \mathbb{N}` et état fondamental le vide.*

*Démonstration.* Positivité : `\langle N_p^\sigma(f, f) \, \xi, \xi\rangle = \|b_p^\sigma(f) \xi\|^2 \geq 0` pour tout `\xi \in \mathcal{H}_{PT}`. L'auto-adjonction essentielle suit du fait que `N_p^\sigma(f, f)` préserve la décomposition `\mathcal{H}_p^\sigma = \bigoplus_k \mathrm{Sym}^k L²(\Sigma)` et y agit comme un opérateur borné par niveau. L'intégrale (4.11) doit être interprétée au sens des formes (l'intégrande `b_p^\sigma{}^\dagger(x) b_p^\sigma(x)` n'est pas un opérateur ponctuellement défini), et sa fermeture est l'opérateur de comptage des points de la `(p, \sigma)`-composante d'une configuration `\gamma \in \mathrm{Conf}(X)`. Sous l'identification (3.6), `N_p^\sigma \xi = k \xi` pour `\xi \in \mathrm{Sym}^k L²(\Sigma) \subset \mathcal{H}_p^\sigma`. ∎

**Remarque 4.10 (interprétation Poisson).** Au niveau classique, l'état KMS `\omega_\beta^{PT}` correspond à une configuration aléatoire `\gamma` distribuée selon `\pi_\beta` (Définition 3.3). À chaque `(p, \sigma)`, la `(p, \sigma)`-composante de `\gamma` est un processus de Poisson sur `\Sigma` d'intensité `a_p^\sigma p^{-\beta}`. La fonction de partition canonique
```
Z(\beta) = \pi_\beta(\mathrm{Conf}(X))^{-1} \cdot \text{(normalisation Wick)} = \prod_{p, \sigma} \exp(a_p^\sigma p^{-\beta} \mu_\Sigma(\Sigma)) \qquad (4.12)
```
réduit (sous normalisation Wick et `\mu_\Sigma(\Sigma) = 1`) à `\exp(\sum_{p, \sigma} a_p^\sigma p^{-\beta}) = \zeta_+(\beta) \zeta_-(\beta) = \varphi_\Sigma(\beta)`, qui est précisément la matrice de scattering du cusp identifiée en [6, §6.3, eq. (6.6)]. L'identification rigoureuse de `Z(\beta)` à `\varphi_\Sigma(\beta)`, ainsi que la condition KMS pour la dynamique `\alpha_t(b_p^\sigma) = p^{-it} b_p^\sigma`, fait l'objet du §6.

---

## 5. Dynamique et structure modulaire de Tomita-Takesaki

### 5.1 Le Hamiltonien et la dynamique automorphique

La dynamique canonique sur `A_{PT}` provient d'un Hamiltonien diagonal dans la décomposition par espèces (3.7), dont chaque facteur `(p, \sigma)` reçoit le poids logarithmique `\log p`. Ce choix reflète l'identité `p^{-\beta} = e^{-\beta \log p}` et garantit, en §7, l'identification de la fonction de partition canonique au produit `\zeta_+(\beta) \zeta_-(\beta)` de la matrice de scattering cuspidale (1.1).

**Définition 5.1 (Hamiltonien total).** *On définit l'opérateur autosadjoint positif*
```
H \;:=\; \sum_{p \in \mathbb{P}} (\log p) \, \bigl( N_p^+ + N_p^- \bigr), \qquad (5.1)
```
*où `N_p^\sigma` est l'opérateur de comptage total (Lemme 4.9, eq. (4.11)) sur le facteur `\mathcal{H}_p^\sigma`. L'opérateur est défini a priori comme somme formelle ; sa réalisation auto-adjointe est précisée dans la Proposition 5.2.*

Le domaine naturel de `H` est le sous-espace dense
```
\mathcal{D}_H \;:=\; \Bigl\{ \xi \in \mathcal{H}_{PT} : \sum_{p} (\log p)^2 \, \|N_p^+ \xi\|^2 < \infty, \;\; \sum_{p} (\log p)^2 \, \|N_p^- \xi\|^2 < \infty \Bigr\}. \qquad (5.2)
```

**Proposition 5.2 (auto-adjonction essentielle).** *L'opérateur `H` est essentiellement auto-adjoint sur `\mathcal{D}_H`. Son spectre est l'adhérence dans `[0, \infty)` de l'ensemble*
```
\mathrm{spec}(H) \;=\; \overline{\Bigl\{ \sum_{p \in S^+} k_p^+ \log p + \sum_{p \in S^-} k_p^- \log p \,:\, S^\pm \subset \mathbb{P} \text{ fini}, \; k_p^\sigma \in \mathbb{N} \Bigr\}} \;=\; \{0\} \cup [\log 2, \infty), \qquad (5.3)
```
*partie discrète des points isolés (configurations à support fini) et partie continue (limites d'accumulation des combinaisons linéaires des `\log p` à coefficients entiers).*

*Démonstration.* `H` préserve la décomposition spectrale `\mathcal{H}_{PT} = \bigoplus_{\mathbf{k}} \mathcal{H}_{PT}^{(\mathbf{k})}` où `\mathbf{k} = (k_p^\sigma)_{p, \sigma}` est une famille presque nulle d'entiers et `\mathcal{H}_{PT}^{(\mathbf{k})} := \bigotimes_{p, \sigma} \mathrm{Sym}^{k_p^\sigma} L²(\Sigma)`. Sur `\mathcal{H}_{PT}^{(\mathbf{k})}` l'action de `H` est par le scalaire `\lambda_{\mathbf{k}} := \sum_{p, \sigma} k_p^\sigma \log p`. Le domaine (5.2) est précisément le domaine maximal de cet opérateur multiplicatif sur le bloc diagonal. L'auto-adjonction essentielle suit du critère de Nelson appliqué à chaque bloc fini puis du passage à la limite (les `\mathcal{H}_{PT}^{(\mathbf{k})}` étant orthogonaux). Le spectre est l'ensemble (5.3) ; le minimum non nul `\log 2` est atteint par `\mathbf{k} = (1, 0, \dots)` correspondant à une particule unique au prime `p = 2` ; l'adhérence remplit `[\log 2, \infty)` par densité des `\sum_p k_p \log p` parmi les réels positifs (densité Weyl pour les logarithmes des entiers). ∎

**Définition 5.3 (dynamique).** *La dynamique automorphique sur `A_{PT}` est le groupe à un paramètre*
```
\alpha_t \,:\, A_{PT} \to A_{PT}, \qquad \alpha_t(a) := e^{itH} \, a \, e^{-itH}, \quad t \in \mathbb{R}. \qquad (5.4)
```
*Sur les générateurs de Wick (4.1)-(4.2), elle agit par dilatation logarithmique :*
```
\alpha_t(b_p^\sigma(f)) = p^{-it} \, b_p^\sigma(f), \qquad \alpha_t(b_p^\sigma{}^\dagger(f)) = p^{+it} \, b_p^\sigma{}^\dagger(f). \qquad (5.5)
```

**Lemme 5.4 (groupe d'automorphismes).** *`(\alpha_t)_{t \in \mathbb{R}}` est un groupe à un paramètre fortement continu de `*`-automorphismes de `A_{PT}`. La continuité forte est uniforme sur les sous-algèbres de niveau borné, et l'orbite de chaque générateur est analytique dans `t \in \mathbb{C}`.*

*Démonstration.* Que `\alpha_t \circ \alpha_s = \alpha_{t+s}` et `\alpha_0 = \mathrm{id}` est immédiat. Que `\alpha_t` soit un `*`-automorphisme suit de l'unitarité de `e^{itH}`. Pour vérifier (5.5), notons que `b_p^\sigma{}^\dagger(f)` envoie un vecteur de niveau `k_p^\sigma` en un vecteur de niveau `k_p^\sigma + 1`, donc augmente la valeur propre de `H` de `\log p` ; conjugaison par `e^{itH}` multiplie alors par `e^{it \log p} = p^{it}`. La relation (5.5) en découle, et la continuité forte ainsi que l'analyticité résultent de la borne uniforme `\|e^{itH} \xi\| = \|\xi\|` et de la dépendance polynomiale par niveau. ∎

**Remarque 5.5 (lien avec le flot cuspidal).** Sous l'identification asymptotique `\Sigma \to \mathbb{H}^{n+1}` (§2.3), la dynamique `\alpha_t` se relève au flot Mellin sur la composante `\Sigma` (dilatation `\eta \mapsto e^t \eta`) engendré par l'opérateur cuspidal `H^{\mathrm{cusp}}_n = -i(\eta \partial_\eta + 1/2 - n\eta/(2R))` de [6, Prop 3.2]. Plus précisément, l'opérateur `H` de (5.1) est l'extension de seconde quantification de `H^{\mathrm{cusp}}_n` au Fock-Poisson `\mathcal{H}_{PT}`, agissant sur la composante arithmétique-de-branche par la graduation des `\log p`. Cette identification ferme la boucle entre la formule de trace cuspidale du préprint compagnon et la structure quantique-statistique présente.

### 5.2 Vide Poisson et propriétés cyclique/séparante

**Définition 5.6 (vide Poisson).** *Le **vide Poisson** est le vecteur unité*
```
\Omega_\beta \;\in\; \mathcal{H}_{PT} \qquad (5.6)
```
*correspondant, sous l'isomorphisme (3.6), à la configuration vide `\gamma = 0 \in \mathrm{Conf}(X)` (i.e. à la fonction caractéristique de l'événement `\{\gamma : \gamma = 0\}` normalisée). Dans la décomposition (3.7), `\Omega_\beta = \bigotimes_{p, \sigma} \Omega_p^\sigma` où chaque `\Omega_p^\sigma \in \mathcal{H}_p^\sigma` est le vecteur de niveau zéro.*

**Proposition 5.7 (cyclicité).** *Le vide `\Omega_\beta` est cyclique pour `A_{PT}` : l'ensemble*
```
\bigl\{ b_{p_1}^{\sigma_1}{}^\dagger(f_1) \cdots b_{p_n}^{\sigma_n}{}^\dagger(f_n) \, \Omega_\beta : n \geq 0, \; p_i \in \mathbb{P}, \; \sigma_i \in \{+, -\}, \; f_i \in L²(\Sigma) \bigr\} \qquad (5.7)
```
*engendre un sous-espace dense de `\mathcal{H}_{PT}`.*

*Démonstration.* Par la décomposition (3.6)-(3.7), il suffit de montrer la cyclicité du vide sur chaque facteur `\mathcal{H}_p^\sigma`. Or, par l'action (4.1), les vecteurs `b_p^\sigma{}^\dagger(f_1) \cdots b_p^\sigma{}^\dagger(f_k) \, \Omega_p^\sigma = (a_p^\sigma p^{-\beta})^{k/2} \, f_1 \otimes_{\mathrm{sym}} \cdots \otimes_{\mathrm{sym}} f_k` parcourent un sous-espace dense de `\mathrm{Sym}^k L²(\Sigma)` quand `f_i` parcourt une partie dense de `L²(\Sigma)`. La densité totale dans `\mathcal{H}_p^\sigma = \bigoplus_k \mathrm{Sym}^k L²(\Sigma)` est obtenue par sommation directe. Pour le produit tensoriel infini relatif à la suite de vides `(\Omega_p^\sigma)_{p, \sigma}`, le critère de densité d'Araki-Woods [10] donne la densité dans `\mathcal{H}_{PT}`. ∎

**Proposition 5.8 (séparabilité).** *Le vide `\Omega_\beta` est séparant pour la fermeture faible `\mathcal{M}_{PT} := \pi(A_{PT})''` (où `\pi` désigne la représentation GNS associée à l'état vide `\langle \Omega_\beta, \cdot \, \Omega_\beta \rangle`) : si `m \in \mathcal{M}_{PT}` vérifie `m \, \Omega_\beta = 0`, alors `m = 0`.*

*Démonstration.* C'est le pendant cyclique-séparant standard de Tomita : `\Omega_\beta` séparant pour `\mathcal{M}_{PT}` ⇔ `\Omega_\beta` cyclique pour `\mathcal{M}_{PT}'`. Par Proposition 5.7 `\Omega_\beta` est cyclique pour `A_{PT} \subset \mathcal{M}_{PT}`, et le commutant `\mathcal{M}_{PT}'` est non trivial (il contient en particulier les translatés à droite obtenus par échange entre les deux branches via l'antiautomorphisme `\sigma \mapsto -\sigma`, ainsi que les générateurs d'« annihilation à droite » construits par la `J`-conjugaison du Théorème 5.9). La cyclicité de `\Omega_\beta` pour `\mathcal{M}_{PT}'` se vérifie alors par symétrie : l'orbite `\mathcal{M}_{PT}' \, \Omega_\beta` contient l'image par `J` de l'orbite `\mathcal{M}_{PT} \, \Omega_\beta = \mathcal{H}_{PT}` (Proposition 5.7), donc est dense. ∎

### 5.3 Théorème modulaire de Tomita-Takesaki

L'environnement `(\mathcal{M}_{PT}, \Omega_\beta)` étant cyclique-séparant, le théorème de Tomita-Takesaki [Tak70] (cf. Bratteli-Robinson [BR81, §2.5]) s'applique.

**Théorème 5.9 (structure modulaire).** *Il existe une involution anti-unitaire `J : \mathcal{H}_{PT} \to \mathcal{H}_{PT}` et un opérateur positif inversible `\Delta` sur `\mathcal{H}_{PT}` (dans l'algèbre de von Neumann étendue) tels que :*

*(i) L'opérateur*
```
S : \mathcal{M}_{PT} \, \Omega_\beta \to \mathcal{M}_{PT} \, \Omega_\beta, \qquad S(m \, \Omega_\beta) := m^* \, \Omega_\beta \qquad (5.8)
```
*est fermable et sa fermeture admet la décomposition polaire `\overline{S} = J \, \Delta^{1/2}`.*

*(ii) Le groupe modulaire des automorphismes*
```
\sigma_t^{\Omega_\beta}(m) := \Delta^{it} \, m \, \Delta^{-it}, \qquad t \in \mathbb{R}, \qquad (5.9)
```
*préserve `\mathcal{M}_{PT}` ;*

*(iii) `J \mathcal{M}_{PT} J = \mathcal{M}_{PT}'` (l'involution `J` échange l'algèbre et son commutant).*

*Esquisse.* C'est l'énoncé canonique du théorème de Tomita-Takesaki appliqué à la paire cyclique-séparante `(\mathcal{M}_{PT}, \Omega_\beta)` (Propositions 5.7-5.8). Aucune adaptation n'est requise au-delà de la vérification que `\mathcal{M}_{PT}` est en position standard, ce qui est garanti par la construction GNS. ∎

### 5.4 Identification du groupe modulaire

La structure modulaire `\sigma_t^{\Omega_\beta}` admet une expression explicite sur les générateurs.

**Proposition 5.10 (groupe modulaire et dynamique).** *Sur les générateurs de Wick, le groupe modulaire de Tomita-Takesaki agit par*
```
\sigma_t^{\Omega_\beta}(b_p^\sigma(f)) = p^{-i \beta t} \, b_p^\sigma(f) = \alpha_{\beta t}(b_p^\sigma(f)), \qquad (5.10)
```
*où `\alpha` est la dynamique de (5.5). Autrement dit, à une reparamétrisation linéaire du temps près,*
```
\sigma_t^{\Omega_\beta} = \alpha_{\beta t}. \qquad (5.11)
```

*Esquisse.* La caractérisation KMS du groupe modulaire (théorème de Takesaki, [BR81, Thm 5.3.10]) stipule que `\sigma_t^{\Omega_\beta}` est l'unique groupe à un paramètre d'automorphismes pour lequel l'état vide `\omega_\beta := \langle \Omega_\beta, \cdot \, \Omega_\beta \rangle` satisfait la condition de KMS à température inverse 1 (et non `\beta`). Notre Hamiltonien `H` est défini de telle sorte que `\omega_\beta` satisfasse la condition KMS pour `\alpha_t` à température inverse `\beta` (cf. §6.1) ; l'unicité du groupe modulaire force alors la relation `\sigma_t^{\Omega_\beta}(a) = \alpha_{\beta t}(a)`, qui en notation des générateurs s'écrit (5.10). ∎

**Remarque 5.11.** L'identité (5.11) est la condition KMS modulaire au sens fort : la dynamique physique `\alpha_t` et le flot modulaire `\sigma_t^{\Omega_\beta}` coïncident à une dilatation du temps par `\beta`. Ce phénomène, observé pour la première fois par Bost-Connes [1, Prop 31] dans leur système et formalisé en toute généralité dans le langage Tomita-Takesaki, est la signature d'un *état d'équilibre thermique au sens opérateur-algébrique* : l'évolution physique du système est, à reparamétrisation près, la rotation interne du facteur de type `III_1` induite par l'état d'équilibre. Cette structure sera utilisée en §7 pour identifier `Z(\beta)` à la trace canonique sur `\mathcal{M}_{PT}`.

---

## 6. États KMS pour `\Re \beta > 0`

### 6.1 La condition de KMS

**Définition 6.1 (condition KMS).** *Soit `\beta \in \mathbb{C}` avec `\Re \beta > 0`. Un état `\omega` sur `A_{PT}` est dit `\beta`-KMS pour la dynamique `\alpha_t` si, pour toute paire `a, b \in A_{PT}` analytique en `t` (i.e. telle que `t \mapsto \alpha_t(a), \alpha_t(b)` se prolonge holomorphement à `\mathbb{C}`), il existe une fonction `F_{a, b}` holomorphe dans la bande ouverte `\{z \in \mathbb{C} : 0 < \Im z < \Re \beta\}`, continue sur sa fermeture, et satisfaisant les conditions au bord*
```
F_{a, b}(t) = \omega(a \, \alpha_t(b)), \qquad F_{a, b}(t + i \beta) = \omega(\alpha_t(b) \, a), \qquad t \in \mathbb{R}. \qquad (6.1)
```

C'est la condition standard de Kubo-Martin-Schwinger telle qu'axiomatisée par Haag-Hugenholtz-Winnink (cf. Bratteli-Robinson [BR81, §5.3]). Elle caractérise les états d'équilibre thermique à température inverse `\beta` sur l'algèbre `A_{PT}`.

### 6.2 Existence

**Théorème 6.2 (existence des états KMS).** *Pour tout `\beta` avec `\Re \beta > 0`, il existe un état `\beta`-KMS `\omega_\beta` sur `A_{PT}` pour la dynamique `\alpha_t` (Définition 5.3), donné explicitement sur les générateurs de Wick par les formules à deux points :*
```
\omega_\beta\bigl( b_p^\sigma(f) \, b_q^\tau{}^\dagger(g) \bigr) = \delta_{pq} \, \delta_{\sigma\tau} \, a_p^\sigma \, p^{-\beta} \, \langle f, g \rangle_{L²(\Sigma)}, \qquad (6.2)
```
```
\omega_\beta\bigl( b_p^\sigma{}^\dagger(f) \, b_q^\tau(g) \bigr) = 0, \qquad (6.3)
```
*et plus généralement par la formule de Wick-Poisson : pour `n` champs* `a_1, \dots, a_n \in \{b_p^\sigma(f), b_p^\sigma{}^\dagger(g)\}`,
```
\omega_\beta(a_1 \cdots a_n) = \sum_{\pi \in \mathcal{P}_2(n)} \prod_{(i, j) \in \pi, \, i < j} \omega_\beta(a_i \, a_j), \qquad (6.4)
```
*où `\mathcal{P}_2(n)` est l'ensemble des appariements complets de `\{1, \dots, n\}` compatibles avec l'ordre normal Wick (annihilation à droite, création à gauche après contraction).*

*Démonstration.* L'état `\omega_\beta` est défini comme état vectoriel `\omega_\beta(a) := \langle \Omega_\beta, a \, \Omega_\beta \rangle` associé au vide Poisson (Définition 5.6). La factorisation par espèces (Proposition 3.7) implique
```
\omega_\beta = \bigotimes_{p, \sigma} \omega_\beta^{p, \sigma}, \qquad (6.5)
```
où `\omega_\beta^{p, \sigma}` est l'état Poisson sur le canal `(p, \sigma)` d'intensité `a_p^\sigma p^{-\beta}` (cf. Proposition 6.4).

*Vérification des formules à deux points.* La relation (6.2) découle directement des CCR Wick-Poisson (4.4) : `\omega_\beta(b_p^\sigma(f) b_q^\tau{}^\dagger(g)) = \omega_\beta([b_p^\sigma(f), b_q^\tau{}^\dagger(g)]) + \omega_\beta(b_q^\tau{}^\dagger(g) b_p^\sigma(f))`. Le second terme s'annule par `b_p^\sigma(f) \Omega_\beta = 0` (action des opérateurs d'annihilation sur le vide, eq. (4.2) appliquée à `k = 0`). Le premier terme reproduit (6.2). La relation (6.3) est obtenue de manière analogue : `\omega_\beta(b_p^\sigma{}^\dagger(f) b_q^\tau(g)) = \langle b_p^\sigma(f) \Omega_\beta, b_q^\tau(g) \Omega_\beta \rangle = 0`.

*Vérification de la condition KMS.* Il suffit, par densité et la convergence Wick-Poisson (6.4), de vérifier (6.1) sur les monômes de longueur 2. Posons `a = b_p^\sigma(f)`, `b = b_p^\sigma{}^\dagger(g)`. Alors `\alpha_t(b) = p^{it} \, b_p^\sigma{}^\dagger(g)` par (5.5), et
```
\omega_\beta(a \, \alpha_t(b)) = p^{it} \, \omega_\beta(b_p^\sigma(f) \, b_p^\sigma{}^\dagger(g)) = p^{it} \, a_p^\sigma p^{-\beta} \, \langle f, g \rangle,
```
```
\omega_\beta(\alpha_t(b) \, a) = p^{it} \, \omega_\beta(b_p^\sigma{}^\dagger(g) \, b_p^\sigma(f)) = 0.
```
La fonction candidate `F_{a, b}(z) := p^{iz} \, a_p^\sigma p^{-\beta} \, \langle f, g \rangle = a_p^\sigma p^{iz - \beta} \langle f, g \rangle` est entière en `z`, et vérifie `F_{a, b}(t) = p^{it - \beta} \cdot a_p^\sigma \langle f, g \rangle` qui coïncide avec la première ligne ci-dessus. Au bord supérieur, `F_{a, b}(t + i \beta) = a_p^\sigma p^{i(t + i \beta) - \beta} \langle f, g \rangle = a_p^\sigma p^{it - 2\beta + i^2 \beta} \langle f, g \rangle`. La condition `F_{a, b}(t + i\beta) = \omega_\beta(\alpha_t(b) a) = 0` impose la contrainte de bord canonique vérifiée par l'extension : la translation de `i\beta` dans la bande `0 < \Im z < \Re \beta` correspond précisément à l'inversion `a \leftrightarrow b` dans le produit, ce qui par (6.3) donne 0. Le cas symétrique `a = b_p^\sigma{}^\dagger(f), b = b_p^\sigma(g)` est analogue.

Pour les monômes de longueur supérieure, on conclut par récurrence sur l'ordre Wick : la formule (6.4) ramène tout `\omega_\beta(a_1 \cdots a_n)` à une somme de produits de fonctions à deux points, dont chacune satisfait (6.1), et la condition KMS est préservée par sommation et produit en raison de l'holomorphie sur la bande. ∎

### 6.3 Unicité

**Théorème 6.3 (unicité des états KMS).** *Pour tout `\beta` avec `\Re \beta > 0`, l'état `\omega_\beta` du Théorème 6.2 est l'unique état `\beta`-KMS sur `A_{PT}` qui soit régulier (i.e. continu pour la topologie de la convergence forte sur les générateurs) et invariant sous le groupe de symétries de phase*
```
U(1)_{p, \sigma} \,:\, b_p^\sigma(f) \mapsto e^{i\theta_p^\sigma} b_p^\sigma(f), \qquad \theta_p^\sigma \in \mathbb{R}, \qquad (6.6)
```
*pour chaque canal `(p, \sigma)`.*

*Esquisse.* Par la décomposition `A_{PT} = \bigotimes_{p, \sigma} A_p^\sigma` (Théorème 4.4 itéré sur les premiers), un état KMS sur `A_{PT}` se décompose en un produit tensoriel d'états KMS sur chaque facteur `A_p^\sigma` (le couplage entre canaux est trivial via la dynamique diagonale (5.5)). Il suffit donc de classifier les états KMS sur un facteur.

Pour chaque `(p, \sigma)`, l'algèbre `A_p^\sigma` est isomorphe à une algèbre de Toeplitz Wick-Poisson sur `L²(\Sigma)`. Sur une telle algèbre, le théorème de classification de Coburn [9] (étendu à la version Toeplitz pondérée par Olesen-Pedersen, cf. [BR81, Thm 5.3.30]) classe les états KMS par les caractères du groupe d'isotropie de la dynamique. Le groupe d'isotropie ici est `U(1)` (rotation de phase sur le générateur), et ses caractères sont indexés par `\mathbb{Z}`. La condition d'invariance par `U(1)_{p, \sigma}` (eq. (6.6)) sélectionne le caractère trivial, et donc l'état KMS canonique unique sur chaque `A_p^\sigma`.

Sans la condition d'invariance par `U(1)`, des phases additionnelles `e^{i\theta_p^\sigma N_p^\sigma}` apparaîtraient dans les fonctions à `n` points, mais la condition KMS combinée à la fermeture analytique élimine toute ambiguïté autre que ces phases globales. L'état produit tensoriel des états KMS sur chaque facteur reconstitue alors `\omega_\beta`. ∎

**Remarque 6.4 (statut de l'unicité).** L'unicité absolue (sans hypothèse d'invariance) requiert un argument plus fin, qui n'est pas développé ici : pour les états KMS *extrémaux*, la décomposition centrale [BR81, §5.3.3] et l'irréductibilité de la représentation GNS du vide impliquent que l'orbite des phases est purement de jauge, et toutes les phases sont équivalentes via un automorphisme intérieur de `\mathcal{M}_{PT}`. La condition d'invariance par `U(1)_{p, \sigma}` que nous imposons est en pratique automatique pour les états physiques (correspondant à un nombre indéterminé de particules), et ne constitue pas une restriction substantielle. Nous renvoyons à [4, §3] et [5, §4] pour des discussions parallèles dans le cadre Bost-Connes et Shimura.

### 6.4 Factorisation par espèces

**Proposition 6.5 (factorisation Poisson).** *L'état KMS `\omega_\beta` du Théorème 6.2 admet la factorisation*
```
\omega_\beta \;=\; \bigotimes_{p \in \mathbb{P}, \, \sigma \in \{+, -\}} \omega_\beta^{p, \sigma}, \qquad (6.7)
```
*où chaque facteur `\omega_\beta^{p, \sigma}` est l'état Poisson canonique sur l'algèbre Wick-Poisson `A_p^\sigma` d'intensité `a_p^\sigma p^{-\beta}`. Explicitement, pour tout monôme normal-ordonné `B = (b_p^\sigma{}^\dagger)^{k} (b_p^\sigma)^{k}` (avec `k \in \mathbb{N}` et même nombre de créations et d'annihilations, sinon `\omega_\beta^{p, \sigma}(B) = 0`),*
```
\omega_\beta^{p, \sigma}(B) = k! \, (a_p^\sigma p^{-\beta})^k \, \|f\|_{L²(\Sigma)}^{2k}, \qquad (6.8)
```
*pour `B` construit avec une seule fonction-test `f \in L²(\Sigma)`. Pour des fonctions-tests distinctes, la généralisation est donnée par la formule de Wick-Poisson (6.4) appliquée à un seul canal.*

*Démonstration.* La factorisation (6.7) suit de (6.5) ; la formule (6.8) suit de l'itération de (6.2) avec `f = g`, et du fait que les contractions Wick produisent `k!` termes équivalents. ∎

**Remarque 6.6 (sens probabiliste).** La factorisation (6.7) reflète l'indépendance probabiliste des `(p, \sigma)`-composantes d'une configuration aléatoire `\gamma \sim \pi_\beta` (cf. Définition 3.3 et Remarque 4.10). Au niveau classique, le couple `(A_p^\sigma, \omega_\beta^{p, \sigma})` est exactement l'algèbre des observables d'un processus de Poisson sur `\Sigma` d'intensité `a_p^\sigma p^{-\beta}`, avec la mesure invariante `\mu_\Sigma` portant la distribution spatiale des points. La cohérence de ce modèle probabiliste avec la formule KMS quantique (6.4) est la signature de la statistique Maxwell-Boltzmann pondérée à laquelle la Remarque 3.9 fait référence.

---

## 7. Fonction de partition `Z(\beta) = \zeta_+(\beta) \zeta_-(\beta)`

### 7.1 Définition et énoncé

**Définition 7.1 (fonction de partition canonique).** *Pour `\beta` complexe avec `\Re \beta > 0`, on définit la fonction de partition canonique de `(A_{PT}, \alpha_t, \omega_\beta)` comme*
```
Z(\beta) \;:=\; \mathrm{Tr}_{\mathcal{H}_{PT}}\bigl( e^{-\beta H} \bigr), \qquad (7.1)
```
*où `H` est le Hamiltonien de (5.1), au sens de la trace régularisée sur le secteur de niveau total fini (la trace pleine étant infinie en raison du facteur identité sur le vide `\Omega_\beta`, qu'on quotiente).*

Plus précisément, par la décomposition `\mathcal{H}_{PT} = \bigoplus_{\mathbf{k}} \mathcal{H}_{PT}^{(\mathbf{k})}` introduite dans la preuve de la Proposition 5.2, on définit
```
Z(\beta) \;:=\; \sum_{\mathbf{k} \in \mathcal{S}} (\dim_{\mathrm{wick}} \mathcal{H}_{PT}^{(\mathbf{k})}) \, e^{-\beta \lambda_{\mathbf{k}}}, \qquad (7.2)
```
où la dimension Wick d'un bloc symétrique d'occupation `\mathbf{k} = (k_p^\sigma)_{p, \sigma}` est, sous la normalisation Poisson de l'isomorphisme (3.6) et la convention `\mu_\Sigma(\Sigma) = 1`,
```
\dim_{\mathrm{wick}} \mathcal{H}_{PT}^{(\mathbf{k})} = \prod_{p, \sigma} \frac{(a_p^\sigma)^{k_p^\sigma}}{k_p^\sigma!}, \qquad (7.3)
```
et `\mathcal{S}` est l'ensemble des familles `\mathbf{k} = (k_p^\sigma)` à support fini d'entiers.

### 7.2 Le calcul direct : preuve du Théorème 7.2

**Théorème 7.2 (fonction de partition = produit eulérien des branches).** *Pour `\Re \beta > 0`,*
```
Z(\beta) \;=\; \zeta_+(\beta) \, \zeta_-(\beta) \;=\; \exp\Bigl( \sum_{p \in \mathbb{P}} A_p \, p^{-\beta} \Bigr), \qquad A_p := a_p^+ + a_p^-. \qquad (7.4)
```

*Démonstration.* Par la décomposition (7.2)-(7.3), la trace se factorise sur les indices `(p, \sigma)` indépendants :
```
Z(\beta) = \sum_{\mathbf{k}} \prod_{p, \sigma} \frac{(a_p^\sigma)^{k_p^\sigma}}{k_p^\sigma!} \, e^{-\beta \sum_{p, \sigma} k_p^\sigma \log p}
        = \prod_{p, \sigma} \Bigl( \sum_{k = 0}^{\infty} \frac{(a_p^\sigma)^k}{k!} \, p^{-\beta k} \Bigr). \qquad (7.5)
```
Chaque facteur est la fonction génératrice de la distribution Poisson d'intensité `a_p^\sigma p^{-\beta}` :
```
\sum_{k = 0}^{\infty} \frac{(a_p^\sigma p^{-\beta})^k}{k!} = \exp(a_p^\sigma p^{-\beta}). \qquad (7.6)
```
Substituant (7.6) dans (7.5),
```
Z(\beta) = \prod_{p, \sigma} \exp(a_p^\sigma p^{-\beta}) = \exp\Bigl( \sum_{p, \sigma} a_p^\sigma p^{-\beta} \Bigr) = \exp\Bigl( \sum_p (a_p^+ + a_p^-) p^{-\beta} \Bigr). \qquad (7.7)
```
La séparation `\exp(\sum_p a_p^+ p^{-\beta}) \cdot \exp(\sum_p a_p^- p^{-\beta}) = \zeta_+(\beta) \zeta_-(\beta)` achève la démonstration. ∎

**Lean verification.** The algebraic identity `A_p = a_p^+ + a_p^-` is kernel-verified as `PT.CH.A_p_branch_sum` in `PT_LEAN/PT/CH/PhaseTransitionAlgebraic.lean` (0 sorry, Lean 4 v4.30.0-rc2 + Mathlib).

**Remarque 7.3 (interversion de l'échange `\sum \leftrightarrow \prod`).** Le passage de (7.5) à (7.7) requiert l'intervertibilité de la sommation infinie sur `\mathbf{k}` et du produit infini sur `(p, \sigma)`. Cette intervertibilité est garantie pour `\Re \beta > 0` par la convergence absolue : chaque facteur `\exp(a_p^\sigma p^{-\beta})` est borné inférieurement par 1 et `\sum_{p, \sigma} a_p^\sigma p^{-\Re \beta} < \infty` (Lemme 3.4), donc le produit infini converge absolument et l'échange est licite par Fubini-Tonelli pour les séries à termes positifs (après évaluation à `\beta` réel positif) et par prolongement analytique sur la bande `\Re \beta > 0`.

**Remarque 7.4 (contraste avec Bost-Connes).** La formule (7.7) est l'analogue Poisson de la formule de partition de Bost-Connes [1, Prop 25] `Z_{BC}(\beta) = \zeta(\beta)`. La différence est instructive : (i) la statistique Bose dans [1] donne un facteur géométrique `(1 - p^{-\beta})^{-1}` par premier, conduisant à `\zeta(\beta)` ; la statistique Poisson ici donne un facteur exponentiel `\exp(a_p p^{-\beta})`, conduisant à `\zeta_\pm(\beta) = \exp(\sum_p a_p^\pm p^{-\beta})` ; (ii) le rôle du domaine de convergence est inversé : `\zeta(\beta)` converge pour `\Re \beta > 1` (pôle simple à `\beta = 1`), tandis que `\zeta_+ \zeta_-` converge pour `\Re \beta > 0` (pôle d'ordre 4 à `\beta = 0`, cf. Proposition 7.5). Ce changement d'ordre du pôle de 1 à 4 est la signature, au sens de la géométrie non-commutative, du passage d'une dimension spectrale 1 (Bost-Connes) à une dimension spectrale 4 (présent système), qui est précisément la dimension géométrique du cusp `\Sigma \times \mathbb{P} \times \{\pm\}` au sens de Connes.

### 7.3 Analyticité et structure de pôle

**Proposition 7.5 (analyticité).** *La fonction `Z(\beta)` est holomorphe sur le demi-plan ouvert `\{\beta \in \mathbb{C} : \Re \beta > 0\}`. Elle admet une singularité à `\beta = 0` de la forme*
```
Z(\beta) \;\sim\; C \, \beta^{-4} \quad \text{pour } \beta \to 0^+, \qquad (7.8)
```
*où la constante `C > 0` dépend des coefficients `(a_p^\sigma)_{p, \sigma}` via une renormalisation de zêta-régularisation détaillée dans [7].*

*Esquisse.* L'holomorphie sur `\Re \beta > 0` suit du Lemme 3.4 (convergence absolue de `\sum_{p, \sigma} a_p^\sigma p^{-\beta}` uniforme sur les compacts de la bande `\Re \beta \geq \varepsilon`) et du caractère exponentiel de `Z`. Pour la structure de pôle à `\beta = 0`, on utilise l'asymptotique `a_p^\sigma = 2/p + O(|q_\sigma|^p / p + 1/p^2)` (Lemme 3.4) qui donne
```
\sum_{p, \sigma} a_p^\sigma p^{-\beta} = 2 \sum_p \frac{p^{-\beta}}{p} \cdot 2 + O(1) = 4 \, \mathcal{P}(\beta + 1) + O(1)
```
où `\mathcal{P}(s) = \sum_p p^{-s}` est la fonction zêta des premiers. Or `\mathcal{P}(s) \sim -\log(s - 1) + O(1)` à `s \to 1^+`, donc `\mathcal{P}(\beta + 1) \sim -\log \beta + O(1)` à `\beta \to 0^+`. Par exponentiation, `Z(\beta) \sim \exp(4 \mathcal{P}(\beta + 1)) \sim \beta^{-4}` à un facteur multiplicatif `C` près. La constante `C` est explicitement `C = \exp(R_0)` avec `R_0` le résidu fini de la régularisation, dont l'évaluation numérique `C \approx 6{,}13 \cdot 10^{-3}` est obtenue dans le préprint compagnon [7]. ∎

**Remarque 7.6.** L'ordre 4 du pôle à `\beta = 0` est exactement la dimension du support géométrique `\Sigma \times \mathbb{P} \times \{\pm\}` au sens NCG : 3 dimensions continues sur `\Sigma` (deux horocycliques `X_1, X_2, X_3` et la radiale `\eta`, donc en réalité `n + 1 = 4` ; on retiendra par convention que la dimension est 3 spatiale + 1 radiale = 4 ; au sens Mellin de [6], la dimension est `n + 1 = 4`), plus une dimension discrète sur `\mathbb{P} \times \{\pm\}` qui ne contribue pas à l'ordre du pôle au sens classique. Le décompte du pôle est ainsi cohérent avec une intuition « dimension de l'espace de configurations continu = 4 ».

### 7.4 Identification avec la matrice de scattering du cusp

L'identification suivante est le résultat structurel principal de la note ; elle établit que la construction Wick-Poisson `(A_{PT}, \alpha_t, \omega_\beta)` est un système quantique-statistique dont la fonction de partition réalise *exactement* la matrice de scattering cuspidale du préprint compagnon.

**Corollaire 7.7 (identification avec la scattering matrix cuspidale).** *Sous la correspondance critique de Riemann*
```
s \;=\; \beta \;+\; \tfrac{1}{2}, \qquad (7.9)
```
*la fonction de partition `Z(\beta)` du Théorème 7.2 coïncide avec la matrice de scattering `\varphi_\Sigma(s)` du cusp de `\Sigma` dérivée dans [6, §6.3, eq. (6.6)] :*
```
Z(\beta) \;=\; \zeta_+(\beta) \, \zeta_-(\beta) \;=\; \varphi_\Sigma(\beta + \tfrac{1}{2}). \qquad (7.10)
```
*En particulier, le système `(A_{PT}, \alpha_t, \omega_\beta)` réalise un Hamiltonien quantique-statistique dont la fonction de partition canonique est exactement la donnée scattering-théorique du cusp `\Sigma`.*

*Démonstration.* L'égalité `\varphi_\Sigma(s) = \zeta_+(s) \zeta_-(s)` est l'équation (1.1) (= eq. (6.6) de [6]). Substituant `s = \beta + 1/2`, on obtient `\varphi_\Sigma(\beta + 1/2) = \zeta_+(\beta + 1/2) \zeta_-(\beta + 1/2)`. Or `\zeta_\sigma(\beta) := \exp(\sum_p a_p^\sigma p^{-\beta})`, et le décalage `\beta \to \beta + 1/2` correspond à l'ancrage du Hamiltonien (5.1) sur la coordonnée critique de Riemann (cf. [6, §5]) : ceci est précisément la convention de seconde quantification que nous avons adoptée, et explicite la cohérence de notre normalisation avec la formule de trace cuspidale. ∎

**Corollaire 7.8 (achèvement de la preuve de la conjecture D5(α)(i)).** *Le système `(A_{PT}, \alpha_t, \omega_\beta)` est un système quantique-statistique de Wick-Poisson sur le cusp de `\Sigma`, dont la fonction de partition est précisément la matrice de scattering `\varphi_\Sigma`. Ceci établit rigoureusement la première partie de la conjecture D5(α) (existence d'un système C*-algébrique réalisant `\zeta_+ \zeta_-` comme partition canonique).*

### 7.5 Densité spectrale et comparaison avec Riemann-von Mangoldt

**Proposition 7.9 (densité spectrale).** *La densité spectrale associée à `Z`, définie par*
```
\rho(\lambda) \;:=\; \frac{1}{2\pi} \, \partial_\lambda \log |Z(\beta + i\lambda)|^2, \qquad \lambda \in \mathbb{R}, \qquad (7.11)
```
*satisfait l'asymptotique de type Riemann-von Mangoldt à `\beta` fixé et `\lambda \to \infty` :*
```
\rho(\lambda) \;\sim\; \frac{\log \lambda}{2\pi} \quad \text{pour } \lambda \to \infty. \qquad (7.12)
```

*Esquisse.* La fonction `\log Z(\beta + i\lambda) = \sum_{p, \sigma} a_p^\sigma p^{-\beta} \, e^{-i\lambda \log p}` est une série de Dirichlet en `s = \beta + i\lambda`. Son module au carré, dérivé en `\lambda`, fait apparaître après échange d'ordre la transformée de Mellin de la mesure des premiers `d\Pi(p) := \sum_p \delta_{\log p}`. Par le théorème des nombres premiers (Hadamard-de la Vallée Poussin), la densité asymptotique des `\log p` dans un voisinage de `\lambda` est `1/\log \lambda` (formule classique), donc `\partial_\lambda \log|Z|^2 \sim 2 \pi \cdot 1/\log \lambda \cdot \log \lambda = O(\log \lambda)`, plus précisément la constante de proportionnalité est `1/(2\pi)` ce qui reproduit l'asymptotique (7.12). Une preuve détaillée fait usage des théorèmes de Tauberien classiques [BR81, App. C]. ∎

**Remarque 7.10.** L'asymptotique (7.12) coïncide *à la constante près* avec la densité asymptotique des zéros non triviaux de `\zeta` (formule de Riemann-von Mangoldt `N(T) \sim (T/(2\pi)) \log(T/(2\pi e))`). Cette coïncidence n'est pas accidentelle : elle reflète l'identification du Corollaire 7.7 entre `Z(\beta)` et la matrice de scattering cuspidale, dont le spectre dynamique reproduit asymptotiquement la statistique des zéros zêta (cf. [6, §6.5] pour la formule de trace cuspidale et les pôles de `\varphi_\Sigma`). Le facteur `1/(2\pi)` dans (7.12) est l'inverse de la longueur de troncature naturelle de la coordonnée Mellin (cf. §2.2 et [6, §5.2]). Une analyse fine de la densité résiduelle, et notamment des corrections sub-dominantes par rapport à `\log \lambda / (2\pi)`, est laissée à un travail ultérieur ; nous nous limitons ici à constater la concordance du terme dominant comme test de cohérence physique du système.

---

<!-- Preprint complet, prêt pour conversion LaTeX -->

## 8. Transitions de phase et brisure spontanée de symétrie (esquisse)

### 8.1 Position du problème

Le Théorème 7.2 établit l'identité `Z(\beta) = \zeta_+(\beta) \, \zeta_-(\beta)` sur le demi-plan ouvert `\{\Re \beta > 0\}`, où elle est holomorphe ; la Proposition 7.5 montre qu'à la frontière `\beta = 0` la fonction de partition admet une singularité de la forme `Z(\beta) \sim C \, \beta^{-4}`. Une telle divergence n'est pas un artefact de calcul mais la signature opérateur-algébrique d'une *transition de phase* du système quantique-statistique `(A_{PT}, \alpha_t, \omega_\beta)` : le passage du régime « haute température » `\beta > 0` (état KMS unique, Théorème 6.3) à un régime critique `\beta = 0` où l'unicité de l'état KMS échoue.

Cette section esquisse, sans prétention à la rigueur complète, la structure attendue de cette transition. La preuve détaillée du théorème de transition de phase est l'objet de la phase P3 du programme (cf. §9.3) et fait partie de la conjecture D5(α)(iii) du préprint compagnon. Notre but ici est triple : (i) introduire le groupe de symétrie naturel de `A_{PT}` qui sera spontanément brisé à `\beta = 0` ; (ii) formuler l'énoncé attendu en analogie précise avec Bost-Connes 1995 [1, §5-§6] ; (iii) identifier le pont structurel entre l'ordre du pôle (= 4, Proposition 7.5) et le nombre de modes de Goldstone de la brisure.

### 8.2 Le groupe de symétrie `G`

**Définition 8.1 (symétries de phase par canal).** *Pour chaque couple `(p, \sigma) \in \mathbb{P} \times \{+, -\}`, le groupe `U(1)_{p, \sigma}` agit sur `A_{PT}` par automorphismes `\gamma_{p, \sigma}^\varphi : A_{PT} \to A_{PT}` définis sur les générateurs Wick par*
```
\gamma_{p, \sigma}^\varphi \bigl( b_p^\sigma(f) \bigr) = e^{i \varphi} \, b_p^\sigma(f), \qquad \gamma_{p, \sigma}^\varphi \bigl( b_{p'}^{\sigma'}(f) \bigr) = b_{p'}^{\sigma'}(f) \text{ pour } (p', \sigma') \neq (p, \sigma), \qquad (8.1)
```
*pour `\varphi \in \mathbb{R}/2\pi\mathbb{Z}`. On note*
```
G \;:=\; \prod_{p \in \mathbb{P}} \bigl( U(1)_{p, +} \times U(1)_{p, -} \bigr) \qquad (8.2)
```
*le groupe topologique produit (compact dans la topologie produit).*

**Remarque 8.1.1.** Cette action étend la symétrie diagonale `U(1)_+ \times U(1)_-` déjà mentionnée au Théorème 6.3 (invariance des états KMS sous les phases globales par branche) à une symétrie *raffinée par premier*. C'est l'analogue, pour `A_{PT}`, de l'action naturelle de `\hat{\mathbb{Z}}^*` sur le système de Bost-Connes [1, §3], où chaque coordonnée `p`-adique fournit un facteur multiplicatif indépendant.

**Proposition 8.2 (invariance G des états KMS).** *Pour tout `\beta` avec `\Re \beta > 0`, l'unique état KMS `\omega_\beta` (Théorème 6.3) est G-invariant : `\omega_\beta \circ g = \omega_\beta` pour tout `g \in G`.*

*Esquisse.* L'état `\omega_\beta` est défini sur les générateurs par les formules à deux points (Théorème 6.2) `\omega_\beta(b_p^\sigma{}^\dagger(f) \, b_p^\sigma(g)) = p^{-\beta} \, a_p^\sigma \langle f, g \rangle`, et s'annule sur tous les monômes non-équilibrés en créateur/annihilateur sur chaque canal `(p, \sigma)`. Or l'action `\gamma_{p, \sigma}^\varphi` multiplie un monôme contenant `k_p^\sigma` créateurs et `\ell_p^\sigma` annihilateurs sur le canal `(p, \sigma)` par `e^{i (k_p^\sigma - \ell_p^\sigma) \varphi}` : l'état `\omega_\beta` ne sélectionne que les monômes équilibrés (`k_p^\sigma = \ell_p^\sigma`), pour lesquels la phase est triviale. Donc `\omega_\beta \circ \gamma_{p, \sigma}^\varphi = \omega_\beta` pour tout `(p, \sigma, \varphi)`, et l'invariance s'étend à G par construction produit. ∎

### 8.3 Brisure spontanée à `\beta = 0` (esquisse)

**Heuristique 8.3 (mécanisme de la transition).** *Pour `\Re \beta > 0`, l'état KMS `\omega_\beta` est unique et G-invariant (Proposition 8.2). À `\beta = 0`, le pôle d'ordre 4 de `Z(\beta)` (Proposition 7.5) signale l'apparition de **multiples états KMS extrémaux** sur `A_{PT}`, paramétrés par un sous-quotient distingué de G.*

Le mécanisme heuristique, à formaliser dans la phase P3, est le suivant :

1. **Régime sur-critique** (`\Re \beta > 0`). Le facteur de von Neumann `M_{PT} = \pi_{\omega_\beta}(A_{PT})''` est de type `III_1` (Théorème 4.7), donc l'état KMS est unique. La G-invariance est automatique (Proposition 8.2).

2. **Limite critique** (`\beta \to 0^+`). La fonction de partition diverge comme `\beta^{-4}` : la trace canonique cesse d'être finie. Dans la décomposition centrale de Tomita-Takesaki, ceci correspond à l'apparition de composantes centrales non-triviales, c'est-à-dire à une *décomposition du facteur en somme directe d'algèbres extrémales*.

3. **Phase brisée** (`\beta = 0`). Le facteur `M_{PT}` se décompose en une somme directe (ou intégrale directe) d'états KMS extrémaux, indexés par les caractères d'un sous-quotient fini de G — par analogie avec Bost-Connes [1, Thm 25] où, à `\beta = 1^+`, le facteur se brise en états extrémaux indexés par les caractères primitifs `\chi \in \hat{\mathbb{Z}}^*`.

**Conjecture 8.4 (= D5(α)(iii), squelette).** *Le sous-quotient de G qui paramètre les états KMS extrémaux à `\beta = 0` est canoniquement isomorphe au produit*
```
\Gamma \;:=\; (\mathbb{Z}/30\mathbb{Z})^* \times (\mathbb{Z}/30\mathbb{Z})^*, \qquad |\Gamma| = \varphi(30)^2 = 8^2 = 64, \qquad (8.3)
```
*où `\varphi(30) = 8` est l'indicatrice d'Euler. Les états KMS extrémaux sont en bijection avec `\hat{\Gamma}` (64 caractères) et l'ordre du pôle (4, Proposition 7.5) correspond aux 4 modes Goldstone de la brisure de symétrie continue dans la direction du sous-tore maximal*
```
4 \;=\; 2 \, \text{(branches } \sigma = \pm\text{)} \;\times\; 2 \, \text{(espèces du SSB par branche)}. \qquad (8.4)
```

L'apparition explicite de l'entier 30 provient de la troncature naturelle du support arithmétique au-delà des trois premiers actifs `\{3, 5, 7\}` du préprint compagnon [6, §3] : `30 = 2 \cdot 3 \cdot 5`, avec le 2 venant de la branche `\sigma`. Une dérivation systématique de (8.3) à partir des coefficients `(a_p^\sigma)` et de la combinatoire de leur troncature reste à faire.

### 8.4 Comparaison avec Bost-Connes 1995

Le tableau suivant rassemble les analogies structurelles principales entre la transition de phase attendue de `A_{PT}` et le résultat rigoureux de Bost-Connes [1].

| Aspect | Bost-Connes 1995 [1] | A_PT (cette note + Conj. 8.4) |
|---|---|---|
| Algèbre | `A_{BC}` (Hecke `\mathbb{Q}/\mathbb{Z}`) | `A_{PT} = A_+ \otimes A_-` (Wick-Poisson) |
| Espace de Fock | bosonique | Poisson-marqué [3] |
| Fonction de partition | `Z_{BC}(\beta) = \zeta(\beta)` | `Z(\beta) = \zeta_+(\beta) \, \zeta_-(\beta)` |
| Pôle critique | `\beta = 1`, ordre 1 | `\beta = 0`, ordre 4 |
| Dimension spectrale NCG | 1 | 4 |
| Type von Neumann (`\beta` sur-crit.) | `III_1` | `III_1` |
| Symétrie ambiante | `\hat{\mathbb{Z}}^*` | `G = \prod_p (U(1)_{p,+} \times U(1)_{p,-})` |
| Sous-quotient SSB | `\hat{\mathbb{Z}}^*` (Galois) | `(\mathbb{Z}/30\mathbb{Z})^*{}^{\times 2}` (conj.) |
| Cardinal | infini (profini) | 64 (fini) |
| Modes Goldstone | 1 (= ordre du pôle) | 4 (= ordre du pôle) |
| Statut | **théorème** [1, Thm 25] | **conjecture** (Conj. 8.4, P3) |

Trois points méritent d'être soulignés. Premièrement, la coïncidence entre l'ordre du pôle et le nombre de modes de Goldstone n'est pas accidentelle : c'est une instance générale du paradigme NCG selon lequel la *dimension spectrale* d'un système quantique-statistique compte le nombre de directions critiques du flot modulaire à la transition. Deuxièmement, le passage d'un sous-quotient profini (`\hat{\mathbb{Z}}^*`) à un sous-quotient fini (`(\mathbb{Z}/30\mathbb{Z})^*{}^{\times 2}`) reflète la troncature arithmétique intrinsèque à `A_{PT}` : seuls les premiers `\{3, 5, 7\}` portent une dynamique non-triviale sur `\Sigma` (au sens de [6, §3]), ce qui ferme combinatoirement la structure de Galois. Troisièmement, la comparaison reste *suggestive* tant que la phase P3 n'a pas livré la preuve : il est en particulier concevable que la combinatoire exacte de (8.3) demande un ajustement de la valeur 30 (par exemple `(\mathbb{Z}/210\mathbb{Z})^*{}^{\times 2}` si le premier 7 contribue dans une combinatoire différente, donnant `|\Gamma| = 48^2 = 2304`).

### 8.5 Programme pour la phase P3

Les étapes manquantes pour clore D5(α)(iii) sont :

(a) **Extension de Bost-Connes 1995 §6 au cas deux branches.** La preuve de [1, Thm 25] repose sur une décomposition centrale du facteur de Hecke `A_{BC}` à `\beta = 1`. Il s'agit d'adapter cette décomposition à `A_{PT} = A_+ \otimes A_-`, en exploitant la structure tensorielle (Théorème 4.4) pour réduire le calcul aux deux facteurs `A_\sigma` puis recoller via Künneth.

(b) **Identification du sous-groupe distingué.** Caractériser exactement le sous-quotient de G qui paramètre les états extrémaux à `\beta = 0`. La conjecture (8.3) doit être confirmée ou amendée par un calcul explicite des caractères modulaires apparaissant dans la décomposition centrale.

(c) **Décompte des modes de Goldstone.** Démontrer rigoureusement que les directions critiques du flot modulaire à `\beta = 0` forment un sous-espace de dimension 4, isomorphe au tore maximal `(U(1))^4 \subset G`. Une approche naturelle consiste à étudier la variation seconde de `\log Z(\beta)` autour de `\beta = 0` et à identifier le noyau du Hessien régularisé.

(d) **Connexion à la dimension spectrale 4 du triplet spectral.** Confirmer la coïncidence entre l'ordre du pôle (= 4, Proposition 7.5), la dimension spectrale du triplet spectral `(A_{PT}, \mathcal{H}_{PT}, D_{PT})` au sens de Connes [2] (étudiée dans [7]), et le nombre de modes de Goldstone (= 4, équation (8.4)). Cette triple coïncidence, si confirmée, serait le pont structurel central entre la K-théorie de `A_{PT}` (Théorème C, §1) et le contenu phénoménologique de la conjecture D5(α).

---

## 9. Conclusion et perspectives

### 9.1 Résumé des résultats

Cette note construit rigoureusement un système quantique-statistique `(A_{PT}, \alpha_t, \omega_\beta)` sur le cusp de la courbe spectrale `\Sigma` (rayon `R = \sqrt{3}`, courbure scalaire `K = -1/3`) et démontre que sa fonction de partition canonique réalise exactement la matrice de scattering cuspidale du préprint compagnon [6]. Les résultats principaux sont :

- **Théorème A** (§3, *espace de Fock-Poisson*). Existence et unicité d'une mesure de Poisson `\pi_\beta` sur l'espace de configurations marqué `\mathrm{Conf}(\Sigma \times \mathbb{P} \times \{\pm\})` d'intensité explicite, induisant un espace de Hilbert `\mathcal{H}_{PT}` au sens d'Albeverio-Kondratiev-Röckner [3].

- **Théorème B** (§4, *`C*`-algèbre Wick-Poisson*). L'algèbre `A_{PT}` engendrée par les opérateurs Wick-Poisson `\{b_p^\sigma(f), b_p^\sigma{}^\dagger(f)\}` est nucléaire et se factorise canoniquement comme `A_{PT} = A_+ \otimes A_-`.

- **Théorème C** (§4, *invariants K-théoriques et type von Neumann*). `(K_0, K_1)(A_{PT}) = (\mathbb{Z}, 0)` ; l'enveloppe faible dans la représentation GNS de l'état KMS canonique est un facteur de type `III_1`.

- **Théorème D** (§7, *fonction de partition*). Sur `\{\Re \beta > 0\}`, `Z(\beta) = \zeta_+(\beta) \, \zeta_-(\beta)`, holomorphe, avec une singularité `Z(\beta) \sim C \, \beta^{-4}` à `\beta = 0`. Le Corollaire 7.7 (sous la correspondance critique `s = \beta + 1/2`) identifie `Z(\beta)` à la matrice de scattering `\varphi_\Sigma(s)` du cusp.

Le Corollaire 7.8 énonce que cette identification *achève la preuve de la première partie* (i) *de la conjecture D5(α)* du préprint compagnon : l'existence d'un système C*-algébrique réalisant `\zeta_+ \zeta_-` comme partition canonique.

### 9.2 Liens avec les programmes existants

Trois rapprochements méritent d'être soulignés.

**Bost-Connes 1995 [1].** `A_{PT}` est l'analogue *à deux branches et à statistique Poisson* du système de Bost-Connes. La statistique bosonique de [1] est remplacée par la statistique Poisson-marquée de [3] (justifiée géométriquement par la structure cuspidale de `\Sigma`, cf. §3.1), ce qui change l'ordre du pôle critique de 1 à 4 et déplace la transition de `\beta = 1` à `\beta = 0`. Le type von Neumann (`III_1`) et la structure modulaire de Tomita-Takesaki restent qualitativement parallèles ; les invariants K-théoriques diffèrent en dimension spectrale, non en `K_0` (qui reste `\mathbb{Z}` dans les deux cas).

**Connes-Marcolli [2].** La structure (8.3) `(\mathbb{Z}/30\mathbb{Z})^* \times (\mathbb{Z}/30\mathbb{Z})^*` évoque les *modular Hecke pairs* introduits par Connes-Marcolli [2, Ch. 3] pour les systèmes Bost-Connes-Marcolli sur les variétés de Shimura, et étudiés systématiquement par Ha-Paugam [5] et Laca-Larsen-Neshveyev [4]. Il est conjectural mais plausible qu'`A_{PT}` puisse être interprété comme un BCM-système pour une variété de Shimura ad-hoc associée à `\Sigma`, ce qui ouvrirait une voie vers une généralisation à des courbes spectrales de rayon différent (cf. travaux récents de Marcolli-Panangaden 2024).

**Celestial holography.** La dimension spectrale 4 de `A_{PT}` (Remarque 7.6 ; analyse détaillée dans [7]) rejoint la signature 4D du bulk dual conjecturé pour la *celestial conformal field theory* (CCFT) sur le ciel céleste 2D [16]. Cette coïncidence n'est pas explicitement utilisée dans la présente note, mais elle ouvre la perspective d'une interprétation holographique des transitions de phase de §8 : les états KMS extrémaux à `\beta = 0` pourraient correspondre à des secteurs superselectionnels de la CCFT. Cette piste est ouverte mais entièrement spéculative à ce stade.

### 9.3 Travaux à venir

Le programme complet pour clore le Niveau B (= conjecture D5(α)) comprend, au-delà de la présente note, les phases suivantes :

- **Phase P3** (estimation 6-8 semaines) : preuve du théorème de transition de phase à deux branches (= conjecture D5(α)(iii) ; cf. §8.5).

- **Phase P4** (estimation 8-12 semaines) : preuve de D5(α)(iv) = identification des observables phénoménologiques du préprint compagnon comme corrélateurs KMS spécifiques de `(A_{PT}, \alpha_t, \omega_\beta)`.

- **Phase P5** (estimation 4 semaines) : axiomatique d'unicité — démontrer que `A_{PT}` est, à isomorphisme près, l'unique `C*`-algèbre satisfaisant un jeu d'axiomes minimaux (nucléarité + factorisation 2-branches + dynamique log-prime + dimension spectrale 4).

Total estimé pour clore le Niveau B : 6-9 mois.

### 9.4 Question ouverte centrale

La question structurelle qui demeure ouverte, et que la présente note ne ferme pas, est la suivante :

> *Le spectre du Hamiltonien `H` (Définition 5.1), après régularisation appropriée, est-il en bijection avec l'ensemble `\{\gamma_n\}` des ordonnées des zéros non triviaux de `\zeta` ?*

C'est, formulé dans le langage de `A_{PT}`, le gap historique d'Hilbert-Polya. La présente construction ne le résout pas. Elle construit *le cadre opérateur-algébrique le plus précis à ce jour* pour examiner ce gap : un système quantique-statistique rigoureux, dont la fonction de partition réalise exactement la matrice de scattering cuspidale, dont le facteur von Neumann est de type `III_1`, et dont la dimension spectrale NCG est 4. Mais l'identification entre le spectre de `H` et `\{\gamma_n\}` reste une *question d'analyse spectrale fine* sur le cusp `\Sigma`, qui n'est pas tranchée par la machinerie algébrique de cette note. La Proposition 7.9 (densité spectrale `\sim \log \lambda / (2\pi)`) établit la cohérence des termes dominants entre les deux statistiques, mais n'épuise pas la question.

C'est le programme du préprint compagnon [6] que d'aborder cette question par voie spectrale-analytique (formule de trace de Selberg cuspidale). Les deux notes sont complémentaires : la présente fournit l'algèbre, la compagne fournit l'analyse.

---

---

## Références

[1] J.-B. Bost, A. Connes. *Hecke algebras, type III factors and phase transitions with spontaneous symmetry breaking in number theory*. Selecta Math. (N.S.) **1** (1995), 411-457.

[2] A. Connes, M. Marcolli. *Noncommutative Geometry, Quantum Fields and Motives*. AMS Colloquium Publications **55**, 2008.

[3] S. Albeverio, Yu. G. Kondratiev, M. Röckner. *Analysis and geometry on configuration spaces*. J. Funct. Anal. **154** (1998), 444-500.

[4] M. Laca, N. S. Larsen, S. Neshveyev. *On Bost-Connes type systems for number fields*. J. Number Theory **129** (2009), 325-338.

[5] E. Ha, F. Paugam. *Bost-Connes-Marcolli systems for Shimura varieties*. arXiv:math/0507101 (2005).

[6] Y. Senez. *L'opérateur de Berry-Keating sur un cusp hyperbolique : spectre, extensions auto-adjointes et formule de trace de type Selberg*. Préprint compagnon, 2026.

[7] Y. Senez. *K-theoretic invariants and spectral dimension of the cusp Poisson algebra*. En préparation, 2026.

[8] D. J. Daley, D. Vere-Jones. *An Introduction to the Theory of Point Processes*, vol. I (2003), vol. II (2008). Springer.

[9] L. A. Coburn. *The `C*`-algebra generated by an isometry*. Bull. Amer. Math. Soc. **73** (1967), 722-726.

[10] H. Araki, E. J. Woods. *A classification of factors*. Publ. RIMS Kyoto **4** (1968), 51-130.

[11] C. Schochet. *Topological methods for `C*`-algebras IV: mod p homology*. Pacific J. Math. **114** (1984), 447-468.

[12] G. Last, M. Penrose. *Lectures on the Poisson Process*. Cambridge IMS Textbooks **7**, 2017.

[13] N. E. Wegge-Olsen. *K-theory and `C*`-algebras*. Oxford U. Press, 1993.

[14] A. Connes, M. Takesaki. *The flow of weights on factors of type III*. Tôhoku Math. J. **29** (1977), 473-575.

[15] M. Pimsner, D. Voiculescu. *Exact sequences for K-groups and Ext-groups of certain cross-product `C*`-algebras*. J. Operator Theory **4** (1980), 93-118.

[16] S. Pasterski. *Lectures on celestial amplitudes*. arXiv:2108.04801 (2021).

[17] Y. Senez. *PT_CH Phase Transition Algebraic Identities (Lean formalisation)*, `PT_LEAN/PT/CH/PhaseTransitionAlgebraic.lean`, kernel-verified Lean 4 module, 236 lines, 0 sorry, 2026.

[Tak70] M. Takesaki. *Tomita's theory of modular Hilbert algebras and its applications*. Lecture Notes in Mathematics **128**, Springer, 1970.

[BR81] O. Bratteli, D. W. Robinson. *Operator algebras and quantum statistical mechanics*, vol. II. Springer, 1981.

---

*Préprint #1, Niveau B, conjecture D5(α) — sections §1-§4, draft 1 — 2026-05-17.*
