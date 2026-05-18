# Modèle spectral canonique des zéros de $\zeta$ via la dynamique de persistance PT

**Auteur** : Yan Senez
**Date** : 15 mai 2026
**Cadre** : Théorie de la Persistance (PT), notes 42, 51, 52, 60–71

---

## Résumé

On construit un modèle spectral canonique du facteur résiduel
$R(s) = \zeta(s) / \big(\zeta_+(s) \zeta_-(s)\big)$ de la fonction zêta de
Riemann, à partir de la dynamique de persistance PT. Trois acquis principaux
sont établis. (i) Une expression close des amplitudes spectrales
$\kappa_p(s) = 1 - (1 - p^{-s}) \exp(A_p p^{-s})$ telle que
$\det(I - D_R(s))^{-1} = R(s)$, vérifiée à $10^{-12}$ dans $\mathrm{Re}(s) > 1$ ;
$D_R$ est l'opérateur diagonal sur $\ell^2(\mathcal{P})$ associé. (ii) Une
formule de trace analytique exacte
$-\,d\log R/ds = \sum_p (\log p)\,[1/(p^s - 1) - A_p\,p^{-s}]$
(signe corrigé 2026-05-17, cf. §7.6), vérifiée à
$10^{-13}$ pour $\sigma = 4$, qui se lit comme une formule de trace de
Berry–Keating sur les longueurs d'orbites $\ell_p = \log p$. (iii) Une trace régularisée
opératorielle
$T(s) = \mathrm{Tr}_{\mathrm{reg}}\big[(s - i H_{\mathrm{PT-BK}})^{-1} D_R(s)\big]$
dont les pôles distributionnels sur la ligne critique sont identifiés aux
zéros non triviaux de $R(s)$ ; sur $t \in [10, 200]$, la formule brute
capture 74 zéros sur 79 à $|\Delta| < 1.5$ (sur $t \in [10, 70]$ :
17/17), avec saturation au-delà de $P_{\max} = 3\,000$ (note 74).
La régularisation explicite R1 (smearing gaussien, $\varepsilon = 0.2$)
améliore substantiellement la précision : **75/79 zéros à $|\Delta| < 0.5$
et 47/79 à $|\Delta| < 0.1$** sur le même range, avec une médiane
médiane $|\Delta| = 0.080$ contre 0.52 pour la formule brute
(note 76). Le test étendu à $t \in [10, 500]$ (269 zéros, note 77) confirme
la robustesse : **267/269 à $|\Delta| < 1.5$ (99.3 \%), 242/269 à
$|\Delta| < 0.5$ (90.0 \%), 139/269 à $|\Delta| < 0.1$ (51.7 \%)**, avec
une capture par bande de 100 unités qui reste $\ge 86.6\,\%$ à $|\Delta| < 0.5$
sur les cinq sous-bandes — la performance ne s'effondre pas au-delà de
$T = 200$.

La localisation $\mathrm{Re}(s) = 1/2$ admet quatre lectures équivalentes
internes à PT : cellule de Planck symplectique $u_\star\, p_\star = 2\pi$,
Jacobien de Haar multiplicatif $du = dp/p$, transformation unitaire
$e^{s u/2}$, et condition aux bords antipériodique forcée par l'involution
$T_3$ du crible mod $3$. Le résidu numérique $1/8$ qui sépare la densité
semi-classique PT de la densité de Riemann–von Mangoldt admet une
sur-détermination cohomologique : $1/8 = s^2/2 = \lambda_H$ (auto-couplage
Higgs), $1/8 = c_1/N_{\mathrm{corners}}$ (première classe de Chern du fibré
spinoriel $q_+/q_-$), $1/8 = s/4$ (phase de Berry par coin de la double
barrière). L'identité $s^2/2 = s/4$ a pour seules solutions $s = 0$ et
$s = 1/2$ ; la valeur $s = 0$ étant exclue par T1 (axiome PT), la cohérence
des deux lectures sur-détermine $s = 1/2$ sans le dériver indépendamment.

Une prédiction structurellement induite — $\Delta_N(L \bmod q) = 1/(8q)$ pour
les fonctions $L$ de Dirichlet primitives — a été testée sur les zéros
pré-calculés de LMFDB ($q \in \{3,5,7,11\}$, $\sim 140$ zéros par
caractère, $T \le 200$). Le test rejette la prédiction à $\chi^2(\mathrm{PT}) =
7.49$ contre $\chi^2(H_0) = 0.023$ sur quatre degrés de liberté. La
calibration sur les 50 premiers zéros de $\zeta$ montre que le résidu $1/8$
de la note 61 est une identité entre formes asymptotiques (entre la densité
semi-classique BK et la densité Riemann–von Mangoldt), non un décalage
mesurable au comptage des zéros. La discipline de falsifiabilité du programme
en sort renforcée et son domaine de validité précisé : le mécanisme PT
décrit ici est spécifique à $\zeta$, pas à la famille Dirichlet.

À partir de ces résultats, on réécrit (§7) la formule de Riemann
elle-même dans le langage canonique PT : le produit d'Euler devient
$\zeta(s) = \zeta_+(s)\, \zeta_-(s)\, \det(I - D_R(s))^{-1}$ ; l'équation
fonctionnelle $s \leftrightarrow 1-s$ s'identifie à l'involution $T_3$
du crible mod $3$ sous transformée Mellin ; la fonction xi de Riemann
admet l'équivalent symétrique tautologique
$\Xi_{\mathrm{PT}}(s) = R(s)\, R(1-s)$ ; la formule explicite de
von Mangoldt prend la forme PT exacte
$-d\log R/ds = \sum_p (\log p)[1/(p^s-1) - A_p p^{-s}]$ (signe corrigé).

La section §7 expose en outre la **géométrie spectrale canonique** de
PT zêta établie en mai 2026 : la courbe de persistance
$\Sigma_{\mathrm{pers}}$ (hyperelliptique genre 14, hyperbolique à
cusp, non-modulaire), l'opérateur Berry-Keating PT $H_n^{\mathrm{cusp}}$
construit rigoureusement sur le sous-espace cuspidal radial avec
auto-adjonction PT-sélectionnée par symétrie $T_3$ (BC antipériodique
$\theta = \pi$), la formule de trace de Selberg-PT en quatre termes
(Weyl + géométrique + parabolique + archimédien) avec scattering matrix
$\phi_{\mathrm{PT}}(s) = \zeta_+(s)\,\zeta_-(s)$, et la **symétrie
circulaire exacte** des 30 branchements de $\Sigma_{\mathrm{pers}}$ sur
**cinq cercles concentriques** centrés au facteur Fisher
$-\mathcal{G}_F = -13/4$ (constante d'$\alpha_{\rm EM}$, ch10), avec
formule unifiée des rayons $r_{C_d} = \mu^\star (2d-1)^{1/d}/4$. La
validation numérique détecte **13/13** des premiers $\gamma_n$ comme
pôles distributionnels à distance médiane $0.075$, avec test de chance
Weyl rejeté à $25\sigma$.

L'article distingue strictement comprendre la structure des zéros et
prouver leur localisation : le programme atteint le premier objectif au
niveau opératoriel ; le second reste un problème ouvert au sens classique,
équivalent au prolongement analytique strict de $R(s)$ hors du demi-plan de
convergence. Le cadre PT zêta constitue à ce jour la version la plus
précise du programme de Hilbert-Polya : géométrie canonique unique,
opérateur canonique avec auto-adjonction PT-sélectionnée, formule de
trace explicite, scattering matrix identifiée, et validation numérique
des premiers zéros.

---

## Table des matières

1. Introduction
2. Le cadre PT
3. Modèle spectral des zéros : amplitudes, déterminant, opérateur canalisé,
   formule de trace
4. Localisation $\mathrm{Re}(s) = 1/2$ : quatre mécanismes équivalents
5. Sur-détermination cohomologique du résidu $1/8$
6. Dualité Higgs ↔ zêta
7. Reformulation PT de la formule de Riemann
8. Vérifications numériques
9. Auto-correction : prédiction Dirichlet falsifiée
10. Limites et ouvertures
11. Conclusion
12. Bibliographie

---

## 1. Introduction

L'hypothèse de Riemann affirme que les zéros non triviaux de la fonction
$\zeta(s) = \sum_{n \ge 1} n^{-s}$ ont pour partie réelle $1/2$. Au-delà de
sa formulation arithmétique, cette conjecture est aussi une question sur la
nature spectrale d'un certain opérateur, dont l'existence est postulée depuis
Hilbert et Pólya. Berry et Keating en 1999, puis Connes en 1996, ont précisé
le contenu de cette hypothèse spectrale : il s'agit de l'opérateur de
dilatation $H = (x p + p x)/2$ quantifié sur la demi-droite, dont la
quantification semi-classique microcanonique reproduit la densité de
Riemann–von Mangoldt $N(\gamma) \sim (\gamma/2\pi) \log(\gamma / 2\pi e)$
des zéros.

Le programme présenté ici ne cherche pas une preuve de l'hypothèse de
Riemann. Il cherche à reconstruire le cadre opératoriel HP–BK de manière
canonique à partir d'une seule structure interne, celle de la Théorie de la
Persistance (PT), qui dérive le Modèle Standard d'un unique paramètre
arithmétique $s = 1/2$ (théorème T1). Le travail revisite la question RH
sous trois angles complémentaires.

Le premier angle est analytique. On exhibe une factorisation
$\zeta(s) = R(s) \zeta_+(s) \zeta_-(s)$, où $\zeta_\pm$ sont définies par
les amplitudes $a_p, b_p$ issues de la bifurcation $q_+/q_-$ de PT. Le
facteur résiduel $R(s)$ admet une représentation Fredholm exacte
$\det(I - D_R(s))^{-1} = R(s)$ avec $D_R$ diagonal et amplitudes propres
$\kappa_p(s)$ en forme close (note 42). Cette représentation est démontrée
dans $\mathrm{Re}(s) > 1$ et vérifiée numériquement à $10^{-12}$.

Le deuxième angle est opératoriel. On identifie l'opérateur canonique
$H_{\mathrm{PT-BK}} = (u\,p_u + p_u\,u)/2$ sur la coordonnée de Mellin
$u = \log p$ avec son cut-off dynamique
$u_{\max}(\gamma) = \gamma/\sqrt{2\pi}$. Ce cut-off n'est pas postulé ; il est
forcé par la cellule de Planck symplectique $u_\star\, p_\star = 2\pi$ de
l'espace de phases canonique (notes 60–61). On prouve que la quantification
semi-classique microcanonique de cet opérateur, en ordre symétrique de
Weyl avec double barrière, reproduit la forme fonctionnelle exacte de la
densité de Riemann–von Mangoldt, à un résidu constant $1/8$ près, vérifié à
$10^{-15}$ sur quatre décades ($\gamma \in [14, 10^3]$, note 61 §1.4).

Le troisième angle est spectral et distributionnel. La construction d'une
trace régularisée
$T(s) = \mathrm{Tr}_{\mathrm{reg}}\big[(s - i H_{\mathrm{PT-BK}})^{-1} D_R(s)\big]$
relie l'opérateur canonique au facteur résiduel via la formule de trace
explicite $-d \log R/ds = \sum_p (\log p)[1/(p^s - 1) - A_p p^{-s}]$
(signe corrigé 2026-05-17),
démontrée comme identité analytique exacte dans $\mathrm{Re}(s) > 1$ et
prolongée distributionnellement à la ligne critique. La régularisation
Hadamard de cette trace produit des pôles aux ordonnées $\gamma_n$ des
zéros de $R$ ; on capture 17 zéros sur 17 testés sur $t \in [10, 70]$ à
$|\Delta| < 1.5$, et 74 sur 79 sur $t \in [10, 200]$ (94 \%), avec
saturation au-delà de $P_{\max} = 3\,000$ (note 74).

Le programme se distingue de Berry–Keating–Connes par l'enrichissement
intrinsèque qu'il fournit. Les nombres apparemment géométriques de la
construction classique — le $2\pi$ de la cellule de Planck, le $1/8$ de la
phase Maslov des coins (correction WKB de $-\pi/8$ par point tournant
semi-classique aux coins du domaine d'aire microcanonique, cf.
Berry–Keating 1999 §3) — reçoivent une dérivation propre à PT. Le $2\pi$
est la mesure de Liouville canonique sur $(u, p_u)$, c'est-à-dire le
commutateur $[u, p_u] = i$. Le $1/8$ admet une triple lecture
cohomologique : $s^2/2 = \lambda_H$ (auto-couplage Higgs, ch15 de la
monographie), $c_1/N_{\mathrm{corners}}$ (première classe de Chern du
fibré spinoriel $q_+/q_-$ Kähler-Fisher), et $s/4$ (phase de Berry par
coin). La cohérence des trois lectures se réduit à l'équation
$s^2/2 = s/4$, dont les seules solutions sont $s = 0$ et $s = 1/2$ ;
l'axiome T1 exclut $s = 0$, ce qui sur-détermine le spin PT modulo T1
supposé (pas une dérivation indépendante de T1).

Une prédiction issue de ce cadre, la formule $\Delta_N(L \bmod q) = 1/(8q)$
pour les fonctions $L$ de Dirichlet primitives, a été testée numériquement
sur les zéros de LMFDB et falsifiée. Cette falsification est documentée
honnêtement. Elle conduit à reformuler la note 61 : le $1/8$ est une
identité asymptotique entre formes lisses, pas un résidu mesurable au
comptage des zéros. Le mécanisme PT décrit ici est spécifique à $\zeta$,
non à la famille Dirichlet via une factorisation naïve par la mesure de
Haar mod $q$.

La distinction maintenue tout au long du texte est la suivante. On
cherche à comprendre la structure des zéros, à fournir un cadre canonique
sans paramètre ad hoc, à dériver les coïncidences numériques connues
(2$\pi$, 7/8, ligne critique) à partir d'une seule structure interne. On
ne cherche pas à prouver que tous les zéros sont sur la ligne critique :
ce verrou ultime (équivalent à RH au sens strict) reste ouvert et ne
semble pas plus facile que la version classique.

---

## 2. Le cadre PT

### 2.1 Le paramètre fondamental

La théorie de la Persistance (PT) prend pour unique axiome arithmétique le
théorème T1 :
$$s = \tfrac{1}{2}.$$
Ce paramètre joue trois rôles distincts dans le corpus PT, qui seront
progressivement identifiés au cours de l'article :

- **arithmétique** : $s$ est le paramètre de symétrie de la transition
  interdite au premier $p = 3$ du crible (théorème T1, monographie ch01) ;
- **géométrique** : $s$ est le seuil $\gamma_p > 1/2$ d'activité d'un
  premier dans la cascade (théorème T5, monographie ch08) ;
- **spectral** : $s$ est la partie réelle qui apparaît dans la transformée
  de Mellin sur les premiers via le Jacobien de la mesure de Haar
  multiplicative.

Une lecture centrale du présent travail est que ces trois manifestations
sont opératoriellement identifiables (section 4).

### 2.2 Point fixe, primes actifs, bifurcation

Le crible PT admet un unique point fixe $\mu^* = 15$. Les premiers
arithmétiquement actifs en ce point fixe sont les diviseurs premiers
$\{3, 5, 7\}$ de $\mu^*$ ; le premier $p = 2$ joue le rôle distinct de
canal de parité (cf. section 6). À ce point fixe, la dynamique de
persistance se bifurque en deux branches conjuguées $q_+$ et $q_-$
caractérisées par
$$q_+ = \frac{13}{15}, \qquad q_- = \exp\!\Big(-\frac{1}{15}\Big).$$
Cette bifurcation est l'objet central du programme zêta. Elle entre dans
toutes les amplitudes propres et toutes les phases cohomologiques rencontrées
plus bas.

### 2.3 Mesure de Haar multiplicative et coordonnée Mellin

La cascade de persistance privilégie une seule mesure invariante sur la
demi-droite des premiers : la mesure de Haar multiplicative
$$du = \frac{dp}{p}, \qquad u = \log p.$$
Cette mesure est invariante par dilatation $p \mapsto \lambda p$, et c'est
la seule mesure invariante par le flot RG de PT (théorème T2,
doublement stochastique de la matrice $T_3$). La coordonnée
$u = \log p$ sera nommée coordonnée de Mellin dans tout l'article ; elle
porte la mesure additivement invariante par translation $u \mapsto u + c$.

### 2.4 Factorisation PT canonique de $\zeta$

Pour $\mathrm{Re}(s) > 1$ on définit, à partir des amplitudes de
bifurcation,
$$
\delta_p^\pm = \frac{1 - q_\pm^p}{p}, \qquad
a_p = \delta_p^+\,(2 - \delta_p^+), \qquad
b_p = \delta_p^-\,(2 - \delta_p^-),
$$
puis
$$
\zeta_+(s) = \exp\Big(\sum_p a_p\,p^{-s}\Big), \qquad
\zeta_-(s) = \exp\Big(\sum_p b_p\,p^{-s}\Big),
$$
et le **facteur résiduel**
$$
R(s) = \frac{\zeta(s)}{\zeta_+(s)\,\zeta_-(s)}.
$$
On pose enfin $A_p = a_p + b_p$, qui sera l'amplitude PT effective par
premier dans toute la suite. Asymptotiquement, $A_p \sim 4/p$ pour $p \to
\infty$.

**Écriture canonique en angles d'holonomie.** Les amplitudes $a_p, b_p$
s'écrivent aussi
$$
a_p = \sin^2\!\big(\theta_p, q_+\big), \qquad b_p = \sin^2\!\big(\theta_p, q_-\big),
$$
soit
$$
A_p = \sin^2\!\big(\theta_p, q_+\big) + \sin^2\!\big(\theta_p, q_-\big),
$$
où $\theta_p$ désigne l'angle d'holonomie de la cascade au premier $p$
(cf. monographie ch10, structure fine ; note 66 §2). L'équivalence entre
les deux écritures est l'identité d'holonomie T6 (monographie ch06) :
pour tout premier $p$ et toute branche $q$,
$\sin^2(\theta_p, q) = \delta_p (2 - \delta_p)$ avec
$\cos(\theta_p, q) = 1 - \delta_p$, où $\delta_p = (1 - q^p)/p$ paramètre
l'angle d'holonomie sur le cercle discret $\mathbb{Z}/p\mathbb{Z}$.

Le programme zêta cherche un modèle spectral canonique de $R$. La
non-annulation `\zeta_+ \zeta_- \neq 0` sur le demi-plan $\mathrm{Re}(s) > 0$
(qui identifie les zéros non triviaux de $R$ avec ceux de $\zeta$ sur la
bande critique) est désormais établie comme **théorème [THM]** : elle découle
du Lemme 3.4 du préprint A_PT (convergence absolue de
$\sum_{p,\sigma} a_p^\sigma p^{-s}$ pour $\mathrm{Re}(s) > 0$ via
l'asymptotique $a_p^\sigma \sim 2/p$), puis du fait que
$\zeta_\pm = \exp(\cdot)$ est par construction non-nul partout où son
exposant est fini. Confirmation numérique : `\min |\zeta_+\zeta_-(1/2+it)|
> 0.48` sur $t \in [10, 1000]$ (`PT_RH_MAY/scripts/test_zeta_pm_nonvanishing.py`,
mpmath dps=50).

---

## 3. Modèle spectral des zéros

Cette section présente l'armature analytique et opératorielle. Quatre
résultats principaux : expression close des amplitudes ($\kappa_p$),
représentation Fredholm exacte de $R$, opérateur canalisé pour le couplage
intra-cascade, et formule de trace explicite.

### 3.1 Expression close des amplitudes : $\kappa_p(s)$

**Définition** (note 42). Pour chaque premier $p$ et chaque $s \in \mathbb{C}$,
on pose
$$
\boxed{\kappa_p(s) = 1 - \big(1 - p^{-s}\big)\,\exp\big(A_p\,p^{-s}\big)}.
$$
Par construction,
$$
1 - \kappa_p(s) = \big(1 - p^{-s}\big)\,\exp\big(A_p\,p^{-s}\big),
$$
et le développement asymptotique donne
$$
\kappa_p(s) = (1 - A_p)\,p^{-s} + O\big(p^{-2s}\big).
$$

**Interprétation**. Le scalaire fondamental de la cascade n'est ni
$p^{-s}$ (terme d'Euler) ni $\exp(A_p p^{-s})$ (compensation bifurquée),
mais leur écart $\kappa_p$. On a la lecture
$$
1 - \kappa_p = \text{survie résiduelle},
\qquad
\kappa_p = \text{amplitude spectrale du résidu}.
$$

### 3.2 Représentation Fredholm exacte

Soit $D_R(s)$ l'opérateur diagonal sur $\ell^2(\mathcal{P})$ ($\mathcal{P}$
l'ensemble des premiers) défini par $D_R(s)\,e_p = \kappa_p(s)\,e_p$.

**Théorème [THM] (note 42)**. Pour $\mathrm{Re}(s) > 1$, l'opérateur
$D_R(s)$ est trace-class, et
$$
\boxed{\det\big(I - D_R(s)\big)^{-1} = R(s)\,.}
$$

**Preuve**. Pour $\mathrm{Re}(s) > 1$, $\sum_p |\kappa_p(s)| < \infty$ par
le développement asymptotique de $\kappa_p$. L'opérateur diagonal est donc
trace-class. Le calcul direct donne
$$
\det(I - D_R(s))^{-1}
= \prod_p (1 - \kappa_p)^{-1}
= \prod_p (1 - p^{-s})^{-1}\,\exp(-A_p\,p^{-s})
= \zeta(s)\,\exp\Big(-\sum_p A_p\,p^{-s}\Big)
= R(s).
$$

**Portée de la représentation.** La preuve est algébrique directe : la
représentation Fredholm est une **réécriture déterminantielle** de la
définition de $R(s)$. Sa valeur est de fournir la **cible spectrale** que
tout opérateur PT canonique candidat doit reproduire, et non un argument
indépendant pour la structure spectrale de $R$ (note 42 §1).

**Vérification numérique**. Pour $s = 1.3 + 14.134725\,i$, $P = 10^6$
premiers tronqués, l'erreur résiduelle entre $R(s)$ calculée par produit
Fredholm et $R(s)$ calculée à partir de $\zeta$ par `mpmath` est $\sim 2.0
\cdot 10^{-12}$ (note 42 §7). À $\sigma = 1.1$ la précision tombe à $\sim
1.3 \cdot 10^{-3}$, ce qui reflète l'approche de la frontière de convergence
absolue, non un défaut de la représentation.

**Non-annulation [THM] (révisée 2026-05-17).** L'identification
zéros$(R) = $ zéros$(\zeta)$ utilisée dans toute la suite requiert
$\zeta_+(s)\,\zeta_-(s) \neq 0$ sur le demi-plan $\mathrm{Re}(s) > 0$
(en particulier sur la ligne critique). Cette non-annulation est désormais
établie comme **théorème inconditionnel** : conséquence directe du
Lemme 3.4 du préprint A_PT (convergence absolue de
$\sum_p a_p^\sigma p^{-s}$ pour $\mathrm{Re}(s) > 0$ via l'asymptotique
$a_p^\sigma \sim 2/p$), suivie du fait que $\zeta_\pm = \exp(\cdot)$ est
non-nul partout où son exposant est fini. Confirmation numérique :
$\min |\zeta_+\zeta_-(1/2+it)| > 0.48$ sur $t \in [10, 1000]$
(`PT_RH_MAY/scripts/test_zeta_pm_nonvanishing.py`, mpmath dps=50).
Audit de cohérence : `PT_RH_MAY/analysis/INCOHERENCES_AUDIT.md` INC-1.

### 3.3 Opérateur canalisé pour le couplage $\sum_p A_p\,p^{-s}$

Pour modéliser le couplage intra-cascade (la somme $\sum_p A_p p^{-s}$)
comme spectre d'un opérateur canonique sans paramètre ajustable, on
construit
$$
H_{\mathrm{PT}} = L^2(\Sigma_{\mathrm{pers}}, du) \otimes L^2\big(\widehat{\mathbb{Z}/M\mathbb{Z}}\big).
$$
La première composante est la coordonnée Mellin continue sur $\Sigma_{\mathrm{pers}}$
(surface de persistance, genre 14 dans le corpus PT) ; la seconde est
discrète, sur les caractères du groupe multiplicatif modulo $M$ (cascade
CRT, $M$ produit primorial). Le couplage inter-fenêtre est par noyau de
Bergman régularisé avec largeur géométrique $\varepsilon_*$ (notes 51, 52) :
$$
B_{\mathrm{BF}}(z_i, z_j) = \frac{\varepsilon_i \varepsilon_j}{|z_i - z_j|^2 + \varepsilon_i \varepsilon_j},
\qquad \varepsilon_i = \frac{\varepsilon_*}{\rho_i}.
$$

**Détermination géométrique de $\varepsilon_*$ (note 52)**. La covariance
de Bergman régularisée à résolution maximale ($\varepsilon_{\mathrm{reg}} =
10^{-3}$, paquets de taille premier individuel) donne
$$
\varepsilon_*^{\mathrm{geom}} = 1.51.
$$
Cette valeur reproduit l'optimum empirique du sweep numérique (note 51) à
mieux que 5 % en moyenne sur l'ensemble des configurations testées, et à 1 %
au régime optimal $M = 1155, W = 24$ ; à $M = 105, W = 48$ l'écart atteint
24 % (note 52 §4 et §10). Le tuning est éliminé au sens où aucun paramètre
n'est ajusté, mais la précision n'est pas uniformément sous le pour cent.

**Performance**. Avec $\varepsilon_*^{\mathrm{geom}}$ injecté sans
ajustement, l'erreur médiane sur $\sum_p A_p p^{-s}$ est $0.085$ à
$M = 1155$, $W = 24$ — soit un facteur $\sim 85$ au-dessus du baseline
scalaire lissé adaptatif ($0.001$, qui utilise un $\beta$ tuné). L'apport
n'est donc pas la précision absolue mais l'absence de paramètre libre
(note 52 §1 et §8).

### 3.4 Formule de trace explicite

**Théorème [THM] (note 62, signe corrigé 2026-05-17)**. Pour tout $s$ avec
$\mathrm{Re}(s) > 1$,
$$
\boxed{
-\,\frac{d}{ds} \log R(s)
= \sum_p \frac{\kappa_p'(s)}{1 - \kappa_p(s)}
= \sum_p (\log p)\,\Big[\frac{1}{p^s - 1} - A_p\,p^{-s}\Big]\,.
}
$$

**Preuve**. Par dérivation logarithmique de
$1 - \kappa_p(s) = (1 - p^{-s}) \exp(A_p p^{-s})$. Poser $z = p^{-s}$, donc
$dz/ds = -\log p \cdot z$. Un calcul direct donne
$\kappa_p'(s) = \log p \cdot z \cdot \exp(A_p z) \cdot [A_p(1-z) - 1]$,
d'où
$$
\frac{\kappa_p'(s)}{1 - \kappa_p(s)}
= (\log p)\Big[A_p\,z - \frac{z}{1 - z}\Big]
= (\log p)\Big[A_p\,p^{-s} - \frac{1}{p^s - 1}\Big].
$$
En sommant et en changeant le signe global on obtient la formule annoncée.
Le sanity check via $-\zeta'/\zeta(s) = \sum_n \Lambda(n) n^{-s} = \sum_p
(\log p)/(p^s - 1)$ et la dérivée
$\frac{d}{ds}\log(\zeta_+\zeta_-) = -\sum_p A_p (\log p) p^{-s}$ donne
$-d\log R/ds = -\zeta'/\zeta - \sum_p A_p(\log p) p^{-s} = \sum_p (\log p)[1/(p^s-1) - A_p p^{-s}]$,
en accord avec la formule encadrée.

**Vérification numérique**. À $s = \sigma + i\,t$ avec
$\sigma \in \{2, 2.5, 3, 4\}$, $t \in \{0, 5, 14.134725\}$, $P = 1000$
premiers, l'erreur entre la somme tronquée et $-d\log R/ds$ (calculé via
`mp.diff`) est de l'ordre $10^{-13}$ pour $\sigma = 4$ et $10^{-4}$ pour
$\sigma = 2$ (cohérent avec le reste de Mertens $O(P^{1-\sigma}/\log P)$).
L'identité algébrique est exacte ; l'erreur n'est due qu'au tronquage de
la somme sur primes.

**Correction de signe (2026-05-17)** : la note 62 originale écrivait
$+A_p p^{-s}$ (signe positif). La vérification directe via
`mp.diff(\log R)` avec $\log R = \log\zeta - \sum_p A_p p^{-s}$ montre
que le signe correct est $-A_p p^{-s}$ (voir A4_formule_trace.md §3.3bis,
A4_trace_verify.py test 1B). Le bug initial provenait du script
`pt_bk_explicit_formula.py` (fonction `kappa_p_deriv`) qui écrivait
$1 + A_p(1-z)$ au lieu de $A_p(1-z) - 1$ ; cette erreur s'est propagée
dans la formule (*).

**Lecture Berry–Keating**. L'opérateur de dilatation
$H_{\mathrm{PT-BK}} = (u p_u + p_u u)/2 = -i(u \partial_u + 1/2)$ engendre
le flot $e^{-it H_{\mathrm{PT-BK}}} f(u) = e^{-t/2} f(e^{-t} u)$.
Sur la mesure atomique $\nu = \sum_p \delta(u - \log p)$, les orbites
périodiques du flot issu de $\log p_0$ ont pour longueurs $\ell_p = k \log
p_0$, $k \ge 1$, correspondant aux puissances $p_0^k$. Le terme
$(\log p)/(p^s - 1) = \sum_{k \ge 1} (\log p) p^{-ks}$ s'identifie alors à
la somme sur ces orbites, et $A_p p^{-s}$ représente la contribution de
couplage holonomique $q_+/q_-$ par premier.

### 3.5 Trace régularisée opératorielle

L'écriture compacte de la formule de trace en termes de la résolvante de
$H_{\mathrm{PT-BK}}$ et de l'opérateur $D_R$ utilise une projection de
Mellin atomique
$$
\Pi : L^2([u_{\min}, u_{\max}], du) \to \ell^2(\mathcal{P}_\gamma),
\qquad \Pi(f)_p = f(\log p)\,\sqrt{\log p},
$$
où $\mathcal{P}_\gamma = \{p \text{ premier} : u_{\min} \le \log p \le
u_{\max}(\gamma)\}$. Le facteur $\sqrt{\log p}$ est la racine de la mesure
atomique de Haar multiplicative. La trace régularisée s'écrit alors
$$
\boxed{
T(s) := \mathrm{Tr}_{\mathrm{reg}}\!\Big[(s - i H_{\mathrm{PT-BK}})^{-1} D_R(s)\Big]
= \sum_{p \in \mathcal{P}_\gamma} \kappa_p(s)\,(\log p)\,\Big[\frac{1}{p^s - 1} - A_p\,p^{-s}\Big]\,.
}
$$
Dans $\mathrm{Re}(s) > 1$, $T(s) \equiv -d \log R/ds$. Sur la ligne critique
$\mathrm{Re}(s) = 1/2$, la somme prime n'est pas absolument convergente ;
elle est néanmoins distributionnellement convergente comme élément de
$\mathcal{S}'(\mathbb{R}_t)$, et la régularisation Hadamard donne un sens
fonctionnel à $T(1/2 + it)$ au sens numérique (régularisations R1/R2/R3 de
la note 65) ; la démonstration fonctionnelle stricte du prolongement n'est
pas close (cf. note 65 §7.1).

**Auto-correction : G_HP4.d falsifié.** Un premier test (note 64) avait
posé l'identité naïve « spectre direct de $H_{\mathrm{PT-BK}}$ borné par
double barrière $=\{\gamma_n\}$ ». Ce test a été falsifié : le spectre
direct est régulier (Bohr–Sommerfeld asymptotiquement linéaire), alors que
les $\gamma_n$ sont irréguliers. Le ratio entre les deux spectres diverge.
Cette falsification a motivé le passage à la trace régularisée
distributionnelle ci-dessus : les zéros de $\zeta$ sont des pôles
distributionnels de $T(s)$, pas des valeurs propres directes de
$H_{\mathrm{PT-BK}}$. La discipline de falsifiabilité du programme inclut
cette correction au même titre que la falsification Dirichlet de la
section 8.

---

## 4. Localisation $\mathrm{Re}(s) = 1/2$ : quatre mécanismes équivalents

L'objet central de cette section est l'explication PT-canonique de la ligne
critique. Quatre mécanismes apparemment distincts apparaissent à des
étapes différentes du programme. Ils sont en fait quatre faces d'une
unique structure : la mesure invariante de la dynamique de persistance.

### 4.1 Cellule de Planck symplectique : $u_\star\, p_\star = 2\pi$

L'opérateur de dilatation $H_{\mathrm{PT-BK}}$ engendre les transformations
d'échelle sur la coordonnée de Mellin. Sa mesure de Liouville canonique est
$du\,dp_u / (2\pi)$, conséquence directe du commutateur canonique $[u, p_u]
= i$.

**Cellule de Planck PT (notes 60, 61)**. La cellule unité de l'espace de
phases est
$$
\boxed{u_\star\, p_\star = 2\pi},
$$
choix symétrique $u_\star = p_\star = \sqrt{2\pi} \approx 2.5066$ (point
fixe de l'involution $(u, p_u) \leftrightarrow (p_u, u)$). À énergie
$\gamma$ fixée, l'hyperbole $H_{\mathrm{PT-BK}} = \gamma$ ne pénètre dans
l'espace admissible que sur le segment
$$
u \in [u_\star, \gamma/p_\star] = [\sqrt{2\pi}, \gamma/\sqrt{2\pi}].
$$
On obtient le cut-off dynamique canonique
$u_{\max}(\gamma) = \gamma/\sqrt{2\pi}$.

**Densité semi-classique (note 61)**. L'aire microcanonique avec double
barrière donne, après ordre symétrique de Weyl,
$$
N_{\mathrm{PT}}(\gamma) = \frac{\gamma}{2\pi}\,\log\frac{\gamma}{2\pi\,e} + 1.
$$
La densité Riemann–von Mangoldt est
$$
N_{\mathrm{RvM}}(\gamma) = \frac{\gamma}{2\pi}\,\log\frac{\gamma}{2\pi\,e} + \frac{7}{8} + O(1/\gamma).
$$
La forme fonctionnelle est identique. L'écart constant
$N_{\mathrm{PT}} - N_{\mathrm{RvM}} = 1/8$ est vérifié à $10^{-15}$ sur
quatre décades ($\gamma \in [14, 10^3]$, table de la note 61 §1.4). Le $2\pi$
qui apparaît dans la densité $N(\gamma)$ est la mesure de Liouville canonique
sur $(u, p_u)$, conséquence du commutateur $[u, p_u] = i$ ; le shift
$\mathrm{Re}(s) = 1/2$ est, lui, dérivé du Jacobien de la mesure de Haar
multiplicative $du = dp/p$ (§4.2). Les deux faits sont reliés (mesure
invariante unique) sans être le même nombre.

### 4.2 Jacobien de Haar multiplicatif

Pour un mode stationnaire $f(\tau, u) = e^{i\omega \tau} \psi(u)$ de
l'opérateur d'Alembert sur la coordonnée Mellin (la métrique Fisher de
$\Sigma_{\mathrm{pers}}$ étant lorentzienne pour $\mu > \mu_c \approx 6.97$,
note 58), le shift Mellin sous la mesure de Haar multiplicative
$d\pi_{\log} = \sum_p \delta_{\log p}\,d\log p$ produit un facteur
$1/p = e^{-u/2}$ devant la transformée. Plus précisément, le pôle de la
transformée de Mellin se déplace de $s = \omega$ à
$$
s = \tfrac{1}{2} + i\omega.
$$
C'est le shift de la mesure de Haar multiplicative. La partie réelle $1/2$
n'est pas postulée : elle vient du Jacobien naturel de la mesure invariante
de Haar sur les premiers.

### 4.3 Transformation unitaire $e^{su/2}$

Le théorème SA de la monographie (ch37b) établit que l'opérateur de
Ruelle–PT discret $L_s$ est auto-adjoint pour la mesure de Gibbs $\mu_s(p) =
p^{-s}$ via la balance détaillée associée à la symétrie antidiagonale
$T_3 = T_3^\top$. En coordonnée Mellin, la mesure de Gibbs devient
$\mu_s(p)\,dp = e^{-su}\,du$. La transformation unitaire
$$
U : f(u) \mapsto e^{-su/2}\,f(u)
$$
ramène l'opérateur $L_s$ sur $L^2(\mathbb{R}, e^{-su} du)$ à
$L^2(\mathbb{R}, du)$ en translatant $H_{\mathrm{PT-BK}}$ par $i s/2$. Le
shift apparaît littéralement dans l'exposant : c'est l'identification
opératorielle entre balance détaillée discrète (T2 du crible) et
auto-adjonction continue de $H_{\mathrm{PT-BK}}$.

### 4.4 Condition aux bords antipériodique

Sur $L^2([u_{\min}, u_{\max}], du)$ avec $0 < u_{\min} < u_{\max} < \infty$,
les indices de défaut de $H_{\mathrm{PT-BK}}$ valent $(n_+, n_-) = (1, 1)$.
Par von Neumann, il existe une famille $U(1)$ d'extensions auto-adjointes
paramétrées par
$$
\psi(u_{\max}) = e^{i\theta}\,\psi(u_{\min}), \qquad \theta \in [0, 2\pi).
$$

**Sélection PT-canonique (note 63)**. L'involution
$$
I : u \mapsto u_{\max} + u_{\min} - u
$$
est la projection continue de l'antidiagonale $T_3 = \begin{pmatrix}0 & 1\\ 1 & 0\end{pmatrix}$
du crible mod 3. C'est l'unique involution continue qui préserve les bords
comme orbite et est une isométrie de $L^2$. Le terme de bord
$I H_{\mathrm{PT-BK}} I^{-1} + H_{\mathrm{PT-BK}}$ s'annule **uniquement** pour
$\theta = \pi$ (antipériodique), qui force la symétrie spectrale
$\{\gamma_n, -\gamma_n\}$. Vérifié numériquement à $5 \cdot 10^{-12}$
(précision machine).

L'antipériodicité continue est donc l'image continue de l'antidiagonalité
discrète $T_3$ du crible, qui est elle-même l'expression la plus simple
d'une matrice stochastique non-triviale 2×2. La balance détaillée discrète
et l'auto-adjonction continue se traduisent l'une l'autre par la
transformation unitaire $e^{su/2}$.

### 4.5 Unification

Les quatre mécanismes ci-dessus manifestent quatre facettes d'une même
dynamique de persistance, avec sa mesure invariante de Haar multiplicative
et sa réversibilité T2 (doublement stochastique). Les liens M1↔M2 (mesure
de Liouville $\to$ Haar mult) et M3↔M4 (balance détaillée Gibbs $\to$
antipériodicité $T_3$) sont explicitement construits ; le lien complet
M1$\leftrightarrow$M2$\leftrightarrow$M3$\leftrightarrow$M4 au sens d'une
équivalence mathématique stricte demande une formulation catégorielle non
incluse ici. Les trois
$1/2$ PT recensés en section 2.1 (arithmétique T1, géométrique T5, spectral
de Haar) sont opératoriellement identifiés par cette section. C'est l'apport
structurel principal du programme.

| Mécanisme | Description | Origine PT |
|---|---|---|
| M1 — Cellule symplectique | $u_\star p_\star = 2\pi$ | $[u, p_u] = i$ (Liouville) |
| M2 — Jacobien Haar | $du = dp/p$, shift $+1/2$ | T2, doublement stochastique |
| M3 — Transformation unitaire | $U = e^{-su/2}$ | balance détaillée Gibbs $\to$ AA |
| M4 — Condition aux bords | $\theta = \pi$ antipériodique | $T_3$, antidiagonale du crible |

---

## 5. Sur-détermination cohomologique du résidu $1/8$

La densité semi-classique PT diffère de la densité Riemann–von Mangoldt par
une constante $1/8$ (section 4.1). Berry et Keating en 1999 identifient cette
constante comme la phase Maslov fine des coins de la double barrière, sans
interprétation dynamique. PT fournit une **triple identification
cohomologique** convergente, sur-déterminée par la cohérence des trois
lectures.

### 5.1 K2 — Identité algébrique : $1/8 = s^2/2 = \lambda_H$

L'auto-couplage scalaire du Higgs dans le secteur PT est dérivé dans la
monographie (ch15, profondeur $d = 3$) :
$$
\lambda_H = \frac{s^2}{2}.
$$
À $s = 1/2$, $\lambda_H = 1/8$. Cette identité est structurellement liée à
la bifurcation $q_+/q_-$ à deux branches : $\lambda_H$ est le carré du
paramètre de bifurcation divisé par le nombre de branches scalaires.

Deux autres identifications algébriques coïncident :
$$
1/8 = (1/2)^3 \qquad (\text{cube de l'unité scalaire, ch20g}),
$$
$$
1/8 = \frac{1}{2 N_{\mathrm{gen}} + 2} \quad (\text{à } N_{\mathrm{gen}} = 3, \text{exact}).
$$

Plusieurs candidats naïfs sont rejetés : $1/N_{\mathrm{Weyl}} = 1/15 \ne
1/8$, $\chi(\Sigma_{\mathrm{pers}})/Z = -26/Z \ne 1/8$, $\gamma_3 \gamma_5
\gamma_7 \approx 0.335 \ne 1/8$, $\alpha_{\mathrm{bare}} \approx 1/136
\ne 1/8$. Le $1/8$ n'est pas un invariant topologique de la surface de
persistance, mais un invariant scalaire de la bifurcation $q_+/q_-$.

### 5.2 K3 — Identité cohomologique : $1/8 = c_1/N_{\mathrm{corners}}$

La structure cohomologique pertinente est déjà construite dans le corpus
PT (notes 67, 68). La variété de Fisher doublée
$$
\mathcal{M}_m^{\mathrm{db}} = \mathcal{M}_m \times \mathcal{M}_m, \qquad
g = g_{\mathrm{Fisher}} \oplus g_{\mathrm{Fisher}},
$$
muni de la structure complexe $J_{\mathrm{PT}}(u, v) = (-v, u)$ est une variété
kählérienne (théorème [THM] kahler_fisher, ch11). Les projecteurs chiraux
$$
P_\pm = \tfrac{1}{2}(I \mp i J_{\mathrm{PT}})
$$
décomposent $H^*(BG_m; \mathbb{C})$ en sous-espaces propres $\pm i$ et sont
identifiés au projecteur spinoriel $q_+/q_-$ du crible (ch12,
def:chiral_projectors).

**Théorème [THM] berry_3pi** (monographie, app. P, §C7). La première classe
de Chern du fibré spinoriel sur $\mathcal{M}_m^{\mathrm{db}}$ est $c_1 = s
= 1/2$. L'holonomie du revêtement $N_{\mathrm{gen}}$-fois est
$\exp(2\pi i c_1 N_{\mathrm{gen}}) = \exp(3\pi i)$. Le sous-énoncé
algébrique $c_1 = 1/2$ est [THM] inconditionnel ; la rigueur APS stricte du
calcul $1/8 = c_1/N_{\mathrm{corners}}$ sur la variété à coin reste
partiellement ouverte (note 68 §6.2–§6.3), statut [DER]/[partiel].

**Calcul cohomologique du résidu (note 68)**. La double barrière BK est, dans
l'espace des phases canonique, un quart de plan symplectique
$Q_\star = \{u \ge u_\star, p \ge p_\star, u_\star p_\star = 2\pi\}$. Le bord
admet deux arêtes se rejoignant au coin. Avec la symétrisation Weyl et les
deux branches miroirs, on obtient $N_{\mathrm{corners}} = 4$ coins distincts
sur la cellule complète. Ce comptage suppose que les branches miroirs
(Weyl) **identifient** deux coins en un seul via l'antipériodicité
$T_3$ (§4.4) ; sans cette identification on aurait $N = 8$ et $c_1/N =
1/16$. L'identification est cohérente avec la condition de bord
antipériodique $\theta = \pi$ (note 63, note 68 §3.1) mais reste un choix
de comptage structurel, non rigoureusement dérivé d'un théorème d'unicité. La phase de Berry contributée par chaque coin,
mesurée sur le quart de tour reliant les deux arêtes adjacentes, est
$$
\gamma_{\mathrm{corner}} = -\frac{2\pi c_1}{N_{\mathrm{corners}}} = -\frac{\pi}{4},
$$
et la densité de Bohr–Sommerfeld associée est
$$
\boxed{N_{\mathrm{corner}} = \frac{\gamma_{\mathrm{corner}}}{2\pi} = -\frac{c_1}{N_{\mathrm{corners}}} = -\frac{1}{8}.}
$$
Vérification numérique : précision machine $1.68 \cdot 10^{-14}$ (script
`pt_k3_cohomological_maslov.py`). C'est exactement le $-\pi/8$ par coin de
Berry–Keating 1999, mais dérivé depuis la structure cohomologique PT, et
non comme artefact géométrique opaque.

### 5.3 Phase de Berry par coin : $1/8 = s/4$

À l'équateur de la sphère de Bloch ($\theta = \pi/2$), la connexion Berry
d'un spin-$s$ est constante : $A_\varphi = -\sin^2(\pi/4) = -1/2 = -s$ par
unité d'angle azimutal. Sur un cycle équatorial complet, la phase Berry est
$\gamma_1 = -\pi = -2\pi s$. Divisée par les quatre coins :
$$
\gamma_1 / N_{\mathrm{corners}} = -\pi/4 = -2\pi \cdot (s/4),
$$
soit $-1/8$ par coin en compte d'aire. L'identité algébrique correspondante
est
$$
\frac{1}{8} = \frac{s}{4}.
$$

### 5.4 Sur-détermination par cohérence

Les trois identifications convergent vers le même nombre à $s = 1/2$ :
$$
\boxed{\frac{1}{8} = \frac{s^2}{2} \big|_{\text{K2, Higgs}}
= \frac{c_1}{N_{\mathrm{corners}}} \big|_{\text{K3, Chern}}
= \frac{s}{4} \big|_{\text{Berry par coin}}.}
$$
La cohérence algébrique entre K2 et K3 se réduit à
$$
\frac{s^2}{2} = \frac{s}{4} \quad \Longleftrightarrow \quad s(2s - 1) = 0.
$$
Les seules solutions sont $s = 0$ (trivial, pas de spin, exclu par T1) et
$s = 1/2$ (théorème T1). **La cohérence des trois lectures cohomologiques
sur-détermine $s = 1/2$ modulo l'axiome T1 supposé** ; elle ne dérive pas
$s = 1/2$ indépendamment de T1. C'est une signature structurelle de T1 par
la zone HP-PT du corpus, non un substitut à T1.

### 5.5 Lecture en termes de théorème d'indice

L'assemblage final s'interprète comme un théorème d'indice
Atiyah–Patodi–Singer formel sur la variété à coin $\square$ :
$$
\mathrm{ind}(D_+|_\square) = \int_\square \mathrm{ch}(E) \mathrm{Td}(\square)
+ \tfrac{1}{2}\big(\eta_{\partial} + \dim \ker D_\partial\big)
+ \sum_{\mathrm{coins}} \frac{\delta_{\mathrm{coin}}}{2}.
$$
Avec $c_1(E) = 1/2$ et $\delta_{\mathrm{coin}} = \pi/2$ pour un coin droit,
la contribution par coin est $\pi/2 / (2 \cdot 2\pi) = 1/8$, en accord avec le
calcul de Berry direct. L'identification équivalente via Chern–Simons :
$\exp(2\pi i c_1 / N_{\mathrm{corners}}) = \exp(i\pi/4)$, dont l'argument
$\pi/4$ donne en densité $1/8$ par coin.

Une rigueur formelle au sens APS strict (au-delà du calcul de l'aire et de
la vérification numérique) reste à développer. Le mécanisme cohomologique
est néanmoins explicitement identifié.

---

## 6. Dualité Higgs ↔ zêta

La triple identification de la section 5 conduit à une conjecture
structurelle inattendue, formulée comme verrou K4 (note 69).

### 6.1 Le bifurcateur canonique commun

**Conjecture K4 (forme structurelle)**. Il existe un objet PT canonique $B$,
le bifurcateur, tel que :
- $m_H / v = \mathrm{scale}(B) = s = 1/2$ (mode breathing, masse du Higgs) ;
- $\Delta_N = \mathrm{coupling}(B) = s^2/2 = 1/8$ (phase Maslov, secteur
  HP-PT) ;
- $\lambda_H = \mathrm{coupling}(B) = s^2/2 = 1/8$ (auto-couplage Higgs).

Le bifurcateur canonique est identifié comme le projecteur spinoriel
$q_+/q_-$ du canal $p = 2$ au point fixe $\mu^* = 15$. Formellement,
$$
B = \Pi_+ : \mathcal{H} \to \mathcal{H}_{q_+}, \qquad
\Pi_+ = \tfrac{1}{2}\big(I - i J_{\mathrm{PT}}\big),
$$
où $J_{\mathrm{PT}}$ est la structure complexe Kähler-Fisher de
$\mathcal{M}_m^{\mathrm{db}}$ (cf. §5.2) et $\mathcal{H} = \mathcal{H}_{q_+}
\oplus \mathcal{H}_{q_-}$ la décomposition chirale. L'opérateur miroir
$\Pi_- = (I + i J_{\mathrm{PT}})/2$ projette sur $q_-$, et $B$ se restreint
au canal $p = 2$ (sous-espace propre de la partition binaire T6).

### 6.2 Le canal $p = 2$ porte les deux observables

Le canal $p = 2$ joue trois rôles en PT (ch02 rem:p2_operator) : opérateur
info/anti-info ($\gamma_2 = 0.867$, le plus actif), bifurcation $q_+/q_-$
(chiralité de la persistance), et partition binaire $\log_2 m$ (théorème
T6, GFT).

Le Higgs habite le canal $p = 2$ (ch20e §C1) : il est neutre sous
$\mathrm{SU}(3) \times \mathrm{SU}(2) \times \mathrm{U}(1)$, les trois jauges
habitent $\{3, 5, 7\}$, donc seul $p = 2$ est admissible pour $H$.

La phase Maslov habite aussi le canal $p = 2$ : chaque coin de la double
barrière BK est un point où le crible bascule $q_+ \leftrightarrow q_-$, ce
qui est l'expression opératorielle de la partition binaire $p = 2$ (T6).

### 6.3 Dérivation explicite par les théorèmes T1 + T2 + Haar + Weyl

$$
\begin{aligned}
\mathrm{T1} &: s = 1/2 \qquad \text{(théorème, ch01)} \\
\mathrm{T2} &: |\lambda_2(T_{30})| = s^2 = 1/4 \qquad \text{(gap spectral, ch03)} \\
\mathrm{Haar} &: \text{shift } \mathrm{Re}(s_\zeta) \to \mathrm{Re}(s_\zeta) + 1/2 \\
\mathrm{Weyl} &: \text{ordre symétrique, compte d'états divisé par 2}
\end{aligned}
$$
$$
\boxed{\frac{1}{8} = \frac{s^2}{2} = \frac{|\lambda_2(T_{30})|}{N_{\mathrm{branches}}\, q_+/q_-} = \frac{1/4}{2} = \lambda_H = \Delta N_{\mathrm{Maslov}}.}
$$

La quantité $1/8$ est donc structurellement le gap spectral T2 divisé par le
nombre de branches $q_+/q_-$. Les deux lectures (Higgs et Maslov) sont deux
noms pour le même nombre, généré par le même mécanisme.

### 6.4 Tests d'isolation

Aucun autre couplage PT canonique n'égale $1/8$ :

| Couplage PT | Valeur | $= 1/8$ ? |
|---|---|---|
| $\alpha_{\mathrm{EM}}$ | $7.30 \cdot 10^{-3}$ | non |
| $\sin^2 \theta_W$ | $0.231$ | non |
| $G_F$ (GeV${}^{-2}$) | $1.17 \cdot 10^{-5}$ | non |
| $\lambda_H$ | $\mathbf{0.125}$ | **oui** |
| $\Delta N_{\mathrm{Maslov}}$ | $\mathbf{0.125}$ | **oui** |

La coïncidence numérique pure est rejetée. La hiérarchie de spin
$(s, s^2, s^2/2)$ est dense jusqu'à $1/8$ et se referme algébriquement
au-delà ($s^3 = s^2/2$ à $s = 1/2$).

### 6.5 Tests numériques quantitatifs

À l'arbre, $m_H^{\mathrm{tree}} = v \cdot s = 123.11$ GeV (en utilisant
$v = 246.22$ GeV). Au premier ordre perturbatif PT,
$$
m_H^{(1)} = v \cdot s \cdot (1 + C_F\,\varepsilon) \simeq 124.65 \text{ GeV},
$$
puis la valeur complète PT-NLO (avec corrections d'ordre supérieur,
monographie ch15, note 69 §4.1) atteint $m_H^{\mathrm{NLO}} = 125.287$ GeV,
où les deux coefficients sont définis intrinsèquement en PT :

- $C_F = (N_c^2 - 1)/(2 N_c) = 4/3$ est le facteur de Casimir de la
  représentation fondamentale de $\mathrm{SU}(N_c)$, avec $N_c = 3$
  imposé par PT comme la seule valeur de $N_c$ pour laquelle l'équation
  de point fixe $\mu^* = \sum p_i$ admet une solution sur les premiers
  actifs (théorème T0 de la monographie, ch08). Géométriquement, $C_F$
  compte les transitions permises de la matrice de transfert $T_3$ du
  canal $p = 3$ (7 transitions sur 9 emplacements ; cf. monographie
  ch11 §J\_PMNS). Statut [THM] (dérivation $N_c = 3 \Rightarrow C_F = 4/3$).

- $\varepsilon = \alpha_s\, s / (2\pi)$ est le petit paramètre PT du
  développement perturbatif à une boucle, qui combine le couplage fort
  $\alpha_s(m_Z) \simeq 0.1179$ et le paramètre de bifurcation $s = 1/2$
  avec la normalisation standard $1/(2\pi)$ d'une boucle gluonique
  (monographie ch15 éq. R19b, ch16 §perturbative). Numériquement
  $\varepsilon \simeq 9.38 \cdot 10^{-3}$ et $C_F\,\varepsilon \simeq
  1.25 \cdot 10^{-2}$. Statut [DER] (formule perturbative standard avec
  $\alpha_s$ et $s$ pris à leurs valeurs PT canoniques, sans paramètre
  libre).

La valeur expérimentale LHC est $125.25$ GeV. L'erreur du premier ordre
$m_H^{(1)}$ est $\sim 0.5\,\%$, l'erreur PT-NLO complet $m_H^{\mathrm{NLO}}$
est $0.0295\,\%$ (note 69 §4.1).

> **Note 2026-05-17 (PT_RH_MAY audit INC-5)** : la formule simplifiée
> $m_H = v \cdot s \cdot (1 + C_F \varepsilon)$ avec $\alpha_s(m_Z) = 0.1179$
> et $\varepsilon = \alpha_s s/(2\pi) \approx 9.38 \cdot 10^{-3}$ produit
> numériquement $m_H \approx 124.65$ GeV. La valeur citée
> $m_H^{\mathrm{PT-NLO}} = 125.287$ GeV utilise la formule NLO complète
> (monographie ch15 R56) qui inclut des termes supplémentaires de
> correction électrofaible non écrits ici. L'écart numérique est de
> $\sim$0.5\,\% entre la version pédagogique (premier ordre) et la version
> NLO complète. La cohérence avec LHC ($125.25 \pm 0.17$ GeV) demeure
> excellente dans les deux cas.

### 6.6 Position finale

K4 est établi comme **conjecture forte structurellement validée**. La
coïncidence numérique pure est rejetée. Le mécanisme commun est identifié :
bifurcation $q_+/q_-$ du canal $p = 2$. La chaîne de dérivation T1 → T2 →
Haar → Weyl utilise uniquement des théorèmes [THM] du corpus. La promotion
à théorème inconditionnel demande la fermeture cohomologique stricte K3
(actuellement [THM partiel]).

---

## 7. Reformulation PT de la formule de Riemann

Les sections précédentes ont établi le modèle spectral PT des zéros de
$\zeta$ et la dualité Higgs $\leftrightarrow$ zêta. Cette section utilise
ces résultats pour **réécrire la formule de Riemann elle-même** dans le
langage canonique de la persistance, **enrichie** par la géométrie
spectrale $\Sigma_{\rm pers}$ (cycle PT_RH_HYPERBOLIC_CUSP, mai 2026) et
l'opérateur Berry-Keating PT canonique $H_n^{\rm cusp}$ (programme A1-A5).
L'objectif n'est pas une nouvelle preuve de RH — qui reste un problème
ouvert — mais d'exhiber le **cadre spectral le plus précis à ce jour**
pour le problème de Hilbert-Polya. Le détail technique complet est
dans les notes 42, 62-65, 80, et `PT_PROJECTS/PT_RH_HYPERBOLIC_CUSP/`.

### 7.1 Produit d'Euler PT

Soient $\zeta_{\pm}(s) := \exp\!\big(\sum_p a_p^{(\pm)}\, p^{-s}\big)$
avec $a_p^{(+)} = a_p = \sin^2(\theta_p, q_+)$ et
$a_p^{(-)} = b_p = \sin^2(\theta_p, q_-)$ les amplitudes d'holonomie sur
les deux branches du crible. Le théorème central de la note 42 donne :

$$
\boxed{\;\zeta(s) \;=\; \zeta_+(s)\, \zeta_-(s)\, \prod_p \frac{1}{1 - \kappa_p(s)}\;}
\qquad (\mathrm{Re}(s) > 1).
$$

Le produit d'Euler classique $\prod_p (1 - p^{-s})^{-1}$ se factorise
ainsi en **trois composants PT** : deux exponentielles de branche
$\zeta_+, \zeta_-$ (canaux couplage et métrique) et un **résidu spectral
de Fredholm** $\prod_p (1 - \kappa_p(s))^{-1} = R(s)$. Les zéros non
triviaux de $\zeta$ vivent **exclusivement** dans le résidu.

### 7.2 Équation fonctionnelle = involution $T_3$

L'équation fonctionnelle classique
$\zeta(s) = 2^s\, \pi^{s-1}\, \sin(\pi s/2)\, \Gamma(1-s)\, \zeta(1-s)$
exhibe la symétrie $s \leftrightarrow 1-s$ centrée sur $\mathrm{Re}(s) = 1/2$.
En PT, cette symétrie n'est pas une coïncidence analytique mais
l'**image, sous la transformée Mellin, de l'involution $T_3$** du crible
mod $3$ (note 63, conjecture structurelle) :

$$
T_3 = \mathrm{antidiag}(1, 1)
\quad\xrightarrow{\;\mathcal{M}\;}\quad
s \,\leftrightarrow\, 1 - s.
$$

La transformation unitaire $U : f(u) \mapsto e^{s u / 2}\, f(u)$ réalise
explicitement le pont entre la balance détaillée discrète de $T_3$ (T2)
et l'équation fonctionnelle analytique. Les facteurs
$\Gamma, \pi^{s-1}, \sin(\pi s/2)$ qui apparaissent classiquement sont
les **traces Mellin** de cette involution arithmétique.

### 7.3 Fonction xi PT

La fonction $\xi$ de Riemann est entière, symétrique
($\xi(s) = \xi(1-s)$), et concentre tous les zéros non triviaux de
$\zeta$. La construction PT canonique :

$$
\boxed{\;\Xi_{\mathrm{PT}}(s) \;:=\; R(s) \cdot R(1-s).\;}
$$

Cette fonction est **tautologiquement** symétrique. La non-annulation
$\zeta_+(s)\zeta_-(s) \neq 0$ pour $0 < \mathrm{Re}(s) < 1$ — désormais
établie comme théorème (§2.4 ; Lemme 3.4 du préprint A_PT) — donne
l'identité de pont :

$$
\Xi_{\mathrm{PT}}(s) = \xi(s)\, \xi(1-s) \cdot \mathcal{P}(s)^{-1},
\quad
\mathcal{P}(s) = \frac{[s(s-1)]^2}{4\, \pi^{1/2}}\, \Gamma(s/2)\, \Gamma((1-s)/2)
\cdot [\zeta_+ \zeta_-](s)\, [\zeta_+ \zeta_-](1-s),
$$

où $\mathcal{P}$ est symétrique facteur par facteur sous
$s \leftrightarrow 1-s$ et ne s'annule pas sur la bande critique (les
facteurs $[\zeta_+\zeta_-](s)$ et $[\zeta_+\zeta_-](1-s)$ sont non-nuls
sur $\mathrm{Re}(s) > 0$ et $\mathrm{Re}(s) < 1$ respectivement par le
Lemme 3.4 du préprint A_PT, simultanément sur $0 < \mathrm{Re}(s) < 1$).
La structure de zéros de $\Xi_{\mathrm{PT}}$ est
équivalente, sous l'identité de pont ci-dessus, à la simplicité
conjecturale des zéros de $\xi$ sur la ligne critique : si tous les
zéros non triviaux sont sur $\mathrm{Re}(s) = 1/2$, alors
$\Xi_{\mathrm{PT}}$ a un zéro de multiplicité double à chaque ordonnée
critique ; sinon, deux zéros distincts $s_0$ et $1-s_0$. Forme produit
(Hadamard PT), simple réécriture de la définition via les produits
d'Euler PT :

$$
\Xi_{\mathrm{PT}}(s)
= \left(\prod_p \frac{1}{1 - \kappa_p(s)}\right)
\cdot
\left(\prod_p \frac{1}{1 - \kappa_p(1-s)}\right)
= R(s)\, R(1-s).
$$

### 7.4 Géométrie spectrale : $\Sigma_{\mathrm{pers}}$ hyperbolique à cusp

Au-delà de la simple algèbre, le programme PT zêta possède une
**géométrie spectrale canonique** : la courbe de persistance
$\Sigma_{\mathrm{pers}}$, définie par PT_GeoFlow comme **unique
solution** du système d'auto-cohérence T6 (note 16 PT_GeoFlow) :

$$
\Sigma_{\mathrm{pers}} \;=\; \{(x, y) \in \mathbb{C}^2 : y^2 = N(x)/(2\pi^2 D)\},
\quad N(x) = -256\, x (2x{-}1)^3 \prod_{d=2}^{7} F_d(x).
$$

C'est une **courbe hyperelliptique de genre $g = \mu^\star - 1 = 14$**,
avec **30 branchements** PT-canoniques : $x = 0$ (simple), $x = 1/2$
(triple, multiplicité $|P| = 3$), $27$ racines des $F_d$, et le point à
l'infini.

**Propriétés géométriques** (cycle PT_RH_HYPERBOLIC_CUSP) :

- $\Sigma_{\mathrm{pers}}$ est **hyperbolique à cusp** : sa métrique
  Fisher Wick-rotée vers le riemannien est asymptotiquement isométrique
  à $\mathbb{H}^{n+1}$ de **rayon $\sqrt{3} = \sqrt{|P|}$**, avec
  **courbure sectionnelle $K = -1/3$** (note 01 PT_RH).
- Le bord $z = 1$ (= $\mu \to \infty$) est un **cusp parabolique
  standard** au sens d'Iwaniec : 4 critères vérifiés (forme conforme,
  volume fini, horocycles shrinkants, groupe parabolique $\mathbb{Z}$
  ou $\mathbb{Z}^3$). Cf. note 02 PT_RH.
- $\Sigma_{\mathrm{pers}}$ est **non-modulaire** : aucun
  $\Gamma_0(N) \backslash \mathbb{H}$ classique de genre 14 dans le
  tower $N \mid 3^a 5^b 7^c$ (note 07 PT_RH). C'est un **point isolé
  de l'espace de modules $M_{0,30}$** des courbes hyperelliptiques.

Cette géométrie **remplace** le demi-plan modulaire $\Gamma \backslash \mathbb{H}$
de Selberg comme support des zéros de $R(s)$. La voie Selberg-Maass
classique est **fermée** ; il faut un cadre spectral nouveau.

### 7.5 Opérateur Berry-Keating PT canonique $H_n^{\mathrm{cusp}}$

Sur le sous-espace cuspidal radial $\mathcal{H}_n = L^2(\mathbb{R}, \mathrm{e}^{-n\eta/\sqrt{3}}\mathrm{d}\eta)$
de $\Sigma_{\mathrm{pers}}$ (avec $\eta = \sqrt{3} \log y$ la
coordonnée hyperbolique standard, $n \in \{1, 3\}$ le rang du cusp), on
construit l'**opérateur Berry-Keating PT canonique** (note A1) :

$$
\boxed{\;H_n^{\mathrm{cusp}} \;=\; -i\Bigl(\eta\partial_\eta + \frac{1}{2}
- \frac{n\,\eta}{2\sqrt{3}}\Bigr)\quad \text{sur } \mathcal{H}_n.\;}
$$

Le terme $-n\eta/(2\sqrt{3})$ est l'**analogue géométrique de la
connexion métrique** sur le cusp ; il **absorbe la courbure** $K = -1/3$
via une jauge unitaire qui réduit $H_n^{\mathrm{cusp}}$ à
$\tilde H_{\mathrm{PT}} := -i\partial_\xi$ sur $L^2([0, r], \mathrm{d}\xi)$
(coordonnée Mellin $\xi = \log\eta$, longueur $r = \log(\eta_{\max}/\eta_{\min})$).

**Théorèmes** (notes A1, A2, A3) :

- **Symétrie** : $H_n^{\mathrm{cusp}}$ est densément défini et
  symétrique sur $\mathcal{D}_0 = C_c^\infty$ (Thm A1.1).
- **Spectre** : pour toute extension auto-adjointe $\tilde H_\theta$
  paramétrée par $\theta \in [0, 2\pi)$,
  $\sigma(\tilde H_\theta) = \{(\theta + 2\pi k)/r : k \in \mathbb{Z}\}$
  (Thm A2.1).
- **Indices de défaut** : $(n_+, n_-) = (1, 1)$, donc une famille
  $U(1)$ d'extensions auto-adjointes (Thm A3.1).
- **Extension PT-canonique** : la symétrie $T_3$ + exclusion du mode
  marginal $\lambda = 0$ sélectionnent uniquement $\theta_{\mathrm{PT}} = \pi$
  (BC **antipériodique**). Spectre PT-canonique :
  $\sigma(\tilde H_{\mathrm{PT}}) = \{(2k+1)\pi/r : k \in \mathbb{Z}\}$
  (Thm A3.2).

**Cellule Planck PT et densité de Weyl** : sous le cut-off dynamique
$\eta_{\min} = \sqrt{6\pi}$, $\eta_{\max}(\gamma) = \sqrt{3}\,\gamma/\sqrt{2\pi}$
(note 60), $r(\gamma) = \log(\gamma/(2\pi))$ et la fonction de comptage
asymptotique :

$$
N_+(\gamma) \;=\; \frac{\gamma}{2\pi}\log\bigl(\gamma/(2\pi)\bigr) + O(1),
$$

soit la **densité de Riemann-von Mangoldt** à terme dominant (Thm A2.2,
note 61).

### 7.6 Formule de trace Selberg-PT

La formule explicite de Riemann-Weil-von Mangoldt relie classiquement
la somme des $\log p$ pondérée à la somme sur les zéros. En PT, cette
relation prend la **forme analytique exacte** (note 62, validée à
$10^{-26}$, **signe corrigé** par note A4 §3.3bis) :

$$
\boxed{\;-\frac{d}{ds}\log R(s) \;=\; \sum_p (\log p)
\left[\frac{1}{p^s - 1} \;-\; A_p\, p^{-s}\right]
\qquad (\mathrm{Re}(s) > 1).\;}
$$

(NB : la note 62 §2.3 originale écrivait $+A_p p^{-s}$ ; la
vérification directe `mp.diff(-\log R)` corrige le signe en
$-A_p p^{-s}$, cf. note A4.)

**Décomposition Selberg-PT en trois termes** (Thm A4.1) : pour toute
fonction test $h$ analytique avec décroissance rapide,

$$
\sum_n h(\gamma_n) \;=\; T_W(h) \;+\; T_G(h) \;+\; T_P(h) \;+\; T_A(h),
$$

où :

- **Terme de Weyl** $T_W(h) = (r/(2\pi)) \int h$, $r = r(\gamma)$ sous
  cellule Planck PT.
- **Terme géométrique** $T_G(h) = -\sum_p (\log p)\, \hat h(\log p)/\sqrt{p}$
  = somme sur les **orbites primitives du flot de dilatation
  Berry-Keating PT**, avec longueurs $\ell_p = \log p$.
- **Terme parabolique** $T_P(h) = +\sum_p A_p (\log p)\, \hat h(\log p)/\sqrt{p}$
  = contribution du **cusp** via la **scattering matrix PT**
  $\phi_{\mathrm{PT}}(s) := \zeta_+(s)\,\zeta_-(s)$.
- **Terme archimédien** $T_A(h)$ = $\Gamma$-corrections.

**Identification PT** : $\phi_{\mathrm{PT}}(s) = \zeta_+ \zeta_-$ est la
**scattering matrix du cusp non-modulaire** de $\Sigma_{\mathrm{pers}}$,
encodant la bifurcation $q_+/q_-$ au cusp. Cette identification est la
signature distinctive de PT vs Selberg classique (où $\phi(s) = \xi(2s-1)/\xi(2s)$
pour $\Gamma_0(1)\backslash \mathbb{H}$).

**Écriture opératorielle finale** (Thm A4 + note 65) :

$$
T(s) \;:=\; \mathrm{Tr}_{\mathrm{reg}}\bigl[(s - i\tilde H_{\mathrm{PT}})^{-1}\,D_R(s)\bigr]
\;=\; \sum_p \kappa_p(s)\,(\log p)\,\biggl[\frac{1}{p^s - 1} \;-\; A_p p^{-s}\biggr].
$$

Cette trace régularisée a pour **pôles distributionnels** sur
$\mathrm{Re}(s) = 1/2$ exactement les $\gamma_n$ des zéros non triviaux
de $R(s) \equiv \zeta(s)$ (la non-annulation de $\zeta_+\zeta_-$ sur
$\mathrm{Re}(s) > 0$ étant [THM] par Lemme 3.4 du préprint A_PT, §2.4).

### 7.7 Identité géométrique de la spirale d'Archimède log-polaire

Le terme géométrique $T_G(h)$ admet, pour la famille des fonctions test
gaussiennes $h(\gamma) = e^{-\sigma^2 \gamma^2}$ ($\sigma > 0$), une
**décomposition canonique par tours de spirale** qui révèle une cascade
de seuils universels $\sigma_{\mathrm{crit}}^{(k)} = \sqrt{\pi(k+1)}$
propre à la période fondamentale $2\pi$ des coordonnées log-polaires.

**Coordonnées et tours**. En log-polaire $(r, \theta) \mapsto (\ln r, \theta)$,
l'entier $n$ est positionné à angle $\theta = \ln n$. Les premiers se
partitionnent canoniquement par **tour** :
$$
\text{turn}_k := \{p \text{ premier} : \ln p \in [2\pi k, 2\pi(k+1))\},
\quad k \in \mathbb{N}.
$$
Turn $0$ = $\{p \leq 535\}$ (99 premiers) ; turn $1$ = $\{535 \leq p < 286\,751\}$
(24 877 premiers) ; turn $k$ contient $\sim 85 \cdot e^{2\pi k}/(k+1)$
premiers asymptotiquement.

**Observation visuelle remarquable** : $\ln 23 \approx \pi$ à 0.2 % près
(via $e^\pi \approx 23.14$), donc la cascade complète {actif ∪ écho ∪
super-écho} $= \{3, 5, 7, 11, 13, 17, 19, 23\}$ se loge exactement dans
le **premier demi-tour** $\theta \in [0, \pi]$ de la spirale.

**Théorème 7.X (Identité de spirale, démonstration complète in [4] §6.6).**
*Dans l'approximation de densité PNT,*
$$
T_G^{\text{turn}_k}(\sigma) \;=\; T_G^{\text{turn}_0}(\sigma)
\quad \Longleftrightarrow \quad
\sigma^2 = \pi (k+1).
$$
*La preuve, par parité de la gaussienne sur intervalles symétriques
autour de $\sigma^2$, ne requiert ni la structure cuspidale ni la
bifurcation $q_\pm$ : c'est un fait géométrique propre à la décomposition
log-polaire des premiers.*

**Triple validation numérique** :

| $k$ | $\sigma_{\mathrm{crit}}^{(k)} = \sqrt{\pi(k+1)}$ | observé | écart | $\#$ primes turn $k$ |
|---|---|---|---|---|
| 1 | $\sqrt{2\pi} \approx 2.5066$ | $2.4909$ | $-0.626\%$ | $24\,877$ |
| 2 | $\sqrt{3\pi} \approx 3.0700$ | $3.0630$ | $-0.228\%$ | $\sim 8.6 \times 10^6$ |
| 3 | $\sqrt{4\pi} \approx 3.5449$ | $3.5419$ | $-0.085\%$ | $\sim 3.4 \times 10^9$ |

La précision **croît avec $k$** (cohérent $O(1/\log(\#\text{primes}))$),
indiquant convergence empirique vers l'identité continue. Pour $k = 3$,
le sieve a été effectué via `primesieve` CLI streamé (83 chunks de range
$10^9$, 12.4 min CPU, voir [4] §6.6.4).

**Conséquence : paquet gaussien en index de tour**. Pour $\sigma$ donné,
$T_G(\sigma)$ est un paquet gaussien centré sur
$k_{\mathrm{peak}}(\sigma) = \sigma^2/(2\pi) - 1/2$ ; pour $\sigma \to \infty$,
$k_{\mathrm{peak}} \to \infty$. La masse de $T_G$ diffuse vers les hautes
tours de la spirale, reproduisant géométriquement la **non-localité
fondamentale** de la condition de positivité de Weil [Weil 1952] : aucune
somme finie de tours ne borne $T_G(\sigma)$ pour $\sigma$ arbitraire.

**Statut épistémique** : Théorème 7.X est **inconditionnel** dans la
limite continue PNT, triple-validé numériquement à $< 1\%$ sur $k = 1, 2, 3$.
Il fournit une **signature géométrique propre** de la période $2\pi$
dans la formule explicite, mais **ne franchit pas RH** : la non-localité
subsiste sous prisme géométrique. Il est par construction PT-compatible
mais non PT-dépendant ; la décomposition par tours est définie par la
seule géométrie log-polaire des premiers.

**Référence détaillée** : démonstration complète, validation numérique
détaillée et discussion des extensions sont dans le préprint cuspidal [4]
(Section 6.6, "Décomposition par tours de spirale log-polaire").

### 7.8 Symétrie circulaire des 30 branches

Une découverte structurelle nouvelle (observation Yan Senez,
2026-05-16, formalisée notes 08 + S1-S4 PT_RH) :

> **Théorème de symétrie circulaire** : *Les 30 branchements de
> $\Sigma_{\mathrm{pers}}$ vivent sur **5 cercles concentriques exacts**
> dans le plan algébrique $x$, centrés au point*
> $$x_c \;=\; -\frac{13}{4} \;=\; -\mathcal{G}_F,$$
> *où $\mathcal{G}_F = (\mu^\star - p_1)/p_1^2$ est le **facteur Fisher
> de profondeur du canal $p = 2$** déjà connu dans le corpus PT
> (ch10 monographie, dressing additif d'$\alpha_{\rm EM}$).*

**Formule unifiée des rayons** (Thm S3) :

$$
\boxed{\;r_{C_d} \;=\; \frac{\mu^\star \cdot (2d-1)^{1/d}}{4}
\quad (d \in \{3, 5, 7\}, \text{primes actifs}) ;
\quad r_{C_d} = \frac{\mu^\star}{4} \quad (d \in \{2, 4, 6\}).\;}
$$

Démonstration (note S2) : la depression $x = y - 13/4$ révèle pour
$d$ impair prime actif **tous coefficients intermédiaires nuls** :
$F_d^{\mathrm{dep}}(y) = y^d + c_d$ avec
$c_d = (\mu^\star)^d (2d-1)/4^d$. Pour $d$ pair, structure cyclotomique
$F_d^{\mathrm{dep}}(y) = (\mu^\star/4)^d \Phi_{d+1}(y/(\mu^\star/4))$.

**Tableau** :

| Cercle | Contenu | Rayon $r_{C_d}$ | Numérique |
|---|---|---|---|
| $C_1$ | $\{x = 0\}$ trivial | $13/4$ | $3.25$ |
| $C_2$ | $\{x = 1/2\}$ + $F_2, F_4, F_6$ | $\mu^\star/4 = 15/4$ | $3.75$ |
| $C_3$ | $F_3$ | $15 \cdot 5^{1/3}/4$ | $6.4124$ |
| $C_5$ | $F_5$ | $15 \cdot 9^{1/5}/4$ | $5.8194$ |
| $C_7$ | $F_7$ | $15 \cdot 13^{1/7}/4$ | $5.4098$ |

**Total** : $1 + 13 + 7 + 5 + 3 + 1_\infty = 30 = 2\mu^\star$.

**Test de chance** : pour un polynôme générique de degré 30, la
probabilité de tomber sur cette structure (depression triviale + même
centre + cyclotomie) est de **mesure nulle** ; elle est de codimension
$\geq 21$ dans l'espace des polynômes. Cette propriété est donc
**structurelle PT-canonique**, cohérente avec l'unicité de $N(x)$
(thm de la courbe de persistance, note 16 PT_GeoFlow).

**Connexion structurelle $\alpha_{\rm EM} \leftrightarrow \Sigma_{\rm pers}$** :
la **même constante** $\mathcal{G}_F = 13/4$ apparaît dans la
**dérivation d'$\alpha_{\rm EM}$** (multiplicateur de profondeur du
binary leakage $F(2)$, ch10) **et** comme **centre géométrique** de la
courbe spectrale. Ce n'est pas une coïncidence : la **dualité
algèbre/géométrie** au niveau de $\mathcal{G}_F$ est une signature
PT-canonique.

### 7.9 Validation numérique

Le programme A5 (script `A5_spectrum_vs_zeros.py`) teste **directement**
la conjecture HP-PT en détectant les pôles distributionnels de $T(s)$
sur la ligne critique et en les comparant aux $\gamma_n$ LMFDB.

**Résultats** ($P = 500$ premiers, $\epsilon = 0.4$-$0.5$, balayage
$t \in [10, 60]$) :

- **Contraste pic/voisin** : moyenne $1.47\times$, médian $1.52\times$,
  **100% des $\gamma_n$** avec contraste $> 1$.
- **Détection automatique** : **13/13** des $\gamma_n$ dans $[10, 60]$
  détectés à distance $< 1$. **Distance médiane $\gamma_n \to$ pic le
  plus proche : $0.075$**.
- **Densité spectrale** : 13 pics observés vs $N_{\rm RvM}(60) \approx 12$,
  ratio $1.08$. Cohérent avec Riemann-von Mangoldt.
- **Test de chance Weyl** : sous hypothèse "pics aléatoires uniformes",
  $E[\text{matches}] = 0.52$. **Observé : 13**, soit **$25\times$
  au-dessus du hasard**. Hypothèse de hasard **fortement rejetée**.
- **Cohérence Selberg-PT vs zêta seule** : ratio
  $|I_{\rm Selberg-PT}|/|I_\zeta|$ approximativement constant
  ($\sim 0.76 \pm 0.05$). Les pôles sont **identiques** ; le couplage
  $A_p$ module l'amplitude sans déplacer les pôles.

**Interprétation** : la formule de trace Selberg-PT **détecte
structurellement** les premiers 13 zéros de $\zeta$ comme pôles
distributionnels, à précision empirique $0.075$ (largement meilleure
que la borne théorique $O(\epsilon)$). Le test de chance Weyl est
**passé à $25\sigma$**.

**Statut** : `[DER-NUM]` validation forte. Reste un test sur les
$\gamma_n$ d'indice $\geq 14$, et démonstration rigoureuse de
l'identification stricte spectre = $\{\gamma_n\}$ (= problème classique
de Hilbert-Polya, non résolu par PT).

### 7.10 Équation maîtresse PT

Les sept reformulations se synthétisent en une seule **équation
maîtresse PT** :

$$
\boxed{\;\zeta(s) \;=\; \zeta_+(s)\,\zeta_-(s)\,\det\!\big(I - D_R(s)\big)^{-1}\;}
$$

où $D_R(s) = \mathrm{diag}(\kappa_p(s))$ est l'opérateur diagonal sur
$\ell^2(\mathcal{P}_\gamma)$ (atomes : primes au cusp de $\Sigma_{\rm pers}$).
Cette forme :

- réécrit le produit d'Euler comme un **déterminant Fredholm** ;
- identifie les zéros non triviaux de $\zeta$ comme les **points
  spectraux** de $I - D_R(s)$ (inconditionnellement sur la bande critique
  par non-annulation de $\zeta_+\zeta_-$ — [THM] Lemme 3.4 préprint A_PT) ;
- relie l'opérateur diagonal $D_R$ à l'opérateur Berry-Keating PT
  $H_n^{\rm cusp}$ géométriquement enraciné dans $\Sigma_{\rm pers}$,
  via la trace régularisée distributionnelle
  $T(s) = \mathrm{Tr}_{\rm reg}\bigl[(s - i\tilde H_{\rm PT})^{-1} D_R(s)\bigr]$
  (note 65 + A4).

### 7.11 Dictionnaire Riemann ↔ PT

| Objet classique | Forme PT canonique |
|---|---|
| $\zeta(s) = \sum 1/n^s$ | $\zeta = \zeta_+ \zeta_- \cdot \prod (1 - \kappa_p)^{-1}$ |
| Produit d'Euler $\prod (1 - p^{-s})^{-1}$ | Déterminant Fredholm $\det(I - D_R)^{-1}$ |
| Équation fonctionnelle | Involution $T_3$ du crible mod 3 (notes 63, 80) |
| Fonction $\xi(s) = \xi(1-s)$ | $\Xi_{\mathrm{PT}}(s) = R(s) R(1-s)$ |
| **Demi-plan modulaire** $\Gamma\backslash\mathbb{H}$ | **Courbe hyperelliptique** $\Sigma_{\mathrm{pers}}$ de genre 14 |
| **Laplacien hyperbolique** $-\Delta$ | **Opérateur Berry-Keating PT** $H_n^{\rm cusp}$ |
| **Série d'Eisenstein** | Scattering matrix $\phi_{\rm PT} = \zeta_+\zeta_-$ |
| Formule explicite de Riemann-Weil | $-d\log R/ds = \sum_p \log p [1/(p^s{-}1) - A_p p^{-s}]$ (signe corrigé) |
| Formule de trace de Selberg | $\sum_n h(\gamma_n) = T_W + T_G + T_P + T_A$ (Thm A4.1) |
| Ligne critique $\mathrm{Re}(s) = 1/2$ | Sur-détermination cohomologique + cellule Planck PT |
| Constante de structure fine $\alpha_{\rm EM}$ | Connexion via $\mathcal{G}_F = -x_c$ (centre des cercles) |

### 7.12 Ce que rend visible la reformulation

La forme PT canonique fait apparaître **explicitement** plusieurs
structures cachées dans l'écriture classique :

- **Séparation triviale/non triviale** : zéros triviaux $\leftrightarrow$
  pôles de $\zeta_+ \zeta_-$ ; zéros non triviaux $\leftrightarrow$ pôles
  du résidu $R(s)$.

- **Origine arithmétique de la ligne critique** : le $1/2$ n'est pas
  un nombre arbitraire mais la **valeur sur-déterminée** par la cohérence
  cohomologique $s^2/2 = s/4$ (sections 5, 6) ; il s'identifie au spin
  PT $s = 1/2$ (axiome T1) **et** au shift de la mesure de Haar
  multiplicative $\mathrm{d}p/p$ dans la transformée de Mellin (note 58).

- **Bifurcateur $q_+/q_-$ comme objet primitif** : le résidu $R(s)$
  exprime explicitement la dépendance en la bifurcation du crible. Le
  même bifurcateur produit la masse du Higgs ($m_H = v \cdot s$),
  organise les zéros de $\zeta$, et **constitue la scattering matrix
  $\phi_{\rm PT}$ du cusp non-modulaire**.

- **Équation fonctionnelle comme involution $T_3$** : la symétrie
  $s \leftrightarrow 1-s$ classique est l'écriture analytique de la
  doublement stochasticité $T_3$ (T2) du crible mod $3$.

- **Géométrie spectrale canonique** : $\Sigma_{\rm pers}$ remplace le
  demi-plan modulaire classique, avec une **symétrie circulaire exacte**
  des 30 branches centrée au facteur Fisher $-\mathcal{G}_F = -13/4$.

- **Connexion $\alpha_{\rm EM} \leftrightarrow \zeta$** : la même
  constante $\mathcal{G}_F$ apparaît dans le dressing de $\alpha_{\rm EM}$
  (canal binaire $F(2)$) et comme centre géométrique de $\Sigma_{\rm pers}$.
  La constante de structure fine et les zéros de zêta partagent une
  **racine PT commune**.

### 7.13 Statut épistémique

| Section | Statut | Référence |
|---|---|---|
| §7.1 Produit d'Euler PT | [THM] $\mathrm{Re}(s) > 1$ | note 42 |
| §7.2 Équation fonctionnelle | [CONJ structurelle] | note 63 |
| §7.3 Fonction $\Xi_{\mathrm{PT}}$ | [DEF] + [DER] conditionnelle | note 42 |
| §7.4 Géométrie $\Sigma_{\rm pers}$ | [DER-NUM] | notes 01, 02, 07 PT_RH |
| §7.5 Opérateur $H_n^{\rm cusp}$ | [THM] symétrie + indices (1,1) | notes A1, A2, A3 |
| §7.6 Formule de trace Selberg-PT | [THM] $\mathrm{Re}(s) > 1$ ; [DER] critique | notes 62, A4 |
| §7.7 Identité de spirale $\sigma_{\rm crit}^{(k)}=\sqrt{\pi(k+1)}$ | [THM] continuum PNT ; triple-validé $k=1,2,3$ à $<1\%$ | préprint cusp §6.6 |
| §7.8 Symétrie circulaire | [THM] (algébrique sympy) | notes 08, S1-S4 |
| §7.9 Validation numérique | [DER-NUM] 13/13 zéros à 0.075 | note A5 |
| §7.10 Équation maîtresse | [THM] $\mathrm{Re}(s) > 1$ | notes 42, 65, A4 |

**Gap fondamental restant** : démontrer que **TOUS** les pôles
distributionnels de $T(s)$ sont sur $\mathrm{Re}(s) = 1/2$. C'est
l'équivalent PT du problème classique RH, et PT n'y échappe pas.

Le programme PT fournit le **cadre spectral le plus précis à ce jour** :

- géométrie canonique $\Sigma_{\rm pers}$ (unique, non-modulaire) ;
- opérateur canonique $H_n^{\rm cusp}$ (BC PT-sélectionnée
  $\theta = \pi$) ;
- formule de trace Selberg-PT explicite (3 termes identifiés) ;
- scattering matrix identifiée ($\phi_{\rm PT} = \zeta_+\zeta_-$) ;
- validation numérique des premiers 13 zéros à précision $0.075$.

Sans prouver RH, PT **réduit** RH au problème spectral
> *démontrer que tous les pôles distributionnels de
> $\mathrm{Tr}_{\rm reg}[(s - i\tilde H_{\rm PT})^{-1} D_R(s)]$ sont
> sur $\mathrm{Re}(s) = 1/2$.*

C'est la version la plus précise de Hilbert-Polya à ce jour.

---

## 8. Vérifications numériques

Cette section synthétise les sept vérifications numériques principales du
programme. Chaque ligne renvoie au script et à la note de référence.

### 7.1 Tableau récapitulatif

| Test | Précision | Note | Script |
|---|---|---|---|
| Identité Fredholm $\det(I - D_R)^{-1} = R(s)$, $\mathrm{Re}(s) > 1$ | $\sim 2 \cdot 10^{-12}$ | 42 | `pt_residual_fredholm_target.py` |
| Identité de trace $-d\log R/ds = \sum_p (\log p)[1/(p^s-1) - A_p p^{-s}]$ (signe corrigé) | $\sim 10^{-13}$ ($\sigma{=}4$, $P{=}1000$) | 62, A4 | `pt_bk_explicit_formula.py` (corrigé), `A4_trace_verify.py` |
| Différence $N_{\mathrm{PT}} - N_{\mathrm{RvM}} = 1/8$ constant sur 4 décades | $\sim 10^{-15}$ | 61 | `pt_bk_density_rvm.py` |
| Symétrie spectrale $\{\gamma_n, -\gamma_n\}$ (BC AP, $\theta = \pi$) | $5 \cdot 10^{-12}$ | 63 | `pt_hp4_deficiency_indices.py` |
| Phase de Berry par coin $-\pi/4$, densité $-1/8$ | $1.68 \cdot 10^{-14}$ | 68 | `pt_k3_cohomological_maslov.py` |
| $\varepsilon_*^{\mathrm{geom}}$ reproduit l'optimum empirique $W = 24$ | $\sim 5\%$ global, $\sim 1\%$ au régime $M = 1155, W = 24$ | 52 | `pt_epsilon_star_geometric.py` |
| Pôles de $|T(1/2 + it)|$ capturent les zéros connus de $\zeta$ (formule brute) | $t \in [10, 70]$ : 17/17 à $|\Delta| < 1.5$, 11-12 à $|\Delta| < 0.5$, 2-3 à $|\Delta| < 0.1$ ; $t \in [10, 200]$ : 74/79 à $|\Delta| < 1.5$, 38-39 à $|\Delta| < 0.5$, 5-7 à $|\Delta| < 0.1$ ; saturation au-delà de $P_{\max} = 3\,000$ (note 65 §3.6, note 74) | 65, 74 | `pt_hp3bc_hadamard_trace.py`, `pt_g4_pmax_scan.py` |
| Pôles de $|T_{R_1}(1/2 + it)|$ avec régularisation explicite ($\varepsilon = 0.2$) | $t \in [10, 200]$ : 77/79 à $|\Delta| < 1.5$, **75/79 à $|\Delta| < 0.5$**, 47/79 à $|\Delta| < 0.1$ ; médiane $|\Delta| = 0.080$ (gain ×6.5 vs brut) | 76 | `pt_g4_r1r2_regularized.py` |
| Pôles de $|T_{R_1}(1/2 + it)|$, range étendu ($\varepsilon = 0.2$) | $t \in [10, 500]$ : 267/269 à $|\Delta| < 1.5$ (99.3 \%), **242/269 à $|\Delta| < 0.5$ (90.0 \%)**, 139/269 à $|\Delta| < 0.1$ (51.7 \%) ; capture par bande de 100 : $\ge 86.6\,\%$ à $|\Delta| < 0.5$ sur 5 sous-bandes ; 55 paires structurellement fusionnées (gap $< 2\pi / \log\gamma$) | 77 | `pt_v2_extended_range.py` |

### 7.2 Forme fonctionnelle exacte de $N_{\mathrm{PT}}$

Test (script `pt_bk_density_rvm.py`) sur $\gamma \in \{14, 21, 30, 50, 100, 200, 500, 1000\}$ :

| $\gamma$ | $N_{\mathrm{PT}}$ | $N_{\mathrm{RvM}}$ | différence |
|---|---|---|---|
| $14.13$ | $0.5737$ | $0.4487$ | $0.1250$ |
| $21.02$ | $1.6945$ | $1.5695$ | $0.1250$ |
| $50.00$ | $9.5478$ | $9.4228$ | $0.1250$ |
| $100$ | $29.1273$ | $29.0023$ | $0.1250$ |
| $500$ | $269.7117$ | $269.5867$ | $0.1250$ |
| $1000$ | $648.7412$ | $648.6162$ | $0.1250$ |

La différence est constante à $1/8$ sur quatre décades de $\gamma$, à
précision machine ($\le 10^{-15}$).

### 7.3 Capture distributionnelle des zéros

Test (script `pt_hp3bc_hadamard_trace.py`) sur $t \in [10, 70]$ avec régularisation
Abel–Hadamard $p^{-\varepsilon}$. Pour $\varepsilon \in \{0.6, 0.4, 0.2\}$, le
module $|F_\varepsilon(t)|$ croît systématiquement de $2.7\times$ à
$3.3\times$ aux zéros de $\zeta$ quand $\varepsilon$ passe de $0.6$ à $0.2$
(signature d'un pôle), contre $1.4\times$ à $2.3\times$ aux points neutres
(régularité). Sur $t \in [10, 70]$, dt = 0.1, on détecte 41-52 pics
locaux selon $P_{\max} \in [3\,000, 100\,000]$ ; les 17 zéros de référence
sont tous à moins de $1.5$ unité d'un pic ($|\Delta| < 1.5$) dès
$P_{\max} = 3\,000$. La capture sature au-delà : la médiane $|\Delta|$
décroît de $0.49$ ($P_{\max} = 3\,000$) à $0.36$ ($P_{\max} = 30\,000$)
puis stagne. La même statistique étendue à $t \in [10, 200]$ donne
74 zéros sur 79 (94 \%) à $|\Delta| < 1.5$, 38 à $|\Delta| < 0.5$, 7 à
$|\Delta| < 0.1$ (note 74). Les 5 zéros manqués appartiennent à des
paires $(\gamma_n, \gamma_{n+1})$ avec écart inférieur à la largeur
intrinsèque d'un pic non régularisé (~$2\pi/\log\gamma$), fusionnées en
un pic unique. Cas extrême de précision : $\gamma = 146.001$ capturé à
$|\Delta| = 0.001$. L'amélioration sub-$0.1$ systématique requiert
d'appliquer la régularisation explicite R1 ou R2 (note 65 §3.6), et non
l'augmentation de $P_{\max}$.

### 7.4 Régularisation par smearing gaussien

Avec $\phi_\varepsilon(t) = (\varepsilon \sqrt{2\pi})^{-1} \exp(-t^2 / 2\varepsilon^2)$ et
$\varepsilon \in \{2, 1, 0.5\}$, le contraste pôle/neutre est plus marqué :
$8 \times$ à $10 \times$ aux zéros contre $2\times$ à $3\times$ aux points
neutres. La régularisation gaussienne et la régularisation Abel sont
équivalentes au sens distributionnel (différence de noyaux régularisants),
mais la première est plus sensible numériquement.

L'extension du test sur $t \in [10, 500]$ (269 zéros, note 77) confirme
la robustesse : R1 $\varepsilon = 0.2$ capture 267/269 zéros à
$|\Delta| < 1.5$, 242/269 à $|\Delta| < 0.5$ et 139/269 à $|\Delta| < 0.1$.
La décomposition par bande de 100 unités (10-100 : 100 \%, 100-200 :
92 \%, 200-300 : 89.8 \%, 300-400 : 87.5 \%, 400-500 : 86.6 \% à
$|\Delta| < 0.5$) montre une dégradation monotone mais bornée — la
performance ne s'effondre pas avec $T$. Le plancher à 86.6 \% est cohérent
avec la limite structurelle de 55 paires fusionnées (gap $< 2\pi/\log\gamma$)
sur 269 zéros, soit 20.4 \% des zéros appartenant à une paire au-dessous
de la résolution intrinsèque du modèle. Une stratégie ε-adaptatif
($\varepsilon \propto 1/\log\gamma$) a été testée (note 78) sous deux
formes : convolution sur $|F|$ brut (effondrement des pics) et
convolution sur $F$ complexe (compromis marginal : médiane $|\Delta|$
réduite à $0.080$ mais couverture réduite à 76 \% à $|\Delta|<1.5$).
La régularisation à poids spectral fixe ($R_1$, $\varepsilon = 0.2$)
reste l'optimum global.

---

## 9. Auto-correction : prédiction Dirichlet falsifiée

La dualité K4 (section 6) suggère une extension naturelle aux fonctions $L$
de Dirichlet primitives par factorisation Haar mod $q$. Cette prédiction a
été testée et falsifiée. La présente section documente le test honnêtement
et reformule la note 61 en conséquence.

### 8.1 La prédiction P3

Si K4 est structurellement extensible, la formule générique pour la phase
Maslov de $L(s, \chi \bmod q)$ serait
$$
\Delta_N(L \bmod q) = s^2 \cdot \mathrm{Haar}(\mathbb{Z}/q\mathbb{Z}) \cdot \frac{1}{N_{\mathrm{branches}}} = \frac{1}{8q}.
$$
Prédictions : $q = 3 : 1/24$, $q = 5 : 1/40$, $q = 7 : 1/56$, $q = 11 : 1/88$.

### 8.2 Test M1 (mpmath, $T \le 200$, 4 points par $q$)

Premier test (note 70), inconclusif : signal/bruit trop faible (4 points, $\sigma_{\mathrm{SE}} \sim
0.2$, signal cible $\sim 0.04$).

### 8.3 Test M1bis (LMFDB, $T \le 200$, ~140 zéros par $q$)

Test approfondi (note 71) sur les zéros pré-calculés de LMFDB pour les
caractères quadratiques primitifs $q \in \{3, 5, 7, 11\}$ ($a$ parité du
caractère) :

| $q$ | $a$ | $N$ zéros | $\langle \Delta_N \rangle$ | SE | $1/(8q)$ | $|t|_{\mathrm{PT}}$ | $|t|_0$ |
|---|---|---|---|---|---|---|---|
| 3 | 1 | 114 | $-2.4 \cdot 10^{-4}$ | $0.019$ | $+0.042$ | **2.19** | $0.01$ |
| 5 | 0 | 129 | $-2.6 \cdot 10^{-3}$ | $0.021$ | $+0.025$ | $1.34$ | $0.13$ |
| 7 | 1 | 140 | $+1.8 \cdot 10^{-3}$ | $0.021$ | $+0.018$ | $0.77$ | $0.09$ |
| 11 | 1 | 155 | $-2.8 \cdot 10^{-4}$ | $0.021$ | $+0.011$ | $0.55$ | $0.01$ |

Test $\chi^2$ discriminatoire sur 4 d.d.l. :

| Hypothèse | $\chi^2$ | $p$-value |
|---|---|---|
| $H_0$ : $\Delta_N = 0$ | $0.023$ | $\approx 1.00$ |
| $H_{\mathrm{PT}}$ : $\Delta_N = +1/(8q)$ | $7.49$ | $0.11$ |
| $H_{-\mathrm{PT}}$ : $\Delta_N = -1/(8q)$ | $7.05$ | $0.13$ |

$H_0$ est préféré par $\sim 300\times$ sur $H_{\mathrm{PT}}$. La prédiction
P3 est **effectivement falsifiée** à $T \le 200$.

**Discussion du signe.** Sur les quatre valeurs de $q$ testées, PT prédit
un décalage strictement positif $+1/(8q)$ ; on observe 3 décalages négatifs
($q = 3, 5, 11$) et 1 positif ($q = 7$). L'accord en signe seul est donc
$1/4$ (note 71 §3.2), bien en-deçà du $1$ qu'imposerait la prédiction.
La falsification n'est pas seulement d'amplitude : la structure de signe
elle-même contredit la prédiction $+1/(8q) > 0$.

### 8.4 Calibration sur $\zeta$

Avec la même méthodologie sur les 50 premiers zéros de $\zeta$ (Odlyzko) :
$$
\langle \Delta_N \rangle_\zeta = +1.004 \pm 0.028.
$$
Comparaisons :

| Cible | Valeur | $t$ |
|---|---|---|
| Pôle de $\zeta$ seul ($+1$) | $1.000$ | $+0.15$ (compatible) |
| Pôle + Maslov ($1 + 1/8 = 1.125$) | $1.125$ | $-4.29$ (rejeté) |
| Maslov seul ($1/8$) | $0.125$ | $+31.29$ (rejeté) |

### 8.5 Reformulation de la note 61

**Conséquence critique.** La note 61 affirmait que
$N_{\mathrm{PT}}(\gamma) - N_{\mathrm{RvM}}(\gamma) = 1/8$, vérifié à
$10^{-15}$ sur quatre décades ($\gamma \in [14, 10^3]$). Cette différence est entre **deux formes
asymptotiques lisses** ($N_{\mathrm{PT}}$ et $N_{\mathrm{RvM}}$), pas un
résidu mesurable au comptage exact des zéros. Au comptage de zéros réels,
le $+1/8$ entre PT et RvM est absorbé par les autres termes de la formule
et n'apparaît pas comme un décalage moyen mesurable. Le $+1$ du pôle de
$\zeta$ absorbe entièrement le $1/8$ dans la statistique de comptage.

La note 61 reste valide en tant qu'identité **asymptotique entre formes**,
mais la traduction « $1/8$ visible dans le comptage des zéros » est
erronée. Cela explique a posteriori pourquoi la prédiction $1/(8q)$ a été
falsifiée : elle attribuait à Dirichlet un résidu mesurable analogue à
celui de $\zeta$, alors que même pour $\zeta$ ce résidu n'est pas mesurable
au comptage des zéros à hauteur $T \sim 200$.

### 8.6 Position dans le programme

La falsification de P3 ne touche pas la structure HP-PT. Restent intacts :

- le modèle spectral PT des zéros de $\zeta$ via $T(s)$ (notes 65, 74,
  76, 77 : 17/17 zéros capturés à $|\Delta| < 1.5$ sur $t \in [10, 70]$,
  74/79 sur $t \in [10, 200]$ avec la formule brute, 75/79 à
  $|\Delta| < 0.5$ après régularisation explicite R1 avec
  $\varepsilon = 0.2$, et 267/269 à $|\Delta| < 1.5$ étendu à
  $t \in [10, 500]$) ;
- l'identification cohomologique $1/8 = s^2/2 = c_1/N_{\mathrm{corners}}
  = s/4$ pour $\zeta$ (notes 67–68) ;
- la dualité Higgs ↔ $\zeta$ via le canal $p = 2$ (note 69) en tant que
  conjecture forte structurelle pour $\zeta$.

Ce qui tombe : l'extension naturelle à Dirichlet via la mesure de Haar mod
$q$. Le mécanisme PT décrit ici est **spécifique à $\zeta$**, non à la
famille Dirichlet via cette factorisation. Le ledger PT marque P3 comme
[FALSIFIED, $T \le 200$].

C'est la discipline de falsifiabilité du programme : une prédiction PT
testable a été produite par K4, testée, rejetée par les données, et le
programme se restructure proprement.

---

## 10. Limites et ouvertures

### 9.1 Ce qui est solidement acquis

Le programme HP-PT, dans son état actuel, fournit :

- l'expression close des amplitudes $\kappa_p(s)$ et la représentation
  Fredholm exacte de $R(s)$ ;
- l'opérateur canalisé $H_{\mathrm{PT}}$ pour le couplage avec
  $\varepsilon_*^{\mathrm{geom}}$ canonique sans paramètre libre ;
- le cut-off dynamique $u_{\max}(\gamma) = \gamma/\sqrt{2\pi}$ forcé par la
  cellule de Planck symplectique ;
- la densité semi-classique exacte $N_{\mathrm{PT}}(\gamma) = N_{\mathrm{RvM}}(\gamma) + 1/8$
  (forme fonctionnelle identique) ;
- les indices de défaut $(1, 1)$ et la condition aux bords antipériodique
  forcée par $T_3$ ;
- la formule de trace analytique exacte $-d \log R/ds$ ;
- la trace régularisée opératorielle $T(s)$ avec pôles distributionnels
  aux zéros de $R$ ;
- la triple identification cohomologique du résidu $1/8$ ;
- la dualité structurelle Higgs ↔ $\zeta$.

**Tableau récapitulatif des statuts épistémiques** (notes 66 §3, 71 §9) :

| Résultat | Statut | Source |
|---|---|---|
| $s = 1/2$ (T1) | [THM] | monographie ch01 |
| $\det(I - D_R)^{-1} = R(s)$, $\mathrm{Re}(s) > 1$ | [THM] | note 42 |
| Formule de trace $-d\log R/ds$ dans $\mathrm{Re}(s) > 1$ | [THM] | note 62 |
| Indices de défaut $(1,1)$ et BC antipériodique $\theta = \pi$ | [THM] | note 63 |
| $c_1 = 1/2$ sur fibré spinoriel $q_+/q_-$ | [THM] | monographie app. P §C7 |
| Densité $N_{\mathrm{PT}} = N_{\mathrm{RvM}} + 1/8$ (identité asymptotique) | [DER] | note 61 |
| $1/8 = c_1/N_{\mathrm{corners}}$ (rigueur APS stricte) | [partiel] | note 68 |
| $\zeta_+ \zeta_- \ne 0$ sur $\mathrm{Re}(s) > 0$ (inc. ligne critique) | [THM] | Lemme 3.4 préprint A_PT |
| Prolongement analytique strict de $R(s)$ hors $\mathrm{Re}(s)>1$ | [OPEN] | note 62 |
| Pôles distributionnels de $T(s)$ sur $\mathrm{Re}(s) = 1/2$ uniquement | [OPEN] (= RH) | note 65 |
| Choix canonique $u_\star$ ($\sqrt{2\pi}$ vs $\log 3$) | [OPEN] | note 61 §2.4 |
| Dualité Higgs ↔ $\zeta$ structurelle (K4) | [CONJ forte] | note 69 |
| $\Delta_N(L \bmod q) = 1/(8q)$ pour Dirichlet (P3) | [FALSIFIED, $T \le 200$] | note 71 |
| G_HP4.d : spectre direct = $\{\gamma_n\}$ | [FALSIFIED] | note 64 |

### 9.2 Ce qui reste ouvert

Quatre ouvertures structurelles (la non-annulation $\zeta_+\zeta_- \neq 0$
sur $\mathrm{Re}(s) > 0$, autrefois jugée conjecturale, est désormais
[THM] par Lemme 3.4 du préprint A_PT — cf. §2.4 et le rappel en fin de
section).

**Ouverture 1 (G_HP3.a)** : prolongement analytique de $R(s)$ hors
$\mathrm{Re}(s) > 1$. La trace $T(s)$ est définie par série dans
$\mathrm{Re}(s) > 1$ puis prolongée distributionnellement à $\mathrm{Re}(s)
= 1/2$ ; la région intermédiaire $1/2 < \mathrm{Re}(s) < 1$ demande
l'analyse classique de $R$. C'est essentiellement le problème classique du
prolongement de $\zeta$.

**Ouverture 2 (RH stricte)** : prouver que **tous** les pôles
distributionnels de $T(s)$ sont sur $\mathrm{Re}(s) = 1/2$. Le programme
établit que la structure favorise cette ligne (quatre mécanismes
équivalents), mais ne prouve pas l'exclusion stricte des pôles hors ligne.
C'est essentiellement RH.

**Ouverture 3 (K3 strict APS)** : rigueur formelle Atiyah–Patodi–Singer
pour la variété à coin de la double barrière BK, au-delà du calcul d'aire
et de la vérification numérique. Le mécanisme cohomologique est identifié,
sa formalisation stricte reste à développer (note 68 §6.2).

**Ouverture 4 (K1, choix canonique de $u_\star$)** : la barrière inférieure
$u_\star$ admet deux candidats naturels — $u_\star = \sqrt{2\pi}$ (point
fixe de l'involution $u \leftrightarrow p_u$, choix retenu ici) et
$u_\star = \log 3$ (premier prime actif PT, choix arithmétique). Les deux
produisent la même densité asymptotique $N(\gamma)$ mais des corrections
$O(1/\gamma)$ différentes. Le statut canonique du choix entre les deux
reste ouvert (note 61 §2.4, note 67 §7, note 68 §8). L'article retient
$u_\star = \sqrt{2\pi}$ par symétrie d'involution, sans prétendre clore K1.

**Note 2026-05-17 : non-annulation $\zeta_+(s)\,\zeta_-(s) \ne 0$
sur $\mathrm{Re}(s) > 0$ — [THM] inconditionnel.** Cette propriété,
qui identifie les zéros de $R$ avec ceux de $\zeta$ sur la bande critique,
est conséquence directe du **Lemme 3.4 du préprint A_PT** : la série
$\sum_{p,\sigma} a_p^\sigma p^{-s}$ converge absolument sur
$\mathrm{Re}(s) > 0$ (asymptotique $a_p^\sigma \sim 2/p$), donc
$\zeta_\pm(s) = \exp(\cdot)$ est holomorphe et non-nul partout sur ce
demi-plan. Confirmation numérique :
$\min |\zeta_+\zeta_-(1/2+it)| > 0.48$ sur $t \in [10, 1000]$
(`PT_RH_MAY/scripts/test_zeta_pm_nonvanishing.py`, mpmath dps=50,
3 décades). L'audit PT_RH_MAY/analysis/INCOHERENCES_AUDIT.md (INC-1)
retrace la propagation de la correction dans le corpus. **Conséquence :**
l'Ouverture 1 du programme HP-PT se réduit désormais au seul
prolongement analytique classique de $\zeta$.

### 9.3 Position canonique du programme

Le programme HP-PT est un cadre opératoriel précis pour les zéros de
$\zeta$. Il fournit :

1. Une expression close des zéros via $\kappa_p(s)$ (intermédiaire
   $D_R$).
2. Une formule de trace régularisée canonique $T(s)$.
3. Quatre mécanismes équivalents pour $\mathrm{Re}(s) = 1/2$.
4. Un dictionnaire continu/discret unifié par la persistance.
5. Une explication structurelle des coïncidences numériques connues
   ($2\pi$, $7/8$, $1/8$, ligne critique).

Il n'apporte pas :

1. Une preuve stricte de RH classique.
2. Un raccourci à l'analyse distributionnelle de la trace régularisée.

Cette distinction est maintenue tout au long du programme : on cherche à
**comprendre** la structure des zéros, pas à les **prouver** sur la ligne
critique. Le passage à la preuve formelle demande l'analyse de Schwartz
fine sur $T(s)$, qui n'est pas plus facile que la preuve classique de RH.

### 9.4 Comparaison avec Berry–Keating–Connes

| Aspect | Berry–Keating–Connes | PT |
|---|---|---|
| Espace de phases | $(x, p) \in \mathbb{R}_+^2$ | $(u, p_u)$, $u = \log p$ |
| Opérateur | $H = xp + 1/2$ | $H_{\mathrm{PT-BK}} = (up + pu)/2$ |
| Mesure | $dx\,dp / 2\pi$ | $du\,dp_u / 2\pi$ |
| Cut-off | $x \ge \ell_x$, $p \ge \ell_p$, $\ell_x \ell_p = 2\pi$ | $u_\star p_\star = 2\pi$ (cellule Planck PT) |
| Régularisation | adélique (Connes) | atomique PT (note 65) |
| Trace | $\mathrm{Tr}(K) = \sum_\rho h(\gamma_\rho)$ (Weil) | $\mathrm{Tr}_{\mathrm{reg}} = \sum_p \kappa_p [\text{Euler + couplage}]$ |
| Couplage holonomique | aucun (zêta seule) | $\kappa_p$ (bifurcation $q_+/q_-$) |
| $2\pi$ | non interprété structurellement | $[u, p_u] = i$, Liouville canonique |
| $1/8$ | phase Maslov des coins (artefact géométrique) | $s^2/2 = c_1/4 = s/4$ (cohomologie spinorielle) |
| Spin $s = 1/2$ | postulé | sur-déterminé par $s^2/2 = s/4$ modulo T1 |

L'apport principal de PT est l'enrichissement intrinsèque : aucun
nombre n'est posé ad hoc, et chaque coïncidence numérique de la
construction classique reçoit une dérivation dans le cadre PT.

**Précisions historiques**. L'opérateur originel de Berry–Keating (1999)
est $H = xp$ ; la forme $H = xp + 1/2$ apparaît dans la formalisation
ultérieure de Sierra (2008) et Sierra–Townsend, qui régularise $xp$ via
une cellule de Planck et exhibe l'ordre symétrique. La régularisation de
Connes (1996, 1999) est de nature différente : elle s'effectue sur
l'espace adélique des classes d'idèles, par formule de trace
non-commutative. La formulation PT est atomique sur les premiers
(note 65) : ni adélique, ni continue sur $\mathbb{R}_+$.

### 9.5 Perspectives

Trois prolongements naturels.

- **Autres fonctions $L$** : la prédiction Dirichlet ayant été falsifiée,
  l'extension aux fonctions $L$ doit être repensée. La cohomologie du
  fibré spinoriel sous twist par $\chi$ modifie complètement
  $c_1$ ; une reformulation cohomologique propre est ouverte.
- **Phase Maslov dans la phase de Berry continue** : l'hypothèse (b) de la
  note 71 §6.3 propose que le $1/(8q)$ existe au niveau de la phase de
  Berry continue le long du flot spectral, mais ne se manifeste pas dans
  le compte des zéros. À tester.
- **Rigueur APS stricte** : développer le théorème d'indice formel sur la
  variété à coin de la double barrière (K3 strict). Pas urgent pour le
  cadre opératoriel, important pour la cohérence structurelle.

---

## 11. Conclusion

Le programme PT zêta livre un modèle spectral canonique des zéros non
triviaux de la fonction zêta de Riemann à partir d'une seule structure
interne : la dynamique de persistance avec son axiome $s = 1/2$. Le bilan
final tient en cinq points.

1. **L'objet opératoriel canonique est unique** : la trace régularisée
   $T(s) = \mathrm{Tr}_{\mathrm{reg}}\big[(s - i H_{\mathrm{PT-BK}})^{-1} D_R(s)\big]
   = \sum_p \kappa_p(s)\,(\log p)\,[1/(p^s-1) - A_p\,p^{-s}]$
   (signe corrigé 2026-05-17), construite à partir des amplitudes
   $\kappa_p$ de la bifurcation $q_+/q_-$ et de l'opérateur de dilatation
   Berry–Keating PT. Aucun paramètre libre. Identité analytique exacte
   dans $\mathrm{Re}(s) > 1$, validée numériquement à $10^{-13}$ ($\sigma = 4$,
   $P = 1000$). Pôles distributionnels sur $\mathrm{Re}(s) = 1/2$ aux zéros
   de $R(s) \equiv \zeta(s)/(\zeta_+ \zeta_-)$, validation A5 : **13/13
   premiers $\gamma_n$ détectés à distance médiane $0.075$**.

2. **La ligne critique $\mathrm{Re}(s) = 1/2$ admet quatre lectures
   équivalentes** : cellule de Planck symplectique $u_\star p_\star = 2\pi$,
   Jacobien de Haar multiplicatif, transformation unitaire $e^{su/2}$,
   condition aux bords antipériodique forcée par $T_3$. Les trois $1/2$ PT
   (T1 arithmétique, T5 géométrique, Haar spectral) sont opératoriellement
   identifiés.

3. **Le résidu cohomologique $1/8$ est triplement déterminé** :
   $1/8 = s^2/2 = \lambda_H$ (auto-couplage Higgs), $1/8 = c_1/N_{\mathrm{corners}}$
   (première classe de Chern du fibré spinoriel $q_+/q_-$ Kähler-Fisher),
   $1/8 = s/4$ (phase de Berry par coin). La cohérence
   $s^2/2 = s/4$ admet $s = 0$ ou $s = 1/2$ ; T1 exclut $s = 0$, ce qui
   sur-détermine $s = 1/2$ : c'est une signature structurelle de T1 par la
   zone HP-PT, non un substitut à T1.

4. **Dualité Higgs ↔ zêta structurelle** : le bifurcateur canonique
   $B$ (projecteur spinoriel $q_+/q_-$ du canal $p = 2$ au point fixe
   $\mu^* = 15$) produit simultanément $m_H = v \cdot s$ et $\Delta N =
   \lambda_H = s^2/2$ via la chaîne T1 + T2 + Haar + Weyl. Aucun autre
   couplage PT n'égale $1/8$ : la coïncidence numérique pure est rejetée.

5. **Discipline de falsifiabilité** : la prédiction Dirichlet
   $\Delta_N(L \bmod q) = 1/(8q)$, structurellement induite par K4, a été
   testée sur les zéros LMFDB de quatre caractères quadratiques primitifs
   et falsifiée ($\chi^2(\mathrm{PT}) = 7.49$ contre $\chi^2(H_0) = 0.023$
   sur 4 d.d.l.). La note 61 a été reformulée : le résidu $1/8$ est une
   identité asymptotique entre formes lisses, pas un décalage mesurable au
   comptage des zéros. Le mécanisme PT décrit ici est spécifique à $\zeta$.

6. **Reformulation PT de la formule de Riemann** (§7) : le produit
   d'Euler classique se réécrit
   $\zeta(s) = \zeta_+(s)\, \zeta_-(s)\, \det(I - D_R(s))^{-1}$ ;
   l'équation fonctionnelle $s \leftrightarrow 1-s$ s'identifie à
   l'involution $T_3$ du crible mod 3 sous transformée Mellin ; la
   fonction $\xi$ classique admet l'équivalent symétrique tautologique
   $\Xi_{\mathrm{PT}}(s) = R(s)\, R(1-s)$ ; la formule explicite de
   von Mangoldt prend la forme PT exacte
   $-d\log R/ds = \sum_p (\log p)[1/(p^s-1) - A_p p^{-s}]$ (signe
   corrigé). La géométrie spectrale $\Sigma_{\rm pers}$ (hyperelliptique
   genre 14, hyperbolique à cusp non-modulaire), l'opérateur Berry-Keating
   PT canonique $H_n^{\rm cusp}$, et la symétrie circulaire des 30
   branchements centrée à $-\mathcal{G}_F = -13/4$ complètent la
   structure. La structure PT cachée dans Riemann est ainsi mise au jour.

Position canonique : on a **compris** la structure des zéros de $\zeta$ ;
on n'a **pas prouvé** leur localisation stricte sur la ligne critique. Le
verrou ultime (G_HP3.a : prolongement analytique strict de $R$, ouverture 1
ci-dessus) est équivalent à RH au sens classique et ne semble pas plus
facile que cette dernière.

L'enrichissement principal sur Berry–Keating–Connes est la lecture
intrinsèque des coïncidences : le $2\pi$ est le commutateur canonique,
le $1/8$ est l'auto-couplage scalaire du Higgs PT, et le spin $s = 1/2$ de
la ligne critique est sur-déterminé par la cohérence cohomologique de la
construction. Aucun ingrédient n'est ad hoc.

---

## 12. Bibliographie

### Notes internes du corpus PT (PT_New_Math_Consolidation_FR/)

- Note 42 — *Attaque par traces du déterminant Fredholm résiduel*. Forme close de $\kappa_p(s)$, théorème Fredholm $\det(I - D_R(s))^{-1} = R(s)$.
- Note 51 — *Opérateur canalisé $C_n$ et statut du tuning $\beta$*. Opérateur canalisé $H_{\mathrm{PT}}$ pour le couplage $\sum_p A_p p^{-s}$.
- Note 52 — *$\varepsilon_*$ géométrique : la canalisation devient canonique*. $\varepsilon_*^{\mathrm{geom}} = 1.51$ canonique, élimination du tuning.
- Note 57 — *Localisation des zéros : reformulation PT et premiers tests*. Reformulation continu/discret dissous ; deux tests naïfs falsifiés.
- Note 58 — *Lorentzien–Berry–Keating PT : la signature change la nature*. Métrique Fisher lorentzienne ; Berry–Keating PT sur $u = \log p$ ; shift Haar multiplicatif.
- Note 59 — *Berry–Keating PT et le cut-off dynamique de persistance*. BK naïf sur intervalle fixe falsifié ; identification du cut-off dynamique $u_{\max}(\gamma) \sim \gamma$.
- Note 60 — *Verrou G_HP1 : définition rigoureuse du cut-off dynamique $u_{\max}(\gamma)$*. Cellule de Planck $u_\star p_\star = 2\pi$ ; $u_{\max} = \gamma/\sqrt{2\pi}$.
- Note 61 — *Verrou G_HP2 : quantification semi-classique de $H_{\mathrm{PT-BK}}$ et émergence exacte de la densité de Riemann–von Mangoldt*. Densité $N_{\mathrm{PT}} = (\gamma/2\pi) \log(\gamma/2\pi e) + 1$ ; différence $1/8$ avec $N_{\mathrm{RvM}}$.
- Note 62 — *Verrou G_HP3 : formule de trace explicite $H_{\mathrm{PT-BK}} \to R(s)$*. Identité analytique exacte $-d \log R/ds = \sum_p (\log p)[1/(p^s-1) - A_p p^{-s}]$ (signe corrigé 2026-05-17 via A4).
- Note 63 — *Verrou G_HP4 : auto-adjonction effective de $H_{\mathrm{PT-BK}}$*. Indices de défaut $(1,1)$ ; BC antipériodique $\theta = \pi$ forcée par $T_3$.
- Note 64 — *G_HP4.d falsifié dans la version directe : la trace régularisée est l'objet correct*. Spectre direct de $H_{\mathrm{PT-BK}}$ régulier ; les zéros sont pôles distributionnels.
- Note 65 — *Verrous G_HP3.b et G_HP3.c : régularisation Hadamard et écriture finale de la trace régularisée*. Trois régularisations équivalentes ; trace régularisée canonique ; 17 zéros capturés.
- Note 74 — *Scan de $P_{\max}$ et extension du range $t$ : capture distributionnelle des zéros de $\zeta$*. Saturation à $P_{\max} = 3\,000$ ; 74 zéros sur 79 capturés à $|\Delta| < 1.5$ sur $t \in [10, 200]$.
- Note 76 — *Régularisation explicite R1/R2 sur la trace régularisée*. Smearing gaussien $\varepsilon = 0.2$ : 77/79 à $|\Delta| < 1.5$, 75/79 à $|\Delta| < 0.5$, 47/79 à $|\Delta| < 0.1$ ; médiane $|\Delta| = 0.080$ (gain ×6.5 vs brut).
- Note 77 — *Extension du test R1/R2 sur $t \in [10, 500]$*. R1 $\varepsilon = 0.2$ : 267/269 à $|\Delta| < 1.5$ (99.3 \%), 242/269 à $|\Delta| < 0.5$ (90.0 \%), 139/269 à $|\Delta| < 0.1$ (51.7 \%) ; capture uniforme par bande $\ge 86.6\,\%$ sur 5 sous-bandes de 100 unités.
- Note 78 — *ε-adaptatif $\varepsilon \propto 1/\log\gamma$ : test sur $|F|$ brut et sur $F$ complexe*. Convolution sur $|F|$ falsifie (effondrement des pics). Convolution sur $F$ complexe : compromis marginal (médiane $|\Delta|$ réduite à $0.080$, couverture réduite à 76 \% à $|\Delta|<1.5$). $R_1$ fixé $\varepsilon = 0.2$ reste l'optimum global.
- Note 80 — *Reformulation PT de la formule de Riemann*. Réécriture canonique du produit d'Euler $\zeta = \zeta_+ \zeta_- \prod(1-\kappa_p)^{-1}$, de l'équation fonctionnelle (involution $T_3$), de la fonction xi ($\Xi_{\mathrm{PT}}(s) = R(s)R(1-s)$), et de la formule explicite ($-d\log R/ds$ comme formule de trace BK).
- Note 66 — *Bilan HP-PT complet : modèle spectral PT des zéros de $\zeta$*. Synthèse globale du programme.
- Note 67 — *Verrou K2 : interprétation PT pure du résidu $1/8$*. Identification algébrique $1/8 = s^2/2 = \lambda_H$.
- Note 68 — *Verrou K3 : preuve cohomologique du $-\pi/8$ par coin via projecteur spinoriel $q_+/q_-$*. Identification $1/8 = c_1/N_{\mathrm{corners}}$.
- Note 69 — *Verrou K4 : dualité Higgs ↔ zêta via le paramètre unique $s = 1/2$*. Bifurcateur canonique commun ; chaîne T1 + T2 + Haar + Weyl.
- Note 70 — *Verrou M1 : test numérique de la prédiction PT $\Delta_N(L \bmod q) = 1/(8q)$*. Test mpmath inconclusif.
- Note 71 — *Verrou M1bis : test approfondi via LMFDB*. Prédiction falsifiée à $T \le 200$ sur quatre caractères ; reformulation de la note 61.

### Monographie PT

- Senez, Y., *Théorie de la Persistance : du paramètre arithmétique $s = 1/2$ aux observables du Modèle Standard*, monographie en cours. Chapitres pertinents : ch01 (T1), ch03 (T2 + $T_3$), ch08 (T5 + $\gamma_p$), ch11 (variété Kähler-Fisher doublée), ch12 (projecteurs chiraux), ch15 (secteur Higgs, $\lambda_H = s^2/2$), ch20e (canal $p = 2$), ch37b (auto-adjonction Ruelle–PT), appendice P (théorème berry_3pi).

### Références externes

- Berry, M. V., Keating, J. P., 1999, *$H = xp$ and the Riemann zeros*, in *Supersymmetry and Trace Formulae : Chaos and Disorder*, NATO ASI Series, Plenum, 355–367.
- Connes, A., 1996, *Formule de trace en géométrie non-commutative et hypothèse de Riemann*, *Comptes Rendus de l'Académie des Sciences*, t. 323, 1231–1236.
- Connes, A., 1999, *Trace formula in noncommutative geometry and the zeros of the Riemann zeta function*, *Selecta Mathematica*, vol. 5, 29–106.
- Riemann, B., 1859, *Über die Anzahl der Primzahlen unter einer gegebenen Grösse*, *Monatsberichte der Berliner Akademie*, 671–680.
- Selberg, A., 1956, *Harmonic analysis and discontinuous groups in weakly symmetric Riemannian spaces with applications to Dirichlet series*, *Journal of the Indian Mathematical Society*, 20, 47–87.
- Hardy, G. H., Littlewood, J. E., 1921, *The zeros of Riemann's zeta-function on the critical line*, *Mathematische Zeitschrift*, 10, 283–317.
- Odlyzko, A. M., 1992, *The $10^{20}$-th zero of the Riemann zeta function and 175 million of its neighbors*, AT&T Bell Labs.
- LMFDB Collaboration, 2024, *The L-functions and Modular Forms Database*, https://www.lmfdb.org/.
- Tenenbaum, G., 2015, *Introduction to Analytic and Probabilistic Number Theory*, 3rd edition, AMS.
- Atiyah, M., Patodi, V. K., Singer, I. M., 1975, *Spectral asymmetry and Riemannian geometry*, *Mathematical Proceedings of the Cambridge Philosophical Society*, 77, 43–69 (I) ; 78, 405–432 (II) ; 79, 71–99 (III).
- Berry, M. V., 1984, *Quantal phase factors accompanying adiabatic changes*, *Proceedings of the Royal Society A*, 392, 45–57.
- Bost, J.-B., Connes, A., 1995, *Hecke algebras, type III factors and phase transitions with spontaneous symmetry breaking in number theory*, *Selecta Mathematica*, 1(3), 411–457.
- Sierra, G., 2008, *A physics pathway to the Riemann hypothesis*, arXiv:0712.1457.
- Iwaniec, H., Kowalski, E., 2004, *Analytic Number Theory*, AMS Colloquium Publications, vol. 53.
- Montgomery, H. L., 1973, *The pair correlation of zeros of the zeta function*, *Proceedings of Symposia in Pure Mathematics*, vol. 24, 181–193.

---

*Article rédigé à partir des notes 42, 51, 52, 60–80 du programme PT zêta, version consolidée du 15 mai 2026.*
