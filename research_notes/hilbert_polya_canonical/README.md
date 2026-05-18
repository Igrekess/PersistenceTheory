# PT_RH_HP_CANONICAL — Reformulation canonique de Hilbert-Polya via PT

**Auteur** : Yan Senez
**Date de création** : 2026-05-16
**Statut** : OUTLINE + DRAFT en cours

---

## Thèse

À partir d'un **unique théorème de théorie des nombres élémentaire** (T1 : transitions interdites mod 3 ⟹ $s = 1/2$) et d'**outils d'analyse complexe standards** (fonction zêta des premiers $P(s)$, théorèmes spectraux von Neumann, classes de Chern, phase de Berry), on construit canoniquement l'opérateur de Berry-Keating-PT $H_{\rm PT-BK}$ sur la coordonnée Mellin $u = \log p$. Tous les ingrédients usuellement postulés ad hoc dans le programme Hilbert-Polya (cellule de Planck $2\pi$, condition aux bords antipériodique, scattering matrix, résidu $1/8$) sont **dérivés** dans ce cadre. La conjecture de Riemann se ramène alors à une **unique** identification spectrale précise.

## Originalité

| Ingrédient Hilbert-Polya | Berry-Keating 1999 | Connes 1999 (adélique) | **PT (cet article)** |
|---|---|---|---|
| Espace de phases | $(x, p) \in \mathbb{R}_+^2$ | adélique | $(u, p_u)$, $u = \log p$ atomique |
| Opérateur | $H = xp$ (heuristique) | scaling adélique | $H_{\rm PT-BK} = -i(u\partial_u + 1/2)$ |
| Cellule de Planck $2\pi$ | postulée | inhérente au formalisme | **dérivée** ($[u, p_u] = i$ Liouville) |
| Cut-off $u_{\max}(\gamma)$ | postulé | régularisation adélique | **dérivé** ($\gamma/\sqrt{2\pi}$ par cellule) |
| Condition aux bords | non spécifiée | non applicable | $\theta = \pi$ **forcée par $T_3$** |
| Scattering matrix | non identifiée | adélique abstraite | $\phi_{\rm PT} = \zeta_+ \zeta_-$ **explicite** |
| Résidu $1/8$ phase Maslov | géométrique opaque | non interprété | $s^2/2 = c_1/N_{\rm corners} = s/4$ **cohomologique** |
| Spin $s = 1/2$ | postulé | postulé | **dérivé** depuis T1 arithmétique |

**Apport principal** : aucun ingrédient n'est postulé ad hoc. Tous sont conséquences de T1 (théorème classique) + outils classiques.

## Structure (~ 25-30 pages)

- **§1 Introduction** : Le programme Hilbert-Polya 1920 → Berry-Keating 1999 → Connes 1999. Ce qui manque : les constantes sont postulées. Annonce.
- **§2 Prérequis classiques** : T1 [THM], Conj A [THM nouveau], conventions $q_+, q_-$, fonctions $\zeta_\pm$.
- **§3 L'opérateur canonique $H_{\rm PT-BK}$** : construction, indices de défaut $(1,1)$, BC $\theta = \pi$ forcée par $T_3$.
- **§4 Quatre mécanismes équivalents pour Re(s) = 1/2** : M1 cellule Planck, M2 Haar, M3 transformation unitaire, M4 BC antipériodique. Convergence vers la mesure invariante unique.
- **§5 Sur-détermination cohomologique du résidu $1/8$** : $1/8 = s^2/2 = c_1/N_{\rm corners} = s/4$. Théorèmes classiques (Chern, Berry, Atiyah-Patodi-Singer formel).
- **§6 Construction du résidu et de la trace régularisée** : $R(s) = \zeta/(\zeta_+ \zeta_-)$, Fredholm, formule de trace explicite, trace régularisée $T(s)$.
- **§7 Reformulation rigoureuse de RH** : **théorème principal**.
- **§8 Comparaison avec Berry-Keating et Connes**.
- **§9 Validation numérique sommaire** (renvoi à PT_RH_VALIDATION pour détails).
- **§10 Limites et ouvertures** : honnêteté épistémique, voies fermées (R50), questions ouvertes.
- **§11 Conclusion**.
- **§12 Bibliographie**.

## Théorème principal (§7)

> **Théorème (Reformulation PT de RH)** : *Soient $T_3 = \mathrm{antidiag}(1,1)$ et soit $H_{\rm PT-BK}$ l'unique extension auto-adjointe sur $L^2([u_{\min}, u_{\max}], du)$ avec $u_\star p_\star = 2\pi$ (cellule de Planck), BC antipériodique $\theta = \pi$ (PT-canonique), forcée par l'antidiagonalité $T_3$. Soit $D_R(s) = \mathrm{diag}(\kappa_p(s))$ avec $\kappa_p(s) = 1 - (1 - p^{-s})\exp(A_p p^{-s})$, $A_p = \sin^2(\theta_p, q_+) + \sin^2(\theta_p, q_-)$.*
>
> *Alors la conjecture de Riemann est équivalente à :*
> $$\mathrm{Tr}_{\rm reg}\bigl[(s - i H_{\rm PT-BK})^{-1} D_R(s)\bigr] \text{ n'a de pôles distributionnels que sur } \mathrm{Re}(s) = 1/2.$$

## Ce que cet article N'EST PAS

- **Pas une preuve** de RH (le théorème ci-dessus est une **équivalence**, pas une preuve)
- Pas une réduction à un problème "strictement plus simple" (le verrou résiduel est équivalent à RH)
- Pas une généralisation aux L-Dirichlet (cf. PT_RH_DIRICHLET_NEG pour la falsification)
- Pas un résultat de validation numérique (cf. PT_RH_VALIDATION)
- Pas une exploration de la dualité Higgs ↔ ζ (cf. PT_HIGGS_ZETA)

## Public cible

- Théoriciens analytiques des nombres (RH, formules explicites)
- Harmonistes (théorie spectrale d'opérateurs non-bornés sur $L^2$ d'une variété)
- Mathématiciens physiciens (Berry-Keating, Hilbert-Polya operationnel)
- Communauté de théorie quantique des nombres

## Cibles éditoriales

- *Selecta Mathematica* (mathématiques pures, formulations canoniques)
- *Comm. Math. Phys.* (physique mathématique, BK et HP)
- *J. Number Theory* (théorie analytique des nombres)

## Liens

- Article principal en cours : `main_fr.tex` (à venir)
- Spinoffs : `PT_HIGGS_ZETA/`, `PT_RH_VALIDATION/`, `PT_RH_DIRICHLET_NEG/`
- Source corpus : `article_PT_zeta_FR.md` (à refondre), notes 42-80 du programme PT zêta
- Cycle parent : `PT_PROJECTS/PT_RH_HYPERBOLIC_CUSP/`
