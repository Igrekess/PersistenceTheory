# Audit monographie + site PT — 2026-04-30

## Portée

Audit effectué sur :

- `PT_MONOGRAPHY` : sources LaTeX FR/EN, PDFs, logs de compilation, ledgers, appendices de scripts.
- `PUBLIC/PersistenceTheory` : PDFs publics, README, site Astro FR/EN, pages de résultats, pages physiques et mathématiques.

Objectif : vérifier la cohérence structurelle de la PT publiée, la synchronisation FR/EN, la cohérence site/monographie, et repérer les points où un statut épistémique ou une nomenclature pourraient être trop forts ou obsolètes.

## Résultat global

La base est structurellement solide : les deux monographies compilées sont présentes, les PDFs publics ont les mêmes nombres de pages que les PDFs maîtres, le site compile, les pages FR/EN du site sont alignées, et les grands compteurs publics sont cohérents avec l’état actuel de la monographie.

L’audit ne révèle pas de contradiction centrale immédiate dans la vitrine actuelle. En revanche, il révèle plusieurs points de finition importants avant une publication “propre” : références LaTeX non résolues dans la version française, traces de vocabulaire `ghost/fantôme` dans la monographie, ledger de statut partiellement daté, et quelques statuts du site qui devaient être nuancés.

## Points vérifiés

### PDFs et publication

- `PT_MONOGRAPHY/main.pdf` : 885 pages.
- `PT_MONOGRAPHY/main_fr.pdf` : 913 pages.
- `PUBLIC/PersistenceTheory/TheTheoryOfPersistence.pdf` : 885 pages.
- `PUBLIC/PersistenceTheory/TheorieDeLaPersistance_FR.pdf` : 913 pages.

Conclusion : les PDFs publics correspondent bien aux PDFs maîtres actuels.

### Structure FR/EN de la monographie

Les deux fichiers `main.tex` et `main_fr.tex` incluent chacun 80 fichiers actifs, avec parité de structure après normalisation des dossiers `frontmatter`/`frontmatter_fr`, `chapters`/`chapters_fr`, `appendices`/`appendices_fr`.

Conclusion : la charpente FR/EN est synchronisée. Les divergences restantes sont dans le contenu et les références internes, pas dans la table active des chapitres.

### Site FR/EN

- Pages Astro FR et EN : alignées.
- Contenus `theorems` : 17 FR / 17 EN.
- Contenus `essays` : 8 FR / 8 EN.
- Build Astro : 232 pages générées, 0 erreur.
- Pagefind : index bilingue FR/EN généré.

Conclusion : le site est techniquement publiable. Les warnings restants sont des hints Astro/TypeScript non bloquants.

### Scripts publics

Le dépôt public expose maintenant une sélection plus propre des scripts compagnons, en particulier pour la gravité quantique/Kerr : les scripts exploratoires initiaux ont été retirés du set public, et les scripts conservés correspondent mieux à la version SOTA démonstrative.

État confirmé précédemment :

- `run_all.py ch13` : 15/15 PASS.
- `run_all.py --list` : 73 entrées visibles, avec 56 scripts canoniques décrits dans l’appendice F.

Conclusion : la partie scripts publics est cohérente avec l’idée “ne publier que ce qui démontre la version actuelle”.

## Corrections effectuées pendant l’audit

### README public

Correction du dernier résidu visible `ghost` :

- avant : `PT↔QED ghost dictionary`
- après : `PT↔QED echo-sector dictionary`

### Site : statuts épistémiques

Corrections FR/EN :

- Relativité : `Fisher devient espace-temps` n’est plus marqué comme `THM` pur, mais comme `BRIDGE/DER`.
- Relativité : `Signature et SO(3,1)` devient `THM/VAL`, ce qui reflète mieux la combinaison preuve structurelle + vérification.
- Temps : `Pas un métronome extérieur` devient `DER`, plus prudent que `THM`.
- Gravité quantique dans la page des statuts : `DER / COND` devient `DER / VAL / OUVERT`, avec mention de la calibration observationnelle encore active.

Ces corrections renforcent la crédibilité du site : elles évitent de transformer trop vite des ponts physiques en théorèmes purs.

## Points ouverts importants

### 1. Références LaTeX non résolues

La version anglaise ne liste pas de références individuelles non résolues dans le scan, même si le log final signale encore des undefined references. La version française contient en revanche une liste explicite de labels non résolus, notamment :

- `prop:g5_curvature`
- `chap:rh_fisher_projective`
- `chap:rh_echo_product`
- `chap:rh_basin_analytical`
- `chap:rh_adelic_closure`
- `prop:pnt_power_saving_circular`
- `eq:ch37b_a_p_asymp`
- `sec:residual_analysis`
- `sec:rh_chain`
- `def:ch37b_hppt`
- `chap:profile_decomposition`
- `sec:T1_forbidden`
- `chap:forbidden`
- `chap:classB_extensions`
- `chap:classB_integration`
- `chap:T3`
- `thm:LF_metric`
- `eq:clw_system`
- `eq:sign_wa`
- `thm:bianchi_backreaction`

Diagnostic : beaucoup de ces références semblent venir de sections RH/cosmologie/dark-energy qui pointent vers des labels anciens, des chapitres non inclus, ou des labels renommés. Ce n’est pas nécessairement une erreur mathématique, mais c’est une erreur de finition documentaire.

Priorité : haute avant diffusion large de la monographie française.

### 2. Nomenclature `ghost/fantôme`

Le site public ne contient plus de résidu visible `ghost/fantôme` après correction du README.

La monographie active contient encore de nombreuses occurrences de `ghost` / `fantôme`, notamment dans le glossaire, les appendices QFT/QG, et certains chapitres physiques. Certaines occurrences peuvent être légitimes dans le sens QFT standard, mais beaucoup semblent désigner l’ancien vocabulaire PT remplacé par `écho`.

Diagnostic : migration de nomenclature incomplète.

Recommandation : ne pas faire un remplacement global automatique. Il faut une passe contrôlée :

- garder les labels techniques si nécessaire pour ne pas casser les références ;
- remplacer la prose visible par `écho`, `secteur inactif`, `queue inactive`, ou `super-écho` selon le cas ;
- ajouter éventuellement une note de transition dans le glossaire : ancien terme `ghost`, terme PT officiel `echo/écho`.

Priorité : haute pour cohérence conceptuelle.

### 3. Ledger de statut daté

Les fichiers :

- `PT_MONOGRAPHY/frontmatter/status_ledger.tex`
- `PT_MONOGRAPHY/frontmatter_fr/status_ledger.tex`

contiennent encore des commentaires de tête du type `Persistence Theory v6` et `29-chapter table`, alors que la monographie active est beaucoup plus large.

Diagnostic : le ledger n’est pas au niveau des chapitres actuels. Il reste utile comme garde-fou, mais il ne doit pas être considéré comme source finale sans comparaison avec les chapitres.

Priorité : moyenne à haute. Le contenu de statut doit être ré-audité après correction des références.

### 4. Versioning interne encore visible

La monographie contient encore quelques traces de versioning interne :

- `As of v6.1` / `Depuis v6.1` dans le chapitre d’audit.
- `scripts_v7` dans les appendices de scripts.
- mentions historiques autour de `PT_RH v1.0` dans certaines sections RH.

Diagnostic : ce n’est pas bloquant si c’est historique ou technique, mais cela entre en tension avec la consigne de présenter PTC/PTP et les résultats dans leur état courant, sans versioning de développement.

Priorité : moyenne.

### 5. Appendix M / rapports d’exécution

L’appendice M reste un instantané de rapport d’exécution. Il peut être utile comme trace, mais il doit être clairement distingué de l’appendice F canonique, sinon le lecteur peut croire que deux états de scripts coexistent.

Priorité : moyenne.

## Solidité des bases PT

À ce stade, l’audit donne une lecture positive mais prudente.

Ce qui est fort :

- Le cœur GFT / capacité-information est présenté de façon stable.
- La distinction entre théorèmes, identités, dérivations, validations et prédictions est bien présente.
- Le site a commencé à refléter correctement la logique “continu natif, discret comme points persistants”.
- Les sections chimie, EI, EA, périodicité et scripts publics sont cohérentes avec la stratégie de démonstration : formules continues, points remarquables, validation numérique.
- La QG est maintenant mieux présentée comme chaîne dérivée/validée avec une frontière observationnelle encore ouverte, ce qui est plus crédible que de la vendre comme entièrement close.

Ce qui doit rester sous contrôle :

- Ne pas promouvoir les ponts physiques en théorèmes purs.
- Ne pas laisser le ledger obsolète contredire les chapitres.
- Ne pas laisser `ghost/fantôme` brouiller la nomenclature officielle `écho`.
- Fermer les références LaTeX non résolues, surtout dans la version française.
- Garder les scripts publics centrés sur les preuves SOTA plutôt que sur les phases exploratoires.

Conclusion : la PT a une présentation de plus en plus robuste, et le site est déjà une bonne vitrine conceptuelle. Mais pour une version “publication sérieuse”, il faut encore une passe de propreté monographique : références, nomenclature, ledger, et suppression du versioning interne visible.

## Prochaine passe recommandée

1. Corriger les références LaTeX non résolues dans `main_fr.log`, puis reconstruire `main_fr.pdf`.
2. Lancer une migration contrôlée `ghost/fantôme` vers `écho`, avec exceptions explicites pour les usages QFT standard.
3. Mettre à jour les ledgers FR/EN pour qu’ils reflètent la monographie active.
4. Nettoyer le versioning interne visible lorsque ce n’est pas historiquement nécessaire.
5. Relancer le build complet monographie + site, puis refaire un scan de cohérence.

