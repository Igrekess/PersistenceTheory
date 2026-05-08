export type MathPageStatus = 'identity' | 'theorem' | 'derivation' | 'exploration' | 'tool';

export interface MathPageContent {
  title: string;
  shortTitle: string;
  kicker: string;
  description: string;
  status: MathPageStatus;
  l1: string[];
  l2: string[];
  l3: string[];
  formulas: string[];
  takeaways: string[];
  source: string;
}

export interface MathPage {
  slug: string;
  visual:
    | 'survivor-wheel'
    | 'survivor-gaps'
    | 'prime-channels'
    | 'sieve-dynamics'
    | 'gft-budget'
    | 'discrete-continuum'
    | 'crt-phase'
    | 'anomalous-dimensions'
    | 'euler-zeta'
    | 'prime-spiral'
    | 'hash-diffusion'
    | 'compression'
    | 'zkp'
    | 'theorem-atlas'
    | 'calculator';
  repositories?: Array<{ name: string; url?: string; status: 'published' | 'planned' }>;
  scripts?: Array<{
    file: string;
    fileEn?: string;
    title: string;
    titleEn?: string;
    summary: string;
    summaryEn?: string;
    expected?: string;
    expectedEn?: string;
  }>;
  fr: MathPageContent;
  en: MathPageContent;
}

export const statusLabels = {
  fr: {
    identity: 'identité',
    theorem: 'théorème',
    derivation: 'dérivation',
    exploration: 'exploration',
    tool: 'outil',
  },
  en: {
    identity: 'identity',
    theorem: 'theorem',
    derivation: 'derivation',
    exploration: 'exploration',
    tool: 'tool',
  },
} as const;

export const mathPages: MathPage[] = [
  {
    slug: 'mecanique-survivants',
    visual: 'survivor-wheel',
    scripts: [
      {
        file: 'survivor_mod30_demo.py',
        title: 'Survivants modulo 30',
        titleEn: 'Survivors modulo 30',
        summary: 'Calcule les résidus survivants modulo 2·3·5 et leurs gaps circulaires.',
        summaryEn: 'Computes the surviving residues modulo 2·3·5 and their circular gaps.',
        expected: 'survivants = [1, 7, 11, 13, 17, 19, 23, 29]',
        expectedEn: 'survivors = [1, 7, 11, 13, 17, 19, 23, 29]',
      },
    ],
    fr: {
      title: 'Mécanique des survivants',
      shortTitle: 'Survivants',
      kicker: 'Mathématique PT / persistance',
      description: 'Comment une mécanique continue de contraintes fait apparaître des points remarquables de persistance.',
      status: 'derivation',
      l1: [
        'La PT ne part pas d’une collection de points isolés. Elle part d’un champ de possibilités, puis regarde ce qui reste reconnaissable quand des contraintes agissent.',
        'Un survivant est un point qui continue à porter une distinction après filtrage. C’est exactement le sens intuitif de la persistance : une forme ne devient mathématiquement intéressante que si elle résiste à la dissipation.',
      ],
      l2: [
        'Dans le crible, les survivants sont les résidus qui ne tombent dans aucun canal éliminé. Dans une loi de probabilité, les survivants sont les régions où la distribution reste distinguable de l’uniforme.',
        'La mécanique générale est donc : espace continu de possibilités, contrainte, perte entropique, puis points de persistance. Les nombres premiers, les gaps, les phases cycliques et les canaux actifs sont des cas particuliers de cette même lecture.',
      ],
      l3: [
        'Techniquement, la PT lit un survivant comme un point où le budget GFT ne se dissipe pas entièrement dans $H(P)$, mais conserve une composante $D_{KL}(P\\|U_m)$.',
        'Le discret n’est donc pas introduit comme ontologie première : il est le lieu où le continu possède des points stationnaires, des seuils ou des résidus invariants sous les contraintes admissibles.',
      ],
      formulas: [
        String.raw`H_{\max}(m)=D_{KL}(P\|U_m)+H(P)`,
        String.raw`\text{survivant}=\text{point où }D_{KL}\text{ reste non nul sous contrainte}`,
      ],
      takeaways: [
        'Le discret est lu comme point remarquable du continu.',
        'Un survivant n’est pas une exception : c’est une trace stable.',
        'Cette page sert de porte d’entrée mathématique générale.',
      ],
      source: 'Monographie : ch01_sieve, ch04_gft, ch24_scope.',
    },
    en: {
      title: 'Mechanics of survivors',
      shortTitle: 'Survivors',
      kicker: 'PT mathematics / persistence',
      description: 'How a continuous mechanics of constraints makes remarkable persistence points appear.',
      status: 'derivation',
      l1: [
        'PT does not start from isolated points. It starts from a field of possibilities, then asks what remains recognizable when constraints act.',
        'A survivor is a point that still carries a distinction after filtering. This is the intuitive meaning of persistence: a form becomes mathematically meaningful when it resists dissipation.',
      ],
      l2: [
        'In the sieve, survivors are residues that fall into no eliminated channel. In a probability law, survivors are regions where the distribution remains distinguishable from uniformity.',
        'The general mechanics is therefore: continuous space of possibilities, constraint, entropic loss, then persistence points. Primes, gaps, cyclic phases, and active channels are special cases of the same reading.',
      ],
      l3: [
        'Technically, PT reads a survivor as a point where the GFT budget is not entirely dissipated into $H(P)$, but keeps a $D_{KL}(P\\|U_m)$ component.',
        'The discrete layer is not introduced as the first ontology: it is the place where the continuum has stationary points, thresholds, or invariant residues under admissible constraints.',
      ],
      formulas: [
        String.raw`H_{\max}(m)=D_{KL}(P\|U_m)+H(P)`,
        String.raw`\text{survivor}=\text{point where }D_{KL}\text{ remains nonzero under constraint}`,
      ],
      takeaways: [
        'The discrete layer is read as a remarkable point of the continuum.',
        'A survivor is not an exception: it is a stable trace.',
        'This page is the general mathematical entry point.',
      ],
      source: 'Monograph: ch01_sieve, ch04_gft, ch24_scope.',
    },
  },
  {
    slug: 'gaps-survivants',
    visual: 'survivor-gaps',
    scripts: [
      {
        file: 'survivor_mod30_demo.py',
        title: 'Gaps de survivants',
        titleEn: 'Survivor gaps',
        summary: 'Reproduit les gaps 6,4,2,4,2,4,6,2 du crible modulo 30.',
        summaryEn: 'Reproduces the 6,4,2,4,2,4,6,2 gaps of the sieve modulo 30.',
        expected: 'gaps = [6, 4, 2, 4, 2, 4, 6, 2]',
        expectedEn: 'gaps = [6, 4, 2, 4, 2, 4, 6, 2]',
      },
    ],
    fr: {
      title: 'Gaps premiers et gaps de survivants',
      shortTitle: 'Gaps',
      kicker: 'Mathématique PT / crible',
      description: 'Lire les écarts premiers comme un cas limite des écarts entre survivants du crible.',
      status: 'theorem',
      l1: [
        'Un gap premier est l’écart entre deux nombres premiers voisins. En PT, on commence plus simplement : un gap est l’écart entre deux positions qui ont survécu au crible.',
        'Quand on augmente la profondeur du crible jusqu’à la fenêtre racine, les survivants restants sont exactement 1 et les nombres premiers. Le problème des gaps premiers devient donc un problème de mécanique des survivants.',
      ],
      l2: [
        'Modulo un primoriel $M_A$, les résidus copremiers à $M_A$ forment une suite circulaire. Les différences successives donnent les gaps du crible.',
        'La réduction exacte aux premiers se produit à $y=\\lfloor\\sqrt{x}\\rfloor$ : tout composite $n\\le x$ possède un facteur premier $\\le\\sqrt{x}$ et ne peut donc plus survivre.',
      ],
      l3: [
        'Le cœur exact est la relation $S(x;\\lfloor\\sqrt{x}\\rfloor)=1+\\pi(x)-\\pi(\\lfloor\\sqrt{x}\\rfloor)$, où $S(x;y)$ compte les entiers non divisibles par les premiers $\\le y$.',
        'La partie encore ouverte n’est pas la réduction aux premiers, mais la loi fermée fine qui prédirait chaque écart $p_{n+1}-p_n$ sans connaître les premiers voisins.',
      ],
      formulas: [
        String.raw`\varphi(M_A)=M_A\prod_{p\in A}(1-1/p)`,
        String.raw`S(x;\lfloor\sqrt{x}\rfloor)=1+\pi(x)-\pi(\lfloor\sqrt{x}\rfloor)`,
      ],
      takeaways: [
        'Les gaps de crible sont exacts.',
        'La fenêtre racine réduit les survivants aux premiers.',
        'La loi fine des gaps premiers reste une frontière de recherche.',
      ],
      source: 'Monographie : ch01_sieve, ch07_convergence, partie Riemann.',
    },
    en: {
      title: 'Prime gaps and survivor gaps',
      shortTitle: 'Gaps',
      kicker: 'PT mathematics / sieve',
      description: 'Reading prime gaps as a limiting case of gaps between sieve survivors.',
      status: 'theorem',
      l1: [
        'A prime gap is the distance between two neighboring primes. PT starts more simply: a gap is the distance between two positions that survived the sieve.',
        'When the sieve depth reaches the square-root window, the remaining survivors are exactly 1 and the primes. Prime gaps become a survivor-mechanics problem.',
      ],
      l2: [
        'Modulo a primorial $M_A$, residues coprime to $M_A$ form a circular sequence. Successive differences give the sieve gaps.',
        'Exact reduction to primes occurs at $y=\\lfloor\\sqrt{x}\\rfloor$: every composite $n\\le x$ has a prime factor $\\le\\sqrt{x}$ and cannot survive.',
      ],
      l3: [
        'The exact core is $S(x;\\lfloor\\sqrt{x}\\rfloor)=1+\\pi(x)-\\pi(\\lfloor\\sqrt{x}\\rfloor)$, where $S(x;y)$ counts integers divisible by no prime $\\le y$.',
        'The open part is not the reduction to primes, but the fine closed law that would predict each gap $p_{n+1}-p_n$ without already knowing the neighboring primes.',
      ],
      formulas: [
        String.raw`\varphi(M_A)=M_A\prod_{p\in A}(1-1/p)`,
        String.raw`S(x;\lfloor\sqrt{x}\rfloor)=1+\pi(x)-\pi(\lfloor\sqrt{x}\rfloor)`,
      ],
      takeaways: [
        'Sieve gaps are exact.',
        'The square-root window reduces survivors to primes.',
        'The fine law of prime gaps remains a research frontier.',
      ],
      source: 'Monograph: ch01_sieve, ch07_convergence, Riemann part.',
    },
  },
  {
    slug: 'pourquoi-nombres-premiers',
    visual: 'prime-channels',
    scripts: [
      {
        file: 'crt_prime_channels_demo.py',
        title: 'Canaux premiers irréductibles',
        titleEn: 'Irreducible prime channels',
        summary: 'Montre que les résidus modulo 30 se factorisent en canaux 2, 3, 5.',
        summaryEn: 'Shows that residues modulo 30 factor into the 2, 3, and 5 channels.',
      },
    ],
    fr: {
      title: 'Pourquoi les nombres premiers ?',
      shortTitle: 'Pourquoi les premiers',
      kicker: 'Mathématique PT / irréductibilité',
      description: 'Pourquoi les premiers apparaissent comme canaux irréductibles de persistance.',
      status: 'theorem',
      l1: [
        'En PT, un premier est un point discret de résonance du crible : il indexe un mode de persistance. La phase cyclique associée à ce premier donne une amplitude, et cette amplitude détermine son statut : actif, frontière, écho ou super-écho.',
        'Les nombres premiers sont les briques multiplicatives des entiers. En PT, ce n’est pas seulement une propriété arithmétique : c’est une propriété de canal.',
        'Un nombre composé mélange plusieurs contraintes. Un premier porte une contrainte irréductible : il ouvre un cycle propre que rien de plus petit ne peut factoriser.',
      ],
      l2: [
        'Le crible retire les multiples d’un canal $p$. Si $p$ est composé, ce canal est déjà expliqué par des canaux plus petits. Seul un premier ajoute une nouvelle obstruction indépendante.',
        'C’est pourquoi la PT lit les premiers comme les directions élémentaires de persistance multiplicative : chaque $p$ introduit un tore cyclique $\\mathbb{Z}/p\\mathbb{Z}$ non réductible.',
      ],
      l3: [
        'Le théorème chinois des restes factorise $\\mathbb{Z}/m\\mathbb{Z}$ en somme de composantes primaires quand les facteurs sont copremiers. Cette factorisation rend les premiers naturels, pas décoratifs.',
        'La lecture PT est que les premiers sont les canaux de phase cyclique minimaux compatibles avec la décomposition CRT et la conservation GFT.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/(pq)\mathbb{Z}\cong\mathbb{Z}/p\mathbb{Z}\oplus\mathbb{Z}/q\mathbb{Z}\quad(p,q\ \text{copremiers})`,
        String.raw`p\ \text{premier}\Rightarrow \mathbb{Z}/p\mathbb{Z}\ \text{canal irréductible}`,
      ],
      takeaways: [
        'Un premier indexe un mode discret de persistance.',
        'Un premier ajoute une contrainte indépendante.',
        'Un composé hérite de contraintes déjà présentes.',
        'La PT transforme les premiers en canaux de persistance.',
      ],
      source: 'Monographie : ch01_sieve, ch05_geometry, ch06_holonomy.',
    },
    en: {
      title: 'Why prime numbers?',
      shortTitle: 'Why primes',
      kicker: 'PT mathematics / irreducibility',
      description: 'Why primes appear as irreducible channels of persistence.',
      status: 'theorem',
      l1: [
        'In PT, a prime is a discrete resonance point of the sieve: it indexes a persistence mode. The cyclic phase associated with that prime gives an amplitude, and that amplitude determines its status: active, boundary, echo, or super-echo.',
        'Prime numbers are the multiplicative building blocks of the integers. In PT, this is not only arithmetic: it is a channel property.',
        'A composite number mixes several constraints. A prime carries an irreducible constraint: it opens its own cycle, which nothing smaller can factor.',
      ],
      l2: [
        'The sieve removes multiples of a channel $p$. If $p$ is composite, that channel is already explained by smaller channels. Only a prime adds a new independent obstruction.',
        'PT therefore reads primes as elementary directions of multiplicative persistence: each $p$ introduces an irreducible cyclic torus $\\mathbb{Z}/p\\mathbb{Z}$.',
      ],
      l3: [
        'The Chinese Remainder Theorem factorizes $\\mathbb{Z}/m\\mathbb{Z}$ into primary components when the factors are coprime. This makes primes structural, not decorative.',
        'In PT, primes are the minimal cyclic-phase channels compatible with CRT decomposition and GFT conservation.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/(pq)\mathbb{Z}\cong\mathbb{Z}/p\mathbb{Z}\oplus\mathbb{Z}/q\mathbb{Z}\quad(p,q\ \text{coprime})`,
        String.raw`p\ \text{prime}\Rightarrow \mathbb{Z}/p\mathbb{Z}\ \text{irreducible channel}`,
      ],
      takeaways: [
        'A prime indexes a discrete persistence mode.',
        'A prime adds an independent constraint.',
        'A composite inherits constraints already present.',
        'PT turns primes into persistence channels.',
      ],
      source: 'Monograph: ch01_sieve, ch05_geometry, ch06_holonomy.',
    },
  },
  {
    slug: 'crible-comme-dynamique',
    visual: 'sieve-dynamics',
    scripts: [
      {
        file: 'sieve_dynamics_demo.py',
        title: 'Dynamique du crible',
        titleEn: 'Sieve dynamics',
        summary: 'Suit la densité de survivants quand on ajoute les contraintes 2, 3, 5, 7.',
        summaryEn: 'Tracks survivor density as the 2, 3, 5, and 7 constraints are added.',
      },
    ],
    fr: {
      title: 'Le crible comme dynamique',
      shortTitle: 'Crible dynamique',
      kicker: 'Mathématique PT / filtration',
      description: 'Lire le crible non comme un simple algorithme, mais comme une dynamique de filtration.',
      status: 'theorem',
      l1: [
        'Le crible d’Ératosthène semble être une méthode pour barrer des multiples. En PT, il devient un laboratoire minimal : une contrainte agit, une partie disparaît, une autre persiste.',
        'C’est important parce qu’une théorie de la persistance a besoin d’un objet simple où voir naître la différence entre bruit, perte, résidu et structure.',
      ],
      l2: [
        'Chaque nouveau premier modifie l’espace des survivants. La densité baisse, les gaps se recomposent, mais les résidus qui restent gardent une structure transportable.',
        'Le crible est donc la lecture discrète d’une mécanique continue sous contrainte : appliquer une contrainte, voir quelles traces restent stables, recommencer.',
      ],
      l3: [
        'La dynamique exacte s’écrit par récurrence de type Legendre/Buchstab : $\\Phi(x,a)=\\Phi(x,a-1)-\\Phi(\\lfloor x/p_a\\rfloor,a-1)$.',
        'PT y ajoute la lecture informationnelle : chaque étape redistribue le budget entre entropie et persistance, au sens GFT.',
      ],
      formulas: [
        String.raw`\Phi(x,a)=\Phi(x,a-1)-\Phi(\lfloor x/p_a\rfloor,a-1)`,
        String.raw`\rho_A=\prod_{p\in A}(1-1/p)`,
      ],
      takeaways: [
        'Le crible est une dynamique de perte et de survivance.',
        'Les gaps se transportent quand la profondeur change.',
        'Le même schéma réapparaît dans les ponts physiques.',
      ],
      source: 'Monographie : ch01_sieve, ch07_convergence.',
    },
    en: {
      title: 'The sieve as a dynamics',
      shortTitle: 'Dynamic sieve',
      kicker: 'PT mathematics / filtration',
      description: 'Reading the sieve not as a mere algorithm, but as a filtration dynamics.',
      status: 'theorem',
      l1: [
        'The sieve of Eratosthenes looks like a method for crossing out multiples. In PT, it becomes a minimal laboratory: a constraint acts, one part disappears, another persists.',
        'This matters because a theory of persistence needs a simple object where one can see the difference between noise, loss, residue, and structure.',
      ],
      l2: [
        'Each new prime modifies the survivor space. Density falls, gaps recombine, but the remaining residues keep a transportable structure.',
        'The sieve is therefore the discrete reading of a continuous mechanics under constraint: apply a constraint, see which traces remain stable, repeat.',
      ],
      l3: [
        'The exact dynamics is written through a Legendre/Buchstab-type recurrence: $\\Phi(x,a)=\\Phi(x,a-1)-\\Phi(\\lfloor x/p_a\\rfloor,a-1)$.',
        'PT adds the informational reading: each step redistributes the budget between entropy and persistence in the GFT sense.',
      ],
      formulas: [
        String.raw`\Phi(x,a)=\Phi(x,a-1)-\Phi(\lfloor x/p_a\rfloor,a-1)`,
        String.raw`\rho_A=\prod_{p\in A}(1-1/p)`,
      ],
      takeaways: [
        'The sieve is a dynamics of loss and survival.',
        'Gaps are transported as depth changes.',
        'The same scheme reappears in physical bridges.',
      ],
      source: 'Monograph: ch01_sieve, ch07_convergence.',
    },
  },
  {
    slug: 'gft-premier-principe-mathematique',
    visual: 'gft-budget',
    scripts: [
      {
        file: 'gft_identity_demo.py',
        title: 'Identité GFT',
        titleEn: 'GFT identity',
        summary: 'Vérifie numériquement log2(m)=D_KL+H sur plusieurs distributions.',
        summaryEn: 'Numerically checks log2(m)=D_KL+H on several distributions.',
        expected: 'résidu < 1e-12',
        expectedEn: 'residual < 1e-12',
      },
    ],
    fr: {
      title: 'GFT comme premier principe mathématique',
      shortTitle: 'GFT mathématique',
      kicker: 'Mathématique PT / identité',
      description: 'Comprendre $\\log_2(m)=D_{KL}+H$ comme conservation exacte du budget informationnel.',
      status: 'identity',
      l1: [
        'Imagine une boîte avec $m$ possibilités. Pour identifier une possibilité, il faut un certain budget de distinctions. La GFT dit que ce budget ne disparaît jamais.',
        'Il se partage en deux parts : ce qui reste structuré, donc informatif, et ce qui s’étale comme bruit ou incertitude. La somme des deux redonne toujours le budget total.',
      ],
      l2: [
        'La formule $\\log_2(m)=D_{KL}(P\\|U_m)+H(P)$ est une identité algébrique pour toute distribution $P$ sur $m$ états.',
        'Sa force PT ne vient pas du fait qu’elle serait difficile à prouver, mais du rôle qu’elle joue : elle devient le premier principe qui interdit de créer ou perdre arbitrairement de la persistance.',
      ],
      l3: [
        'Avec $U_m$ uniforme, $D_{KL}(P\\|U_m)=\\sum_i p_i\\log_2(mp_i)$ et $H(P)=-\\sum_i p_i\\log_2 p_i$. En les additionnant, les termes en $p_i\\log_2 p_i$ s’annulent.',
        'La monographie classe GFT-ID comme identité exacte. Les extensions structurelles, elles, doivent garder leur statut propre : dérivation, pont ou validation selon le cas.',
      ],
      formulas: [
        String.raw`\log_2(m)=D_{KL}(P\|U_m)+H(P)`,
        String.raw`D_{KL}(P\|U_m)=\sum_i p_i\log_2(mp_i)`,
        String.raw`H(P)=-\sum_i p_i\log_2 p_i`,
      ],
      takeaways: [
        'GFT est une identité, pas un fit.',
        'Elle sépare structure persistante et entropie.',
        'Elle fournit le premier principe de la persistance.',
      ],
      source: 'Monographie : ch04_gft, ch14_thermodynamics, notation GFT.',
    },
    en: {
      title: 'GFT as a mathematical first principle',
      shortTitle: 'Mathematical GFT',
      kicker: 'PT mathematics / identity',
      description: 'Understanding $\\log_2(m)=D_{KL}+H$ as exact conservation of the information budget.',
      status: 'identity',
      l1: [
        'Imagine a box with $m$ possibilities. To identify one possibility, a certain budget of distinctions is needed. GFT says this budget never disappears.',
        'It splits into two parts: what remains structured and informative, and what spreads as noise or uncertainty. The sum always gives back the total budget.',
      ],
      l2: [
        'The formula $\\log_2(m)=D_{KL}(P\\|U_m)+H(P)$ is an algebraic identity for any distribution $P$ on $m$ states.',
        'Its PT strength does not come from being hard to prove, but from the role it plays: it becomes the first principle forbidding arbitrary creation or loss of persistence.',
      ],
      l3: [
        'With $U_m$ uniform, $D_{KL}(P\\|U_m)=\\sum_i p_i\\log_2(mp_i)$ and $H(P)=-\\sum_i p_i\\log_2 p_i$. Adding them cancels the $p_i\\log_2 p_i$ terms.',
        'The monograph classifies GFT-ID as an exact identity. Structural extensions keep their own status: derivation, bridge, or validation depending on the case.',
      ],
      formulas: [
        String.raw`\log_2(m)=D_{KL}(P\|U_m)+H(P)`,
        String.raw`D_{KL}(P\|U_m)=\sum_i p_i\log_2(mp_i)`,
        String.raw`H(P)=-\sum_i p_i\log_2 p_i`,
      ],
      takeaways: [
        'GFT is an identity, not a fit.',
        'It separates persistent structure and entropy.',
        'It provides the first principle of persistence.',
      ],
      source: 'Monograph: ch04_gft, ch14_thermodynamics, GFT notation.',
    },
  },
  {
    slug: 'pont-discret-continu',
    visual: 'discrete-continuum',
    scripts: [
      {
        file: 'cyclic_phase_demo.py',
        title: 'Résidu discret, phase continue',
        titleEn: 'Discrete residue, continuous phase',
        summary: 'Transforme des canaux Z/pZ en phases angulaires continues.',
        summaryEn: 'Turns Z/pZ channels into continuous angular phases.',
      },
    ],
    fr: {
      title: 'Pont discret-continu',
      shortTitle: 'Discret-continu',
      kicker: 'Mathématique PT / dissolution',
      description: 'Pourquoi la PT ne dit pas simplement que le continu émerge du discret.',
      status: 'derivation',
      l1: [
        'La lecture la plus juste est subtile : en PT, le discret n’est pas le matériau brut dont le continu serait fabriqué. Le discret est l’endroit où une mécanique continue marque des points persistants.',
        'C’est comme une onde stationnaire : l’onde est continue, mais certains nœuds ou ventres deviennent remarquables. Ces points ne sont pas ajoutés à l’onde ; ils sont révélés par elle.',
      ],
      l2: [
        'Le crible donne des résidus discrets, mais les phases, les angles, les métriques de Fisher et les dimensions anomales sont des objets continus.',
        'La PT dissout donc la frontière : elle ne choisit pas entre discret et continu, elle montre comment les points discrets sont des points de persistance d’une mécanique continue.',
      ],
      l3: [
        'La monographie formule ce passage via CRT, phase cyclique, géométrie de Fisher et limite métrique. La structure discrète $\\mathbb{Z}/p\\mathbb{Z}$ possède déjà des caractères continus via son dual de Pontryagin.',
        'La phrase prudente est : la PT donne une co-description où le discret est le spectre remarquable du continu contraint, pas une approximation grossière du continu.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/p\mathbb{Z}\longrightarrow \widehat{\mathbb{Z}/p\mathbb{Z}}\longrightarrow \theta_p`,
        String.raw`g^F_{ij}=\mathbb{E}[\partial_i\log P\,\partial_j\log P]`,
      ],
      takeaways: [
        'Le discret est point de persistance du continu.',
        'La phase cyclique relie résidu et angle.',
        'La métrique de Fisher donne la lecture géométrique continue.',
      ],
      source: 'Monographie : préface, ch05_geometry, ch06_holonomy, ch13_relativity, ch24_scope.',
    },
    en: {
      title: 'Discrete-continuous bridge',
      shortTitle: 'Discrete-continuous',
      kicker: 'PT mathematics / dissolution',
      description: 'Why PT does not simply say that the continuum emerges from the discrete.',
      status: 'derivation',
      l1: [
        'The correct reading is subtle: in PT, the discrete layer is not the raw material from which the continuum is made. The discrete layer is where a continuous mechanics marks persistent points.',
        'It is like a standing wave: the wave is continuous, yet some nodes or antinodes become remarkable. Those points are not added to the wave; they are revealed by it.',
      ],
      l2: [
        'The sieve gives discrete residues, but phases, angles, Fisher metrics, and anomalous dimensions are continuous objects.',
        'PT therefore dissolves the boundary: it does not choose between discrete and continuous; it shows how discrete points are persistence points of a continuous mechanics.',
      ],
      l3: [
        'The monograph formulates this passage through CRT, cyclic phase, Fisher geometry, and metric limit. The discrete structure $\\mathbb{Z}/p\\mathbb{Z}$ already has continuous characters through Pontryagin duality.',
        'The careful statement is: PT gives a co-description where the discrete layer is the remarkable spectrum of the constrained continuum, not a crude approximation to it.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/p\mathbb{Z}\longrightarrow \widehat{\mathbb{Z}/p\mathbb{Z}}\longrightarrow \theta_p`,
        String.raw`g^F_{ij}=\mathbb{E}[\partial_i\log P\,\partial_j\log P]`,
      ],
      takeaways: [
        'The discrete layer is a persistence point of the continuum.',
        'Cyclic phase connects residue and angle.',
        'The Fisher metric gives the continuous geometric reading.',
      ],
      source: 'Monograph: preface, ch05_geometry, ch06_holonomy, ch13_relativity, ch24_scope.',
    },
  },
  {
    slug: 'crt-holonomie-phase-cyclique',
    visual: 'crt-phase',
    scripts: [
      {
        file: 'cyclic_phase_demo.py',
        title: 'Phase cyclique',
        titleEn: 'Cyclic phase',
        summary: 'Calcule δp et sin²θp pour quelques canaux premiers.',
        summaryEn: 'Computes delta_p and sin^2(theta_p) for a few prime channels.',
      },
    ],
    fr: {
      title: 'CRT, holonomie et phase cyclique',
      shortTitle: 'CRT et holonomie',
      kicker: 'Mathématique PT / phase',
      description: 'Comment CRT et phase cyclique forcent les produits de canaux.',
      status: 'theorem',
      l1: [
        'Quand plusieurs cycles indépendants coexistent, on peut lire chacun séparément puis recomposer l’ensemble. C’est l’intuition du théorème chinois des restes.',
        'La PT ajoute que chaque cycle porte une phase. La persistance d’un canal n’est donc pas seulement une case qui reste ouverte : c’est une orientation cyclique stable.',
      ],
      l2: [
        'CRT décompose un module composite en canaux premiers copremiers. Pontryagin transforme cette décomposition en produit de caractères.',
        'C’est ce mécanisme qui rend naturelles les formes produits de la PT, notamment dans les ponts où les canaux actifs $3,5,7$ contribuent multiplicativement.',
      ],
      l3: [
        'La phase cyclique est encodée par $\\theta_p$ avec $\\sin^2\\theta_p=\\delta_p(2-\\delta_p)$. Elle transforme une profondeur de canal en amplitude persistante.',
        'La monographie traite la forme produit BA5 comme verrouillée par CRT + Pontryagin sous les axiomes de pont concernés.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/(p_1p_2p_3)\mathbb{Z}\cong\mathbb{Z}/p_1\mathbb{Z}\oplus\mathbb{Z}/p_2\mathbb{Z}\oplus\mathbb{Z}/p_3\mathbb{Z}`,
        String.raw`\sin^2\theta_p=\delta_p(2-\delta_p)`,
        String.raw`\alpha_{\rm sieve}=\prod_{p\in\{3,5,7\}}\sin^2\theta_p`,
      ],
      takeaways: [
        'CRT donne l’indépendance des canaux.',
        'Pontryagin transforme l’indépendance en produit.',
        'La phase cyclique porte l’amplitude persistante.',
      ],
      source: 'Monographie : ch05_geometry, ch06_holonomy, ch09_bridge, BA5.',
    },
    en: {
      title: 'CRT, holonomy, and cyclic phase',
      shortTitle: 'CRT and holonomy',
      kicker: 'PT mathematics / phase',
      description: 'How CRT and cyclic phase force channel products.',
      status: 'theorem',
      l1: [
        'When several independent cycles coexist, each can be read separately and then recomposed. This is the intuition behind the Chinese Remainder Theorem.',
        'PT adds that each cycle carries a phase. A persistent channel is not only an open slot: it is a stable cyclic orientation.',
      ],
      l2: [
        'CRT decomposes a composite module into coprime prime channels. Pontryagin duality turns this decomposition into a product of characters.',
        'This mechanism makes PT product forms natural, especially in bridges where the active channels $3,5,7$ contribute multiplicatively.',
      ],
      l3: [
        'Cyclic phase is encoded by $\\theta_p$ with $\\sin^2\\theta_p=\\delta_p(2-\\delta_p)$. It turns channel depth into persistent amplitude.',
        'The monograph treats the BA5 product form as locked by CRT + Pontryagin under the relevant bridge axioms.',
      ],
      formulas: [
        String.raw`\mathbb{Z}/(p_1p_2p_3)\mathbb{Z}\cong\mathbb{Z}/p_1\mathbb{Z}\oplus\mathbb{Z}/p_2\mathbb{Z}\oplus\mathbb{Z}/p_3\mathbb{Z}`,
        String.raw`\sin^2\theta_p=\delta_p(2-\delta_p)`,
        String.raw`\alpha_{\rm sieve}=\prod_{p\in\{3,5,7\}}\sin^2\theta_p`,
      ],
      takeaways: [
        'CRT gives channel independence.',
        'Pontryagin turns independence into a product.',
        'Cyclic phase carries persistent amplitude.',
      ],
      source: 'Monograph: ch05_geometry, ch06_holonomy, ch09_bridge, BA5.',
    },
  },
  {
    slug: 'dimensions-anomales',
    visual: 'anomalous-dimensions',
    scripts: [
      {
        file: 'gamma_threshold_demo.py',
        title: 'Seuil γp = 1/2',
        titleEn: 'Threshold gamma_p = 1/2',
        summary: 'Illustre la coupure active/inactive sur 3,5,7 puis 11,13,17.',
        summaryEn: 'Illustrates the active/inactive cutoff on 3,5,7 then 11,13,17.',
      },
    ],
    fr: {
      title: 'Dimensions anomales',
      shortTitle: 'Dimensions anomales',
      kicker: 'Mathématique PT / seuils',
      description: 'Pourquoi $\\gamma_p$ mesure la sensibilité d’un canal et sélectionne les actifs.',
      status: 'derivation',
      l1: [
        'Une dimension anomale n’est pas une dimension spatiale cachée. C’est un exposant de sensibilité : il dit à quelle vitesse un canal réagit quand la profondeur change.',
        'Si cette sensibilité reste au-dessus du seuil $1/2$, le canal peut porter une direction active. Si elle passe sous le seuil, le canal devient écho.',
      ],
      l2: [
        'La PT lit les canaux premiers via une fonction $\\gamma_p$. Les premiers $3,5,7$ restent actifs ; à partir de $11$, la contribution tombe du côté inactif.',
        'Cela donne une explication mathématique compacte au fait que certaines structures se ferment à trois directions actives plutôt qu’à une infinité.',
      ],
      l3: [
        'Dans la monographie, le seuil $\\gamma_p=1/2$ est relié à la condition GFT par canal et à la symétrie fondamentale $s=1/2$.',
        'Le statut précis dépend du niveau : la définition est mathématique, l’identification physique des canaux actifs relève des ponts et dérivations associés.',
      ],
      formulas: [
        String.raw`\gamma_p>\frac12\Rightarrow\text{canal actif}`,
        String.raw`\gamma_p<\frac12\Rightarrow\text{écho inactif}`,
      ],
      takeaways: [
        'Une dimension anomale est un exposant de sensibilité.',
        '$1/2$ est le seuil de persistance active.',
        'La coupure $3,5,7$ puis $11$ devient lisible.',
      ],
      source: 'Monographie : ch06_holonomy, ch08_fixed_point, ch23_audit.',
    },
    en: {
      title: 'Anomalous dimensions',
      shortTitle: 'Anomalous dimensions',
      kicker: 'PT mathematics / thresholds',
      description: 'Why $\\gamma_p$ measures channel sensitivity and selects active channels.',
      status: 'derivation',
      l1: [
        'An anomalous dimension is not a hidden spatial dimension. It is a sensitivity exponent: it says how fast a channel reacts when depth changes.',
        'If that sensitivity stays above the $1/2$ threshold, the channel can carry an active direction. If it falls below the threshold, the channel becomes an echo.',
      ],
      l2: [
        'PT reads prime channels through a function $\\gamma_p$. The primes $3,5,7$ remain active; from $11$ onward, the contribution falls to the inactive side.',
        'This gives a compact mathematical explanation for why some structures close at three active directions rather than infinitely many.',
      ],
      l3: [
        'In the monograph, the threshold $\\gamma_p=1/2$ is tied to the per-channel GFT condition and to the fundamental symmetry $s=1/2$.',
        'The precise status depends on the level: the definition is mathematical, while the physical identification of active channels belongs to the associated bridges and derivations.',
      ],
      formulas: [
        String.raw`\gamma_p>\frac12\Rightarrow\text{active channel}`,
        String.raw`\gamma_p<\frac12\Rightarrow\text{inactive echo}`,
      ],
      takeaways: [
        'An anomalous dimension is a sensitivity exponent.',
        '$1/2$ is the active-persistence threshold.',
        'The $3,5,7$ then $11$ cutoff becomes readable.',
      ],
      source: 'Monograph: ch06_holonomy, ch08_fixed_point, ch23_audit.',
    },
  },
  {
    slug: 'riemann-zeta-lecture-pt',
    visual: 'euler-zeta',
    repositories: [
      { name: 'Igrekess/PT_Riemann', status: 'planned' },
      { name: 'Igrekess/PT_NT', status: 'planned' },
    ],
    fr: {
      title: 'Riemann et zêta en lecture PT',
      shortTitle: 'Riemann',
      kicker: 'Mathématique PT / programme',
      description: 'Présenter la lecture PT de Riemann comme programme de recherche, sans sur-vendre une preuve fermée.',
      status: 'exploration',
      l1: [
        'L’hypothèse de Riemann parle de l’ordre caché dans les nombres premiers. La PT s’y intéresse parce que les premiers sont précisément des survivants du crible.',
        'La question PT devient : la ligne critique $1/2$ est-elle la trace spectrale du même seuil de persistance que celui qui apparaît ailleurs dans la théorie ?',
      ],
      l2: [
        'La monographie contient un programme Riemann : reformuler la question via monotonie de Čencov, DPI, opérateurs spectraux et décompositions liées au crible.',
        'Il faut rester net : cette page ne doit pas annoncer une preuve standard de RH. Elle doit exposer le dictionnaire PT, les résultats partiels et les obstructions encore identifiées.',
      ],
      l3: [
        'Les chapitres Riemann explorent la voie où la condition spectrale $\\lambda\\ge s^2=1/4$ serait équivalente à RH sous hypothèses PT spécifiques.',
        'Le point fort public est la cohérence du seuil $1/2$ ; le point ouvert est l’assemblage complet sans hypothèse auxiliaire excessive.',
      ],
      formulas: [
        String.raw`\Re(s)=\frac12`,
        String.raw`\lambda\ge s^2=\frac14`,
        String.raw`\zeta(s)=0\Rightarrow \Re(s)=\frac12\quad\text{(RH, non revendiquée comme fermée ici)}`,
      ],
      takeaways: [
        'La page doit être attirante mais prudente.',
        'Le seuil $1/2$ relie crible, spectre et persistance.',
        'Le statut est programme de recherche, pas théorème public fermé.',
      ],
      source: 'Monographie : partie Riemann, ch21_RH_cencov_dpi, ch37_RH_bifurcation, ch37b_RH_proofs.',
    },
    en: {
      title: 'Riemann and zeta in PT reading',
      shortTitle: 'Riemann',
      kicker: 'PT mathematics / programme',
      description: 'Presenting the PT reading of Riemann as a research programme without overselling a closed proof.',
      status: 'exploration',
      l1: [
        'The Riemann Hypothesis is about hidden order in the prime numbers. PT cares because primes are precisely sieve survivors.',
        'The PT question becomes: is the critical line $1/2$ the spectral trace of the same persistence threshold that appears elsewhere in the theory?',
      ],
      l2: [
        'The monograph contains a Riemann programme: reformulating the question through Čencov monotonicity, DPI, spectral operators, and sieve-related decompositions.',
        'The page must stay clear: it should not announce a standard proof of RH. It should expose the PT dictionary, partial results, and remaining obstructions.',
      ],
      l3: [
        'The Riemann chapters explore a route where the spectral condition $\\lambda\\ge s^2=1/4$ would be equivalent to RH under specific PT hypotheses.',
        'The public strength is the coherence of the $1/2$ threshold; the open point is complete assembly without excessive auxiliary hypotheses.',
      ],
      formulas: [
        String.raw`\Re(s)=\frac12`,
        String.raw`\lambda\ge s^2=\frac14`,
        String.raw`\zeta(s)=0\Rightarrow \Re(s)=\frac12\quad\text{(RH, not claimed closed here)}`,
      ],
      takeaways: [
        'The page should be attractive but careful.',
        'The $1/2$ threshold links sieve, spectrum, and persistence.',
        'The status is research programme, not public closed theorem.',
      ],
      source: 'Monograph: Riemann part, ch21_RH_cencov_dpi, ch37_RH_bifurcation, ch37b_RH_proofs.',
    },
  },
  {
    slug: 'spirales-premieres',
    visual: 'prime-spiral',
    repositories: [
      { name: 'Igrekess/PT_SPIRALS', status: 'planned' },
    ],
    scripts: [
      {
        file: 'prime_spiral_demo.py',
        title: 'Spirale première légère',
        titleEn: 'Light prime spiral',
        summary: 'Produit les coordonnées de quelques points premiers sur une spirale d’Archimède.',
        summaryEn: 'Produces coordinates for a few prime points on an Archimedean spiral.',
      },
    ],
    fr: {
      title: 'Spirales premières',
      shortTitle: 'Spirales',
      kicker: 'Mathématique PT / visualisation',
      description: 'Utiliser les spirales d’Ulam, Sacks ou Archimède comme visualisation des survivants premiers.',
      status: 'exploration',
      l1: [
        'Les spirales de nombres premiers rendent visible une chose surprenante : les premiers ne ressemblent pas à un bruit uniforme.',
        'La PT peut les utiliser comme vitrine pédagogique : quand le crible agit, les survivants gardent des alignements, des bandes et des traces géométriques.',
      ],
      l2: [
        'Une spirale ne prouve pas une loi. Elle transforme une suite arithmétique en image et révèle des corrélations que l’œil saisit vite.',
        'Pour la section mathématique du site, c’est un excellent pont simple : on voit avant de calculer pourquoi les survivants ont une géométrie.',
      ],
      l3: [
        'Le projet PT_SPIRALS peut fournir des figures comparant Ulam, Sacks et Archimède, avec contrôles par classes modulaires.',
        'Le statut doit rester visuel et exploratoire sauf lorsqu’une statistique modulaire précise est calculée et référencée.',
      ],
      formulas: [
        String.raw`n\mapsto r(n)e^{i\theta(n)}`,
        String.raw`\text{classe modulaire}\rightarrow\text{alignement visuel}`,
      ],
      takeaways: [
        'Très bon support pour expliquer simplement.',
        'Ne pas confondre visualisation et preuve.',
        'Bon pont vers gaps, crible et classes modulaires.',
      ],
      source: 'Dépôt GitHub à publier : Igrekess/PT_SPIRALS ; monographie ch01_sieve.',
    },
    en: {
      title: 'Prime spirals',
      shortTitle: 'Spirals',
      kicker: 'PT mathematics / visualization',
      description: 'Using Ulam, Sacks, or Archimedean spirals as visualizations of prime survivors.',
      status: 'exploration',
      l1: [
        'Prime spirals make something surprising visible: primes do not look like uniform noise.',
        'PT can use them as a pedagogical showcase: when the sieve acts, survivors keep alignments, bands, and geometric traces.',
      ],
      l2: [
        'A spiral does not prove a law. It turns an arithmetic sequence into an image and reveals correlations the eye grasps quickly.',
        'For the mathematics section, it is an excellent plain-language bridge: one sees before calculating why survivors have geometry.',
      ],
      l3: [
        'The PT_SPIRALS project can provide figures comparing Ulam, Sacks, and Archimedean spirals, with checks by modular classes.',
        'The status should remain visual and exploratory unless a precise modular statistic is computed and referenced.',
      ],
      formulas: [
        String.raw`n\mapsto r(n)e^{i\theta(n)}`,
        String.raw`\text{modular class}\rightarrow\text{visual alignment}`,
      ],
      takeaways: [
        'Very strong for public explanation.',
        'Do not confuse visualization with proof.',
        'Good bridge toward gaps, sieve, and modular classes.',
      ],
      source: 'GitHub repository to publish: Igrekess/PT_SPIRALS; monograph ch01_sieve.',
    },
  },
  {
    slug: 'cryptographie-fonctions-sens-unique',
    visual: 'hash-diffusion',
    repositories: [
      { name: 'Igrekess/PT_Cryptographie', status: 'planned' },
      { name: 'Igrekess/sha256-sieve-structure', status: 'planned' },
    ],
    scripts: [
      {
        file: 'hash_kl_demo.py',
        title: 'Mini diagnostic hash/KL',
        titleEn: 'Mini hash/KL diagnostic',
        summary: 'Mesure la divergence à l’uniforme sur des préfixes SHA-256 jouets.',
        summaryEn: 'Measures divergence from uniformity on toy SHA-256 prefixes.',
      },
    ],
    fr: {
      title: 'Cryptographie et fonctions à sens unique',
      shortTitle: 'Cryptographie',
      kicker: 'Mathématique PT / information',
      description: 'Lire l’asymétrie facile/difficile comme perte contrôlée de structure persistante.',
      status: 'exploration',
      l1: [
        'Une fonction à sens unique est facile à calculer dans un sens, mais difficile à inverser. Intuitivement, elle garde le résultat tout en dispersant le chemin.',
        'La PT donne une langue naturelle pour cela : l’aller conserve assez de structure pour produire une sortie, mais l’inverse doit reconstruire une persistance perdue dans l’entropie.',
      ],
      l2: [
        'Dans une fonction de hachage, beaucoup d’entrées mènent à des sorties de taille fixe. Le budget d’information se contracte ; l’inversion demande de retrouver la structure effacée.',
        'La page doit rester une application conceptuelle, sauf pour les modules où les scripts mesurent effectivement des signatures de crible ou de structure.',
      ],
      l3: [
        'Les projets PT_Cryptographie et sha256-sieve-structure peuvent soutenir une page expérimentale : spectres de résidus, classes survivantes, mesures KL, collisions et tests de robustesse.',
        'Le statut canonique doit être prudent : lecture PT de l’asymétrie informationnelle, pas preuve générale de sécurité cryptographique.',
      ],
      formulas: [
        String.raw`f:x\mapsto y\quad\text{facile},\qquad f^{-1}(y)\quad\text{difficile}`,
        String.raw`\Delta D_{KL}\downarrow\Rightarrow\text{structure inverse moins accessible}`,
      ],
      takeaways: [
        'La cryptographie illustre très bien persistance vs entropie.',
        'Ne pas revendiquer de preuve de sécurité générale.',
        'Les scripts peuvent donner des démonstrations expérimentales.',
      ],
      source: 'Dépôts GitHub à publier : Igrekess/PT_Cryptographie, Igrekess/sha256-sieve-structure ; monographie ch04_gft.',
    },
    en: {
      title: 'Cryptography and one-way functions',
      shortTitle: 'Cryptography',
      kicker: 'PT mathematics / information',
      description: 'Reading easy/hard asymmetry as controlled loss of persistent structure.',
      status: 'exploration',
      l1: [
        'A one-way function is easy to compute forward, but hard to invert. Intuitively, it keeps the result while dispersing the path.',
        'PT gives a natural language for this: the forward map preserves enough structure to produce an output, while inversion must reconstruct persistence lost into entropy.',
      ],
      l2: [
        'In a hash function, many inputs lead to fixed-size outputs. The information budget contracts; inversion requires recovering erased structure.',
        'The page should remain a conceptual application unless scripts actually measure sieve or structural signatures.',
      ],
      l3: [
        'The PT_Cryptographie and sha256-sieve-structure projects can support an experimental page: residue spectra, survivor classes, KL measures, collisions, and robustness tests.',
        'The canonical status must remain careful: PT reading of informational asymmetry, not a general proof of cryptographic security.',
      ],
      formulas: [
        String.raw`f:x\mapsto y\quad\text{easy},\qquad f^{-1}(y)\quad\text{hard}`,
        String.raw`\Delta D_{KL}\downarrow\Rightarrow\text{inverse structure less accessible}`,
      ],
      takeaways: [
        'Cryptography illustrates persistence vs entropy very well.',
        'Do not claim a general security proof.',
        'Scripts can provide experimental demonstrations.',
      ],
      source: 'GitHub repositories to publish: Igrekess/PT_Cryptographie, Igrekess/sha256-sieve-structure; monograph ch04_gft.',
    },
  },
  {
    slug: 'compression-information',
    visual: 'compression',
    repositories: [
      { name: 'Igrekess/pt-compress', status: 'planned' },
    ],
    scripts: [
      {
        file: 'compression_gft_demo.py',
        title: 'Compression et GFT',
        titleEn: 'Compression and GFT',
        summary: 'Compare une chaîne redondante et une chaîne pseudo-aléatoire via entropie empirique.',
        summaryEn: 'Compares a redundant string and a pseudo-random string through empirical entropy.',
      },
    ],
    fr: {
      title: 'Compression et information',
      shortTitle: 'Compression',
      kicker: 'Mathématique PT / GFT',
      description: 'Compresser comme extraire ce qui persiste et rejeter ce qui est entropique.',
      status: 'exploration',
      l1: [
        'Compresser un fichier, c’est enlever ce qui se répète ou ce qui n’aide pas à reconstruire l’essentiel. En langage PT : on cherche ce qui persiste.',
        'C’est une des meilleures portes d’entrée vers GFT : l’information structurée se garde, le bruit coûte cher, et le budget total impose une limite.',
      ],
      l2: [
        'Une compression efficace augmente la part exploitable de structure relativement à une représentation brute. Elle ne crée pas d’information ; elle réorganise le budget.',
        'La PT peut présenter la compression comme un cas concret de la partition $D_{KL}+H$ : la structure repérable contre l’uniforme d’un côté, l’entropie irréductible de l’autre.',
      ],
      l3: [
        'Le projet pt-compress peut servir de laboratoire : mesurer entropie, redondance, divergence à l’uniforme, et coût de reconstruction.',
        'Le point mathématique canonique reste GFT. Les performances d’un compresseur particulier relèvent d’une validation expérimentale.',
      ],
      formulas: [
        String.raw`\text{budget brut}=\text{structure compressible}+\text{résidu entropique}`,
        String.raw`\log_2(m)=D_{KL}+H`,
      ],
      takeaways: [
        'Compresser, c’est extraire la persistance.',
        'GFT donne le langage du budget.',
        'Très bon pont pédagogique vers l’information.',
      ],
      source: 'Dépôt GitHub à publier : Igrekess/pt-compress ; monographie ch04_gft, ch_PM.',
    },
    en: {
      title: 'Compression and information',
      shortTitle: 'Compression',
      kicker: 'PT mathematics / GFT',
      description: 'Compression as extracting what persists and rejecting what is entropic.',
      status: 'exploration',
      l1: [
        'Compressing a file means removing what repeats or what does not help reconstruct the essential part. In PT language: one looks for what persists.',
        'This is one of the best entry points into GFT: structured information is kept, noise is costly, and the total budget imposes a limit.',
      ],
      l2: [
        'Efficient compression increases the usable share of structure relative to a raw representation. It does not create information; it reorganizes the budget.',
        'PT can present compression as a concrete case of the $D_{KL}+H$ partition: structure detectable against uniformity on one side, irreducible entropy on the other.',
      ],
      l3: [
        'The pt-compress project can serve as a laboratory: measuring entropy, redundancy, divergence from uniformity, and reconstruction cost.',
        'The canonical mathematical point remains GFT. The performance of a particular compressor is an experimental validation matter.',
      ],
      formulas: [
        String.raw`\text{raw budget}=\text{compressible structure}+\text{entropic residue}`,
        String.raw`\log_2(m)=D_{KL}+H`,
      ],
      takeaways: [
        'Compression extracts persistence.',
        'GFT gives the budget language.',
        'A very good pedagogical bridge to information.',
      ],
      source: 'GitHub repository to publish: Igrekess/pt-compress; monograph ch04_gft, ch_PM.',
    },
  },
  {
    slug: 'zkp-preuves-sans-revelation',
    visual: 'zkp',
    repositories: [
      { name: 'Igrekess/PT_ZKP', status: 'planned' },
    ],
    scripts: [
      {
        file: 'zkp_toy_demo.py',
        title: 'Preuve jouet sans révélation',
        titleEn: 'Toy proof without revealing',
        summary: 'Montre la différence entre propriété vérifiable et témoin caché dans un protocole jouet.',
        summaryEn: 'Shows the difference between a verifiable property and a hidden witness in a toy protocol.',
      },
    ],
    fr: {
      title: 'ZKP : prouver sans révéler',
      shortTitle: 'ZKP',
      kicker: 'Mathématique PT / preuve',
      description: 'Pourquoi les preuves zero-knowledge parlent naturellement de persistance de structure.',
      status: 'exploration',
      l1: [
        'Une preuve zero-knowledge permet de convaincre quelqu’un qu’on sait quelque chose sans révéler ce quelque chose.',
        'Vu par la PT, c’est presque une image parfaite : on transmet une trace persistante suffisante pour prouver, sans transmettre toute la structure.',
      ],
      l2: [
        'Le vérificateur reçoit un invariant, pas l’objet complet. La preuve conserve la persistance de la propriété, mais dissipe l’information inutile sur le témoin.',
        'Cette séparation entre propriété persistante et contenu caché est exactement dans l’esprit GFT.',
      ],
      l3: [
        'Le projet PT_ZKP peut être présenté comme application math-info : protocoles, contraintes, témoins, vérification et budget informationnel.',
        'Statut : application exploratoire de la lecture PT, pas remplacement de la théorie cryptographique standard des ZKP.',
      ],
      formulas: [
        String.raw`\text{témoin caché}\quad w`,
        String.raw`\text{preuve publique}\quad \pi\Rightarrow V(x,\pi)=1`,
        String.raw`\text{persistance de propriété}\ne\text{révélation du témoin}`,
      ],
      takeaways: [
        'Un ZKP transmet une persistance sans tout révéler.',
        'Bon pont entre preuve, information et cryptographie.',
        'À présenter comme application exploratoire.',
      ],
      source: 'Dépôt GitHub à publier : Igrekess/PT_ZKP ; monographie ch04_gft.',
    },
    en: {
      title: 'ZKP: proving without revealing',
      shortTitle: 'ZKP',
      kicker: 'PT mathematics / proof',
      description: 'Why zero-knowledge proofs naturally speak about persistence of structure.',
      status: 'exploration',
      l1: [
        'A zero-knowledge proof lets someone be convinced that you know something without revealing that thing.',
        'From PT, this is almost a perfect image: one transmits a persistent trace sufficient to prove, without transmitting the whole structure.',
      ],
      l2: [
        'The verifier receives an invariant, not the complete object. The proof preserves the property persistence, but dissipates unnecessary information about the witness.',
        'This separation between persistent property and hidden content is exactly in the spirit of GFT.',
      ],
      l3: [
        'The PT_ZKP project can be presented as a math-info application: protocols, constraints, witnesses, verification, and information budget.',
        'Status: exploratory application of the PT reading, not a replacement for standard cryptographic ZKP theory.',
      ],
      formulas: [
        String.raw`\text{hidden witness}\quad w`,
        String.raw`\text{public proof}\quad \pi\Rightarrow V(x,\pi)=1`,
        String.raw`\text{property persistence}\ne\text{witness revelation}`,
      ],
      takeaways: [
        'A ZKP transmits persistence without revealing everything.',
        'Good bridge between proof, information, and cryptography.',
        'To be presented as an exploratory application.',
      ],
      source: 'GitHub repository to publish: Igrekess/PT_ZKP; monograph ch04_gft.',
    },
  },
  {
    slug: 'atlas-theoremes-pt',
    visual: 'theorem-atlas',
    fr: {
      title: 'Atlas des théorèmes PT',
      shortTitle: 'Atlas théorèmes',
      kicker: 'Mathématique PT / carte',
      description: 'Une carte de lecture distinguant identités, théorèmes, ponts, dérivations et validations.',
      status: 'tool',
      l1: [
        'Une théorie ambitieuse devient vite intimidante. L’atlas sert à répondre à une question simple : qu’est-ce qui est prouvé, qu’est-ce qui est dérivé, qu’est-ce qui est testé ?',
        'C’est essentiel pour la crédibilité du site : la PT doit montrer sa force, mais aussi ses statuts avec honnêteté.',
      ],
      l2: [
        'L’atlas relie GFT, T0-T6, BA5, les ponts E/F/G, les dérivations physiques et les validations numériques.',
        'L’objectif n’est pas de remplacer la page Théorèmes, mais d’ajouter une carte narrative : de quelle brique dépend quel résultat ?',
      ],
      l3: [
        'Le niveau technique doit afficher les catégories épistémiques : ID, THM, BRIDGE, DER, VAL, PRED, META.',
        'Un résultat plus précis numériquement ne doit jamais être promu si son statut logique ne le permet pas. L’atlas est un garde-fou contre la confusion.',
      ],
      formulas: [
        String.raw`\text{ID}\rightarrow\text{THM}\rightarrow\text{BRIDGE}\rightarrow\text{DER}\rightarrow\text{VAL/PRED}`,
        String.raw`\text{précision numérique}\ne\text{statut théorématique}`,
      ],
      takeaways: [
        'Clarifie la force logique des résultats.',
        'Protège le site contre le sur-claiming.',
        'Donne une carte de lecture aux nouveaux lecteurs.',
      ],
      source: 'Monographie : frontmatter/status_ledger.tex, NOMENCLATURE_MAP.md, ch23_audit.',
    },
    en: {
      title: 'Atlas of PT theorems',
      shortTitle: 'Theorem atlas',
      kicker: 'PT mathematics / map',
      description: 'A reading map distinguishing identities, theorems, bridges, derivations, and validations.',
      status: 'tool',
      l1: [
        'An ambitious theory quickly becomes intimidating. The atlas answers a simple question: what is proved, what is derived, what is tested?',
        'This is essential for credibility: PT must show its strength, but also its statuses honestly.',
      ],
      l2: [
        'The atlas connects GFT, T0-T6, BA5, bridges E/F/G, physical derivations, and numerical validations.',
        'It does not replace the Theorems page; it adds a narrative map: which result depends on which brick?',
      ],
      l3: [
        'The technical level should display epistemic categories: ID, THM, BRIDGE, DER, VAL, PRED, META.',
        'A numerically more precise result must never be promoted if its logical status does not allow it. The atlas is a guardrail against confusion.',
      ],
      formulas: [
        String.raw`\text{ID}\rightarrow\text{THM}\rightarrow\text{BRIDGE}\rightarrow\text{DER}\rightarrow\text{VAL/PRED}`,
        String.raw`\text{numerical precision}\ne\text{theorem status}`,
      ],
      takeaways: [
        'Clarifies the logical strength of results.',
        'Protects the site against over-claiming.',
        'Gives new readers a map.',
      ],
      source: 'Monograph: frontmatter/status_ledger.tex, NOMENCLATURE_MAP.md, ch23_audit.',
    },
  },
  {
    slug: 'calculateur-persistance',
    visual: 'calculator',
    scripts: [
      {
        file: 'gft_identity_demo.py',
        title: 'Identité GFT',
        titleEn: 'GFT identity',
        summary: 'Même calcul que le widget, mais sous forme de script reproductible.',
        summaryEn: 'The same calculation as the widget, but as a reproducible script.',
        expected: 'résidu < 1e-12',
        expectedEn: 'residual < 1e-12',
      },
    ],
    fr: {
      title: 'Calculateur de persistance',
      shortTitle: 'Calculateur',
      kicker: 'Mathématique PT / outil',
      description: 'Manipuler directement la partition GFT entre entropie et information persistante.',
      status: 'tool',
      l1: [
        'Le calculateur rend GFT tactile : choisis une distribution, et vois le budget total se partager entre structure persistante et incertitude.',
        'Quand la distribution est uniforme, tout est entropie. Quand elle est concentrée, la persistance augmente.',
      ],
      l2: [
        'L’outil calcule $H(P)$, $D_{KL}(P\\|U_m)$ et $\\log_2(m)$ pour une distribution discrète normalisée.',
        'Il montre visuellement que la somme ne bouge pas : déplacer la masse de probabilité ne change pas le budget total, seulement sa répartition.',
      ],
      l3: [
        'Les calculs sont faits en base 2 : $H=-\\sum p_i\\log_2 p_i$ et $D_{KL}=\\sum p_i\\log_2(p_i/(1/m))$. Les zéros sont ignorés dans les sommes logarithmiques.',
        'Cette page illustre l’identité GFT ; elle ne dépend d’aucun ajustement ni d’aucune donnée empirique.',
      ],
      formulas: [
        String.raw`H(P)=-\sum_i p_i\log_2 p_i`,
        String.raw`D_{KL}(P\|U_m)=\sum_i p_i\log_2(mp_i)`,
        String.raw`\log_2(m)=D_{KL}+H`,
      ],
      takeaways: [
        'Un outil simple pour comprendre le principe fondamental.',
        'Uniforme : entropie maximale, persistance nulle.',
        'Concentré : persistance forte, entropie faible.',
      ],
      source: 'Monographie : ch04_gft.',
    },
    en: {
      title: 'Persistence calculator',
      shortTitle: 'Calculator',
      kicker: 'PT mathematics / tool',
      description: 'Directly manipulating the GFT partition between entropy and persistent information.',
      status: 'tool',
      l1: [
        'The calculator makes GFT tactile: choose a distribution and see the total budget split between persistent structure and uncertainty.',
        'When the distribution is uniform, everything is entropy. When it is concentrated, persistence rises.',
      ],
      l2: [
        'The tool computes $H(P)$, $D_{KL}(P\\|U_m)$, and $\\log_2(m)$ for a normalized discrete distribution.',
        'It visually shows that the sum does not move: shifting probability mass does not change the total budget, only its partition.',
      ],
      l3: [
        'Computations use base 2: $H=-\\sum p_i\\log_2 p_i$ and $D_{KL}=\\sum p_i\\log_2(p_i/(1/m))$. Zeros are ignored in logarithmic sums.',
        'This page illustrates the GFT identity; it depends on no fit and no empirical data.',
      ],
      formulas: [
        String.raw`H(P)=-\sum_i p_i\log_2 p_i`,
        String.raw`D_{KL}(P\|U_m)=\sum_i p_i\log_2(mp_i)`,
        String.raw`\log_2(m)=D_{KL}+H`,
      ],
      takeaways: [
        'A simple tool to understand the fundamental principle.',
        'Uniform: maximal entropy, zero persistence.',
        'Concentrated: strong persistence, low entropy.',
      ],
      source: 'Monograph: ch04_gft.',
    },
  },
];

export function getMathPage(slug: string): MathPage | undefined {
  return mathPages.find((page) => page.slug === slug);
}
