/-
PT — Persistence Theory, Lean 4 Formalisation (umbrella module).

Importing `PT` loads the full library. Modules are grouped by topic;
within each group they are listed alphabetically. See `README.md` at
the root of `PT_LEAN/` for the per-module status table, the critical
path T1 → T7 → W7-1, and the build instructions.
-/

-- Sieve (32 modules) — the discrete arithmetic core
import PT.Sieve.AdmissibleRMasterCharacterisation
import PT.Sieve.AdmissibleResiduesArithmetic
import PT.Sieve.AdmissibleResiduesCubeSumMod
import PT.Sieve.AdmissibleResiduesFifthPowerSumMod
import PT.Sieve.AdmissibleResiduesFourthPowerSumMod
import PT.Sieve.AdmissibleResiduesProductMod
import PT.Sieve.AdmissibleResiduesSixthPowerSumMod
import PT.Sieve.AdmissibleResiduesSquareSumMod
import PT.Sieve.AdmissibleSemigroupQuotient
import PT.Sieve.AdmissibleSquareGroupStructure
import PT.Sieve.Bimodality
import PT.Sieve.BimodalityCardinality
import PT.Sieve.BimodalityCharacterFormula
import PT.Sieve.BimodalityT1Projection
import PT.Sieve.CoprimeAdmissibilityProduct
import PT.Sieve.KleinFourEmbedding
import PT.Sieve.LegendreLogParity
import PT.Sieve.N1AtomicUniqueness
import PT.Sieve.N2SelfCoherence
import PT.Sieve.N3aFactorisation
import PT.Sieve.N4DimensionCascade
import PT.Sieve.PrimitiveRootCascade
import PT.Sieve.PrimorialCoprime
import PT.Sieve.SixRough
import PT.Sieve.T1AntidiagOrbits
import PT.Sieve.T1ForbiddenTransitions
import PT.Sieve.T1OrbitsZMod5
import PT.Sieve.T1OrbitsZMod7
import PT.Sieve.T3Antidiagonal
import PT.Sieve.T6aFieldStructure
import PT.Sieve.T6bAxiomsC1C4
import PT.Sieve.TotientCascadeIdentities

-- Stochastic (24 modules) — stationary distribution `s = 1/2`,
-- spectral analysis of `T₃` and its Kronecker products
import PT.Stochastic.SHalf
import PT.Stochastic.T2T3CesaroLimit
import PT.Stochastic.T2T3KroneckerEigenvalues
import PT.Stochastic.T2T3PerronEigenvectorUniqueness
import PT.Stochastic.T2T3SpectralMixing
import PT.Stochastic.T2T3StationaryUniqueness
import PT.Stochastic.T2T3T5KroneckerSpectrum
import PT.Stochastic.T30AntiSector
import PT.Stochastic.T30CharPolyComplete
import PT.Stochastic.T30FullDeterminantIdentity
import PT.Stochastic.T30FullEigenpairCount
import PT.Stochastic.T30FullSpectralAnalysis
import PT.Stochastic.T30PerronAntiCommutator
import PT.Stochastic.T30PerronUniqueness
import PT.Stochastic.T30SpectralMixingExtended
import PT.Stochastic.T30TraceCubeIdentity
import PT.Stochastic.T30TraceDeterminant
import PT.Stochastic.T30TraceFifthPowerIdentity
import PT.Stochastic.T30TraceFormulaExtended
import PT.Stochastic.T30TraceFourthPowerIdentity
import PT.Stochastic.T30TracePowerSequence
import PT.Stochastic.T30TraceSquaredIdentity
import PT.Stochastic.T3CesaroLimit
import PT.Stochastic.T3SpectralDecomposition

-- Conservation (25 modules) — T2 conservation exponent, gap-sum
-- identities, moments of the prime-gap distribution
import PT.Conservation.ConservationActivePrimes
import PT.Conservation.ConservationID
import PT.Conservation.ConservationIDExtensions
import PT.Conservation.ConservationIDPrimorial
import PT.Conservation.CumulativeBoundsExtended
import PT.Conservation.GapBoundedBelow
import PT.Conservation.GapDistributionAllMomentsMaster
import PT.Conservation.GapDistributionFifthMoment
import PT.Conservation.GapDistributionKurtosisExtended
import PT.Conservation.GapDistributionMeanExtended
import PT.Conservation.GapDistributionMomentsSummary
import PT.Conservation.GapDistributionSixthMoment
import PT.Conservation.GapDistributionSkewExtended
import PT.Conservation.GapDistributionVariance
import PT.Conservation.GapDistributionVarianceExtended
import PT.Conservation.GapEntropyBound
import PT.Conservation.GapMaxBound
import PT.Conservation.GapMomentsNExtended
import PT.Conservation.GapParityDecomposition
import PT.Conservation.GapStatisticalMoments
import PT.Conservation.GapSumMasterIdentity
import PT.Conservation.GapVarianceSmallN
import PT.Conservation.PrimeGapMoments
import PT.Conservation.PtPrimeStructuralTheorem
import PT.Conservation.T2Alpha

-- Information (26 modules) — GFT identity, Bekenstein bound, Shannon
-- and KL machinery, Shore-Johnson-Čencov axioms (T6c)
import PT.Information.BekensteinBound
import PT.Information.BekensteinExtensions
import PT.Information.BinaryEntropyMonotonicity
import PT.Information.ConditionalEntropy
import PT.Information.CrossEntropyBound
import PT.Information.EntropyAdditivityCorners
import PT.Information.EntropyBoundsTight
import PT.Information.EntropyDifferenceSymmetric
import PT.Information.EntropyJointProduct
import PT.Information.EntropyMasterFramework
import PT.Information.EntropyMonotonicity
import PT.Information.EntropyOfBinaryDistribution
import PT.Information.EntropyTernaryDistribution
import PT.Information.GFTIdentity
import PT.Information.GFTOnZMod30
import PT.Information.GFTSpecialMValues
import PT.Information.GFTSpecialisations
import PT.Information.InfoTheoreticMaster
import PT.Information.KLAdditivityFromMI
import PT.Information.KLAdditivityProduct
import PT.Information.L0MaxEntropy
import PT.Information.MutualInfoDistributional
import PT.Information.MutualInformationBasic
import PT.Information.RelativeEntropyAdditivity
import PT.Information.ShannonEntropyConcavity
import PT.Information.T6cChencov

-- Holonomy (26 modules) — `sin² θ_p = δ_p (2 − δ_p)`, the
-- active-prime criterion `γ_p > 1/2 ⇔ p ∈ {3, 5, 7}`,
-- and the `α = ∏ sin² θ_p` cascade
import PT.Holonomy.ActivePrimeCriterion
import PT.Holonomy.ActivePrimeMargins
import PT.Holonomy.AlphaGammaRelation
import PT.Holonomy.AlphaInverseCascadeIdentity
import PT.Holonomy.AlphaInversePowerSequence
import PT.Holonomy.AlphaInverseTimesGammaSum
import PT.Holonomy.AlphaPowerSequence
import PT.Holonomy.AlphaTimesGammaSum
import PT.Holonomy.CouplingReconstructionBounds
import PT.Holonomy.CyclicPhaseAlgebraic
import PT.Holonomy.CyclicPhaseIdentity
import PT.Holonomy.CyclicPhaseInversion
import PT.Holonomy.CyclicPhaseTable
import PT.Holonomy.GammaPrimorialProduct
import PT.Holonomy.GammaProduct
import PT.Holonomy.GammaRatio
import PT.Holonomy.GammaSumActive
import PT.Holonomy.GammaSumProduct
import PT.Holonomy.GammaTablesExtended
import PT.Holonomy.HolonomyMasterFramework
import PT.Holonomy.InvAlphaSquaredBracket
import PT.Holonomy.InverseSinSq
import PT.Holonomy.InverseSinSqProduct
import PT.Holonomy.SinSqProductBounds
import PT.Holonomy.SinSqProductChain
import PT.Holonomy.SinSqRatios

-- Fixed point (4 modules) — T7 `μ* = 15`
import PT.FixedPoint.DimensionProtection
import PT.FixedPoint.FixedPointMasterTheorem
import PT.FixedPoint.T7GlobalUniqueness
import PT.FixedPoint.T7MuStar

-- Number theory (1 module) — Mertens M3 scaffold; the PT-needed
-- compactness `‖mertensSum x − log log x‖ < ∞` is kernel-verified
import PT.NumberTheory.T5Mertens

-- EML (7 modules) — the two-branch Sheffer algebra `q±`
import PT.EML.EMLAlgebra
import PT.EML.EMLDepth3
import PT.EML.EMLIdentities
import PT.EML.EMLSheffer3Args
import PT.EML.QParameterMonotonicity
import PT.EML.QPlusQMinusComparison
import PT.EML.QSheffer

-- Analysis (1 module) — W7-1 spiral identity of the Weil explicit
-- formula, forward direction
import PT.Analysis.W7SpiralIdentity

-- Bridge (3 modules) — math/physics dissolution, derivation chain,
-- status graph
import PT.Bridge.MathPhysicsDissolution
import PT.Bridge.PTCascadeDerivationChain
import PT.Bridge.StatusGraphFormalisation

-- CRT decoupling (7 modules) — Theorem 6.1 of the PT-CRT preprint,
-- abstract algebraic + empirical + spectral reduction + ergodic
-- closure, with a concrete sieve dynamical-system instance.
-- The submodule has its own `PT/CrtDecoupling/README.md`.
import PT.CrtDecoupling.Empirical
import PT.CrtDecoupling.Main
import PT.CrtDecoupling.Phase3
import PT.CrtDecoupling.Phase4
import PT.CrtDecoupling.Phase4Instance
import PT.CrtDecoupling.SpectralReduction
import PT.CrtDecoupling.Tensor

-- Top-level / unified (1 module) — meta-theorem assembling the
-- master results into a single citation-ready statement
import PT.PTGrandUnifiedTheorem

/-!
# PT — Persistence Theory, Lean 4 Formalisation

Importing `PT` loads the full library. The PT critical path

```
T1 → T3 → s = 1/2 → T2 → L0 → T7 → pt_T5_mertens_compactness → W7-1
```

is *kernel-verified without `sorry`*: twenty-two foundational theorems
on this path are formally proved. Around 150 secondary modules cover
the surrounding ecosystem.

Across the 179 modules under `PT/`, **178 are entirely sorry-free**.
The single remaining `sorry` lives in
`PT/Information/G3FisherUniqueness.lean` and concerns the classical
uniqueness up to scale of the Fisher metric — a result the monograph
already marks `\leanExternal`. The W7-1 reverse direction
(`J σ k = J σ 0 ⟹ σ² = π(k+1)`) is not yet present as a theorem in
the file (left as future work, not as a `sorry`). Neither is on the
PT critical path.

The module catalog (by topic) and the build instructions live in
`README.md` at the root of `PT_LEAN/`.
-/
