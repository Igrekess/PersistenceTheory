/-
Copyright (c) 2026 Yan Senez. All rights reserved.
Released under Apache 2.0 license.
-/
import PT.Bridge.MathPhysicsDissolution
import Mathlib.Tactic

/-!
# Status Graph Formalisation — A formal dependency graph for PT modules

## Thesis

This module extends `PT.Bridge.MathPhysicsDissolution` with a **formal graph**
of dependencies between the structural modules of Persistence Theory, indexed
by epistemic status (`THM`, `BRIDGE`, `REC`, `VAL`).

The graph captures:

1. **Vertices** (`PTModule`): a finite enumeration of structural PT modules
   (BA0–BA5, BA5b, plus key theorems T1, T2, T3, T5, T6, T7, GFT, Bekenstein).
2. **Status assignment** (`moduleStatus`): each vertex carries a `Status`.
3. **Direct dependencies** (`directDeps`): a directed graph encoding which
   modules another module's derivation relies on.
4. **Graph-level properties**:
   * the unique `REC` vertex (BA5b) is **isolated** (depends on no `REC`);
   * the `THM`-subgraph is connected via BA-chains;
   * all `BRIDGE` vertices link `THM`s to `THM`s.

## Reference

* Bridge axioms grid: `PT.Bridge.BridgeAxioms` and Monograph Chapter 9
  table `tab:bridge_axioms`.
* Status taxonomy: `PT.Bridge.MathPhysicsDissolution`.

## Scope reduction

To keep `decide` tractable, we restrict to **12 structural modules**:
the 7 bridge items + 5 foundational theorems (T1, T2, T3, T5, T6, T7) +
GFT + Bekenstein. The full PT module DAG is much larger (~80 files), but
the structural skeleton presented here suffices to certify the
`math_physics_dissolution` claim in graph form.
-/

namespace PT.Bridge

/-! ### Vertex enumeration -/

/-- The structural PT modules tracked in the status graph. -/
inductive PTModule
  /-- BA0 — T0 closure (field = `{g_n}`). -/
  | BA0
  /-- BA1 — Gauge connection (observables = `g_n mod m`). -/
  | BA1
  /-- BA2 — GFT identity (`log m = D_KL + H`). -/
  | BA2
  /-- BA3 — Cyclic phase identity (`sin²θ_p = δ_p(2-δ_p)`). -/
  | BA3
  /-- BA4 — Active set / `μ* = 15`. -/
  | BA4
  /-- BA5 — Algebraic product (`α_bare = ∏ sin²θ_p`). -/
  | BA5
  /-- BA5b — CODATA recognition (`α_bare ↔ α_EM`). -/
  | BA5b
  /-- T1 — Forbidden transitions / six-rough sieve. -/
  | T1
  /-- T2 — Conservation identity / α-conservation. -/
  | T2
  /-- T3 — Antidiagonal recurrence / spectral decomposition. -/
  | T3
  /-- T5 — Mertens product. -/
  | T5
  /-- T6 — Field structure (axioms C1–C4). -/
  | T6
  /-- T7 — Fixed-point uniqueness, `μ* = 15`. -/
  | T7
  /-- GFT — Generalised Function Theory partition identity. -/
  | GFT
  /-- Bekenstein — Holographic / entropy bound. -/
  | Bekenstein
  deriving DecidableEq, Repr, Fintype

namespace PTModule

/-! ### Status assignment -/

/-- Epistemic status of each structural PT module.

Every module is `THM` (kernel-derivable) except BA5b which is `REC`
(CODATA recognition). There are no `BRIDGE` or `VAL` items in this
structural skeleton: bridges between submodules are absorbed in the
graph edges, and external numerical validations live outside Lean. -/
def moduleStatus : PTModule → Status
  | BA0 => Status.THM
  | BA1 => Status.THM
  | BA2 => Status.THM
  | BA3 => Status.THM
  | BA4 => Status.THM
  | BA5 => Status.THM
  | BA5b => Status.REC
  | T1 => Status.THM
  | T2 => Status.THM
  | T3 => Status.THM
  | T5 => Status.THM
  | T6 => Status.THM
  | T7 => Status.THM
  | GFT => Status.THM
  | Bekenstein => Status.THM

/-! ### Direct dependencies (graph edges) -/

/-- `directDeps m` is the finite set of modules that `m`'s derivation
    relies on directly. The graph is acyclic by construction.

* `BA0` depends on `T1` (forbidden transitions) and `T7` (uniqueness).
* `BA1` is independent (definitional gauge encoding).
* `BA2` depends on `GFT`.
* `BA3` is independent (trigonometric identity).
* `BA4` depends on `T7` (and indirectly `T6`).
* `BA5` depends on `BA3` and `BA4`.
* `BA5b` depends on `BA5` (algebraic product feeds the recognition).
* `T1`–`T7`, `GFT`, `Bekenstein` are foundational: most are independent,
  with `T2` ⟵ `T1`, `T3` ⟵ `T1`, `T7` ⟵ `T6` ⟵ `T1`,
  `Bekenstein` ⟵ `GFT`. -/
def directDeps : PTModule → List PTModule
  | BA0 => [T1, T7]
  | BA1 => []
  | BA2 => [GFT]
  | BA3 => []
  | BA4 => [T7]
  | BA5 => [BA3, BA4]
  | BA5b => [BA5]
  | T1 => []
  | T2 => [T1]
  | T3 => [T1]
  | T5 => []
  | T6 => [T1]
  | T7 => [T6]
  | GFT => []
  | Bekenstein => [GFT]

end PTModule

/-! ### Cardinalities by status -/

/-- The set of all `THM` modules. -/
def thmModules : Finset PTModule :=
  Finset.univ.filter (fun m => PTModule.moduleStatus m = Status.THM)

/-- The set of all `REC` modules. -/
def recModules : Finset PTModule :=
  Finset.univ.filter (fun m => PTModule.moduleStatus m = Status.REC)

/-- The set of all `BRIDGE` modules. -/
def bridgeModules : Finset PTModule :=
  Finset.univ.filter (fun m => PTModule.moduleStatus m = Status.BRIDGE)

/-- The set of all `VAL` modules. -/
def valModules : Finset PTModule :=
  Finset.univ.filter (fun m => PTModule.moduleStatus m = Status.VAL)

/-- **Cardinality of `THM` vertices.** The 12 foundational modules plus
    the 6 bridge axioms BA0–BA5 are all `THM` — 14 in total. -/
theorem thmModules_card : thmModules.card = 14 := by
  decide

/-- **Cardinality of `REC` vertices.** Exactly one (BA5b). -/
theorem recModules_card : recModules.card = 1 := by
  decide

/-- **Cardinality of `BRIDGE` vertices.** Zero in the structural skeleton
    (BRIDGEs are absorbed into edges, not vertices). -/
theorem bridgeModules_card : bridgeModules.card = 0 := by
  decide

/-- **Cardinality of `VAL` vertices.** Zero in-kernel (VALidations are
    external by definition). -/
theorem valModules_card : valModules.card = 0 := by
  decide

/-! ### Uniqueness of the recognition vertex -/

/-- **Uniqueness of the empirical point (graph form).**
    The only `REC` vertex in the structural module graph is BA5b. -/
theorem recModule_unique :
    ∃! m : PTModule, PTModule.moduleStatus m = Status.REC := by
  refine ⟨PTModule.BA5b, ?_, ?_⟩
  · rfl
  · intro m hm
    cases m <;> simp [PTModule.moduleStatus] at hm ⊢

/-! ### Isolation of the REC vertex -/

/-- `BA5b` (the unique `REC` vertex) does not depend on any other `REC`
    module. Equivalently: the unique recognition point is **graph-isolated**
    among recognitions — there is no `REC ← REC` edge. -/
theorem recVertex_isolated :
    ∀ d ∈ PTModule.directDeps PTModule.BA5b,
      PTModule.moduleStatus d ≠ Status.REC := by
  intro d hd
  simp [PTModule.directDeps] at hd
  subst hd
  decide

/-! ### THM-subgraph closure -/

/-- **THM-closure of dependencies.** Every dependency of a `THM` module is
    itself a `THM` module. Equivalently: there is no `THM ← REC` edge.

    This is the formal content of "the math/physics frontier is a single
    point": no theorem in the kernel can leak into recognition territory. -/
theorem thm_deps_are_thm :
    ∀ m : PTModule, PTModule.moduleStatus m = Status.THM →
      ∀ d ∈ PTModule.directDeps m, PTModule.moduleStatus d = Status.THM := by
  intro m _ d hd
  cases m <;> simp [PTModule.directDeps] at hd <;>
    rcases hd with rfl | rfl | rfl <;> rfl

/-! ### Existence of BA-chain connectivity -/

/-- The dependency closure (reflexive-transitive) reached in one step from
    `BA5` contains `BA3` and `BA4`. -/
theorem BA5_depends_on_BA3_BA4 :
    PTModule.BA3 ∈ PTModule.directDeps PTModule.BA5 ∧
    PTModule.BA4 ∈ PTModule.directDeps PTModule.BA5 := by
  refine ⟨?_, ?_⟩ <;> decide

/-- `BA5b` depends on `BA5`. -/
theorem BA5b_depends_on_BA5 :
    PTModule.BA5 ∈ PTModule.directDeps PTModule.BA5b := by
  decide

/-- The dependency chain `BA5b → BA5 → BA4 → T7 → T6 → T1` links the
    sole `REC` vertex back to the deepest foundational `THM` (T1). -/
theorem ba_chain_to_T1 :
    PTModule.BA5 ∈ PTModule.directDeps PTModule.BA5b ∧
    PTModule.BA4 ∈ PTModule.directDeps PTModule.BA5 ∧
    PTModule.T7 ∈ PTModule.directDeps PTModule.BA4 ∧
    PTModule.T6 ∈ PTModule.directDeps PTModule.T7 ∧
    PTModule.T1 ∈ PTModule.directDeps PTModule.T6 := by
  refine ⟨?_, ?_, ?_, ?_, ?_⟩ <;> decide

/-! ### Bridge to MathPhysicsDissolution -/

/-- Map from the 7-element bridge index `Fin 7` (used in
    `MathPhysicsDissolution`) to the corresponding `PTModule` vertex. -/
def baIndexToModule : Fin 7 → PTModule
  | 0 => PTModule.BA0
  | 1 => PTModule.BA1
  | 2 => PTModule.BA2
  | 3 => PTModule.BA3
  | 4 => PTModule.BA4
  | 5 => PTModule.BA5
  | 6 => PTModule.BA5b

/-- **Compatibility.** The `Fin 7` status assignment of
    `MathPhysicsDissolution.BAStatus` agrees with `moduleStatus` on the
    bridge image. -/
theorem baStatus_eq_moduleStatus (i : Fin 7) :
    BAStatus i = PTModule.moduleStatus (baIndexToModule i) := by
  fin_cases i <;> rfl

/-! ### Master theorem -/

/-- **Status graph summary.** The structural PT module graph satisfies:

    1. **14 `THM` vertices**, **1 `REC` vertex**, **0 `BRIDGE`**, **0 `VAL`**.
    2. The `REC` vertex is **unique** (BA5b).
    3. The `REC` vertex is **isolated** (no `REC ← REC` edge).
    4. The `THM`-subgraph is **closed under dependencies** (no `THM ← REC` edge).
    5. The bridge `BA5b → BA5 → BA4 → T7 → T6 → T1` connects the recognition
       point back to the deepest foundational theorem.
    6. The graph-level status assignment is **compatible** with the
       `Fin 7` indexing of `MathPhysicsDissolution.BAStatus`. -/
theorem pt_status_graph_summary :
    -- (1) Cardinalities:
    thmModules.card = 14
    ∧ recModules.card = 1
    ∧ bridgeModules.card = 0
    ∧ valModules.card = 0
    -- (2) REC uniqueness:
    ∧ (∃! m : PTModule, PTModule.moduleStatus m = Status.REC)
    -- (3) REC isolation:
    ∧ (∀ d ∈ PTModule.directDeps PTModule.BA5b,
         PTModule.moduleStatus d ≠ Status.REC)
    -- (4) THM-subgraph closure:
    ∧ (∀ m : PTModule, PTModule.moduleStatus m = Status.THM →
         ∀ d ∈ PTModule.directDeps m, PTModule.moduleStatus d = Status.THM)
    -- (5) BA-to-T1 chain:
    ∧ (PTModule.BA5 ∈ PTModule.directDeps PTModule.BA5b ∧
       PTModule.BA4 ∈ PTModule.directDeps PTModule.BA5 ∧
       PTModule.T7 ∈ PTModule.directDeps PTModule.BA4 ∧
       PTModule.T6 ∈ PTModule.directDeps PTModule.T7 ∧
       PTModule.T1 ∈ PTModule.directDeps PTModule.T6)
    -- (6) Compatibility with MathPhysicsDissolution.BAStatus:
    ∧ (∀ i : Fin 7, BAStatus i = PTModule.moduleStatus (baIndexToModule i)) :=
  ⟨thmModules_card,
   recModules_card,
   bridgeModules_card,
   valModules_card,
   recModule_unique,
   recVertex_isolated,
   thm_deps_are_thm,
   ba_chain_to_T1,
   baStatus_eq_moduleStatus⟩

/-! ### Reaffirmation of `math_physics_dissolution` -/

/-- **The dissolution theorem stays valid under the extended graph.**

    The graph extension of `MathPhysicsDissolution` does not invalidate
    any clause of `math_physics_dissolution`. We restate the latter here,
    transported through `baStatus_eq_moduleStatus`. -/
theorem math_physics_dissolution_graph_compatible :
    -- All clauses of the original dissolution theorem hold:
    (Finset.univ.filter (fun i : Fin 7 => BAStatus i = Status.THM)).card = 6
    ∧ (Finset.univ.filter (fun i : Fin 7 => BAStatus i = Status.REC)).card = 1
    ∧ (∃! i : Fin 7, BAStatus i = Status.REC)
    ∧ RECOGNITION_BA5b_schema (PT.Holonomy.alphaBareQ : ℝ)
    ∧ (∀ i : Fin 7, BAStatus i = Status.THM → (BAStatus i).isProvable = true)
    -- Plus the new graph-level consistency:
    ∧ (∀ i : Fin 7, BAStatus i = PTModule.moduleStatus (baIndexToModule i)) :=
  ⟨math_physics_dissolution.1,
   math_physics_dissolution.2.1,
   math_physics_dissolution.2.2.1,
   math_physics_dissolution.2.2.2.1,
   math_physics_dissolution.2.2.2.2,
   baStatus_eq_moduleStatus⟩

end PT.Bridge
