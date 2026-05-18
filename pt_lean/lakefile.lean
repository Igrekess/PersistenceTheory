import Lake
open Lake DSL

package «PT» where
  -- Use Mathlib's recommended Lean options
  leanOptions := #[
    ⟨`pp.unicode.fun, true⟩,
    ⟨`autoImplicit, false⟩
  ]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "master"

@[default_target]
lean_lib «PT» where
  -- All modules under PT/ are part of the library
  globs := #[.andSubmodules `PT]
