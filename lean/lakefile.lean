import Lake
open Lake DSL

package «BootstrapTowerTransition» where
  srcDir := "."

@[default_target]
lean_lib «BootstrapTowerTransition» where
  roots := #[`BootstrapTowerTransition]
