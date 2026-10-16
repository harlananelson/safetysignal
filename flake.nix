{
  description = "safetysignal — Bayesian GPS (MGPS) signal detection for pharmacovigilance";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
        };

        # Same build faers-mobi does from the GitHub source; here from this checkout.
        safetysignal = pkgs.rPackages.buildRPackage {
          pname = "safetysignal";
          version = "0.2.1";
          src = self;
          propagatedBuildInputs = with pkgs.rPackages; [
            cli dplyr rlang tibble
          ];
        };

        rWithPkgs = pkgs.rWrapper.override {
          packages = with pkgs.rPackages; [
            # Imports (DESCRIPTION)
            cli
            dplyr
            rlang
            tibble

            # Suggests
            ggplot2
            knitr
            rmarkdown
            testthat

            # Dev
            devtools
            roxygen2
            languageserver
          ];
        };

      in {
        packages.default = safetysignal;

        devShells.default = pkgs.mkShell {
          name = "safetysignal";

          packages = [
            rWithPkgs
            pkgs.git
          ];

          env = {
            LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [
              pkgs.stdenv.cc.cc.lib
            ];
          };

          shellHook = ''
            # nix-name: give this shell a name when getpwuid has none for the uid.
            # Pure nix shells on iuhlwolfap1122 otherwise prompt "I have no name!".
            # No LD_PRELOAD. The nixpkgs nss_wrapper is linked to a newer glibc than
            # the host, and preloading it breaks host binaries started from the shell.
            if ! id -un >/dev/null 2>&1; then
            _nn_dir="''${TMPDIR:-/tmp}/nix-name-$$"
            mkdir -p "$_nn_dir"
            _nn_name=$(logname 2>/dev/null || printf '%s' "''${USER:-nixuser}")
            printf '%s\n' '#!/bin/sh' "printf '%s\\n' '$_nn_name'" > "$_nn_dir/whoami"
            chmod +x "$_nn_dir/whoami"
            export PATH="$_nn_dir:''${PATH}"
            export USER="$_nn_name"
            export LOGNAME="$_nn_name"
            export PS1="$_nn_name@\h:\w\\$ "
            unset _nn_dir _nn_name
            fi
            export USER="$(whoami 2>/dev/null || echo unknown)"
            export LANG=C.UTF-8
            export LC_ALL=C.UTF-8

            echo ""
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
            echo "  safetysignal dev environment"
            echo "  R: $(R --version | head -1)"
            echo "  Test:  Rscript -e 'devtools::test()'"
            echo "  Check: Rscript -e 'devtools::check()'"
            echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
          '';
        };
      });
}
