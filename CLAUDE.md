# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

`rms-host-galileo` is the Galileo host plug-in for `oops` (rms-oops). It is a port of
`rms-oops/src/oops/hosts/galileo/` to the newer `oops.Host` / `oops.host` API, part of
splitting the mission hosts out of rms-oops (see the sibling `rms-host-cassini`). SSI is
ported.

## Standards

Tooling and conventions follow **rms-polymath**. `.claude/rules/*.md`, which were copied from
it, are the authoritative coding, testing, documentation, dependency and environment
standards. The repo layout, CI and scripts come from `rms-devenv/repo_template`, adapted as
rms-polymath adapted them, by way of `rms-host-cassini`. Model changes to packaging, CI or
scripts on rms-polymath's `pyproject.toml`, `.flake8`, `.github/workflows/` and `scripts/`,
and keep them in step with rms-host-cassini. Some details in the rules are
polymath-specific, for example the two published `.pyi` stubs and `stubtest`. Ignore those
until this repo has the equivalent.

- Hand-aligned columns (assignments, dict colons, `from x     import`) are deliberate. Never
  run `ruff format` or any other auto-formatter.
- Keep the file banners: a line of 90 `#` characters at the top and bottom of each module,
  with a `# host_galileo/<file>.py` header.

## Current state

- `scripts/run-all-checks.sh` is the single source of truth for the checks: ruff, flake8
  `E12`/`E13`, mypy on `tests/`, pytest, pyroma, pip-audit, the Sphinx build (`-W -n`),
  codespell and PyMarkdown. `ruff format` is disabled. `.github/workflows/run-tests.yml`
  runs the same set. `scripts/setup-venv.sh` builds the `venv/` that the scripts expect;
  set `VENV` to use another.
- The `oops.Host` API is unreleased. CI, `.readthedocs.yaml` and `scripts/setup-venv.sh`
  install rms-oops from its `mrs_260925_host_reorg` branch, and `docs/conf.py` ignores the
  unresolvable `oops.Host` and `oops.host` references. Drop all four once a release has
  them.
- `docs/conf.py` resolves the short class names in docstrings (`FCPath`, `Observation`,
  ...) through intersphinx with `_SHORT_NAMES`. Add a name there when a nitpicky build
  reports a class from another RMS package.
- Coverage must stay at or above 65% (`fail_under` in `pyproject.toml`, and the project and
  patch targets in `codecov.yml`); it is at 100%. The ubuntu / Python 3.13 CI job uploads
  to Codecov and fails if the upload fails. The upload authenticates through GitHub OIDC
  (`use_oidc`), so it needs no `CODECOV_TOKEN` secret.
- Tests never touch SPICE kernels or spicedb. `tests/oops/conftest.py` replaces them with
  fakes for every test, so new tests get this automatically. Unlike Cassini, Galileo
  furnishes a fixed kernel list through the SPICE file cache; the fakes make the cache
  return paths unread and record what `cspyce.furnsh` is given.

## Architecture gotchas

- The oops support lives in the `host_galileo.oops` subpackage: `_Galileo` in
  `_galileo.py` and `GalileoSSI` in `ssi.py`. `GalileoSSI` calls `_register()` at import
  time, and `oops/__init__.py` imports it. `Host.from_file()` dispatches through the
  registered hosts' `_detect_in_*` methods, so the host is invisible until
  `host_galileo.oops` (or its module) is imported.
- `host_galileo/__init__.py` imports the hosts lazily, through a PEP 562 `__getattr__`, so
  that importing `host_galileo` does not import oops. Keep any import of
  `host_galileo.oops` out of `host_galileo/__init__.py`.
- `_Galileo` is never instantiated; all its state lives in class attributes. It keeps the
  monthly CK/SPK machinery that Cassini uses, but every month is marked loaded and
  `load_kernels()` furnishes a fixed list (`_KERNELS`) instead. `tdb_from_sclk()` loads
  only the leap-second and clock kernels, so it works without `initialize()`.
- The metadata of a VICAR image come from its own header: `_vicar_to_pds3()` copies
  each value under its PDS3 keyword name and builds the clock count, `IMAGE_TIME` and
  `FILTER_NAME` from the VICAR fields. No detached `.LBL` is read.
- `IMAGE_TIME` is the middle of the exposure. When it is `UNK`, the start time comes from
  `SPACECRAFT_CLOCK_START_COUNT` plus the shutter delay for the `FRAME_DURATION`
  (`_SHUTTER_DELAY`), one frame earlier for an `EXTENDED` exposure.
- An SSI supplemental index writes `SPACECRAFT_CLOCK_START_COUNT` in SPICE form and the
  index in label form, so the default validation of `oops.Host._read_index_rows()`, which
  compares every shared column, rejects every real pair. `GalileoSSI.from_index()` passes
  `validate=_INDEX_KEYS` to match rows on `VOLUME_ID` and `FILE_SPECIFICATION_NAME`
  instead. Unless one is given, it reads the supplement beside the index
  (`oops.Host.SUPPLEMENT`), which must exist: an index without one raises
  `FileNotFoundError`. Index values win the merge, so the clock count stays in label
  form. `_read_index_rows()` is a classmethod that raises `WrongHostError` for another
  host's index, so call it through `GalileoSSI`.
  `oops.Host.from_index()` reads the rows itself with the default validation, so it still
  fails on any index whose supplement sits beside it.
- The hosts apply the standard overrides through `oops.Host._apply_overrides()`, passing
  them to `_make_snapshot()` as one `overrides` dict, as rms-host-cassini does.
- TARGET_NAME repairs and the star targets that become light sources are in `ssi.py`
  (`_TARGET_NAME_REPAIRS`, `_TARGET_STARS`).
- Running anything needs SPICE kernels: the fixed list is read through
  `spicedb.get_spice_filecache_prefix()`, and the instrument kernel through `spicedb`
  (`SPICE_PATH`, or `$OOPS_RESOURCES/SPICE`, plus `SPICE_SQLITE_DB_NAME`). `spicedb`
  currently ships inside rms-oops.

## Repo etiquette

- Commit subjects are a plain, capitalized, imperative sentence with no type prefix and no
  trailing period. PRs are squash-merged onto `main`.
- Branch names: `<initials>_<YYMMDD>_<topic>`.
- Versions come from `setuptools_scm`; never hand-write a version.
