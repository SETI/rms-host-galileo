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
and keep them in step with rms-host-cassini. The rules' references to the two published
`.pyi` stubs and `stubtest` apply here to the stubs described below.

- Hand-aligned columns (assignments, dict colons, `from x     import`) are deliberate. Never
  run `ruff format` or any other auto-formatter.
- Keep the file banners: a line of 90 `#` characters at the top and bottom of each module,
  with a `# host_galileo/<file>.py` header (`.pyi` for a stub).

## Current state

- `scripts/run-all-checks.sh` is the single source of truth for the checks: ruff, flake8
  `E12`/`E13`, mypy on `tests/`, pytest, pyroma, stubtest, pip-audit, the Sphinx build
  (`-W -n`), codespell and PyMarkdown. `ruff format` is disabled.
  `.github/workflows/run-tests.yml` runs the same set. `scripts/setup-venv.sh` builds the
  `venv/` that the scripts expect; set `VENV` to use another.
- The host needs the `oops.Host` support methods of the rms-oops host contract
  (`_index_rows`, `_set_path_attributes`, the new `_apply_overrides`), which are not
  released yet, so `pyproject.toml` installs rms-oops from that PR branch through a git
  URL. Restore a version floor once rms-oops is released.
- The type information is published in two stubs: `host_galileo/__init__.pyi` declares
  `GalileoSSI` in full, and `host_galileo/oops/__init__.pyi` re-exports it. `stubtest`
  checks them against the runtime in the check script and in CI, so adding, renaming or
  re-signing a public member means updating `host_galileo/__init__.pyi` in the same
  change. The `[tool.mypy] exclude`, the override list in `pyproject.toml` and
  `.stubtest-allowlist` keep it off the stub-less modules (`_galileo.py`,
  `_ssi_labels.py` and `ssi.py`). `from_file` and `from_index` are each published as two
  overloads, one without `select` that returns the narrower type. The two
  `type: ignore[override]` comments are real mismatches with `oops.Host`: `from_file`
  names its first parameter `fileinfo`, neither method takes the `**kwargs` that
  `oops.Host` declares, and the results are typed more narrowly.
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
  return paths unread and record what `cspyce.furnsh` is given. `fake_spice` also copies
  the oops Path and Frame registries, caches, and each class's waypoint and wayframe
  pools for every test, as rms-host-new-horizons does, until oops has one helper for it.
  The writers at the end of the conftest build small synthetic products (a VICAR image
  and its PDS3 label) and an index with its supplement; no test reads archive files, and
  the cut-down files in `test_files` are part of the repository.

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
- What `from_file()` reads follows the rules all four host packages share. Given the
  VICAR image, only that file is read: its header, plus its array unless `astrometry` is
  True. Given the PDS3 label (a name ending in `.lbl`), the label is parsed with the
  strict method (`_PDS3_DEFAULTS`, which the caller's `pds3_options` update), and unless
  `astrometry` is True the image file that `^IMAGE` names is found case-tolerantly
  through oops's `_resolve_pds3_filename` and read. One `VicarImage` read gives both the
  header and the array, laid out by the VICAR header, never by the record in `^IMAGE`
  (52 GO_0019 labels put it one record early). The header is read with the loose VICAR
  reader (`strict=False`).
- Every label and image read of `from_file()` goes through `_read_product()`, which turns
  each failure to read or interpret either file into `PdsHostError`, `VicarHostError` or
  `HostFormatError` naming the file. Interpretation (`_ssi_labels` and `_interpret()`) is
  plain Python with no SPICE call, so the `KeyError`, `TypeError` and `ValueError` it
  catches come from the files; the clock conversion runs later, in `_make_snapshot()`,
  so a cspyce error keeps its SPICE code. It also checks that the array has the
  dimensions the label gives it. Index rows get no checks.
- `_ssi_labels.py` turns each source into metadata items and merges them. Label items
  are `Pds3Label.as_dict()` without None values and PDS3 null texts
  (`oops.host.is_pds3_null()`), with `CUT_OUT_WINDOW` and `TRUTH_WINDOW` from their
  `_list` keys, and without an `ORBIT_NUMBER` outside 0-35. Header items are
  `VicarLabel.as_dict(unique=False)` without -999 and null texts, with each converted
  value only under its PDS3 name (`_VICAR_PDS3_PAIRS`, the clock count, the SCET/ERT
  times as datetimes unless the year is -32768, RATE, GAIN, FILTER, and the FIBE/MOFIBE
  exposure type). No `ORBIT_NUMBER` is derived from a header. `merge_items()` keeps the
  source that `_RULES` names for an item both give, or the label for an item without a
  rule. Each rule carries its evidence and cites the label-differences issue
  (SETI/rms-host-galileo#41), which stays open. `label_dict` and `header_dict` are the
  parsers' dictionaries as they come.
- `IMAGE_TIME` is the middle of the exposure. Without it (absent from merged metadata,
  or a null text in an index row), the start comes from the public
  `GalileoSSI.tstart_from_sclk()`: `SPACECRAFT_CLOCK_START_COUNT` plus the shutter delay
  for the `FRAME_DURATION` (`_SHUTTER_DELAY`), one frame earlier for an `EXTENDED`
  exposure before `_PHASE2_FIRST_RIM`. It loads only the leap-second and clock kernels.
- `GalileoSSI.from_index()` reads the index through `oops.Host._index_rows()`, which
  streams both tables and pairs them by position (`MERGE_SUPPLEMENT_BY_FILE = False`),
  checking that each pair names the same file. It passes `validate=()`, because the
  supplement writes `SPACECRAFT_CLOCK_START_COUNT` in SPICE form and the index in label
  form. The supplement is required: `supplement=None` raises `HostIndexError` unless
  `row_dicts` is given, and a missing supplement, default or explicit, raises
  `FileNotFoundError`. Each row is used
  as the table gives it, `_mask` entries included; a `CUT_OUT_WINDOW` of -1s means no
  window. `filepath` is the relative `VOLUME_ID/FILE_SPECIFICATION_NAME`. The method
  returns a generator expression, so its call-time errors are raised at once.
- The file subfields come from `oops.Host._set_path_attributes()`, called after `data` is
  inserted. The standard overrides go through `oops.Host._apply_overrides()`, passed to
  `_make_snapshot()` as one `overrides` dict, and `select` through
  `oops.Host._apply_select()`, as from a tuple holding the one Snapshot.
- No pointing is checked: a gap in the C kernels shows up as the cspyce error when the
  frame is used.
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
