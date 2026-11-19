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
- The minimum rms-oops in `pyproject.toml` is 0.5.0, the first release with
  `oops.OopsSpiceError` (the `oops.Host` API came in 0.4.0). Until 0.5.0 is on PyPI,
  develop against the rms-oops branch that adds it, installed into `venv/`.
- The type information is published in two stubs: `host_galileo/__init__.pyi` declares
  `GalileoSSI` in full, and `host_galileo/oops/__init__.pyi` re-exports it. `stubtest`
  checks them against the runtime in the check script and in CI, so adding, renaming or
  re-signing a public member means updating `host_galileo/__init__.pyi` in the same
  change. The `[tool.mypy] exclude`, the override list in `pyproject.toml` and
  `.stubtest-allowlist` keep it off the stub-less modules (`_galileo.py`,
  `_ssi_labels.py` and `ssi.py`). The two `type: ignore[override]` comments are real
  mismatches with `oops.Host`: `from_file` names its first parameter `fileinfo`, and
  `from_index` returns `list[Snapshot]`.
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
- What `from_file()` reads follows the operator's rules of 2026-11-18, which all four
  host packages share. Given the VICAR image, only that file is read: its header, plus its
  array unless `astrometry` is True. Given the PDS3 label, the label is read, and unless
  `astrometry` is True, so is the whole image file, found from `^IMAGE` in the label's
  directory, in the label's case or the other (`_image_path`). The array always comes
  from the VICAR file's own record layout (`_read_vicar_image`), never from the record
  number in `^IMAGE`: 52 GO_0019 labels put it a record early, and the label of GO_0018
  E4/EUROPA/C0374685452R describes the wrong layout. The holdings cannot be corrected.
- A PDS3 label that the fast parser rejects, for its syntax (`PdsHostError`) or for a
  value it cannot convert such as an hour of 25 (`OopsValueError`), is retried with the
  strict parser and then the loose one (`_parse_fileinfo`); a requested method starts the
  chain there, and if all fail the result is `PdsHostError`. The parsers always run with
  `sources=False`, whatever the caller passes. A parsed `Pds3Label`, `VicarLabel`, or
  `VicarImage` passed in is checked with the detectors (`WrongHostError`). A label with
  no `^IMAGE` raises `PdsHostError`, and an image shorter than its header says (numpy's
  bare `ValueError`) raises `VicarHostError`.
- `_ssi_labels.py` turns each source into metadata items and merges them, all in forms
  `json` can write. `label_json()` makes `label_dict`: `as_dict()` with each time as the
  parser's ISO text (`_fmt`, no trailing `Z`), each set as its `_list`, no None (the
  parser's `END`), and no key the parser adds beside the keywords (the lower-case
  suffixes `_fmt`, `_list`, `_offset`, `_unit`, `_day`, `_sec`, `_unwrap`, `_source`,
  and `objects`; only a pointer's `_offset` holds more than its keyword, and the host
  never uses it). `header_json()` makes `header_dict`: `as_dict(unique=True)` with each
  repeated parameter keyed `"NAME[n]"` (`TASK[0]`, `TASK[1]`), the form all three VICAR
  hosts share. `label_items()` drops PDS3 nulls (`UNK`, `N/A`, `NULL`, any case, and
  lists of nothing else) and null list elements, and drops an `ORBIT_NUMBER` outside the
  tour (0-35): the labels of
  GO_0003-GO_0009 and GO_0016 give a pre-Jupiter frame the first two characters of its
  observation ID read as a 16-bit integer (12613 for `E1`). `header_items()` strips
  string padding, drops fill values (numeric -999, `PICNO` or `SPICE_C_ID` "?", time
  fields out of range), and gives every value it converts only under its PDS3 name
  (design D4): the `_VICAR_PDS3_PAIRS` keywords, `RIM` and `MOD91` as the clock count,
  the `SCET`/`ERT` fields as `IMAGE_TIME`/`EARTH_RECEIVED_TIME`, `RATE` and `GAIN` as
  `FRAME_DURATION` and `GAIN_MODE_ID`. `MOD10`, `MOD8`, `FIBE`, and `MOFIBE` are not
  wholly held by any PDS3 item and keep their names. `ORBIT_NUMBER` comes from an
  observation ID that starts with a tour orbit code (the IDs of the Earth flybys `E1` and
  `E2`, Gaspra `GA`, Ida `ID`, Venus `VE`, optical navigation `OP`, and calibrations `PH`
  name no orbit). `row_items()` cleans an index row the same way (numpy values to Python;
  nulls, -999, and the supplement's `-1` window dropped). It ignores the table reader's
  masks: the indexes declare real values as null constants (`EXPOSURE_DURATION` 0, the
  `ON_CHIP_MOSAIC_FLAG` N), and honoring them drops zero-exposure rows.
- `merge_items()` keeps every item of either source. Where both give an item and the
  values differ, `resolve()` keeps the source that `_RULES` names for that item, the same
  for every image, or the PDS3 label for an item with no rule; there is no per-image test.
  Each rule cites its evidence (survey counts, SPICE checks, DATASET.CAT) and #41, which
  lists every difference, each rule's reason, and the images on which its source is
  wrong. A value missing from the winning source falls to the other; nothing else does.
  `dict` is the merged result; `label_dict` and `header_dict` are present only when read.
  A missing value is never `None` or a placeholder: no key, no subfield.
- `_make_snapshot()` computes everything from the merged items, so a label read with its
  image gets the best values of both. The exposure, a time (`IMAGE_TIME`, or the clock
  count and `FRAME_DURATION`), and `TELEMETRY_FORMAT_ID` are required and raise
  `HostFormatError` naming the item and the file; an index row without a telemetry format
  (one read without its supplement) gets the full frame instead (`_INDEX_DEFAULT_MODE`).
- Missing pointing raises `oops.OopsSpiceError` (`_check_pointing`), built with
  `from_spice_error()` from the cspyce exception (only `LookupError`, `OSError`,
  `RuntimeError`, `ValueError`) at the midtime after any cadence override, naming the
  image file. Only `from_file()` checks, and only when no `frame` is given;
  `from_index()` never checks and never drops a row. `abspath` is removed only for an
  astrometry read of a remote label, whose image was never fetched.
- `IMAGE_TIME` is the middle of the exposure. Without it, the start time comes from
  `SPACECRAFT_CLOCK_START_COUNT` plus the shutter delay for the `FRAME_DURATION`
  (`_SHUTTER_DELAY`), one frame earlier for an `EXTENDED` exposure taken before the
  phase-2 flight software of the Jupiter tour (`_PHASE2_FIRST_RIM`).
- An SSI supplemental index writes `SPACECRAFT_CLOCK_START_COUNT` in SPICE form and the
  index in label form, so the default validation of `oops.Host._read_index_rows()`, which
  compares every shared column, rejects every real pair. `GalileoSSI.from_index()` passes
  `validate=_INDEX_KEYS` to match rows on `VOLUME_ID` and `FILE_SPECIFICATION_NAME`
  instead. By default (`supplement=oops.Host.SUPPLEMENT`, as in `oops.Host.from_index`)
  it reads the supplement beside the index, which must exist: an index without one
  raises `FileNotFoundError`, and `supplement=None` raises `HostIndexError` rather than
  meaning "search". It takes `label=` like the base, refuses any `select` but None
  (`OopsValueError`, since each row is one observation), and refuses any other keyword
  (`OopsTypeError`). Index values win the merge, so the clock count stays in label
  form. A row's `FILE_SPECIFICATION_NAME` names the label; the Snapshot's `filepath`
  names the image beside it, with `.IMG` (or `.img`) in place of `.LBL`, as a full path
  or URL: `_volumes_directory()` finds the volume directories from the index's own
  location (beside the volume of a volume's `INDEX` directory, `volumes` in place of
  `metadata` in an RMS Node metadata tree, else the index's directory), and the row's
  `VOLUME_ID` and file specification complete the path (design D5). `abspath` is set
  only for a local file or one already in the cache.
  `_read_index_rows()` is a classmethod that raises `WrongHostError` for another host's
  index, so call it through `GalileoSSI`.
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
