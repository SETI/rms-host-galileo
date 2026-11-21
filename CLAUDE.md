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
  neither method takes the `**kwargs` that `oops.Host` declares.
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
  return paths unread and record what `cspyce.furnsh` is given. The `file_cache`
  fixture gives every test its own empty default `FileCache`, so no test depends on what
  the machine's shared cache holds. `fake_spice` copies the oops Path and Frame
  registries, caches, and each class's waypoint and wayframe pools for every test, as
  rms-host-new-horizons does, until oops has one helper for it (the issue is cited in
  the fixture's comment).
- `test_files/metadata_versions/` holds the GO_0021 index label of every version of the
  metadata tree that the holdings keep (v1.0 to v4 and the current one), for the tests
  of `ratios_rotated()`.

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
  E4/EUROPA/C0374685452R describes the wrong layout. Those product labels are part of
  the holdings, which stay as they are.
  Whenever a file path is known, `header_dict` comes from a `VicarLabel` read from the
  file (`_read_vicar_header`), never from a `VicarImage`, which replaces the file's
  `HOST` and `BHOST` with the reading machine's; a `VicarImage` passed in with a file
  path has its header read again from the file. A `VicarImage` passed in without one
  (as `VicarImage.from_file()` makes it) gives only the header it holds,
  `fileinfo.label`, so its `HOST` and `BHOST` are the reading machine's. This stays
  until rms-vicar keeps the header as the file writes it (the issue is cited in the
  comment in `from_file()`). The array of a `VicarImage` with more than one band raises
  `VicarHostError` (`_array_2d`).
  Parsed content without a file path leaves out the path subfields, and a full read of a
  `Pds3Label` or `VicarLabel` without one raises `HostFormatError` (`_known_path`).
  `from_file()` refuses `select` other than None (`oops.Host._check_select`), as
  `from_index()` does. Neither takes `**kwargs`, so Python's own `TypeError` reports an
  unknown keyword.
- A PDS3 label that the fast parser rejects, for its syntax or for a value it cannot
  convert such as an hour of 25, is retried with the strict parser and then the loose
  one: `_parse_fileinfo` calls `oops.Host._read_fileinfo(..., pds3_fallback=True)`; a
  requested method starts the chain there, and if all fail the result is oops's
  `PdsHostError` ("cannot parse the PDS3 label <path>: method fast: ...; method strict:
  ...; method loose: ..."). The host parses label text as the file writes it, so it
  keeps the default `oops.Host._correct_label_text`. The parsers always run with
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
  hosts share. `label_items()` drops nulls (`oops.host.is_pds3_null()`: `UNK`, `UNKNOWN`,
  `N/A`, `NULL`, `NOT APPLICABLE`, `---`, and the empty string, any case, and lists of
  nothing else) and null list elements, and drops an `ORBIT_NUMBER` outside the tour
  (0-35): the labels of GO_0003-GO_0009 and GO_0016 give a pre-Jupiter frame the first
  two characters of its observation ID read as a 16-bit integer (12613 for `E1`). It
  renames the pre-tour `TIME_FROM_CLOSEST_APPROACH` and `SATELLITE_TIME_FROM_CLST_APR` to
  the tour names (`_LABEL_NAMES`).
  `drop_fills()` drops, from every source once its items have PDS3 names, the values
  that DATASET.CAT documents as not implemented, not applicable, or not available
  (sections C and D, for targets other than a planet or asteroid and for images without
  SPICE data): `_PDS3_FILLS` (`CONE_ANGLE`, `SMEAR_AZIMUTH`, `CENTER_RING_RADIUS`,
  `CENTRAL_BODY_DISTANCE`, `HORIZONTAL_PIXEL_SCALE`, `VERTICAL_PIXEL_SCALE`, and
  `SLANT_DISTANCE` 0, `SOLAR_DISTANCE` 777909100), a negative value of
  `_NEVER_NEGATIVE` (the pixel scales and slant distance, under their PDS3 and VICAR
  names; the index's cut-short negative values are dropped by the same rule until the
  index is corrected), the two times from closest approach
  that begin with `-73T` (`_UNKNOWN_TIME_ITEMS`), and each of `_ZERO_GROUPS` (the
  pointing; `CENTER_LATITUDE`/`CENTER_LONGITUDE`; `INTERCEPT_POINT_LATITUDE`/`LONGITUDE`,
  which the header's `LAT` and `LON` become; the three azimuths) when the source gives
  at least two of its items and all it gives are 0. One 0 of a group is kept.
  `header_items()` strips string padding, drops fill values (numeric -999, and
  `_VICAR_FILLS`: `PICNO` or `SPICE_C_ID` "?", `TCA` "CLOSEST TIME", `SMEAR` 0.1, the
  documented `RAD`, `CONE`, `SMRAZ`, `HSCL`, `VSCL`, `SLRANGE`, and `PLRANGE` 0 and
  `SOLRANGE` 777909100; time fields out of range), and gives every value it converts
  only under its PDS3 name (design D4): the `_VICAR_PDS3_PAIRS` keywords, `RIM` and
  `MOD91` as the clock count, the `SCET`/`ERT` fields as
  `IMAGE_TIME`/`EARTH_RECEIVED_TIME`, `RATE` and `GAIN` as `FRAME_DURATION` and
  `GAIN_MODE_ID`; then `drop_fills()`. `MOD10`, `MOD8`, `FIBE`, and `MOFIBE` are not
  wholly held by any PDS3 item and keep their names. `ORBIT_NUMBER` comes from an
  observation ID that starts with a tour orbit code; the IDs of the earlier phases start
  with a code of the phase instead (the Earth flybys `E1` and `E2`, Gaspra `GA`, Ida
  `ID`, Venus `VE`, optical navigation `OP`, and calibrations `PH`). `row_items()` cleans
  an index row the same way (numpy values to Python; nulls, -999, `ORBIT_NUMBER` outside
  0-35 (the index's -99), every value the table reader masks because the table declares
  it as the column's constant (the `-1` windows and ICT values, the 0 compression ratios)
  except the one real value of each `_REAL_NULL_CONSTANTS` column (`EXPOSURE_DURATION`
  0, `ON_CHIP_MOSAIC_FLAG` N, `PRODUCT_VERSION_ID` 0, so zero-exposure rows keep their
  exposure; another declared value in those columns is dropped), the tour's undeclared
  `SMEAR_MAGNITUDE` 0.1 (`_INDEX_PLACEHOLDERS`), and `drop_fills()`, as for a label).
  Following the declarations keeps a corrected table right. With `rotated_ratios` it
  renames the three compression ratios (`_INDEX_RATIO_COLUMNS`), which every SSI index
  writes in rotated columns. `ratios_rotated()` decides that once per table from the
  order of the column names in the index label alone (MIN, MAX, RATIO in byte order):
  every version of the metadata tree, each made again with a new
  `PRODUCT_CREATION_TIME`, kept that order and the rotated values, and the requested
  correction renames the columns. It gives the tour's `DATA_SET_ID` the labels' version
  1.1 (`_INDEX_VALUE_REPAIRS`) and the shortened `NTV_SAT_TIME_FROM_CLOSEST_APR` the
  labels' name (`_INDEX_NAMES`); both leave a corrected row as it is. Each index
  workaround's comment names its rms-data-projects item.
- `merge_items()` keeps every item of either source. Where both give an item and the
  values differ, `resolve()` keeps the source that `_RULES` names for that item, the same
  for every image, or the PDS3 label for an item with no rule; there is no per-image test.
  Each rule cites its evidence (survey counts, SPICE checks, DATASET.CAT) and the issue
  that lists every difference, each rule's reason, and the images on which its source is
  wrong. Only a value missing from the winning source falls to the other.
  `dict` is the merged result; `label_dict` and `header_dict` are present only when read.
  A missing value is never `None` or a placeholder: no key, no subfield.
- `_make_snapshot()` computes everything from the merged items, so a label read with its
  image gets the best values of both. The exposure, a time (`IMAGE_TIME`, or the clock
  count and `FRAME_DURATION`), and `TELEMETRY_FORMAT_ID` are required and raise
  `HostFormatError` naming the item and the file; an index row without a telemetry format
  (one read without its supplement) gets the full frame instead (`_INDEX_DEFAULT_MODE`).
- Missing pointing goes through `oops.Host._check_spice_frame()`, called before
  `_apply_overrides()` with the cadence and the overrides, only when no `frame` is given.
  It builds the frame and the time as the overrides will: `GLL_SCAN_PLATFORM` is
  relative to J2000, so the C kernels are read at the midtime after any cadence override
  (or the tracker's start or end), moved by a "T" timeshift and again by an "F" one.
  `from_file()` raises `oops.OopsSpiceError` at load, naming the image file. In
  `from_index()`, every
  row gives a Snapshot; a row whose pointing the C kernels lack gets an
  `oops.frame.SpiceErrorFrame`, which raises the same error whenever it is evaluated. A
  clock count that the clock kernel cannot convert raises `OopsSpiceError` from
  `_Galileo.tdb_from_sclk()`. The path subfields come from
  `oops.Host._apply_overrides()`, which leaves the file system as it is. `abspath` is
  passed to it only by a full `from_file()` read that retrieves the image file for its
  array: `_read_vicar_image()` calls `FCPath.retrieve()` and returns the local path with
  the array. An astrometry read, a `VicarImage` passed in, and `from_index()` give none,
  even for a file that is local or already in the cache.
- `IMAGE_TIME` is the middle of the exposure. Without it, the start time comes from
  `SPACECRAFT_CLOCK_START_COUNT` plus the shutter delay for the `FRAME_DURATION`
  (`_SHUTTER_DELAY`), one frame earlier for an `EXTENDED` exposure taken before the
  phase-2 flight software of the Jupiter tour (`_PHASE2_FIRST_RIM`).
- `GalileoSSI.from_index()` reads the index through `oops.Host._index_rows()`, which
  streams both tables (`PdsTable(defer_load=True).iter_dicts`, a few rows at a time),
  pairs, validates, and filters each row as the iteration reaches it; it is the host's
  only reader of tables. `MERGE_SUPPLEMENT_BY_FILE = False`: every SSI metadata index lists
  its files in its supplement's order (a survey of all 23, 29,324 rows), so the rows are
  paired by position, and a row whose files differ raises `HostIndexError` at that row.
  The host reads the Ring-Moon Systems Node metadata indexes and their supplements. The
  supplement writes `SPACECRAFT_CLOCK_START_COUNT` in SPICE
  form and the index in label form, so the default validation, which compares every
  shared column, would reject every real pair until the index is corrected;
  `validate=()` leaves the pairing's comparison of `FILE_SPECIFICATION_NAME` as the only
  check, which is right on a corrected table too.
  By default (`supplement=oops.Host.SUPPLEMENT`) the supplement beside the index is read,
  and it must exist (`FileNotFoundError`); `supplement=None` raises `HostIndexError`
  unless `row_dicts` is given. Index values win the merge, so the clock count stays in
  label form. The index label is read once, by `oops.Host._read_index_label()` (which
  returns `label` if given), for
  `ratios_rotated()`, also with `row_dicts`, so the label must exist then; it is passed
  to `_index_rows()` as `label`. A row's `FILE_SPECIFICATION_NAME` names the label; the
  Snapshot's `filepath` names the image beside it, with `.IMG` (or `.img`) in place of
  `.LBL`, as a full path or URL: `oops.Host._volumes_directory()` finds the volume
  directories from the index's own location (`volumes` in place of `metadata` in an RMS
  Node metadata tree, else the index's directory), and `oops.Host._index_file_path()`
  joins the row's `VOLUME_ID` and file specification to it (design D5).
  `_detect_in_index()` decides from the label's `INSTRUMENT_HOST_NAME` and
  `INSTRUMENT_ID`; for a label without them (the metadata versions before v4) it
  returns None, so that the host is identified from the first row, which both
  `oops.Host.from_index()` and `_index_rows()` read at call time. A table value that `PdsTable` cannot read as
  its column's number is yielded as a string with a UserWarning, not an error of the
  table; the docs say what oops does with malformed rows rather than promising a check.
- `GalileoSSI.from_index()` returns an iterator. It checks its arguments and reads the
  labels (and the first index row for a label that does not name the host, so a
  missing index table then raises at the call), then returns a generator that builds
  each row's Snapshot when the iteration reaches it. It must not be a generator function, which would defer every call-time
  error to the first `next()`.
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
