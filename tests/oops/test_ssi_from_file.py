##########################################################################################
# tests/oops/test_ssi_from_file.py
##########################################################################################
"""Tests for GalileoSSI.from_file(): which files it reads, the metadata it merges, the
path subfields, and the errors for a label or image file it cannot read."""

import datetime as dt
import pathlib
from typing import Any

import numpy as np
import oops
import pytest
import vicar
from conftest import (EXPOSURE_MS, SCLK_COUNT, SSI_PIXELS, fake_scs2e, image_array,
                      write_pds3_label, write_product, write_text, write_vicar_image)
from pdsparser import Pds3Label

from host_galileo.oops.ssi import GalileoSSI

# The trimmed test image, for the cutout window [2, 1, 399, 392]
TRIMMED = image_array()[1:400, 0:392]


def count_vicar_reads(monkeypatch: pytest.MonkeyPatch) -> list[pathlib.Path]:
    """Record the file of each read of an image file by the host; return the record."""

    reads: list[pathlib.Path] = []
    read_vicar = GalileoSSI._read_vicar

    def recording(filepath: Any, *, astrometry: bool) -> Any:
        reads.append(filepath)
        return read_vicar(filepath, astrometry=astrometry)

    monkeypatch.setattr(GalileoSSI, '_read_vicar', staticmethod(recording))
    return reads

##########################################################################################
# What is read
##########################################################################################


def test_label_gives_the_trimmed_array(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path))
    assert type(obs) is oops.observation.Snapshot
    assert np.array_equal(obs.data, TRIMMED)


def test_label_reads_the_image_file_once(tmp_path: pathlib.Path,
                                         monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_product(tmp_path)
    reads = count_vicar_reads(monkeypatch)
    GalileoSSI.from_file(path)
    assert [str(path) for path in reads] == [str(tmp_path / 'C0349632100R.IMG')]


def test_label_array_ignores_the_record_of_the_image_pointer(
        tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path, record=1))
    assert np.array_equal(obs.data, TRIMMED)


def test_label_finds_the_image_file_in_another_case(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path, image_name='c0349632100r.img'))
    assert obs.basename == 'C0349632100R.IMG'


def test_label_reads_the_window_in_label_order(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, CUT_OUT_WINDOW=[1, 2, 300, 200],
                         header={'CUT_OUT_WINDOW': None})
    obs = GalileoSSI.from_file(path)
    assert obs.data.shape == (300, 200)
    assert obs.data[0, 0] == image_array()[0, 1]


def test_full_fov_gives_the_full_array(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path), full_fov=True)
    assert obs.data.shape == (SSI_PIXELS, SSI_PIXELS)
    assert obs.fov is GalileoSSI._FOVS['HIS']


def test_label_astrometry_reads_no_image_file(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path)
    (tmp_path / 'C0349632100R.IMG').unlink()
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert 'data' not in obs.subfields
    assert 'header_dict' not in obs.subfields
    assert obs.filter == 'VIOLET'


def test_image_file_reads_no_label(tmp_path: pathlib.Path) -> None:
    write_product(tmp_path, OBSERVATION_ID='FROM_THE_LABEL')
    obs = GalileoSSI.from_file(tmp_path / 'C0349632100R.IMG')
    assert obs.dict['OBSERVATION_ID'] == 'G1GSGLOBAL02'
    assert 'label_dict' not in obs.subfields


def test_image_file_gives_its_array(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path))
    assert np.array_equal(obs.data, TRIMMED)


def test_image_file_astrometry_reads_only_the_header(
        tmp_path: pathlib.Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = write_vicar_image(tmp_path)

    def fail(*args: Any, **kwargs: Any) -> None:
        raise AssertionError('the VICAR image array was read')

    monkeypatch.setattr(vicar.VicarImage, '__init__', fail)
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert 'data' not in obs.subfields
    assert obs.filter == 'VIOLET'


def test_vicar_image_given_is_not_read_again(tmp_path: pathlib.Path,
                                             monkeypatch: pytest.MonkeyPatch) -> None:
    image = vicar.VicarImage(write_vicar_image(tmp_path), strict=False)
    reads = count_vicar_reads(monkeypatch)
    obs = GalileoSSI.from_file(image)
    assert np.array_equal(obs.data, TRIMMED)
    assert reads == []


def test_vicar_image_given_for_astrometry(tmp_path: pathlib.Path) -> None:
    image = vicar.VicarImage(write_vicar_image(tmp_path), strict=False)
    obs = GalileoSSI.from_file(image, astrometry=True)
    assert 'data' not in obs.subfields
    assert obs.header_dict['EXP'] == EXPOSURE_MS


def test_vicar_header_given_reads_the_image(tmp_path: pathlib.Path) -> None:
    header = vicar.VicarLabel(write_vicar_image(tmp_path), strict=False)
    assert np.array_equal(GalileoSSI.from_file(header).data, TRIMMED)


def test_vicar_header_given_for_astrometry(tmp_path: pathlib.Path) -> None:
    header = vicar.VicarLabel(write_vicar_image(tmp_path), strict=False)
    assert 'data' not in GalileoSSI.from_file(header, astrometry=True).subfields


def test_parsed_label_given(tmp_path: pathlib.Path) -> None:
    label = Pds3Label(write_product(tmp_path), method='strict')
    assert np.array_equal(GalileoSSI.from_file(label).data, TRIMMED)


def test_label_is_parsed_with_the_strict_method(tmp_path: pathlib.Path,
                                                monkeypatch: pytest.MonkeyPatch) -> None:
    methods: list[str] = []
    read_fileinfo = oops.Host._read_fileinfo.__func__

    def recording(cls: type, /, fileinfo: Any, **kwargs: Any) -> Any:
        methods.append(kwargs['pds3_options']['method'])
        return read_fileinfo(cls, fileinfo, **kwargs)

    monkeypatch.setattr(GalileoSSI, '_read_fileinfo', classmethod(recording))
    GalileoSSI.from_file(write_product(tmp_path), astrometry=True)
    GalileoSSI.from_file(write_product(tmp_path), astrometry=True,
                         pds3_options={'method': 'fast'})
    assert methods == ['strict', 'fast']


def test_from_file_initializes_ssi(tmp_path: pathlib.Path) -> None:
    GalileoSSI.from_file(write_product(tmp_path), astrometry=True)
    assert GalileoSSI._initialized


@pytest.mark.parametrize('writer', [write_product, write_vicar_image],
                         ids=['label', 'image'])
def test_host_from_file_dispatches_to_ssi(tmp_path: pathlib.Path, writer: Any) -> None:
    obs = oops.Host.from_file(writer(tmp_path))
    assert np.array_equal(obs.data, TRIMMED)

##########################################################################################
# The metadata
##########################################################################################


def test_label_dict_is_the_parsed_label(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path))
    assert obs.label_dict['IMAGE_TIME_fmt'] == '1996-06-26T08:46:04.890'


def test_header_dict_is_the_parsed_header(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path))
    assert obs.header_dict['EXP'] == EXPOSURE_MS


def test_header_wins_the_exposure(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, EXPOSURE_DURATION=45.83, header={'EXP': 45.8333})
    assert GalileoSSI.from_file(path).dict['EXPOSURE_DURATION'] == pytest.approx(45.8333)


def test_label_wins_the_target(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, TARGET_NAME='J RINGS', header={'TARGET': 'JUPITER'})
    assert GalileoSSI.from_file(path).dict['TARGET_NAME'] == 'J RINGS'


def test_header_gives_an_item_that_the_label_leaves_null(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, TELEMETRY_FORMAT_ID='UNK')
    assert GalileoSSI.from_file(path).dict['TELEMETRY_FORMAT_ID'] == 'HIS'


def test_header_items_take_their_pds3_names(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True)
    assert obs.dict['SPACECRAFT_CLOCK_START_COUNT'] == SCLK_COUNT
    assert obs.dict['IMAGE_TIME'] == dt.datetime(1996, 6, 26, 8, 46, 4, 890000)
    assert 'RIM' not in obs.dict
    assert obs.dict['MOD10'] == 0


def test_label_and_header_give_one_time_type(tmp_path: pathlib.Path) -> None:
    label = GalileoSSI.from_file(write_product(tmp_path), astrometry=True)
    image = GalileoSSI.from_file(tmp_path / 'C0349632100R.IMG', astrometry=True)
    assert type(label.dict['IMAGE_TIME']) is type(image.dict['IMAGE_TIME'])


def test_label_orbit_number_outside_the_tour_is_missing(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, ORBIT_NUMBER=12613)
    assert 'ORBIT_NUMBER' not in GalileoSSI.from_file(path, astrometry=True).dict


def test_image_file_of_an_observation_without_an_orbit_code(
        tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path, PA='OPNAV01A')
    obs = GalileoSSI.from_file(path, astrometry=True)
    assert obs.observation_id == 'OPNAV01A'
    assert 'ORBIT_NUMBER' not in obs.dict


def test_label_image_time_unk_starts_from_the_clock(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, IMAGE_TIME='UNK', header={'SCETYEAR': -32768})
    obs = GalileoSSI.from_file(path)
    assert obs.time[0] == pytest.approx(fake_scs2e(-77, '3496321:0:0:0') + 1.150,
                                        abs=1.e-6)
    assert obs.time_from_sclk is True


def test_image_file_without_a_time_starts_from_the_clock(tmp_path: pathlib.Path) -> None:
    # A phase-1 header: FIBE marks the exposure EXTENDED and RATE 3 gives 30.333 s
    path = write_vicar_image(tmp_path, SCETYEAR=-32768, RIM=2492218, MOFIBE=None,
                             FIBE='0001', RATE=3)
    obs = GalileoSSI.from_file(path, astrometry=True)
    expected = fake_scs2e(-77, '2492218:0:0:0') + 2.817 - 30.333
    assert obs.time[0] == pytest.approx(expected, abs=1.e-6)


def test_zero_exposure_loads(tmp_path: pathlib.Path) -> None:
    assert GalileoSSI.from_file(write_vicar_image(tmp_path, EXP=0.)).texp == 0.

##########################################################################################
# Path subfields
##########################################################################################


@pytest.mark.parametrize(('name', 'expected'), [
    ('filepath', 'C0349632100R.IMG'), ('basename', 'C0349632100R.IMG'),
    ('image_url', 'C0349632100R.IMG'), ('image_abspath', 'C0349632100R.IMG'),
    ('label_url', 'C0349632100R.LBL'), ('label_abspath', 'C0349632100R.LBL'),
])
def test_label_path_subfield(tmp_path: pathlib.Path, name: str, expected: str) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path))
    assert str(getattr(obs, name)).endswith(expected)


def test_label_astrometry_has_no_image_abspath(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_product(tmp_path), astrometry=True)
    assert 'image_abspath' not in obs.subfields
    assert obs.basename == 'C0349632100R.IMG'


@pytest.mark.parametrize('astrometry', [False, True])
def test_image_file_path_subfields(tmp_path: pathlib.Path, astrometry: bool) -> None:
    path = write_vicar_image(tmp_path)
    obs = GalileoSSI.from_file(path, astrometry=astrometry)
    assert obs.image_abspath == path.resolve()
    assert 'label_url' not in obs.subfields

##########################################################################################
# select
##########################################################################################


@pytest.mark.parametrize('select', [None, 0, -1, (0,)], ids=['none', '0', '-1', 'tuple'])
def test_select_gives_the_snapshot(tmp_path: pathlib.Path, select: Any) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True,
                               select=select)
    assert type(obs) is oops.observation.Snapshot


def test_select_slice_gives_a_tuple(tmp_path: pathlib.Path) -> None:
    obs = GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True,
                               select=slice(None))
    assert type(obs) is tuple


def test_select_out_of_range(tmp_path: pathlib.Path) -> None:
    with pytest.raises(oops.OopsIndexError, match='select index out of range: 1'):
        GalileoSSI.from_file(write_vicar_image(tmp_path), astrometry=True, select=1)

##########################################################################################
# Errors
##########################################################################################


def test_unknown_keyword(tmp_path: pathlib.Path) -> None:
    with pytest.raises(TypeError, match="unexpected keyword argument 'method'"):
        GalileoSSI.from_file(write_product(tmp_path),
                             method='strict')


def test_label_of_another_host(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, SPACECRAFT_NAME='CASSINI ORBITER')
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI PDS3 label'):
        GalileoSSI.from_file(path, astrometry=True)


def test_image_file_of_another_host(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path, MISSION='CASSINI')
    with pytest.raises(oops.host.WrongHostError, match='not a Galileo SSI VICAR file'):
        GalileoSSI.from_file(path)


def test_label_syntax_error(tmp_path: pathlib.Path) -> None:
    path = write_text(tmp_path / 'C0349632100R.LBL', ['PDS_VERSION_ID = PDS3',
                                                      'IMAGE_TIME = = 3'])
    with pytest.raises(oops.host.PdsHostError,
                       match=r'cannot parse the PDS3 label .*C0349632100R\.LBL'):
        GalileoSSI.from_file(path)


def test_label_value_that_cannot_be_converted(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, IMAGE_TIME='1996-02-30T00:00:00')
    with pytest.raises(oops.host.PdsHostError,
                       match=r'cannot parse the PDS3 label .*C0349632100R.LBL: day is '
                             r'out of range'):
        GalileoSSI.from_file(path)


def test_label_without_an_image_pointer(tmp_path: pathlib.Path) -> None:
    lines = ['PDS_VERSION_ID = PDS3', 'SPACECRAFT_NAME = "GALILEO ORBITER"',
             'INSTRUMENT_NAME = "SOLID STATE IMAGING SYSTEM"']
    path = write_text(tmp_path / 'C0349632100R.LBL', lines)
    with pytest.raises(oops.host.PdsHostError,
                       match=r'no \^IMAGE pointer in Galileo SSI PDS3 label '
                             '.*C0349632100R.LBL'):
        GalileoSSI.from_file(path)


def test_label_without_its_image_file(tmp_path: pathlib.Path) -> None:
    path = write_pds3_label(tmp_path)
    with pytest.raises(FileNotFoundError):
        GalileoSSI.from_file(path)


def test_truncated_image_file(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path)
    path.write_bytes(path.read_bytes()[:-1000])
    with pytest.raises(oops.host.VicarHostError,
                       match=r'cannot read the Galileo SSI image file '
                             '.*C0349632100R.IMG: cannot reshape'):
        GalileoSSI.from_file(path)


def test_image_file_that_is_not_vicar(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path)
    (tmp_path / 'C0349632100R.IMG').write_bytes(b'not a VICAR file ' * 100)
    with pytest.raises(oops.host.VicarHostError,
                       match=r'cannot read the Galileo SSI image file '
                             '.*C0349632100R.IMG: Missing LBLSIZE'):
        GalileoSSI.from_file(path)


def test_image_file_with_bands(tmp_path: pathlib.Path) -> None:
    path = write_vicar_image(tmp_path, array=np.zeros((2, 10, 10), dtype='uint8'))
    with pytest.raises(oops.host.VicarHostError,
                       match=r'cannot read the Galileo SSI image file '
                             '.*C0349632100R.IMG: VICAR array shape'):
        GalileoSSI.from_file(path)


def test_array_that_does_not_match_its_label(tmp_path: pathlib.Path) -> None:
    path = write_product(tmp_path, lines=400)
    with pytest.raises(oops.host.HostFormatError,
                       match=r'C0349632100R.IMG holds a 800x800 array, where its label '
                             r'.*C0349632100R.LBL gives 400x800'):
        GalileoSSI.from_file(path)


@pytest.mark.parametrize(('changes', 'header', 'message'), [
    ({'EXPOSURE_DURATION': None}, {'EXP': None}, "'EXPOSURE_DURATION'"),
    ({'TELEMETRY_FORMAT_ID': 'XED'}, {}, "'XED'"),
    ({'IMAGE_TIME': None, 'FRAME_DURATION': 5.}, {'SCETYEAR': -32768},
     'unrecognized Galileo SSI FRAME_DURATION: 5.0'),
], ids=['no-exposure', 'unknown-telemetry-format', 'unknown-frame-duration'])
def test_metadata_that_cannot_be_interpreted(tmp_path: pathlib.Path,
                                             changes: dict[str, Any],
                                             header: dict[str, Any],
                                             message: str) -> None:
    path = write_product(tmp_path, header=header, **changes)
    with pytest.raises(oops.host.HostFormatError,
                       match=f'cannot interpret the Galileo SSI product '
                             f'.*C0349632100R.LBL: {message}'):
        GalileoSSI.from_file(path)

##########################################################################################
