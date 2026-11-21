##########################################################################################
# host_galileo/oops/_galileo.py
##########################################################################################

import re

import numpy as np

import julian
import spicedb
import cspyce
import oops

from oops.body import Body
from spicedb   import get_spice_filecache_prefix

# Mission targets, rough time divisions
TIMELINE = [
    {'ET': 0, 'UTC': '1989-10-18T12:00:00.00', 'target': 'VENUS',   'moons': False},
    {'ET': 0, 'UTC': '1990-07-08T12:00:00.00', 'target': 'EARTH',   'moons': True},
    {'ET': 0, 'UTC': '1991-05-29T12:00:00.00', 'target': 'GASPRA',  'moons': False},
    {'ET': 0, 'UTC': '1992-07-08T12:00:00.00', 'target': 'EARTH',   'moons': True},
    {'ET': 0, 'UTC': '1993-04-28T12:00:00.00', 'target': 'IDA',     'moons': True},
    {'ET': 0, 'UTC': '1994-07-16T12:00:00.00', 'target': 'SL9',     'moons': False},
    {'ET': 0, 'UTC': '1994-07-22T12:00:00.00', 'target': 'JUPITER', 'moons': True}
]

for _phase in TIMELINE:
    _phase['ET'] = julian.tdb_from_tai(julian.tai_from_iso(_phase['UTC']))

# The fixed set of kernels furnished by _Galileo.load_kernels(), relative to the SPICE
# file cache prefix
_KERNELS = [
    'General/LSK/naif0012.tls',
    'Galileo/IK/gll36001.ti',
    'Galileo/FK/gll_v0.tf',
    'Galileo/SPK/de421.bsp',
    'Galileo/SPK/de432s.bsp',
    'Galileo/CK/ckc03b_plt.bc',
    'Galileo/CK/ckc09b_plt.bc',
    'Galileo/CK/ckc10b_plt.bc',
    'Galileo/CK/ckc20f_plt.bc',
    'Galileo/CK/ckc21f_plt.bc',
    'Galileo/CK/ckc22f_plt.bc',
    'Galileo/CK/ckc23f_plt.bc',
    'Galileo/CK/ckc30f_plt.bc',
    'Galileo/CK/cke04b_plt.bc',
    'Galileo/CK/cke06b_plt.bc',
    'Galileo/CK/cke11b_plt.bc',
    'Galileo/CK/cke12f_plt.bc',
    'Galileo/CK/cke14f_plt.bc',
    'Galileo/CK/cke15f_plt.bc',
    'Galileo/CK/cke16f_plt.bc',
    'Galileo/CK/cke17f_plt.bc',
    'Galileo/CK/cke18f_plt.bc',
    'Galileo/CK/cke19f_plt.bc',
    'Galileo/CK/cke26f_plt.bc',
    'Galileo/CK/ckg01b_plt.bc',
    'Galileo/CK/ckg02b_plt.bc',
    'Galileo/CK/ckg07b_plt.bc',
    'Galileo/CK/ckg08b_plt.bc',
    'Galileo/CK/ckg28f_plt.bc',
    'Galileo/CK/ckg29f_plt.bc',
    'Galileo/CK/cki24f_plt.bc',
    'Galileo/CK/cki25f_plt.bc',
    'Galileo/CK/cki27f_plt.bc',
    'Galileo/CK/cki31f_plt.bc',
    'Galileo/CK/cki32f_plt.bc',
    'Galileo/CK/ckj0cav3_plt.bc',
    'Galileo/CK/ckj0cduh_plt.bc',
    'Galileo/CK/ckj0cv3_plt.bc',
    'Galileo/CK/ckj0eav3_plt.bc',
    'Galileo/CK/ckj0ebv3_plt.bc',
    'Galileo/CK/ckj0ecv3_plt.bc',
    'Galileo/CK/ckjaap_plt.bc',
    'Galileo/CK/ckjaav3_plt.bc',
    'Galileo/CK/ckjabp_plt.bc',
    'Galileo/CK/ckjabv3_plt.bc',
    'Galileo/CK/gll_plt_pre_1990_v00.bc',
    'Galileo/CK/gll_plt_pre_1991_v00.bc',
    'Galileo/CK/gll_plt_pre_1992_v00.bc',
    'Galileo/CK/gll_plt_pre_1993_v00.bc',
    'Galileo/CK/gll_plt_pre_1994_v00.bc',
    'Galileo/CK/gll_plt_pre_1995_v00.bc',
    'Galileo/CK/gll_plt_pre_1996_v00.bc',
    'Galileo/CK/gll_plt_pre_1997_v00.bc',
    'Galileo/CK/gll_plt_pre_1998_v00.bc',
    'Galileo/CK/gll_plt_pre_1999_v00.bc',
    'Galileo/CK/gll_plt_pre_2000_v00.bc',
    'Galileo/CK/gll_plt_pre_2001_v00.bc',
    'Galileo/SPK/de421.bsp',
    'Galileo/SPK/de432s.bsp',
    'Galileo/SPK/gll_951120_021126_raj2007.bsp',
    'Galileo/SPK/gll_951120_021126_raj2021.bsp',
    'Galileo/SPK/s000131a.bsp',
    'Galileo/SPK/s000615a.bsp',
    'Galileo/SPK/s020128a.bsp',
    'Galileo/SPK/s030916a.bsp',
    'Galileo/SPK/s960730a.bsp',
    'Galileo/SPK/s970311a.bsp',
    'Galileo/SPK/s971125a.bsp',
    'Galileo/SPK/s980326a.bsp',
]


class _Galileo:
    """An instance-free class to hold Galileo-specific parameters.

    Attributes:
        START_TIME (str): The start of the mission as an ISO date.
        STOP_TIME (str): The end of the mission as an ISO date.
        MONTHS (int): The number of equal "months" into which the mission is divided for
            the purpose of loading kernels.
        SPACECRAFT_ID (int): The NAIF ID of the Galileo spacecraft.
        SCLK_KERNEL (str): The spacecraft clock kernel, relative to the SPICE root.
        TDB0 (float): The mission start time in seconds TDB.
        TDB1 (float): The mission stop time in seconds TDB.
        DTDB (float): The duration of one "month" in seconds.
        SLOP (float): The margin, in seconds, by which each month is extended when
            deciding which kernels apply to it.
        CK_LOADED (numpy.ndarray): Boolean array with one flag per month, True if the C
            kernels for that month have been furnished.
        CK_LIST (numpy.ndarray): Object array holding, for each month, the list of
            KernelInfo objects for the C kernels needed within that month.
        CK_DICT (dict[str, KernelInfo]): The furnished C kernels, keyed by filespec.
        SPK_LOADED (numpy.ndarray): Boolean array with one flag per month, True if the SP
            kernels for that month have been furnished.
        SPK_LIST (numpy.ndarray): Object array holding, for each month, the list of
            KernelInfo objects for the SP kernels needed within that month.
        SPK_DICT (dict[str, KernelInfo]): The furnished SP kernels, keyed by filespec.
        loaded_instruments (list[str]): The names of the instruments whose kernels have
            been loaded.
        initialized (bool): True after :meth:`initialize` has been called.
        sclk_loaded (bool): True after the spacecraft clock kernel has been furnished.
    """

    START_TIME = '1989-10-18'
    STOP_TIME  = '2003-09-21'
    MONTHS = 167        # 14 years * 12 months/year - 1 month

    SPACECRAFT_ID = -77
    SCLK_KERNEL = 'Galileo/SCLK/mk00062a.tsc'

    TDB0 = julian.tdb_from_tai(julian.tai_from_iso(START_TIME))
    TDB1 = julian.tdb_from_tai(julian.tai_from_iso(STOP_TIME))
    DTDB = (TDB1 - TDB0) / MONTHS
    SLOP = 43200.

    # We load CK and SPK files on a very rough month-by-month basis. This is simpler than
    # a more granular approach involving detailed calendar calculations. We divide the
    # period October 18, 1989 to September 21, 2003 up into "months" of equal length.
    # Given any TDB, we quickly determine the month within which it falls. Each month is
    # associated with a list of kernels that should be loaded whenever information is
    # needed about any time within that month +/- 12 hours. Galileo does not yet use this
    # mechanism: initialize() leaves the monthly lists empty and load_kernels() furnishes
    # a fixed set of kernels instead.

    CK_LOADED = np.zeros(MONTHS, dtype='bool')      # True if month was loaded
    CK_LIST   = np.empty(MONTHS, dtype='object')    # Kernels needed by month
    CK_DICT   = {}      # Dictionary keyed by filespec returns kernel info
                        # object, but only if loaded.

    SPK_LOADED  = np.zeros(MONTHS, dtype='bool')
    SPK_LIST    = np.empty(MONTHS, dtype='object')
    SPK_DICT    = {}

    loaded_instruments = []

    initialized = False
    sclk_loaded = False

    @staticmethod
    def initialize(planets=None, asof=None, irregulars=True):
        """Initialize the Galileo mission internals.

        After the first call, later calls to this function are ignored.

        Parameters:
            planets (int | collections.abc.Sequence[int], optional): The number of one
                planet, from 1 to 9, or a sequence of such numbers, to pass to
                :meth:`~oops.Body.define_solar_system`. None or 0 means all.
            asof (str, optional): A date, passed to
                :meth:`~oops.Body.define_solar_system` so that the database chooses the
                kernels that existed before it; None to ignore. Where the database does
                not list a kernel of some kind from before the date, it chooses as if no
                date were given. The leap-second kernel that this method loads first
                ignores the date, and so does the fixed set that :meth:`load_kernels`
                furnishes later, which overrides these kernels wherever the two overlap.
            irregulars (bool, optional): True to include the irregular satellites; False
                otherwise.
        """

        if _Galileo.initialized:
            return

        # Make sure the leap seconds have been loaded
        oops.spice.load_leap_seconds()

        # Define some important paths and frames
        Body.define_solar_system(_Galileo.START_TIME, _Galileo.STOP_TIME, asof=asof,
                                 planets=planets, irregulars=irregulars)

        _ = oops.path.SpicePath('GLL', 'JUPITER')

        spicedb.open_db()

        # This means no SPK will ever be loaded; handling is manual
        _Galileo.initialize_kernels([], _Galileo.SPK_LIST)
        _Galileo.SPK_LOADED = np.ones(_Galileo.MONTHS, dtype='bool')

        # This means no CK will ever be loaded; handling is manual
        _Galileo.initialize_kernels([], _Galileo.CK_LIST)
        _Galileo.CK_LOADED = np.ones(_Galileo.MONTHS, dtype='bool')

        spicedb.close_db()

        _Galileo.initialized = True

    @staticmethod
    def reset():
        """Reset the internal parameters.

        Can be useful for debugging.
        """
        _Galileo.loaded_instruments = []

        _Galileo.CK_LOADED = np.zeros(_Galileo.MONTHS, dtype='bool')
        _Galileo.CK_LIST = np.empty(_Galileo.MONTHS, dtype='object')
        _Galileo.CK_DICT = {}

        _Galileo.SPK_LOADED = np.zeros(_Galileo.MONTHS, dtype='bool')
        _Galileo.SPK_LIST = np.empty(_Galileo.MONTHS, dtype='object')
        _Galileo.SPK_DICT = {}

        _Galileo.initialized = False
        _Galileo.sclk_loaded = False

    @staticmethod
    def load_kernels():
        """Furnish the fixed set of SPICE kernels needed by the Galileo mission."""

        _Galileo.load_sclk()

        paths = get_spice_filecache_prefix().retrieve(_KERNELS)
        for path in paths:
            cspyce.furnsh(path)

    @staticmethod
    def load_sclk():
        """Furnish the leap-second and Galileo spacecraft clock kernels if not yet loaded.

        Unlike :meth:`initialize`, this loads nothing else, so a clock conversion costs
        two text kernels rather than the full mission kernel set.
        """

        if _Galileo.sclk_loaded:
            return

        oops.spice.load_leap_seconds()

        path = get_spice_filecache_prefix().retrieve(_Galileo.SCLK_KERNEL)
        cspyce.furnsh(path)

        _Galileo.sclk_loaded = True

    @staticmethod
    def tdb_from_sclk(count, *, filepath=None):
        """Convert a Galileo spacecraft clock count to seconds TDB.

        The kernels needed are furnished on the first call through :meth:`load_sclk`, so
        no call to :meth:`initialize` is required. The result is the start of the frame
        that the count names, which is not the start of the exposure: the shutter opens
        up to 6.5 s later, depending on the frame duration, or one frame earlier for an
        extended exposure.

        Parameters:
            count (str): The clock count as one to four integer fields (RIM, mod-91,
                mod-10 and mod-8), separated by any non-alphanumeric characters, as in
                the label form "RRRRRRRR.MM" or the SPICE form "RRRRRRRR:MM:S:T".
                Missing trailing fields are taken as zero.
            filepath (FCPath, optional): The image file whose count this is, named in an
                error.

        Returns:
            float: The time in seconds TDB.

        Raises:
            OopsValueError: If the count has more than four fields or a field is not an
                unsigned integer.
            OopsSpiceError: If the clock kernel cannot convert the count, such as one
                outside the clock's partitions. Its message is "no time for Galileo clock
                count '<count>' for <file>: <SPICE error code>", and its cause is the
                SPICE exception. A failure of the conversion whose message lacks a SPICE
                error code is raised as it is.
        """

        _Galileo.load_sclk()

        fields = re.split(r'[^0-9A-Za-z]+', count.strip())
        if len(fields) > 4 or not all(field.isdigit() for field in fields):
            raise oops.OopsValueError('invalid Galileo spacecraft clock count: '
                                      f'{count!r}')

        fields += ['0'] * (4 - len(fields))
        try:
            return cspyce.scs2e(_Galileo.SPACECRAFT_ID,
                                ':'.join(str(int(field)) for field in fields))
        except oops.OopsSpiceError.CSPYCE_ERRORS as err:
            if not oops.OopsSpiceError.is_spice_error(err):
                raise
            raise oops.OopsSpiceError.from_spice_error(
                err, description=f'no time for Galileo clock count {count!r}',
                filepath=filepath) from err

    ######################################################################################
    # Initialize the kernel lists
    ######################################################################################

    @staticmethod
    def initialize_kernels(kernels, lists):
        """Initialize the monthly lists of KernelInfo objects.

        After initialization, `lists[m]` holds the KernelInfo objects needed within the
        specified month.

        Parameters:
            kernels (list[KernelInfo]): The KernelInfo objects to distribute among the
                months of the mission.
            lists (numpy.ndarray): Object array with one entry per month; each entry is
                replaced by the list of KernelInfo objects that apply within that month,
                extended by :attr:`SLOP` at each end.
        """
        for i in range(_Galileo.MONTHS):
            lists[i] = []

        for kernel in kernels:

            # Find the range of months applicable, extended by 12 hours
            t0 = cspyce.str2et(kernel.start_time) - _Galileo.SLOP
            t1 = cspyce.str2et(kernel.stop_time)  + _Galileo.SLOP

            m1 = int((t0 - _Galileo.TDB0) // _Galileo.DTDB)
            m2 = int((t1 - _Galileo.TDB0) // _Galileo.DTDB) + 1

            m1 = max(m1, 0)     # ignore time limits outside mission duration
            m2 = min(m2, _Galileo.MONTHS - 1)

            # Add this kernel to each month's list
            for m in range(m1, m2+1):
                lists[m] += [kernel]

    ######################################################################################
    # Routines for managing the loading other kernels
    ######################################################################################

    @staticmethod
    def load_instruments(instruments=None, asof=None):
        """Load the SPICE kernels for the Galileo mission.

        This also defines the basic paths and frames.

        It is generally only to be called once.

        Parameters:
            instruments (list[str], optional): The names of the instruments whose kernels
                to load. On the first call, 'SSI' is always included.
            asof (str, optional): A date or date-time in ISO format, so that the
                database chooses the instrument, frame, and spacecraft clock kernels
                that existed before it; where it does not list one of some kind from
                before the date, it chooses as if no date were given. If None, the most
                recent versions are loaded. The fixed set that :meth:`load_kernels`
                furnishes later overrides these kernels wherever the two overlap.
        """

        # Copy so that the caller's list is never modified
        instruments = [] if instruments is None else list(instruments)

        # Load the default instruments on the first pass
        if _Galileo.loaded_instruments == []:
            instruments += ['SSI']

        # Skip instruments that are already loaded, and duplicates
        instruments = [inst for inst in dict.fromkeys(instruments)
                       if inst not in _Galileo.loaded_instruments]

        # On later calls, return quickly if there's nothing to do
        if instruments == []:
            return

        # Check the formatting of the "as of" date
        if asof is not None:
            (day, sec) = julian.day_sec_from_iso(asof)
            asof = julian.ymdhms_format_from_day_sec(day, sec)

        # Furnish the instrument kernels with their frame and spacecraft clock kernels
        spicedb.open_db()
        _ = spicedb.furnish_inst(_Galileo.SPACECRAFT_ID, inst=instruments, asof=asof)
        spicedb.close_db()

        _Galileo.loaded_instruments += instruments

    ######################################################################################
    # Routines for managing text kernel information
    ######################################################################################

    @staticmethod
    def spice_instrument_kernel(inst, asof=None):
        """A dictionary containing the Instrument Kernel information.

        Also furnishes it for use by the SPICE tools.

        Parameters:
            inst (str | list | tuple): One of "SSI", etc.
            asof (str, optional): An optional date in the past, in ISO date or date-time
                format. If provided, the database chooses the instrument kernel that
                existed before that date, or chooses as if no date were given if none
                did. Otherwise, the most recent kernel is chosen.

        Returns:
            tuple[dict, str]: The dictionary generated by :func:`textkernel.from_file`
            and the name of the kernel.
        """
        if asof is not None:
            (day,sec) = julian.day_sec_from_iso(asof)
            asof = julian.ymdhms_format_from_day_sec(day, sec)

        spicedb.open_db()
        kernel_info = spicedb.select_inst(_Galileo.SPACECRAFT_ID, types='IK', inst=inst,
                                          asof=asof)
        spicedb.furnish_kernels(kernel_info, fast=True)
        spicedb.close_db()

        return (spicedb.as_dict(kernel_info), spicedb.as_names(kernel_info)[0])

    @staticmethod
    def used_kernels(time, inst, return_all_planets=False, ck=True):
        """The kernels associated with a Galileo observation.

        The list covers a selected range of times.

        Parameters:
            time (tuple[float, float]): A (start, stop) tuple of times in seconds TDB.
            inst (str | list | tuple): The instrument name, e.g., 'ssi'.
            return_all_planets (bool, optional): Include kernels for all planets, not just
                the target of the mission phase containing the start time.
            ck (bool, optional): True (default) to include CK (pointing) kernels; False to
                exclude them, e.g. for an observation whose pointing came from a custom
                frame rather than SPICE, where any furnished CK is unrelated to how the
                observation was actually pointed.

        Returns:
            list[str]: The basenames of the furnished kernels that apply to the
            observation.
        """
        # Determine targets based on mission phase, or return_all_planets
        if return_all_planets:
            targets = ['MERCURY', 'VENUS', 'EARTH', 'MARS', 'JUPITER', 'SATURN', 'URANUS',
                       'NEPTUNE']
            moons = [False, False, True, True, True, True, True, True]
        else:
            times = np.array([phase['ET'] for phase in TIMELINE])
            i = np.max(np.where(time[0] >= times))
            targets = [TIMELINE[i]['target']]
            moons = [TIMELINE[i]['moons']]

        # Specify desired body ids
        bodies = []
        for target, moon in zip(targets, moons, strict=True):
            if target == 'SL9':       # TODO find a kernel for SL9
                continue
            body = cspyce.bodn2c(target)
            barycenter = np.trunc(body/100).astype(int)
            bodies += [barycenter, body]

            # Add moons if requested. Body keeps no list of the Moon, and none for the
            # asteroids, so Earth gets the Moon and an asteroid gets nothing.
            if moon:
                if target == 'EARTH':
                    bodies += [301]
                else:
                    bodies += getattr(Body, target + '_MOONS_LOADED', [])

        types = None
        if not ck:
            types = [t for t in spicedb.KERNEL_TYPE_SORT_ORDER if t != 'CK']

        # Return the relevant used kernels
        return spicedb.used_basenames(types=types, time=time, inst=inst,
                                      sc=_Galileo.SPACECRAFT_ID, bodies=bodies)

##########################################################################################
