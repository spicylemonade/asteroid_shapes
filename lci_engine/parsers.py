"""
Data ingestion and parsing modules for ALCDEF lightcurve data and MPCORB orbital elements.

References:
    - ALCDEF format: https://alcdef.org
    - MPCORB format: https://minorplanetcenter.net/iau/info/MPOrbitFormat.html
"""

import re
import gzip
import zipfile
import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple


@dataclass
class LightcurveSession:
    """A single observing session (lightcurve) from ALCDEF."""
    asteroid_number: int
    asteroid_name: str
    mpc_desig: str
    session_date: str
    session_time: str
    observer: str
    mpc_code: str
    telescope: str
    phase_angle: float  # degrees
    pabl: float  # phase angle bisector longitude
    pabb: float  # phase angle bisector latitude
    filter_band: str
    mag_band: str
    differ_mags: bool
    jd: np.ndarray  # Julian dates
    mag: np.ndarray  # magnitudes
    mag_err: np.ndarray  # magnitude uncertainties
    metadata: dict = field(default_factory=dict)


@dataclass
class OrbitalElements:
    """Orbital elements from MPCORB."""
    designation: str  # packed designation
    number: int  # asteroid number (0 if unnumbered)
    name: str
    H: float  # absolute magnitude
    G: float  # slope parameter
    epoch: str  # packed epoch
    M: float  # mean anomaly (deg)
    omega: float  # argument of perihelion (deg)
    Omega: float  # longitude of ascending node (deg)
    i: float  # inclination (deg)
    e: float  # eccentricity
    n: float  # mean daily motion (deg/day)
    a: float  # semi-major axis (AU)
    num_obs: int
    num_opp: int
    arc: str
    rms: float


def parse_alcdef_file(content: str) -> List[LightcurveSession]:
    """Parse an ALCDEF format text file into a list of LightcurveSession objects.

    Each file may contain multiple lightcurve sessions separated by
    STARTMETADATA/ENDMETADATA blocks followed by DATA lines.
    """
    sessions = []
    blocks = content.split('STARTMETADATA')

    for block in blocks[1:]:  # Skip text before first STARTMETADATA
        parts = block.split('ENDMETADATA', 1)
        if len(parts) < 2:
            continue

        metadata_text, data_text = parts

        # Parse metadata
        meta = {}
        for line in metadata_text.strip().split('\n'):
            line = line.strip()
            if '=' in line:
                key, _, val = line.partition('=')
                meta[key.strip()] = val.strip()

        # Parse data points
        jd_list, mag_list, err_list = [], [], []
        for line in data_text.strip().split('\n'):
            line = line.strip()
            if line.startswith('DATA='):
                parts = line[5:].split('|')
                if len(parts) >= 3:
                    try:
                        jd_list.append(float(parts[0]))
                        mag_list.append(float(parts[1]))
                        err_list.append(float(parts[2]))
                    except ValueError:
                        continue

        if not jd_list:
            continue

        # Parse numeric fields safely
        def safe_float(val, default=0.0):
            try:
                return float(val)
            except (ValueError, TypeError):
                return default

        session = LightcurveSession(
            asteroid_number=int(meta.get('OBJECTNUMBER', '0')),
            asteroid_name=meta.get('OBJECTNAME', ''),
            mpc_desig=meta.get('MPCDESIG', ''),
            session_date=meta.get('SESSIONDATE', ''),
            session_time=meta.get('SESSIONTIME', ''),
            observer=meta.get('OBSERVERS', ''),
            mpc_code=meta.get('MPCCODE', ''),
            telescope=meta.get('TELESCOPE', ''),
            phase_angle=safe_float(meta.get('PHASE', '0')),
            pabl=safe_float(meta.get('PABL', '0')),
            pabb=safe_float(meta.get('PABB', '0')),
            filter_band=meta.get('FILTER', ''),
            mag_band=meta.get('MAGBAND', ''),
            differ_mags=meta.get('DIFFERMAGS', 'FALSE').upper() == 'TRUE',
            jd=np.array(jd_list),
            mag=np.array(mag_list),
            mag_err=np.array(err_list),
            metadata=meta,
        )
        sessions.append(session)

    return sessions


def load_alcdef_asteroid(zip_path: str, asteroid_number: int = None,
                         asteroid_name: str = None) -> List[LightcurveSession]:
    """Load ALCDEF data for a specific asteroid from the zip archive.

    Args:
        zip_path: Path to ALCDEF_ALL.zip
        asteroid_number: Asteroid number (e.g., 433 for Eros)
        asteroid_name: Asteroid name (e.g., 'Eros')

    Returns:
        List of LightcurveSession objects
    """
    zf = zipfile.ZipFile(zip_path)

    target_files = []
    for fname in zf.namelist():
        if asteroid_number is not None:
            if fname.startswith(f'ALCDEF_{asteroid_number}_'):
                target_files.append(fname)
        elif asteroid_name is not None:
            if asteroid_name.lower() in fname.lower():
                target_files.append(fname)

    all_sessions = []
    for fname in target_files:
        content = zf.read(fname).decode('utf-8', errors='replace')
        sessions = parse_alcdef_file(content)
        all_sessions.extend(sessions)

    zf.close()
    return all_sessions


def _unpack_mpc_epoch(packed: str) -> str:
    """Convert MPC packed epoch (e.g., 'K25BL') to readable date."""
    if len(packed) < 5:
        return packed

    century_map = {'I': '18', 'J': '19', 'K': '20'}
    char_map = {chr(i): i - 55 for i in range(65, 91)}  # A=10, B=11, ...
    char_map.update({str(i): i for i in range(10)})

    try:
        century = century_map.get(packed[0], '20')
        year = century + packed[1:3]
        month = char_map.get(packed[3], 1)
        day = char_map.get(packed[4], 1)
        return f"{year}-{month:02d}-{day:02d}"
    except (KeyError, ValueError, IndexError):
        return packed


def parse_mpcorb_line(line: str) -> Optional[OrbitalElements]:
    """Parse a single line of MPCORB.DAT into OrbitalElements."""
    if len(line) < 160:
        return None

    try:
        designation = line[0:7].strip()
        H = float(line[8:13].strip()) if line[8:13].strip() else 99.0
        G = float(line[14:19].strip()) if line[14:19].strip() else 0.15
        epoch = line[20:25].strip()
        M = float(line[26:35].strip())
        omega = float(line[37:46].strip())
        Omega = float(line[48:57].strip())
        i = float(line[59:68].strip())
        e = float(line[70:79].strip())
        n = float(line[80:91].strip())
        a = float(line[92:103].strip())

        # Parse number from designation
        num = 0
        try:
            num = int(designation)
        except ValueError:
            pass

        # Parse additional fields
        num_obs = 0
        num_opp = 0
        arc = ''
        rms = 0.0
        name = ''

        try:
            ref = line[107:116].strip()
            num_obs = int(line[117:122].strip()) if line[117:122].strip() else 0
            num_opp = int(line[123:126].strip()) if line[123:126].strip() else 0
            arc = line[127:136].strip()
            rms_str = line[137:141].strip()
            rms = float(rms_str) if rms_str else 0.0
        except (ValueError, IndexError):
            pass

        try:
            name = line[175:194].strip()
        except IndexError:
            name = ''

        return OrbitalElements(
            designation=designation, number=num, name=name,
            H=H, G=G, epoch=epoch, M=M, omega=omega, Omega=Omega,
            i=i, e=e, n=n, a=a,
            num_obs=num_obs, num_opp=num_opp, arc=arc, rms=rms,
        )
    except (ValueError, IndexError):
        return None


class MPCORBDatabase:
    """Interface to the MPCORB orbital elements database."""

    def __init__(self, gz_path: str):
        self.gz_path = gz_path
        self._by_number: Dict[int, OrbitalElements] = {}
        self._by_name: Dict[str, OrbitalElements] = {}
        self._loaded = False

    def load(self, max_records: int = None):
        """Parse MPCORB.DAT.gz into memory."""
        if self._loaded:
            return

        count = 0
        with gzip.open(self.gz_path, 'rt', errors='replace') as f:
            in_data = False
            for line in f:
                if not in_data:
                    if line.startswith('-----'):
                        in_data = True
                    continue

                elem = parse_mpcorb_line(line)
                if elem is not None:
                    if elem.number > 0:
                        self._by_number[elem.number] = elem
                    if elem.name:
                        self._by_name[elem.name.lower()] = elem
                    count += 1

                    if max_records and count >= max_records:
                        break

        self._loaded = True

    def get_by_number(self, number: int) -> Optional[OrbitalElements]:
        """Look up orbital elements by asteroid number."""
        if not self._loaded:
            self.load()
        return self._by_number.get(number)

    def get_by_name(self, name: str) -> Optional[OrbitalElements]:
        """Look up orbital elements by asteroid name."""
        if not self._loaded:
            self.load()
        return self._by_name.get(name.lower())

    @property
    def count(self) -> int:
        if not self._loaded:
            self.load()
        return len(self._by_number) + len(self._by_name) - len(
            set(self._by_number.keys()) &
            {v.number for v in self._by_name.values() if v.number > 0}
        )


def estimate_diameter_from_H(H: float, albedo: float = 0.15) -> float:
    """Estimate asteroid diameter in km from absolute magnitude H and albedo.

    D = 1329 / sqrt(albedo) * 10^(-H/5)

    Reference: Bowell1989
    """
    return 1329.0 / np.sqrt(albedo) * 10.0 ** (-H / 5.0)
