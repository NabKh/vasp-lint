"""Readers for the VASP input files of one calculation directory.

The parsers are deliberately small and dependency-free. They read only what
the rules need and never fail on a file they do not understand: a file that
cannot be parsed is reported as missing.
"""

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

TRUE_WORDS = {"T", "TRUE"}
FALSE_WORDS = {"F", "FALSE"}


@dataclass
class Incar:
    tags: dict = field(default_factory=dict)  # TAG -> raw string value
    lines: dict = field(default_factory=dict)  # TAG -> line number

    def has(self, key):
        return key in self.tags

    def str(self, key, default=None):
        return self.tags.get(key, default)

    def int(self, key, default=None):
        value = self.tags.get(key)
        if value is None:
            return default
        try:
            return int(float(value.split()[0]))
        except (ValueError, IndexError):
            return default

    def float(self, key, default=None):
        value = self.tags.get(key)
        if value is None:
            return default
        try:
            return float(value.split()[0].lower().replace("d", "e"))
        except (ValueError, IndexError):
            return default

    def bool(self, key, default=False):
        value = self.tags.get(key)
        if value is None:
            return default
        word = value.split()[0].strip(".").upper() if value.split() else ""
        if word in TRUE_WORDS:
            return True
        if word in FALSE_WORDS:
            return False
        return default

    def list(self, key):
        """Numbers of an array tag, with the N*value shorthand expanded."""
        value = self.tags.get(key)
        if value is None:
            return None
        out = []
        for token in value.split():
            if "*" in token:
                count, number = token.split("*", 1)
                out += [float(number)] * int(count)
            else:
                out.append(float(token))
        return out


@dataclass
class Poscar:
    species: list  # element symbols, or None for a VASP 4 POSCAR
    counts: list
    lattice: list  # three lattice vectors in Angstrom, scale applied

    @property
    def nions(self):
        return sum(self.counts)


@dataclass
class PotcarEntry:
    symbol: str  # e.g. "Mo_sv"
    element: str  # e.g. "Mo"
    zval: float
    enmax: float
    enmin: float


@dataclass
class Kpoints:
    mode: str  # "mesh", "auto", "line", "explicit", "generalized"
    mesh: list = None  # subdivisions, for mode "mesh"
    length: float = None  # R_k, for mode "auto"
    nkpts: int = None  # number of listed points, for mode "explicit"


def read_incar(path):
    incar = Incar()
    for number, raw in enumerate(Path(path).read_text(errors="replace").splitlines(), 1):
        line = re.split(r"[#!]", raw, maxsplit=1)[0]
        for statement in line.split(";"):
            if "=" not in statement:
                continue
            key, value = statement.split("=", 1)
            key = key.strip().upper()
            if key:
                incar.tags[key] = value.strip()
                incar.lines[key] = number
    return incar


def read_poscar(path):
    lines = Path(path).read_text(errors="replace").splitlines()
    try:
        scale = float(lines[1].split()[0])
        lattice = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
        tokens = lines[5].split()
        if all(t.isdigit() for t in tokens):
            species, counts = None, [int(t) for t in tokens]
        else:
            species = [re.match(r"[A-Za-z]+", t).group(0) for t in tokens]
            species = [s.split("_")[0] for s in species]
            counts = [int(t) for t in lines[6].split()[: len(species)]]
    except (IndexError, ValueError, AttributeError):
        return None
    if scale < 0:  # negative scale is the cell volume
        scale = (-scale / abs(_det(lattice))) ** (1 / 3)
    lattice = [[scale * x for x in row] for row in lattice]
    return Poscar(species=species, counts=counts, lattice=lattice)


def read_potcar(path):
    entries = []
    text = Path(path).read_text(errors="replace")
    for block in re.split(r"End of Dataset", text):
        titel = re.search(r"TITEL\s*=\s*\S+\s+(\S+)", block)
        zval = re.search(r"ZVAL\s*=\s*([-\d.]+)", block)
        enmax = re.search(r"ENMAX\s*=\s*([-\d.]+)", block)
        enmin = re.search(r"ENMIN\s*=\s*([-\d.]+)", block)
        if not (titel and zval and enmax):
            continue
        symbol = titel.group(1)
        element = re.match(r"[A-Z][a-z]?", symbol)
        entries.append(
            PotcarEntry(
                symbol=symbol,
                element=element.group(0) if element else symbol,
                zval=float(zval.group(1)),
                enmax=float(enmax.group(1)),
                enmin=float(enmin.group(1)) if enmin else float("nan"),
            )
        )
    return entries or None


def read_kpoints(path):
    lines = [line.strip() for line in Path(path).read_text(errors="replace").splitlines()]
    try:
        n = int(lines[1].split()[0])
        style = lines[2][:1].upper()
        if style == "L":
            return Kpoints(mode="line")
        if n == 0 and style in "GM":
            return Kpoints(mode="mesh", mesh=[int(x) for x in lines[3].split()[:3]])
        if n == 0 and style == "A":
            return Kpoints(mode="auto", length=float(lines[3].split()[0]))
        if n == 0:
            return Kpoints(mode="generalized")
        return Kpoints(mode="explicit", nkpts=n)
    except (IndexError, ValueError):
        return None


def _det(m):
    return (
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    )


def reciprocal_lengths(lattice):
    """|b_i| with b_i . a_j = delta_ij (no 2 pi), in 1/Angstrom."""
    a1, a2, a3 = lattice
    volume = _det(lattice)

    def cross(u, v):
        return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]

    return [math.dist([0, 0, 0], cross(u, v)) / abs(volume) for u, v in ((a2, a3), (a3, a1), (a1, a2))]


@dataclass
class Calculation:
    """Everything the rules can look at for one directory."""

    path: Path
    incar: Incar
    poscar: Poscar = None
    potcar: list = None
    kpoints: Kpoints = None
    has_kpoints_opt: bool = False
    files: set = field(default_factory=set)

    @property
    def species(self):
        if self.poscar and self.poscar.species:
            return self.poscar.species
        if self.potcar:
            return [p.element for p in self.potcar]
        return None

    def kmesh(self):
        """Subdivisions of the regular mesh VASP will use, or None if unknown."""
        if self.kpoints is not None:
            if self.kpoints.mode == "mesh":
                return self.kpoints.mesh
            if self.kpoints.mode == "auto" and self.poscar:
                b = reciprocal_lengths(self.poscar.lattice)
                return [int(max(1, self.kpoints.length * bi + 0.5)) for bi in b]
            return None
        if self.poscar:  # no KPOINTS file: VASP uses KSPACING (default 0.5)
            spacing = self.incar.float("KSPACING", 0.5)
            b = reciprocal_lengths(self.poscar.lattice)
            return [max(1, math.ceil(bi * 2 * math.pi / spacing)) for bi in b]
        return None


def load(directory):
    directory = Path(directory)
    incar = read_incar(directory / "INCAR")
    images = incar.int("IMAGES", 0)
    # for a nudged elastic band, the structures live in 00, 01, ...
    structure_dir = directory / "00" if images > 0 and not (directory / "POSCAR").exists() else directory

    def maybe(reader, name, where=directory):
        path = where / name
        return reader(path) if path.is_file() else None

    potcar = maybe(read_potcar, "POTCAR")
    if potcar is None and structure_dir != directory:
        potcar = maybe(read_potcar, "POTCAR", structure_dir)
    return Calculation(
        path=directory,
        incar=incar,
        poscar=maybe(read_poscar, "POSCAR", structure_dir),
        potcar=potcar,
        kpoints=maybe(read_kpoints, "KPOINTS"),
        has_kpoints_opt=(directory / "KPOINTS_OPT").is_file(),
        files={p.name for p in directory.iterdir()},
    )
