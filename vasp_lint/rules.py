"""The checks. Every rule links the VASP wiki page its advice comes from."""

import math
from dataclasses import dataclass

from .inputs import read_poscar

WIKI = "https://www.vasp.at/wiki/index.php/"

ERROR, WARNING, INFO = "error", "warning", "info"

D3_ELEMENTS = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn"}
LANTHANIDES = {"La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"}
ACTINIDES = {"Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr"}
F_ELEMENTS = LANTHANIDES | ACTINIDES
TETRAHEDRON = {-4, -5, -14, -15}
RELAX_IBRION = {1, 2, 3}
PHONON_IBRION = {5, 6, 7, 8}
VOLUME_ISIF = {3, 6, 7, 8}


@dataclass
class Rule:
    id: str
    severity: str
    title: str
    url: str
    check: callable


@dataclass
class Finding:
    rule: str
    severity: str
    message: str
    url: str
    file: str = "INCAR"
    line: int = None


RULES = []


def hit(message, tag=None, severity=None, file=None):
    """One finding from a check: the INCAR tag it points at, and optional overrides."""
    return message, tag, severity, file


def rule(id, severity, title, page):
    def register(check):
        RULES.append(Rule(id, severity, title, WIKI + page, check))
        return check

    return register


# ---------------------------------------------------------------- helpers


def nsw(calc):
    return calc.incar.int("NSW", 0)


def ibrion(calc):
    default = -1 if nsw(calc) in (-1, 0) else 0
    return calc.incar.int("IBRION", default)


def ismear(calc):
    return calc.incar.int("ISMEAR", 1)


def isif(calc):
    default = 0 if ibrion(calc) == 0 or calc.incar.bool("LHFCALC") else 2
    return calc.incar.int("ISIF", default)


def noncollinear(calc):
    return calc.incar.bool("LNONCOLLINEAR") or calc.incar.bool("LSORBIT")


def moves_ions(calc):
    return nsw(calc) > 0 and ibrion(calc) in RELAX_IBRION | {0}


def hybrid(calc):
    return calc.incar.bool("LHFCALC")


def metagga(calc):
    return calc.incar.has("METAGGA") and calc.incar.str("METAGGA").upper() not in ("", "NONE", "--")


# ---------------------------------------------------------------- smearing


@rule("VL001", WARNING, "ISMEAR not set", "Smearing_technique")
def ismear_unset(calc):
    if not calc.incar.has("ISMEAR"):
        yield hit(
            "ISMEAR is not set, so VASP uses Methfessel-Paxton (ISMEAR = 1), which can give wrong "
            "results for semiconductors and insulators. Set it explicitly; ISMEAR = 0 with "
            "SIGMA = 0.03-0.1 is the safe choice when you do not know whether the system has a gap."
        )


@rule("VL002", WARNING, "Methfessel-Paxton smearing", "Smearing_technique")
def methfessel_paxton(calc):
    if calc.incar.has("ISMEAR") and ismear(calc) > 0:
        extra = ""
        if ibrion(calc) in PHONON_IBRION:
            extra = " For phonons in a gapped system the error in the frequencies can exceed 20%."
        yield hit(
            f"ISMEAR = {ismear(calc)} (Methfessel-Paxton) is meant for metals. For a semiconductor "
            f"or insulator use ISMEAR = 0 (Gaussian) or the tetrahedron method.{extra}",
            "ISMEAR",
        )


@rule("VL003", WARNING, "Tetrahedron method while ions move", "Smearing_technique")
def tetrahedron_forces(calc):
    if ismear(calc) in (-4, -5) and moves_ions(calc):
        yield hit(
            f"ISMEAR = {ismear(calc)} with ionic steps (NSW > 0): the tetrahedron method is not "
            "variational, so forces and stress can be wrong by 5-10% in metals. It is fine only "
            "if the gap stays open; for metals use ISMEAR = 1 or 2.",
            "ISMEAR",
        )


@rule("VL004", WARNING, "Tetrahedron method with too few k-points", "Smearing_technique")
def tetrahedron_kpoints(calc):
    if ismear(calc) not in TETRAHEDRON:
        return
    if calc.kpoints is not None and calc.kpoints.mode == "explicit":
        n = calc.kpoints.nkpts
    else:
        mesh = calc.kmesh()
        if mesh is None:
            return
        n = math.prod(mesh)
    if n < 4:
        yield hit(
            f"ISMEAR = {ismear(calc)} needs at least 4 k-points to form a tetrahedron; this setup has {n}.",
            "ISMEAR",
        )


# ---------------------------------------------------------------- basis set


@rule("VL005", WARNING, "ENCUT not set", "ENCUT")
def encut_unset(calc):
    if not calc.incar.has("ENCUT"):
        yield hit(
            "ENCUT is not set, so the default (largest ENMAX in the POTCAR) depends on the "
            "species. Set it by hand so that energies of different calculations can be compared."
        )


@rule("VL006", ERROR, "ENCUT below the POTCAR cutoffs", "POTCAR")
def encut_low(calc):
    encut = calc.incar.float("ENCUT")
    if encut is None or not calc.potcar:
        return
    enmax = max(p.enmax for p in calc.potcar)
    enmin = max(p.enmin for p in calc.potcar)
    if encut < enmin:
        yield hit(
            f"ENCUT = {encut:g} eV is below the largest ENMIN ({enmin:g} eV), the minimum viable cutoff.", "ENCUT"
        )
    elif encut < enmax:
        yield hit(
            f"ENCUT = {encut:g} eV is below the largest ENMAX ({enmax:g} eV), the recommended cutoff.",
            "ENCUT",
            WARNING,
        )


@rule("VL007", WARNING, "Cell relaxation with a low ENCUT (Pulay stress)", "ISIF")
def pulay_stress(calc):
    if not (isif(calc) in VOLUME_ISIF and nsw(calc) > 0 and ibrion(calc) in RELAX_IBRION):
        return
    encut = calc.incar.float("ENCUT")
    if encut is None:
        yield hit(
            f"ISIF = {isif(calc)} changes the cell volume, but ENCUT is not set, so it defaults to "
            "ENMAX. Volume relaxations need about 1.3 x ENMAX to limit the Pulay stress.",
            "ISIF",
        )
    elif calc.potcar:
        enmax = max(p.enmax for p in calc.potcar)
        if encut < 1.3 * enmax - 1e-6:
            yield hit(
                f"ISIF = {isif(calc)} changes the cell volume with ENCUT = {encut:g} eV < "
                f"1.3 x ENMAX = {1.3 * enmax:.0f} eV. The basis is not adjusted when the cell "
                "changes; raise ENCUT to limit the Pulay stress.",
                "ENCUT",
            )


@rule("VL008", WARNING, "LREAL = .TRUE. or On", "LREAL")
def lreal_true(calc):
    value = calc.incar.str("LREAL")
    if value is None:
        return
    word = value.split()[0].strip(".").upper()
    if word in ("T", "TRUE", "O", "ON"):
        yield hit(f"LREAL = {value} is not recommended. Use LREAL = Auto (or .FALSE.).", "LREAL")


@rule("VL009", INFO, "Real-space projection in a small cell", "LREAL")
def lreal_small(calc):
    value = calc.incar.str("LREAL")
    if value is None or not calc.poscar:
        return
    word = value.split()[0].strip(".").upper()
    if word in ("A", "AUTO") and calc.poscar.nions < 30:
        yield hit(
            f"LREAL = Auto with {calc.poscar.nions} atoms. Real-space projection is recommended "
            "above about 30 atoms and always adds a small error; for a small cell LREAL = .FALSE. "
            "costs little.",
            "LREAL",
        )


# ---------------------------------------------------------------- electrons and bands


@rule("VL010", ERROR, "Too few bands", "NBANDS")
def nbands_low(calc):
    nbands = calc.incar.int("NBANDS")
    if nbands is None:
        return
    nelect = calc.incar.float("NELECT")
    if nelect is None:
        if not (calc.potcar and calc.poscar) or len(calc.potcar) != len(calc.poscar.counts):
            return
        nelect = sum(p.zval * n for p, n in zip(calc.potcar, calc.poscar.counts))
    occupied = math.ceil(nelect) if noncollinear(calc) else math.ceil(nelect / 2)
    if nbands < occupied + 1:
        yield hit(
            f"NBANDS = {nbands}, but {nelect:g} electrons fill {occupied} bands. VASP needs all "
            f"occupied bands plus at least one empty band (NBANDS >= {occupied + 1}).",
            "NBANDS",
        )


@rule("VL011", ERROR, "Hybrid functional with ALGO = Fast", "LHFCALC")
def hybrid_algo_fast(calc):
    if hybrid(calc) and (calc.incar.str("ALGO", "") or "").upper().startswith("F"):
        yield hit(
            "ALGO = Fast is not properly supported for hybrid functionals, and VASP prints no "
            "warning. Use ALGO = All, Damped or Normal.",
            "ALGO",
        )


@rule("VL012", ERROR, "ICHARG = 11 with a hybrid functional", "Band-structure_calculation_using_hybrid_functionals")
def hybrid_icharg11(calc):
    if hybrid(calc) and calc.incar.int("ICHARG") == 11:
        yield hit(
            "ICHARG = 11 (fixed density) must never be used with a hybrid functional: exact "
            "exchange needs the orbitals on a regular k-mesh. Use KPOINTS_OPT or zero-weight k-points.",
            "ICHARG",
        )


@rule("VL013", WARNING, "ICHARG = 11 without CHGCAR", "ICHARG")
def icharg11_chgcar(calc):
    if calc.incar.int("ICHARG") == 11 and "CHGCAR" not in calc.files:
        yield hit(
            "ICHARG = 11 reads the charge density from CHGCAR, which is not in this directory. "
            "Make sure the job copies it from the self-consistent run.",
            "ICHARG",
        )


# ---------------------------------------------------------------- k-points


@rule("VL014", ERROR, "Line-mode KPOINTS in a self-consistent run", "KPOINTS")
def line_mode_scf(calc):
    if calc.kpoints is None or calc.kpoints.mode != "line":
        return
    if hybrid(calc) or metagga(calc):
        yield hit(
            "KPOINTS is in line mode, but hybrid and meta-GGA functionals always need a regular "
            "mesh. Put the mesh in KPOINTS and the path in KPOINTS_OPT.",
            None,
            ERROR,
            "KPOINTS",
        )
    elif calc.incar.int("ICHARG") != 11:
        yield hit(
            "KPOINTS is in line mode, which is not suitable for a self-consistent run. Set "
            "ICHARG = 11 and read the CHGCAR of a previous run.",
            None,
            ERROR,
            "KPOINTS",
        )


@rule("VL015", ERROR, "KPOINTS_OPT without a regular mesh", "KPOINTS_OPT")
def kpoints_opt_mesh(calc):
    if calc.has_kpoints_opt and calc.kpoints is not None and calc.kpoints.mode == "line":
        yield hit(
            "KPOINTS_OPT is present, so KPOINTS must hold a uniform mesh, not a line-mode path.",
            None,
            ERROR,
            "KPOINTS",
        )


@rule(
    "VL016", WARNING, "Hybrid band structure without HFRCUT = -1", "Band-structure_calculation_using_hybrid_functionals"
)
def hybrid_hfrcut(calc):
    if calc.has_kpoints_opt and hybrid(calc) and calc.incar.int("HFRCUT", 0) != -1:
        yield hit(
            "Hybrid band structure on KPOINTS_OPT: the default treatment of the Coulomb "
            "singularity gives discontinuous bands. HFRCUT = -1 is recommended, and converges "
            "best for systems with a gap.",
            "HFRCUT" if calc.incar.has("HFRCUT") else "LHFCALC",
        )


@rule("VL017", WARNING, "KSPACING ignored", "KSPACING")
def kspacing_ignored(calc):
    if calc.incar.has("KSPACING") and "KPOINTS" in calc.files:
        yield hit(
            "KSPACING is only used when there is no KPOINTS file. This directory has one, so KSPACING is ignored.",
            "KSPACING",
        )


@rule("VL018", INFO, "Fully automatic KPOINTS", "KPOINTS")
def kpoints_auto(calc):
    if calc.kpoints is not None and calc.kpoints.mode == "auto":
        yield hit(
            "The fully automatic KPOINTS mode (length R_k) is deprecated. Prefer KSPACING, or give "
            "the mesh explicitly for production runs.",
            None,
            INFO,
            "KPOINTS",
        )


# ---------------------------------------------------------------- structure files


@rule("VL019", ERROR, "POTCAR does not match POSCAR", "POTCAR")
def potcar_order(calc):
    if not (calc.potcar and calc.poscar):
        return
    potcar = [p.element for p in calc.potcar]
    if len(potcar) != len(calc.poscar.counts):
        yield hit(
            f"POTCAR has {len(potcar)} species ({' '.join(potcar)}) but POSCAR has "
            f"{len(calc.poscar.counts)} species counts.",
            None,
            ERROR,
            "POTCAR",
        )
    elif calc.poscar.species and potcar != calc.poscar.species:
        yield hit(
            f"POTCAR order ({' '.join(potcar)}) differs from POSCAR ({' '.join(calc.poscar.species)}). "
            "The potentials must follow the order of the species in POSCAR.",
            None,
            ERROR,
            "POTCAR",
        )


@rule("VL020", ERROR, "Wrong number of MAGMOM values", "MAGMOM")
def magmom_count(calc):
    if not calc.incar.has("MAGMOM") or not calc.poscar:
        return
    try:
        values = calc.incar.list("MAGMOM")
    except ValueError:
        return
    nions = calc.poscar.nions
    if noncollinear(calc):
        expected = 3 * nions
    elif calc.incar.int("ISPIN", 1) == 2:
        expected = nions
    else:
        return
    if len(values) != expected:
        yield hit(f"MAGMOM has {len(values)} values; this cell with {nions} atoms needs {expected}.", "MAGMOM")


@rule("VL021", ERROR, "DFT+U arrays do not match the species", "LDAUL")
def ldau_lengths(calc):
    if not calc.incar.bool("LDAU") or not calc.poscar:
        return
    nspecies = len(calc.poscar.counts)
    for tag in ("LDAUL", "LDAUU", "LDAUJ"):
        try:
            values = calc.incar.list(tag)
        except ValueError:
            continue
        if values is not None and len(values) != nspecies:
            yield hit(f"{tag} has {len(values)} values but POSCAR has {nspecies} species; give one per species.", tag)


@rule("VL022", WARNING, "DFT+U with LMAXMIX too low", "LMAXMIX")
def ldau_lmaxmix(calc):
    if not calc.incar.bool("LDAU"):
        return
    try:
        ls = calc.incar.list("LDAUL") or []
    except ValueError:
        return
    needed = 6 if 3 in ls else 4 if 2 in ls else 2
    lmaxmix = calc.incar.int("LMAXMIX", 2)
    if lmaxmix >= needed:
        return
    kind = "f" if needed == 6 else "d"
    if calc.incar.int("ICHARG") == 11:
        yield hit(
            f"DFT+U band structure on {kind} states with LMAXMIX = {lmaxmix}: LMAXMIX = {needed} is "
            "strictly required, both here and in the run that wrote CHGCAR.",
            "LMAXMIX" if calc.incar.has("LMAXMIX") else "LDAU",
            ERROR,
        )
    else:
        yield hit(
            f"DFT+U on {kind} states with LMAXMIX = {lmaxmix}: LMAXMIX = {needed} often converges "
            "much faster, and is required if this CHGCAR is used for ICHARG = 11 later.",
            "LMAXMIX" if calc.incar.has("LMAXMIX") else "LDAU",
        )


@rule("VL023", WARNING, "LASPH off for d/f elements", "LASPH")
def lasph(calc):
    if calc.incar.bool("LASPH") or not calc.species:
        return
    method = [
        name
        for name, on in (
            ("DFT+U", calc.incar.bool("LDAU")),
            ("hybrid", hybrid(calc)),
            ("meta-GGA", metagga(calc)),
            ("vdW-DF", calc.incar.bool("LUSE_VDW")),
        )
        if on
    ]
    heavy = [s for s in calc.species if s in D3_ELEMENTS | F_ELEMENTS]
    if method and heavy:
        yield hit(
            f"{', '.join(method)} with {', '.join(heavy)} and LASPH = .FALSE.: the aspherical "
            "gradient corrections inside the PAW spheres are essential here. Set LASPH = .TRUE.",
            "LASPH" if calc.incar.has("LASPH") else None,
        )


# ---------------------------------------------------------------- ionic motion and phonons


@rule("VL024", ERROR, "Molecular dynamics without POTIM", "POTIM")
def md_potim(calc):
    if ibrion(calc) == 0 and nsw(calc) > 0 and not calc.incar.has("POTIM"):
        yield hit("IBRION = 0 (molecular dynamics) needs POTIM, the time step in fs; VASP stops without it.", "IBRION")


@rule("VL025", WARNING, "Molecular dynamics with symmetry", "ISYM")
def md_isym(calc):
    if ibrion(calc) == 0 and nsw(calc) > 0 and calc.incar.int("ISYM", 2) != 0:
        yield hit(
            "ISYM = 0 should be set for molecular dynamics (IBRION = 0).",
            "ISYM" if calc.incar.has("ISYM") else "IBRION",
        )


@rule("VL026", INFO, "Phonons without symmetry", "IBRION")
def phonon_symmetry(calc):
    if ibrion(calc) in (5, 7):
        better = ibrion(calc) + 1
        yield hit(
            f"IBRION = {ibrion(calc)} is deprecated. IBRION = {better} uses symmetry and can be "
            "orders of magnitude cheaper with the same accuracy.",
            "IBRION",
        )


@rule("VL027", WARNING, "Phonons with a loose EDIFF", "IBRION")
def phonon_ediff(calc):
    if ibrion(calc) in PHONON_IBRION and calc.incar.float("EDIFF", 1e-4) > 1e-6:
        yield hit(
            f"Phonons (IBRION = {ibrion(calc)}) with EDIFF = {calc.incar.float('EDIFF', 1e-4):g}: "
            "the default 1E-4 often gives unacceptably large errors; use EDIFF <= 1E-6.",
            "EDIFF" if calc.incar.has("EDIFF") else "IBRION",
        )


@rule("VL028", ERROR, "Dispersion correction with DFPT phonons", "IVDW")
def dfpt_ivdw(calc):
    if ibrion(calc) in (7, 8) and calc.incar.int("IVDW", 0) != 0:
        yield hit(
            f"IVDW = {calc.incar.int('IVDW')} is not available for DFPT phonons (IBRION = 7 or 8). "
            "Use finite differences (IBRION = 6) instead.",
            "IVDW",
        )


# ---------------------------------------------------------------- parallelisation


@rule("VL029", WARNING, "NCORE and NPAR both set", "NCORE")
def ncore_npar(calc):
    if calc.incar.has("NCORE") and calc.incar.has("NPAR"):
        yield hit("Do not set NCORE and NPAR together; NPAR takes precedence. Keep only NCORE.", "NPAR")


# ---------------------------------------------------------------- nudged elastic band


@rule("VL030", ERROR, "Incomplete NEB image directories", "Nudged_elastic_bands")
def neb_images(calc):
    images = calc.incar.int("IMAGES", 0)
    if images <= 0:
        return
    structures, missing = [], []
    for i in range(images + 2):
        folder = calc.path / f"{i:02d}"
        if (folder / "POSCAR").is_file():
            structures.append((folder.name, read_poscar(folder / "POSCAR")))
        else:
            missing.append(folder.name)
    if missing:
        yield hit(
            f"IMAGES = {images} needs {images + 2} folders 00-{images + 1:02d}, each with a POSCAR "
            f"(end points in 00 and {images + 1:02d}). Missing: {', '.join(missing)}.",
            "IMAGES",
        )
    reference = None
    for name, poscar in structures:
        if poscar is None:
            continue
        key = (poscar.species, poscar.counts)
        if reference is None:
            reference = (name, key)
        elif key != reference[1]:
            yield hit(f"{name}/POSCAR has a different element order or count than {reference[0]}/POSCAR.", "IMAGES")


def run(calc, ignore=()):
    """Apply every rule to one calculation and return the findings."""
    findings = []
    for r in RULES:
        if r.id in ignore:
            continue
        for message, tag, severity, file in r.check(calc) or ():
            findings.append(
                Finding(
                    rule=r.id,
                    severity=severity or r.severity,
                    message=message,
                    url=r.url,
                    file=file or "INCAR",
                    line=calc.incar.lines.get(tag) if (file or "INCAR") == "INCAR" else None,
                )
            )
    return findings
