import json

import pytest

from vasp_lint import RULES, lint
from vasp_lint.cli import main

POSCAR_MOS2 = """MoS2
1.0
 3.16 0.00 0.00
-1.58 2.7367 0.00
 0.00 0.00 12.30
Mo S
2 4
Direct
0.333 0.667 0.25
0.667 0.333 0.75
0.333 0.667 0.62
0.333 0.667 0.88
0.667 0.333 0.12
0.667 0.333 0.38
"""

POSCAR_FEO = (
    """FeO
1.0
4.3 0 0
0 4.3 0
0 0 4.3
Fe O
4 4
Direct
"""
    + "0 0 0\n" * 8
)


def potcar(*entries):
    """Minimal POTCAR text: (symbol, zval, enmax, enmin) per species."""
    blocks = []
    for symbol, zval, enmax, enmin in entries:
        blocks.append(
            f"  PAW_PBE {symbol} 01Jan2000\n"
            f"   TITEL  = PAW_PBE {symbol} 01Jan2000\n"
            f"   POMASS =   1.0; ZVAL   =   {zval}    mass and valenz\n"
            f"   ENMAX  =  {enmax}; ENMIN  =  {enmin} eV\n"
            " End of Dataset\n"
        )
    return "".join(blocks)


MOS2_POTCAR = potcar(("Mo_sv", 14.0, 241.6, 181.2), ("S", 6.0, 258.7, 194.0))
GOOD_INCAR = "ENCUT = 450\nISMEAR = 0\nSIGMA = 0.05\n"
MESH = "auto\n0\nGamma\n9 9 2\n"


def make(tmp_path, incar, poscar=POSCAR_MOS2, pot=MOS2_POTCAR, kpoints=MESH, extra=None):
    (tmp_path / "INCAR").write_text(incar)
    if poscar:
        (tmp_path / "POSCAR").write_text(poscar)
    if pot:
        (tmp_path / "POTCAR").write_text(pot)
    if kpoints:
        (tmp_path / "KPOINTS").write_text(kpoints)
    for name, text in (extra or {}).items():
        (tmp_path / name).write_text(text)
    return tmp_path


def ids(directory):
    return {f.rule for f in lint(directory)}


def test_rule_ids_unique_and_linked():
    assert len({r.id for r in RULES}) == len(RULES)
    assert all(r.url.startswith("https://www.vasp.at/wiki/") for r in RULES)


def test_clean_input_has_no_findings(tmp_path):
    assert ids(make(tmp_path, GOOD_INCAR)) == set()


@pytest.mark.parametrize(
    "incar, rule",
    [
        ("ENCUT = 450\n", "VL001"),
        ("ENCUT = 450\nISMEAR = 1\n", "VL002"),
        ("ENCUT = 450\nISMEAR = -5\nNSW = 50\nIBRION = 2\n", "VL003"),
        ("ISMEAR = 0\n", "VL005"),
        ("ENCUT = 150\nISMEAR = 0\n", "VL006"),
        ("ENCUT = 300\nISMEAR = 0\nISIF = 3\nNSW = 50\nIBRION = 2\n", "VL007"),
        ("ISMEAR = 0\nISIF = 3\nNSW = 50\nIBRION = 2\n", "VL007"),
        (GOOD_INCAR + "LREAL = .TRUE.\n", "VL008"),
        (GOOD_INCAR + "LREAL = On\n", "VL008"),
        (GOOD_INCAR + "LREAL = Auto\n", "VL009"),
        (GOOD_INCAR + "NBANDS = 26\n", "VL010"),
        (GOOD_INCAR + "LHFCALC = .TRUE.\nALGO = Fast\n", "VL011"),
        (GOOD_INCAR + "LHFCALC = .TRUE.\nICHARG = 11\n", "VL012"),
        (GOOD_INCAR + "ICHARG = 11\n", "VL013"),
        (GOOD_INCAR + "KSPACING = 0.3\n", "VL017"),
        (GOOD_INCAR + "ISPIN = 2\nMAGMOM = 2*1.0 3*0.0\n", "VL020"),
        (GOOD_INCAR + "LSORBIT = .TRUE.\nMAGMOM = 6*0.0\n", "VL020"),
        (GOOD_INCAR + "LDAU = .TRUE.\nLDAUL = 2\nLDAUU = 3\nLDAUJ = 0\nLMAXMIX = 4\n", "VL021"),
        (GOOD_INCAR + "IBRION = 0\nNSW = 100\nISYM = 0\n", "VL024"),
        (GOOD_INCAR + "IBRION = 0\nNSW = 100\nPOTIM = 1\n", "VL025"),
        (GOOD_INCAR + "IBRION = 5\nEDIFF = 1E-7\n", "VL026"),
        (GOOD_INCAR + "IBRION = 6\n", "VL027"),
        (GOOD_INCAR + "IBRION = 8\nEDIFF = 1E-7\nIVDW = 12\n", "VL028"),
        (GOOD_INCAR + "NCORE = 4\nNPAR = 2\n", "VL029"),
    ],
)
def test_rule_fires(tmp_path, incar, rule):
    assert rule in ids(make(tmp_path, incar))


def test_magmom_count_correct_is_quiet(tmp_path):
    incar = GOOD_INCAR + "ISPIN = 2\nMAGMOM = 2*1.0 4*0.0\n"
    assert "VL020" not in ids(make(tmp_path, incar))


def test_nbands_enough_is_quiet(tmp_path):
    assert "VL010" not in ids(make(tmp_path, GOOD_INCAR + "NBANDS = 27\n"))


def test_nbands_noncollinear_counts_one_electron_per_band(tmp_path):
    incar = GOOD_INCAR + "LNONCOLLINEAR = .TRUE.\nNBANDS = 40\n"
    assert "VL010" in ids(make(tmp_path, incar))


def test_encut_between_enmin_and_enmax_is_warning(tmp_path):
    findings = lint(make(tmp_path, "ENCUT = 220\nISMEAR = 0\n"))
    assert [(f.rule, f.severity) for f in findings] == [("VL006", "warning")]


def test_volume_relax_with_high_encut_is_quiet(tmp_path):
    incar = "ENCUT = 450\nISMEAR = 0\nISIF = 3\nNSW = 50\nIBRION = 2\n"
    assert "VL007" not in ids(make(tmp_path, incar))


def test_tetrahedron_gamma_only(tmp_path):
    incar = "ENCUT = 450\nISMEAR = -5\n"
    assert "VL004" in ids(make(tmp_path, incar, kpoints="g\n0\nGamma\n1 1 1\n"))
    assert "VL004" not in ids(make(tmp_path, incar, kpoints=MESH))


def test_tetrahedron_kspacing_mesh(tmp_path):
    # no KPOINTS: a large KSPACING gives a 1x1x1 mesh for this cell
    assert "VL004" in ids(make(tmp_path, "ENCUT = 450\nISMEAR = -5\nKSPACING = 5\n", kpoints=None))


def test_line_mode_kpoints(tmp_path):
    line = "path\n20\nLine-mode\nreciprocal\n0 0 0 G\n0.5 0 0 M\n"
    assert "VL014" in ids(make(tmp_path, GOOD_INCAR, kpoints=line))
    assert "VL014" not in ids(make(tmp_path, GOOD_INCAR + "ICHARG = 11\n", kpoints=line, extra={"CHGCAR": ""}))
    assert "VL014" in ids(make(tmp_path, GOOD_INCAR + "LHFCALC = .TRUE.\n", kpoints=line))
    assert "VL015" in ids(make(tmp_path, GOOD_INCAR + "ICHARG = 11\n", kpoints=line, extra={"KPOINTS_OPT": line}))


def test_hybrid_kpoints_opt_hfrcut(tmp_path):
    incar = GOOD_INCAR + "LHFCALC = .TRUE.\nHFSCREEN = 0.2\nALGO = All\n"
    line = "path\n20\nLine-mode\nreciprocal\n0 0 0 G\n0.5 0 0 M\n"
    assert "VL016" in ids(make(tmp_path, incar, extra={"KPOINTS_OPT": line}))
    assert "VL016" not in ids(make(tmp_path, incar + "HFRCUT = -1\n", extra={"KPOINTS_OPT": line}))


def test_fully_automatic_kpoints(tmp_path):
    assert "VL018" in ids(make(tmp_path, GOOD_INCAR, kpoints="auto\n0\nAuto\n20\n"))


def test_potcar_order(tmp_path):
    swapped = potcar(("S", 6.0, 258.7, 194.0), ("Mo_sv", 14.0, 241.6, 181.2))
    assert "VL019" in ids(make(tmp_path, GOOD_INCAR, pot=swapped))
    single = potcar(("Mo_sv", 14.0, 241.6, 181.2))
    assert "VL019" in ids(make(tmp_path, GOOD_INCAR, pot=single))


def test_ldau_lmaxmix_and_lasph(tmp_path):
    pot = potcar(("Fe", 8.0, 267.9, 200.9), ("O", 6.0, 400.0, 300.0))
    incar = "ENCUT = 520\nISMEAR = 0\nLDAU = .TRUE.\nLDAUL = 2 -1\nLDAUU = 5.3 0\nLDAUJ = 0 0\n"
    found = ids(make(tmp_path, incar, poscar=POSCAR_FEO, pot=pot))
    assert {"VL022", "VL023"} <= found
    fixed = incar + "LMAXMIX = 4\nLASPH = .TRUE.\n"
    assert not {"VL021", "VL022", "VL023"} & ids(make(tmp_path, fixed, poscar=POSCAR_FEO, pot=pot))


def test_ldau_band_structure_is_error(tmp_path):
    pot = potcar(("Fe", 8.0, 267.9, 200.9), ("O", 6.0, 400.0, 300.0))
    incar = "ENCUT = 520\nISMEAR = 0\nLDAU = .TRUE.\nLDAUL = 2 -1\nLDAUU = 5.3 0\nLDAUJ = 0 0\nICHARG = 11\n"
    findings = lint(make(tmp_path, incar, poscar=POSCAR_FEO, pot=pot, extra={"CHGCAR": ""}))
    assert any(f.rule == "VL022" and f.severity == "error" for f in findings)


def test_neb_folders(tmp_path):
    make(tmp_path, GOOD_INCAR + "IMAGES = 2\nIBRION = 3\nNSW = 50\n", poscar=None)
    for i in ("00", "01", "02"):
        (tmp_path / i).mkdir()
        (tmp_path / i / "POSCAR").write_text(POSCAR_MOS2)
    messages = [f.message for f in lint(tmp_path) if f.rule == "VL030"]
    assert len(messages) == 1 and "Missing: 03." in messages[0]
    (tmp_path / "03").mkdir()
    (tmp_path / "03" / "POSCAR").write_text(POSCAR_FEO)
    messages = [f.message for f in lint(tmp_path) if f.rule == "VL030"]
    assert len(messages) == 1 and "different element order" in messages[0]


def test_incar_comments_and_semicolons(tmp_path):
    incar = "ENCUT = 450 # cutoff\nISMEAR = 0 ; SIGMA = 0.05 ! smearing\nLREAL = .TRUE.\n"
    findings = lint(make(tmp_path, incar))
    assert [(f.rule, f.line) for f in findings] == [("VL008", 3)]


def test_cli_exit_codes_and_json(tmp_path, capsys):
    make(tmp_path, GOOD_INCAR)
    assert main([str(tmp_path)]) == 0
    make(tmp_path, GOOD_INCAR + "LHFCALC = .TRUE.\nICHARG = 11\n")
    assert main([str(tmp_path)]) == 1
    capsys.readouterr()
    assert main([str(tmp_path), "--format", "json", "--ignore", "VL012,VL013"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data[0]["findings"] == []


def test_cli_finds_nested_calculations(tmp_path, capsys):
    for name in ("relax", "bands"):
        (tmp_path / name).mkdir()
        make(tmp_path / name, GOOD_INCAR)
    assert main([str(tmp_path)]) == 0
    assert "2 calculation(s)" in capsys.readouterr().out
