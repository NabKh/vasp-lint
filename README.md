# vasp-lint

[![tests](https://github.com/NabKh/vasp-lint/actions/workflows/tests.yml/badge.svg)](https://github.com/NabKh/vasp-lint/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)

Check VASP inputs (INCAR, POSCAR, POTCAR, KPOINTS) for common mistakes **before** you submit the job.
Every check cites the VASP wiki page it comes from, so you can see why it matters.

It catches the mistakes that waste a day of queue time or, worse, run fine and give wrong numbers:
ICHARG = 11 with a hybrid functional, ALGO = Fast with HSE, MAGMOM with the wrong length,
a POTCAR in a different order from the POSCAR, ENCUT too low for a cell relaxation,
Methfessel-Paxton smearing on a semiconductor, LASPH off for DFT+U on 3d metals, and more.

- No dependencies: plain Python 3.9+, works on login nodes.
- Reads the files only; nothing is uploaded anywhere. POTCARs stay on your machine.
- Works on one folder, a whole project tree, in CI or as a pre-commit hook.

## Install

```bash
pip install git+https://github.com/NabKh/vasp-lint
```

## Use

```bash
vasp-lint                 # every folder with an INCAR below the current one
vasp-lint relax/ bands/   # specific calculations
vasp-lint --list-rules    # all checks with their wiki links
```

Example, for an HSE band-structure INCAR with a few classic mistakes:

```text
INCAR:3: VL011 error: ALGO = Fast is not properly supported for hybrid functionals, and VASP prints no warning. Use ALGO = All, Damped or Normal.
    https://www.vasp.at/wiki/index.php/LHFCALC
INCAR:4: VL012 error: ICHARG = 11 (fixed density) must never be used with a hybrid functional: exact exchange needs the orbitals on a regular k-mesh. Use KPOINTS_OPT or zero-weight k-points.
    https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals
INCAR:6: VL020 error: MAGMOM has 2 values; this cell with 6 atoms needs 6.
    https://www.vasp.at/wiki/index.php/MAGMOM
INCAR: VL001 warning: ISMEAR is not set, so VASP uses Methfessel-Paxton (ISMEAR = 1), which can give wrong results for semiconductors and insulators. Set it explicitly; ISMEAR = 0 with SIGMA = 0.03-0.1 is the safe choice when you do not know whether the system has a gap.
    https://www.vasp.at/wiki/index.php/Smearing_technique
INCAR:7: VL008 warning: LREAL = .TRUE. is not recommended. Use LREAL = Auto (or .FALSE.).
    https://www.vasp.at/wiki/index.php/LREAL
INCAR:4: VL013 warning: ICHARG = 11 reads the charge density from CHGCAR, which is not in this directory. Make sure the job copies it from the self-consistent run.
    https://www.vasp.at/wiki/index.php/ICHARG
.: note: no readable POTCAR; checks that need them were skipped

1 calculation(s): 3 error(s), 3 warning(s), 0 info
```

Options:

| Option | Effect |
|---|---|
| `--ignore VL001,VL009` | skip rules |
| `--no-info` | hide info-level findings |
| `--strict` | exit with 1 on warnings as well as errors |
| `--format json` | machine-readable output |

Exit code: 0 if clean, 1 if there are errors (or warnings with `--strict`), 2 if no INCAR was found.

Checks that need the POSCAR, POTCAR or KPOINTS are skipped when that file is missing, and vasp-lint says so.

### From Python

```python
from vasp_lint import lint

for f in lint("relax"):
    print(f.rule, f.severity, f.line, f.message)
```

### pre-commit

```yaml
repos:
  - repo: https://github.com/NabKh/vasp-lint
    rev: main
    hooks:
      - id: vasp-lint
```

## Rules

| ID | Severity | Check | VASP wiki |
|---|---|---|---|
| VL001 | warning | ISMEAR not set | [Smearing technique](https://www.vasp.at/wiki/index.php/Smearing_technique) |
| VL002 | warning | Methfessel-Paxton smearing | [Smearing technique](https://www.vasp.at/wiki/index.php/Smearing_technique) |
| VL003 | warning | Tetrahedron method while ions move | [Smearing technique](https://www.vasp.at/wiki/index.php/Smearing_technique) |
| VL004 | warning | Tetrahedron method with too few k-points | [Smearing technique](https://www.vasp.at/wiki/index.php/Smearing_technique) |
| VL005 | warning | ENCUT not set | [ENCUT](https://www.vasp.at/wiki/index.php/ENCUT) |
| VL006 | error | ENCUT below the POTCAR cutoffs | [POTCAR](https://www.vasp.at/wiki/index.php/POTCAR) |
| VL007 | warning | Cell relaxation with a low ENCUT (Pulay stress) | [ISIF](https://www.vasp.at/wiki/index.php/ISIF) |
| VL008 | warning | LREAL = .TRUE. or On | [LREAL](https://www.vasp.at/wiki/index.php/LREAL) |
| VL009 | info | Real-space projection in a small cell | [LREAL](https://www.vasp.at/wiki/index.php/LREAL) |
| VL010 | error | Too few bands | [NBANDS](https://www.vasp.at/wiki/index.php/NBANDS) |
| VL011 | error | Hybrid functional with ALGO = Fast | [LHFCALC](https://www.vasp.at/wiki/index.php/LHFCALC) |
| VL012 | error | ICHARG = 11 with a hybrid functional | [Band-structure calculation using hybrid functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals) |
| VL013 | warning | ICHARG = 11 without CHGCAR | [ICHARG](https://www.vasp.at/wiki/index.php/ICHARG) |
| VL014 | error | Line-mode KPOINTS in a self-consistent run | [KPOINTS](https://www.vasp.at/wiki/index.php/KPOINTS) |
| VL015 | error | KPOINTS_OPT without a regular mesh | [KPOINTS OPT](https://www.vasp.at/wiki/index.php/KPOINTS_OPT) |
| VL016 | warning | Hybrid band structure without HFRCUT = -1 | [Band-structure calculation using hybrid functionals](https://www.vasp.at/wiki/index.php/Band-structure_calculation_using_hybrid_functionals) |
| VL017 | warning | KSPACING ignored | [KSPACING](https://www.vasp.at/wiki/index.php/KSPACING) |
| VL018 | info | Fully automatic KPOINTS | [KPOINTS](https://www.vasp.at/wiki/index.php/KPOINTS) |
| VL019 | error | POTCAR does not match POSCAR | [POTCAR](https://www.vasp.at/wiki/index.php/POTCAR) |
| VL020 | error | Wrong number of MAGMOM values | [MAGMOM](https://www.vasp.at/wiki/index.php/MAGMOM) |
| VL021 | error | DFT+U arrays do not match the species | [LDAUL](https://www.vasp.at/wiki/index.php/LDAUL) |
| VL022 | warning | DFT+U with LMAXMIX too low | [LMAXMIX](https://www.vasp.at/wiki/index.php/LMAXMIX) |
| VL023 | warning | LASPH off for d/f elements | [LASPH](https://www.vasp.at/wiki/index.php/LASPH) |
| VL024 | error | Molecular dynamics without POTIM | [POTIM](https://www.vasp.at/wiki/index.php/POTIM) |
| VL025 | warning | Molecular dynamics with symmetry | [ISYM](https://www.vasp.at/wiki/index.php/ISYM) |
| VL026 | info | Phonons without symmetry | [IBRION](https://www.vasp.at/wiki/index.php/IBRION) |
| VL027 | warning | Phonons with a loose EDIFF | [IBRION](https://www.vasp.at/wiki/index.php/IBRION) |
| VL028 | error | Dispersion correction with DFPT phonons | [IVDW](https://www.vasp.at/wiki/index.php/IVDW) |
| VL029 | warning | NCORE and NPAR both set | [NCORE](https://www.vasp.at/wiki/index.php/NCORE) |
| VL030 | error | Incomplete NEB image directories | [Nudged elastic bands](https://www.vasp.at/wiki/index.php/Nudged_elastic_bands) |

The checks follow the recommendations on the VASP wiki. They are guidance, not proof that a calculation is converged:
always converge ENCUT, k-points and smearing for your own system.

## Contributing

Found a mistake that bit you and is not caught? Open an issue with the INCAR (no POTCAR, please; it is licensed)
and the wiki page that describes the problem. Pull requests with a rule and a test are welcome:

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

## License

MIT
