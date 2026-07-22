"""
Convert specex_linelist_desi.txt (NIST-sourced, air wavelengths above 2000A
per NIST's default display convention -- see porting-notes.md "wavelength
offset resolved") to vacuum wavelengths, to test the hypothesis that the
shared ~0.5-0.9A wavelength-residual-vs-truth offset seen throughout this
project is caused by comparing a vacuum-calibrated trace against an air
line list.

Uses the standard IAU/Morton (2000, ApJS 130, 403) air->vacuum conversion
(the same formula used by SDSS, astropy-adjacent pipelines, etc.):
  sigma2 = (1e4 / wave_air_A)^2      [wave in Angstroms, sigma in 1/micron]
  n = 1 + 5.792105e-8/(238.0185e-8 - ... )  -- see air_to_vac() below for the
  exact coefficients as normally quoted (dimensionless correction ~2.7-2.9e-4
  across the optical/NIR).
  wave_vacuum = wave_air * n

Only the wavelength column (whitespace-split field index 1) is touched;
everything else (species name, score, intensity, trailing comments,
'#'-comment lines, blank lines) is preserved byte-for-byte.

Usage:
  python testing/air_to_vacuum_linelist.py
    (reads py/specex/data/specex_linelist_desi.txt, writes
     py/specex/data/specex_linelist_desi_vacuum.txt)
"""
import os


def air_to_vac(wave_air_angstrom):
    """Morton (2000) IAU standard air->vacuum conversion, valid >2000A."""
    s2 = (1e4 / wave_air_angstrom) ** 2
    n = 1 + 5.792105e-2 / (238.0185 - s2) + 1.67917e-3 / (57.362 - s2)
    return wave_air_angstrom * n


def main():
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    in_path = os.path.join(base, "py/specex/data/specex_linelist_desi.txt")
    out_path = os.path.join(base, "py/specex/data/specex_linelist_desi_vacuum.txt")

    n_converted = 0
    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                fout.write(line)
                continue
            parts = stripped.split()
            if len(parts) < 3:
                fout.write(line)
                continue
            try:
                wave_air = float(parts[1])
            except ValueError:
                fout.write(line)
                continue
            wave_vac = air_to_vac(wave_air)
            parts[1] = f"{wave_vac:.4f}"
            fout.write(" ".join(parts) + "\n")
            n_converted += 1

    print(f"Converted {n_converted} lines: {in_path} -> {out_path}")


if __name__ == "__main__":
    main()
