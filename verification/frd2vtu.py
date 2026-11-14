"""Convert CalculiX .frd results to VTK .vtu files for ParaView.

Reads the mesh and every nodal result block (displacement, stress, contact,
error estimate, ...) from an ASCII .frd file and writes, next to it:

  <job>.vtu        the mesh with the results of the last increment, under their
                   plain names (DISP, STRESS, CONTACT, ...), plus von_Mises when
                   a STRESS block is present
  <job>_inc<k>.vtu one file per increment, and
  <job>.pvd        a time series of them    (only when there is more than one
                   increment, e.g. contact or nonlinear steps)

Only numpy is needed. fea_common.run_ccx calls convert() after every successful
CalculiX run, so every script that runs ccx through it gets these files.

Usage (for .frd files from elsewhere, or after editing a deck by hand)
    uv run frd2vtu.py runs/plate_hole/ccx/plate_CPS8_n4.frd
    uv run frd2vtu.py runs/hertz/ccx/*.frd

In ParaView: File > Open the .vtu (or the .pvd to step through the increments),
Apply, colour by STRESS, von_Mises, CONTACT (CPRESS), ... For the deformed
shape: Filters > Warp By Vector, vector DISP, scale factor as needed.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# frd element type -> (VTK cell type, permutation from frd node order to VTK order)
_ID = lambda n: list(range(n))
FRD_TO_VTK = {
    1: (12, _ID(8)),                                            # he8   -> VTK_HEXAHEDRON
    4: (25, _ID(12) + [16, 17, 18, 19, 12, 13, 14, 15]),        # he20  -> VTK_QUADRATIC_HEXAHEDRON
    3: (10, _ID(4)),                                            # te4   -> VTK_TETRA
    7: (5, _ID(3)),                                             # tr3   -> VTK_TRIANGLE
    9: (9, _ID(4)),                                             # qu4   -> VTK_QUAD
    10: (23, _ID(8)),                                           # qu8   -> VTK_QUADRATIC_QUAD
}
# frd he20 lists the four vertical mid-edge nodes (Abaqus 17-20) before the
# top-face mid-edge nodes (Abaqus 13-16); VTK follows the Abaqus order.


def _fields(line: str, start: int, n: int, width: int = 12) -> list[float]:
    return [float(line[start + width * k:start + width * (k + 1)]) for k in range(n)]


def read_frd(path: Path):
    """Nodes, cells, and result blocks [(name, components, {node: values}, time)]."""
    lines = Path(path).read_text(errors="replace").splitlines()
    nodes: dict[int, list[float]] = {}
    cells: list[tuple[int, list[int]]] = []
    results: list[tuple[str, list[str], dict[int, list[float]], float]] = []
    time = 0.0
    i = 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("    2C"):                      # node block
            i += 1
            while not lines[i].startswith(" -3"):
                nodes[int(lines[i][3:13])] = _fields(lines[i], 13, 3)
                i += 1
        elif ln.startswith("    3C"):                    # element block
            i += 1
            while not lines[i].startswith(" -3"):
                etype = int(lines[i][13:18])
                i += 1
                conn = []
                while i < len(lines) and lines[i].startswith(" -2"):
                    body = lines[i][3:]
                    conn += [int(body[k:k + 10]) for k in range(0, len(body.rstrip()), 10)]
                    i += 1
                cells.append((etype, conn))
        elif ln.startswith("  100C"):                    # result-block header: step time
            time = float(ln[12:24])
        elif ln.startswith(" -4"):                        # result block
            name, ncomp = ln[5:13].strip(), int(ln[13:18])
            comps = []
            i += 1
            while lines[i].startswith(" -5"):
                comps.append(lines[i][5:13].strip())
                i += 1
            comps = [c for c in comps if c != "ALL"][:ncomp]
            data = {}
            while not lines[i].startswith(" -3"):
                if lines[i].startswith(" -1"):
                    body = lines[i][13:].rstrip()
                    data[int(lines[i][3:13])] = _fields(lines[i], 13, len(body) // 12)
                i += 1
            results.append((name, comps, data, time))
        i += 1
    return nodes, cells, results


def von_mises(s: np.ndarray) -> np.ndarray:
    sxx, syy, szz, sxy, syz, szx = s.T
    return np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2) + 3 * (sxy**2 + syz**2 + szx**2))


def write_vtu(path: Path, nodes, cells, results) -> None:
    """One .vtu: the mesh plus the given result blocks, each under its own name."""
    ids = sorted(nodes)
    index = {n: k for k, n in enumerate(ids)}
    pts = np.array([nodes[n] for n in ids])
    conn, offsets, types = [], [], []
    for etype, c in cells:
        if etype not in FRD_TO_VTK:
            raise ValueError(f"frd element type {etype} is not supported by this converter")
        vtk_type, perm = FRD_TO_VTK[etype]
        conn += [index[c[p]] for p in perm]
        offsets.append(len(conn))
        types.append(vtk_type)

    arrays = []
    for name, comps, data, *_ in results:
        ncomp = max(len(v) for v in data.values())
        a = np.zeros((len(ids), ncomp))
        for n, v in data.items():
            a[index[n], :len(v)] = v
        arrays.append((name, comps[:ncomp], a))
        if name == "STRESS" and ncomp == 6:
            arrays.append(("von_Mises", [], von_mises(a)[:, None]))

    def block(a: np.ndarray, fmt: str) -> str:
        return "\n".join(" ".join(fmt % v for v in row) for row in np.atleast_2d(a))

    out = ['<?xml version="1.0"?>',
           '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
           "<UnstructuredGrid>",
           f'<Piece NumberOfPoints="{len(ids)}" NumberOfCells="{len(cells)}">',
           "<PointData>"]
    for label, comps, a in arrays:
        names_attr = "".join(f' ComponentName{j}="{c}"' for j, c in enumerate(comps))
        out += [f'<DataArray type="Float64" Name="{label}" NumberOfComponents="{a.shape[1]}"{names_attr}'
                ' format="ascii">', block(a, "%.8g"), "</DataArray>"]
    out += ['<DataArray type="Int32" Name="node_id" format="ascii">', " ".join(map(str, ids)), "</DataArray>",
            "</PointData>", "<Points>",
            '<DataArray type="Float64" NumberOfComponents="3" format="ascii">', block(pts, "%.10g"),
            "</DataArray>", "</Points>", "<Cells>",
            '<DataArray type="Int32" Name="connectivity" format="ascii">', " ".join(map(str, conn)), "</DataArray>",
            '<DataArray type="Int32" Name="offsets" format="ascii">', " ".join(map(str, offsets)), "</DataArray>",
            '<DataArray type="UInt8" Name="types" format="ascii">', " ".join(map(str, types)), "</DataArray>",
            "</Cells>", "</Piece>", "</UnstructuredGrid>", "</VTKFile>", ""]
    path.write_text("\n".join(out))


def convert(frd: Path) -> list[Path]:
    """Write <job>.vtu (last increment) and, for several increments, <job>_inc<k>.vtu
    plus <job>.pvd. Returns the files written."""
    frd = Path(frd)
    nodes, cells, results = read_frd(frd)
    times = sorted({r[3] for r in results})
    last = [r for r in results if r[3] == times[-1]] if times else []
    main = frd.with_suffix(".vtu")
    write_vtu(main, nodes, cells, last)
    written = [main]
    for old in frd.parent.glob(f"{frd.stem}_inc*.vtu"):    # stale series from an earlier run
        old.unlink()
    if len(times) > 1:
        entries = []
        for k, t in enumerate(times, start=1):
            f = frd.with_name(f"{frd.stem}_inc{k}.vtu")
            write_vtu(f, nodes, cells, [r for r in results if r[3] == t])
            entries.append(f'<DataSet timestep="{t:.8g}" part="0" file="{f.name}"/>')
            written.append(f)
        pvd = frd.with_suffix(".pvd")
        pvd.write_text("\n".join(['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1">',
                                  "<Collection>", *entries, "</Collection>", "</VTKFile>", ""]))
        written.append(pvd)
    return written


def main() -> None:
    files = [Path(a) for a in sys.argv[1:]]
    if not files:
        raise SystemExit(__doc__)
    for frd in files:
        written = convert(frd)
        nodes, cells, results = read_frd(frd)
        names = list(dict.fromkeys(r[0] for r in results))
        print(f"{frd} -> {', '.join(p.name for p in written)}  ({len(nodes)} nodes, {len(cells)} elements, "
              f"results: {', '.join(names) or 'none'})")


if __name__ == "__main__":
    main()
