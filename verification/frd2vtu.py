"""Convert CalculiX .frd results to VTK .vtu files for ParaView.

Reads the mesh and every nodal result block (displacement, stress, error
estimate, ...) from an ASCII .frd file and writes an unstructured-grid .vtu
next to it. Adds von Mises stress when a STRESS block is present. Only numpy is
needed.

Usage
    uv run frd2vtu.py runs/plate_hole/ccx/plate_CPS8_n4.frd
    uv run frd2vtu.py runs/plate_hole/ccx/*.frd

Then in ParaView: File > Open the .vtu, Apply, colour by STRESS (or von_Mises).
To see the deformed shape: Filters > Warp By Vector, vector DISP, scale factor
as needed. Results are written once per step; for several steps the arrays are
suffixed _step<n>.
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
    lines = Path(path).read_text(errors="replace").splitlines()
    nodes: dict[int, list[float]] = {}
    cells: list[tuple[int, list[int]]] = []
    results: list[tuple[str, list[str], dict[int, list[float]]]] = []
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
            results.append((name, comps, data))
        i += 1
    return nodes, cells, results


def von_mises(s: np.ndarray) -> np.ndarray:
    sxx, syy, szz, sxy, syz, szx = s.T
    return np.sqrt(0.5 * ((sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2) + 3 * (sxy**2 + syz**2 + szx**2))


def write_vtu(path: Path, nodes, cells, results) -> None:
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

    names = [r[0] for r in results]
    arrays = []
    for k, (name, comps, data) in enumerate(results):
        label = name if names.count(name) == 1 else f"{name}_step{names[:k + 1].count(name)}"
        ncomp = max(len(v) for v in data.values())
        a = np.zeros((len(ids), ncomp))
        for n, v in data.items():
            a[index[n], :len(v)] = v
        arrays.append((label, comps[:ncomp], a))
        if name == "STRESS" and ncomp == 6:
            arrays.append((label.replace("STRESS", "von_Mises"), [], von_mises(a)[:, None]))

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


def main() -> None:
    files = [Path(a) for a in sys.argv[1:]]
    if not files:
        raise SystemExit(__doc__)
    for frd in files:
        nodes, cells, results = read_frd(frd)
        out = frd.with_suffix(".vtu")
        write_vtu(out, nodes, cells, results)
        print(f"{frd} -> {out}  ({len(nodes)} nodes, {len(cells)} elements, "
              f"results: {', '.join(r[0] for r in results) or 'none'})")


if __name__ == "__main__":
    main()
