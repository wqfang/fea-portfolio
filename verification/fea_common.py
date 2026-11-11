"""Shared helpers for the FEA verification suite.

Structured hexahedral meshes, Abaqus-format input decks that both CalculiX (ccx)
and Abaqus accept, solver launching, and .dat parsing for both solvers.
Units are consistent: mm, N, MPa.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# Element families. "order" = polynomial order of the element's shape functions.
ELEMENT_ORDER = {
    "C3D8": 1,   # full integration: locks in bending
    "C3D8R": 1,  # reduced integration: hourglass-controlled
    "C3D8I": 1,  # incompatible modes: fixes bending locking
    "C3D20": 2,  # quadratic, full integration
    "C3D20R": 2, # quadratic, reduced integration
}


@dataclass
class HexMesh:
    """Structured hex mesh of the box [0,L]x[0,H]x[0,B]."""
    nodes: dict[int, tuple[float, float, float]]
    elements: list[list[int]]
    order: int
    grid_id: dict[tuple[int, int, int], int] = field(repr=False)
    shape: tuple[int, int, int]  # (nx, ny, nz) element divisions

    @property
    def n_nodes(self) -> int:
        return len(self.nodes)

    @property
    def n_dof(self) -> int:
        return 3 * len(self.nodes)

    def node_set(self, predicate) -> list[int]:
        return sorted(n for n, xyz in self.nodes.items() if predicate(*xyz))


def box_hex_mesh(L: float, H: float, B: float, nx: int, ny: int, nz: int,
                 order: int = 1) -> HexMesh:
    """Build a structured 8-node (order 1) or 20-node serendipity (order 2) mesh.

    Node numbering and element connectivity follow the Abaqus/CalculiX convention.
    """
    s = 2 if order == 2 else 1  # grid points per element edge
    gx, gy, gz = s * nx + 1, s * ny + 1, s * nz + 1
    grid_id: dict[tuple[int, int, int], int] = {}
    nodes: dict[int, tuple[float, float, float]] = {}
    nid = 0
    for k in range(gz):
        for j in range(gy):
            for i in range(gx):
                # Serendipity elements have no face-centre or body-centre nodes:
                # drop grid points with two or more odd indices.
                if order == 2 and (i % 2) + (j % 2) + (k % 2) >= 2:
                    continue
                nid += 1
                grid_id[(i, j, k)] = nid
                nodes[nid] = (L * i / (gx - 1), H * j / (gy - 1), B * k / (gz - 1))

    elements = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                I, J, K = s * i, s * j, s * k
                if order == 1:
                    c = [(I, J, K), (I + 1, J, K), (I + 1, J + 1, K), (I, J + 1, K),
                         (I, J, K + 1), (I + 1, J, K + 1), (I + 1, J + 1, K + 1), (I, J + 1, K + 1)]
                else:
                    c = [(I, J, K), (I + 2, J, K), (I + 2, J + 2, K), (I, J + 2, K),
                         (I, J, K + 2), (I + 2, J, K + 2), (I + 2, J + 2, K + 2), (I, J + 2, K + 2),
                         # mid-edge nodes 9-12 (bottom face), 13-16 (top face), 17-20 (vertical)
                         (I + 1, J, K), (I + 2, J + 1, K), (I + 1, J + 2, K), (I, J + 1, K),
                         (I + 1, J, K + 2), (I + 2, J + 1, K + 2), (I + 1, J + 2, K + 2), (I, J + 1, K + 2),
                         (I, J, K + 1), (I + 2, J, K + 1), (I + 2, J + 2, K + 1), (I, J + 2, K + 1)]
                elements.append([grid_id[p] for p in c])
    return HexMesh(nodes, elements, order, grid_id, (nx, ny, nz))


def face_traction_loads(mesh: HexMesh, i_face: int, total_force: float) -> dict[int, float]:
    """Consistent nodal forces for a uniform traction on the grid plane i = i_face.

    The traction acts on the y-z face of the box and sums to total_force.
    Linear quad face: 1/4 of the face force at each corner.
    8-node serendipity face: -1/12 at corners, +1/3 at mid-edges.
    """
    nx, ny, nz = mesh.shape
    loads: dict[int, float] = {}
    f_face = total_force / (ny * nz)
    s = 2 if mesh.order == 2 else 1
    for k in range(nz):
        for j in range(ny):
            J, K = s * j, s * k
            if mesh.order == 1:
                pts = [((J, K), 0.25), ((J + 1, K), 0.25), ((J + 1, K + 1), 0.25), ((J, K + 1), 0.25)]
            else:
                corners = [(J, K), (J + 2, K), (J + 2, K + 2), (J, K + 2)]
                mids = [(J + 1, K), (J + 2, K + 1), (J + 1, K + 2), (J, K + 1)]
                pts = [(p, -1.0 / 12.0) for p in corners] + [(p, 1.0 / 3.0) for p in mids]
            for (jj, kk), w in pts:
                nid = mesh.grid_id[(i_face, jj, kk)]
                loads[nid] = loads.get(nid, 0.0) + w * f_face
    return loads


def _id_lines(ids, per_line: int = 10) -> str:
    ids = list(ids)
    return "\n".join(", ".join(str(i) for i in ids[p:p + per_line])
                     for p in range(0, len(ids), per_line))


def write_deck(path: Path, mesh: HexMesh, etype: str, E: float, nu: float,
               fixed_nodes: list[int], cloads: dict[int, tuple[int, float]],
               output_nset: list[int], heading: str) -> None:
    """Write a linear-static deck readable by both ccx and Abaqus/Standard."""
    out = [f"*HEADING\n{heading}", "*NODE, NSET=NALL"]
    out += [f"{n}, {x:.10g}, {y:.10g}, {z:.10g}" for n, (x, y, z) in mesh.nodes.items()]
    out.append(f"*ELEMENT, TYPE={etype}, ELSET=EALL")
    for e, conn in enumerate(mesh.elements, start=1):
        # Abaqus allows at most 16 entries on the first data line of an element
        ids = [e] + conn
        out.append(", ".join(map(str, ids[:16])) + ("," if len(ids) > 16 else ""))
        if len(ids) > 16:
            out.append(", ".join(map(str, ids[16:])))
    out += ["*NSET, NSET=FIXED", _id_lines(fixed_nodes),
            "*NSET, NSET=TIP", _id_lines(output_nset),
            "*MATERIAL, NAME=MAT", "*ELASTIC", f"{E:.10g}, {nu:.10g}",
            "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT",
            "*STEP", "*STATIC", "*BOUNDARY", "FIXED, 1, 3, 0.", "*CLOAD"]
    out += [f"{n}, {dof}, {f:.12g}" for n, (dof, f) in sorted(cloads.items())]
    out += ["*NODE PRINT, NSET=TIP", "U", "*END STEP", ""]
    path.write_text("\n".join(out))


_NUM = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"
_ROW = re.compile(rf"^\s*(\d+)\s+({_NUM})\s+({_NUM})\s+({_NUM})\s*$")


def parse_dat_displacements(dat_path: Path) -> dict[int, np.ndarray]:
    """Read nodal displacements printed by *NODE PRINT from a ccx or Abaqus .dat file.

    Takes the last printed table (final increment) in either format.
    """
    text = Path(dat_path).read_text(errors="replace")
    marker = "N O D E   O U T P U T"  # Abaqus; ccx has no such banner
    if marker in text:
        text = text[text.rindex(marker):]
    elif "displacements" in text:
        text = text[text.rindex("displacements"):]
    disp = {}
    for line in text.splitlines():
        m = _ROW.match(line)
        if m:
            vals = [float(v.replace("D", "E").replace("d", "E")) for v in m.groups()[1:]]
            disp[int(m.group(1))] = np.array(vals)
    if not disp:
        raise ValueError(f"no displacement table found in {dat_path}")
    return disp


def _find_versioned_ccx(ccx: str) -> str:
    """Homebrew installs CalculiX as e.g. ccx_2.23 with no plain `ccx`; a shell
    alias is invisible to subprocess, so search PATH for the newest ccx_<ver>."""
    if ccx != "ccx":
        raise FileNotFoundError(f"CalculiX executable not found: {ccx}")
    found = {}
    for d in os.environ.get("PATH", "").split(os.pathsep):
        for p in Path(d).glob("ccx_[0-9]*") if Path(d).is_dir() else ():
            ver = tuple(int(v) for v in p.name[4:].split(".") if v.isdigit())
            if os.access(p, os.X_OK) and ver:
                found.setdefault(ver, str(p))
    if not found:
        raise FileNotFoundError("CalculiX not found: no `ccx` or `ccx_<version>` on PATH; "
                                "install it or pass --ccx /path/to/ccx")
    return found[max(found)]


def run_ccx(deck: Path, ccx: str = "ccx", threads: int | None = None) -> Path:
    """Run CalculiX on a deck; returns the .dat path. Raises on solver failure."""
    exe = shutil.which(ccx) or _find_versioned_ccx(ccx)
    env = dict(os.environ)
    # Default to 1 thread: multithreaded SPOOLES in ccx 2.23 intermittently returns
    # wrong solutions (observed 0.06 mm vs 1.403 mm on the same deck at 8 threads).
    env["OMP_NUM_THREADS"] = env["CCX_NPROC_EQUATION_SOLVER"] = str(threads or 1)
    job = deck.with_suffix("")
    res = subprocess.run([exe, "-i", job.name], cwd=deck.parent, env=env,
                         capture_output=True, text=True)
    (deck.parent / f"{job.name}.ccx.log").write_text(res.stdout + res.stderr)
    dat = job.with_suffix(".dat")
    if res.returncode != 0 or "*ERROR" in res.stdout or not dat.exists():
        raise RuntimeError(f"ccx failed on {deck.name}; see {job.name}.ccx.log")
    return dat
