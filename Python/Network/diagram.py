"""
Wiring diagram for the Mark 2 cycle, generated straight from ``cycle.yaml``.

No import of ``model.py``, no ``ast`` parse, no solve, no CEA.  The YAML already
carries the topology: pressure ``states`` are graph nodes, ``elements`` are the
edges, and each Volume names the node it sits on.

Emits ``network_diagram.md``: a Mermaid flowchart of the flow path plus a
shaft-side flowchart.

Run:  python diagram.py [--config PATH] [--out PATH]
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import yaml

# kwarg names that carry the upstream / downstream pressure station of an element
UPSTREAM_KEYS = ("upstream_pressure", "upstream_total_pressure")
DOWNSTREAM_KEYS = ("downstream_pressure", "discharge_pressure", "downstream_total_pressure")

SIDE_PREFIXES = [
    ("ox", "Oxidizer"),
    ("fuel", "Fuel"),
    ("jacket", "Regen / TCA fuel"),
    ("regen", "Regen / TCA fuel"),
    ("main_fuel_valve", "Fuel"),
    ("chamber", "Combustion"),
    ("gas_generator", "Combustion"),
    ("turbine", "Combustion"),
]

REF_TAGS = ("p", "m", "t", "s", "d", "param")


# --------------------------------------------------------------------------
# reference helpers
# --------------------------------------------------------------------------
def ref_parts(value):
    """('p', 'ox_manifold') for 'p:ox_manifold.head', ('state', None) for {state: ...}."""
    if isinstance(value, str) and ":" in value:
        tag, rest = value.split(":", 1)
        if tag in REF_TAGS:
            return tag, rest.split(".")[0]
    if isinstance(value, dict) and set(value) == {"state"}:
        return "state", None
    return None, None


def station_ref(value):
    tag, name = ref_parts(value)
    return name if tag == "p" else None


def flow_ref(value):
    tag, name = ref_parts(value)
    return name if tag == "m" else None


# --------------------------------------------------------------------------
# model view of the YAML
# --------------------------------------------------------------------------
class CycleView:
    def __init__(self, config: dict):
        self.name = config.get("name", "cycle")
        self.elements = config.get("elements", [])

        self.stations = [
            spec.get("name", key)
            for key, spec in config.get("states", {}).items()
            if spec.get("kind") == "pressure"
        ]

        self.volume_name: dict[str, str] = {}
        for el in self.elements:
            if el.get("type") == "Volume":
                station = station_ref((el.get("kwargs") or {}).get("pressure"))
                if station:
                    self.volume_name[station] = el["name"]

    def flow_elements(self):
        for el in self.elements:
            if el.get("type") == "Volume":
                continue
            kwargs = el.get("kwargs") or {}
            has_port = any(k in kwargs for k in UPSTREAM_KEYS + DOWNSTREAM_KEYS)
            if el.get("category") == "flow" or has_port:
                yield el

    def shaft_elements(self):
        return [el for el in self.elements if el.get("category") == "shaft"]


# --------------------------------------------------------------------------
# graph assembly
# --------------------------------------------------------------------------
class Graph:
    def __init__(self):
        self.node_label: dict[str, str] = {}
        self.node_kind: dict[str, str] = {}          # 'place' | 'boundary'
        self.edges: list[tuple[str, str, str]] = []


def build_graph(view: CycleView) -> Graph:
    g = Graph()

    def place(key: str) -> str:
        if key not in g.node_label:
            name = view.volume_name.get(key)
            g.node_label[key] = f"{key}<br/>{name}" if name else key
            g.node_kind[key] = "place"
        return key

    def boundary(name: str) -> str:
        if name not in g.node_label:
            g.node_label[name] = name
            g.node_kind[name] = "boundary"
        return name

    for station in view.stations:
        place(station)

    for el in view.flow_elements():
        kwargs = el.get("kwargs") or {}
        name = el["name"]

        upstream = next(
            (station_ref(kwargs[k]) for k in UPSTREAM_KEYS if station_ref(kwargs.get(k))),
            None,
        )
        downstream = None
        downstream_boundary = None
        for key in DOWNSTREAM_KEYS:
            if key not in kwargs:
                continue
            station = station_ref(kwargs[key])
            if station:
                downstream = station
            elif ref_parts(kwargs[key])[0] == "state":
                downstream_boundary = "turbine exit" if "turbine" in name.lower() else f"{name} exit"

        carried = (
            flow_ref(kwargs.get("mass_flow"))
            or flow_ref((el.get("exports") or {}).get("mass_flow"))
            or el.get("type", "")
        )
        label = f"{name}<br/>{carried}"

        src = place(upstream) if upstream else None
        if downstream:
            dst = place(downstream)
        elif downstream_boundary:
            dst = boundary(downstream_boundary)
        else:
            dst = None

        if src is None and dst is not None:
            src = boundary(f"{name} inlet")
        if dst is None and src is not None:
            dst = boundary(f"{name} outlet")
        if src is None and dst is None:
            continue

        g.edges.append((src, dst, label))

    return g


def side_of(key: str) -> str:
    for prefix, name in SIDE_PREFIXES:
        if key.startswith(prefix):
            return name
    return "Other"


# --------------------------------------------------------------------------
# mermaid
# --------------------------------------------------------------------------
def sanitize(key: str) -> str:
    return "n_" + "".join(c if (c.isalnum() or c == "_") else "_" for c in key)


def mermaid_flow(g: Graph) -> str:
    lines = ["flowchart LR"]

    grouped: dict[str, list[str]] = defaultdict(list)
    for key in g.node_label:
        grouped[side_of(key) if g.node_kind[key] == "place" else "Boundaries"].append(key)

    for group, keys in grouped.items():
        lines.append(f'    subgraph {sanitize(group)}["{group}"]')
        for key in keys:
            nid = sanitize(key)
            if g.node_kind[key] == "boundary":
                lines.append(f'        {nid}[["{g.node_label[key]}"]]')
            else:
                lines.append(f'        {nid}(["{g.node_label[key]}"])')
        lines.append("    end")

    for src, dst, label in g.edges:
        lines.append(f'    {sanitize(src)} -->|"{label}"| {sanitize(dst)}')

    return "\n".join(lines)


def mermaid_shaft(view: CycleView) -> str:
    lines = ["flowchart LR"]
    seen: set[str] = set()

    def node(nid: str, label: str, shape: str = "rect") -> None:
        if nid in seen:
            return
        seen.add(nid)
        if shape == "round":
            lines.append(f'    {sanitize(nid)}(["{label}"])')
        else:
            lines.append(f'    {sanitize(nid)}["{label}"]')

    shaft = view.shaft_elements()
    for el in shaft:
        name = el["name"]
        node(name, name)
        for arg, value in (el.get("kwargs") or {}).items():
            tag, ref_name = ref_parts(value)
            if tag not in ("s", "d", "m", "p"):
                continue
            node(ref_name, ref_name, shape="round")
            lines.append(f'    {sanitize(ref_name)} -->|"{arg}"| {sanitize(name)}')

    if not shaft:
        lines.append("    empty[No shaft components found]")
    return "\n".join(lines)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main() -> None:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=here / "cycle.yaml")
    ap.add_argument("--out", type=Path, default=here / "network_diagram.md")
    args = ap.parse_args()

    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    view = CycleView(config)
    g = build_graph(view)

    md = [
        f"# {view.name} cycle network - static wiring diagram",
        "",
        f"Generated by `diagram.py` from `{args.config.name}`. No solve, no CEA.",
        "",
        "## Flow topology",
        "",
        "```mermaid",
        mermaid_flow(g),
        "```",
        "",
        "## Shaft / synchronous link",
        "",
        "```mermaid",
        mermaid_shaft(view),
        "```",
        "",
    ]
    args.out.write_text("\n".join(md), encoding="utf-8")

    print(f"wrote {args.out}")
    print(
        f"  {len(g.node_label)} nodes, {len(g.edges)} edges, "
        f"{len(view.stations)} pressure stations, "
        f"{sum(1 for _ in view.flow_elements())} flow elements"
    )


if __name__ == "__main__":
    main()
