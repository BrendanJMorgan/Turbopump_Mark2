"""
Mark 2 primary-flow cycle model - generic assembler.

The cycle topology is not written here.  It lives in ``cycle.yaml``:

    * ``states``    solver unknowns and fixed boundary States
    * ``elements``  FullFlow flow elements / volumes, addressed by id
    * ``build``     an ordered manifest that interleaves ``physics:`` hooks
                    (named callables in ``physics.py``) with ``element:`` refs

``build_cycle`` walks that manifest: it creates every State, then for each build
step either calls the named physics hook or instantiates the referenced element,
resolving its kwargs from the reference grammar documented in ``cycle.yaml``.

Unknowns and residuals are unchanged from the hand-wired model:

    branch mass flows            - DarcyWeisbach momentum residuals
    node pressures               - Volume mass balances
    pump mass flows              - pump-curve discharge-pressure residuals
    chamber / gas generator P    - Volume mass balances against nozzle / turbine
    fuel shaft speed             - fuel shaft torque balance
    turbine torque               - turbine efficiency residual
"""

from __future__ import annotations

import importlib
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import yaml

from fullflow import (
    Balance,
    ConstantDensityPump,
    DarcyWeisbach,
    DischargeCoefficient,
    GasTurbine,
    IsentropicNozzle,
    Lookup,
    Map,
    Network,
    Rotor,
    State,
    Volume,
)

from components import ChokedFlowParameter, SynchronousLink

DEFAULT_CONFIG = Path(__file__).resolve().parent / "cycle.yaml"

# type name in the YAML -> FullFlow / custom component class
COMPONENT_TYPES: dict[str, type] = {
    "DarcyWeisbach": DarcyWeisbach,
    "ConstantDensityPump": ConstantDensityPump,
    "DischargeCoefficient": DischargeCoefficient,
    "IsentropicNozzle": IsentropicNozzle,
    "GasTurbine": GasTurbine,
    "Rotor": Rotor,
    "Volume": Volume,
    "ChokedFlowParameter": ChokedFlowParameter,
    "SynchronousLink": SynchronousLink,
}

# YAML state kind -> (CycleModel dict attribute, default track group)
_STATE_KINDS = {
    "pressure": ("pressure", "pressure"),
    "flow": ("mass_flow", "mass_flow"),
    "temperature": ("temperature", "temperature"),
    "scalar": ("scalar", "shaft"),
}

# reference tag -> BuildContext container attribute
_REF_TAGS = {
    "p": "pressure",
    "m": "mass_flow",
    "t": "temperature",
    "s": "scalar",
    "d": "derived",
    "param": "params",
}


@dataclass
class CycleModel:
    """Handles into the built network, for post-processing and for the driver."""

    network: Network
    config: dict
    pressure: dict[str, State] = field(default_factory=dict)
    mass_flow: dict[str, State] = field(default_factory=dict)
    temperature: dict[str, State] = field(default_factory=dict)
    scalar: dict[str, State] = field(default_factory=dict)
    shaft: dict[str, State] = field(default_factory=dict)
    performance: dict[str, State] = field(default_factory=dict)
    derived: dict[str, Any] = field(default_factory=dict)
    components: dict[str, Any] = field(default_factory=dict)


def _maybe_number(value):
    """Coerce numeric strings to int/float.

    PyYAML follows YAML 1.1, which only recognises an exponent when it carries an
    explicit sign, so ``2.5e5`` loads as the string ``"2.5e5"``.  Rather than
    litter the YAML with ``+`` signs, numbers are coerced here at the point of
    use.  Non-numeric strings (``hp``, ``O2(L)``) pass through untouched.
    """
    if not isinstance(value, str):
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def load_config(config: str | Path | dict) -> dict:
    if isinstance(config, Mapping):
        data = dict(config)
        data.setdefault("__dir__", str(Path.cwd()))
        return data
    path = Path(config)
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    data["__dir__"] = str(path.resolve().parent)
    return data


class BuildContext:
    """Shared state threaded through every physics hook and element build."""

    def __init__(self, model: CycleModel, config: dict):
        self.model = model
        self.net = model.network
        self.config = config
        self.config_dir = Path(config["__dir__"])

        self.params: dict[str, Any] = {
            **config.get("parameters", {}),
            **config.get("param_blocks", {}),
        }
        self.maps_cfg: dict[str, Any] = config.get("maps", {})
        self.propellants = SimpleNamespace(**config.get("propellants", {}))

        self.pressure = model.pressure
        self.mass_flow = model.mass_flow
        self.temperature = model.temperature
        self.scalar = model.scalar
        self.derived = model.derived
        self.components = model.components
        self.performance = model.performance

        self._element_specs = {e["id"]: e for e in config.get("elements", [])}
        self.derived["ambient_pressure"] = State(self.params["ambient_pressure"])

    # -- helpers used by physics hooks -------------------------------------
    def track(self, path: str, state: State) -> None:
        self.net.track(path, state)

    def set_temperature(self, name: str, state: State, *, track: bool = True) -> State:
        self.temperature[name] = state
        if track:
            self.track(f"temperature/{name}", state)
        return state

    def set_performance(self, name: str, state: State, *, track: bool = True) -> State:
        self.performance[name] = state
        if track:
            self.track(f"performance/{name}", state)
        return state

    def load_map(self, key: str, inputs: dict[str, State],
                 outputs: list[str]) -> Map:
        spec = self.maps_cfg[key]
        return Map.from_hdf5(
            spec.get("name", key),
            self.net,
            self.config_dir / spec["file"],
            spec["group"],
            inputs=inputs,
            outputs=outputs,
            extrapolate=spec.get("extrapolate", False),
        )

    # -- reference resolution -------------------------------------------------
    def resolve(self, spec: Any) -> Any:
        if isinstance(spec, str):
            return self._resolve_ref(spec)
        if isinstance(spec, Mapping):
            if set(spec) == {"sum"}:
                terms = [self.resolve(s) for s in spec["sum"]]
                total = terms[0]
                for term in terms[1:]:
                    total = total + term
                return total
            if set(spec) == {"state"}:
                return State(self.resolve(spec["state"]))
            return {k: self.resolve(v) for k, v in spec.items()}
        if isinstance(spec, list):
            return [self.resolve(s) for s in spec]
        return spec

    def _resolve_ref(self, ref: str):
        tag, sep, rest = ref.partition(":")
        if not sep or tag not in _REF_TAGS:
            return _maybe_number(ref)  # literal: number-like string or e.g. mode: hp
        container = getattr(self, _REF_TAGS[tag])
        parts = rest.split(".")
        try:
            obj = container[parts[0]]
        except KeyError as exc:
            raise KeyError(
                f"unresolved reference {ref!r}: {parts[0]!r} not in {tag!r} namespace "
                f"(available: {sorted(container)})"
            ) from exc
        for attr in parts[1:]:
            obj = getattr(obj, attr)
        return _maybe_number(obj)

    # -- construction ------------------------------------------------------
    def declare_states(self) -> None:
        for state_id, spec in self.config.get("states", {}).items():
            kind = spec["kind"]
            dict_attr, group = _STATE_KINDS[kind]
            container = getattr(self, dict_attr)
            name = spec.get("name", state_id)

            if "value" in spec:
                state = State(self.resolve(spec["value"]))
            else:
                bounds = (
                    tuple(_maybe_number(b) for b in spec["bounds"])
                    if "bounds" in spec
                    else (None, None)
                )
                state = State(_maybe_number(spec["guess"]), bounds=bounds)

            container[name] = state

            track = spec.get("track", f"{group}/{name}")
            if track:
                self.track(track, state)

    def build_element(self, element_id: str):
        spec = self._element_specs[element_id]
        cls = COMPONENT_TYPES[spec["type"]]

        kwargs: dict[str, Any] = {}
        spread = spec.get("spread")
        for block in ([spread] if isinstance(spread, str) else (spread or [])):
            kwargs.update(self.resolve(block))
        for key, value in (spec.get("kwargs") or {}).items():
            kwargs[key] = self.resolve(value)

        obj = cls(spec.get("name", element_id), self.net, **kwargs)
        self.components[spec.get("store", element_id)] = obj

        for attr, target in (spec.get("exports") or {}).items():
            tag, _, key = target.partition(":")
            getattr(self, _REF_TAGS[tag])[key] = getattr(obj, attr)
        return obj


def _call_hook(ref: str, ctx: BuildContext) -> None:
    module_name, sep, func_name = ref.rpartition(":")
    module = importlib.import_module(module_name if sep else "physics")
    getattr(module, func_name)(ctx)


def build_cycle(config: str | Path | dict = DEFAULT_CONFIG) -> CycleModel:
    config = load_config(config)
    net = Network(config.get("name", "cycle"))
    model = CycleModel(network=net, config=config)
    ctx = BuildContext(model, config)

    ctx.declare_states()

    for step in config["build"]:
        (kind, ref), = step.items()
        if kind == "physics":
            _call_hook(ref, ctx)
        elif kind == "element":
            ctx.build_element(ref)
        else:
            raise ValueError(f"unknown build step {kind!r} (expected 'physics' or 'element')")

    return model


if __name__ == "__main__":
    import sys

    m = build_cycle(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CONFIG)
    print(f"built {m.network.name}: "
          f"{len(m.components)} components, "
          f"{len(m.pressure)} pressures, {len(m.mass_flow)} flows")
