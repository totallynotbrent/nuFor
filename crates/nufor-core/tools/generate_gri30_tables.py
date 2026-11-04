#!/usr/bin/env python3
"""Generate the nuFor GRI-Mech 3.0 tables and the Cantera falloff oracle.

usage:

    python3 generate_gri30_tables.py [yaml_path]

the yaml path defaults to the gri30.yaml shipped with Cantera 3.2.0. two
files are written:

    ../src/gri30_generated.rs                   generated tables
    ../data/reference/gri30_falloff_oracle.csv  cantera oracle fixture

both outputs are deterministic: the same inputs give byte-identical
files, because every floating point value is emitted with python's
shortest round-tripping repr and nothing carries a timestamp.

unit conventions, si with m, kmol, s:

    concentrations      kmol/m3 (numerically mol/l)
    gas constant        8314.46261815324 j/(kmol k)
    activation energy   j/kmol (the yaml's cal/mol times 4184)
    bimolecular rate    m3/(kmol s)
    trimolecular rate   m6/(kmol2 s)

the rate constants are copied from cantera's own converted rate objects,
which is the same as applying the factors by hand. a yaml
pre-exponential factor in cm3/(mol s) becomes m3/(kmol s) through
1e-6 m3/cm3 times 1e-3 kmol/mol, that is 1e-3 overall, because the
concentrations move from mol/cm3 to kmol/m3 by a factor 1e3 in both the
concentration and the reciprocal-concentration places. a yaml
pre-exponential factor in cm6/(mol2 s) becomes m6/(kmol2 s) through
1e-12 times 1e-6, that is 1e-6 overall, and a yaml pre-exponential
factor in 1/s is unchanged. a yaml activation energy in cal/mol becomes
j/kmol through 4184.

the third-body concentration follows cantera: the bath is the sum of
every species concentration times its collision efficiency, with a
species missing from the efficiency table contributing at the default
efficiency of one.
"""

import csv
import hashlib
import io
import math
import os
import re
import subprocess
import sys

DEFAULT_YAML = (
    "/home/ebarleb/.local/lib/python3.12/site-packages/cantera/data/gri30.yaml"
)

HERE = os.path.dirname(os.path.abspath(__file__))
CORE_DIR = os.path.dirname(HERE)
RS_OUT = os.path.join(CORE_DIR, "src", "gri30_generated.rs")
CSV_OUT = os.path.join(CORE_DIR, "data", "reference", "gri30_falloff_oracle.csv")

RU_KMOL = 8314.46261815324

TEMPERATURES_K = (600.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0)
PRESSURES_ATM = (0.1, 1.0, 10.0)

MIXTURES = (
    ("H2:0.30;H2O:0.30;O2:0.20;AR:0.20", "hydrogen-steam"),
    ("CH4:0.20;C2H4:0.20;O2:0.20;H2O:0.20;N2:0.20", "hydrocarbon"),
    ("CH2:0.05;N2O:0.05;AR:0.90", "radical-dilute"),
)


def die(message):
    """abort the generator with a diagnostic"""
    sys.stderr.write("generate_gri30_tables.py: %s\n" % message)
    raise SystemExit(1)


def ruslit(value):
    """format a float as a rust literal that round-trips exactly"""
    value = float(value)
    if value != value or value in (float("inf"), float("-inf")):
        die("cannot represent %r in rust" % value)
    text = repr(value)
    if "." not in text and "e" not in text and "E" not in text:
        text += ".0"
    return text


TOKEN_RE = re.compile(r"^\s*(\d+(?:\.\d+)?)?\s*([A-Za-z0-9_()\[\]]+)\s*$")


def parse_side(side):
    """parse one side of an equation into ordered (species, coefficient) pairs"""
    cleaned = side.replace("(+M)", "")
    pairs = []
    for chunk in cleaned.split("+"):
        token = chunk.strip()
        if not token:
            continue
        if token.upper() == "M":
            continue
        match = TOKEN_RE.match(token)
        if not match:
            die("cannot parse side fragment %r of %r" % (chunk, side))
        coefficient = float(match.group(1)) if match.group(1) else 1.0
        pairs.append((match.group(2), coefficient))
    return pairs


def parse_reaction_entry(entry):
    """normalise one yaml reaction entry into a plain dictionary"""
    equation = entry["equation"]
    reversible = False
    reactant_text = ""
    product_text = ""
    if "<=>" in equation:
        reactant_text, product_text = equation.split("<=>", 1)
        reversible = True
    elif "=>" in equation:
        reactant_text, product_text = equation.split("=>", 1)
    else:
        die("cannot parse equation %r" % equation)
    return {
        "equation": equation,
        "reactants": parse_side(reactant_text),
        "products": parse_side(product_text),
        "reversible": reversible,
        "duplicate": bool(entry.get("duplicate", False)),
        "efficiencies": dict(entry.get("efficiencies", {}) or {}),
        "troe": dict(entry["Troe"]) if "Troe" in entry else None,
    }


def sha256_of(path):
    """return the hex sha256 digest of a file"""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 16), b""):
            digest.update(block)
    return digest.hexdigest()


def format_rust(text):
    """run rustfmt over the generated source so the output is stable"""
    completed = subprocess.run(
        ["rustfmt", "--edition", "2021"],
        input=text,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        die("rustfmt failed: %s" % completed.stderr.strip())
    return completed.stdout


def write_if_changed(path, text):
    """write a file atomically, leaving it alone when the content matches"""
    directory = os.path.dirname(path)
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    payload = text.encode("utf-8")
    if os.path.exists(path):
        with open(path, "rb") as handle:
            if handle.read() == payload:
                return
    temporary = path + ".new"
    with open(temporary, "wb") as handle:
        handle.write(payload)
    os.replace(temporary, path)


def branch_rate(branch, temperature):
    """evaluate one cantera arrhenius branch in cantera's si units"""
    return (
        branch.pre_exponential_factor
        * temperature ** branch.temperature_exponent
        * math.exp(-branch.activation_energy / (RU_KMOL * temperature))
    )


def main(argv):
    """generate the tables and the oracle fixture"""
    yaml_path = argv[1] if len(argv) > 1 else DEFAULT_YAML
    if not os.path.exists(yaml_path):
        die("mechanism yaml not found: %s" % yaml_path)

    sys.path.insert(0, "/home/ebarleb/.local/lib/python3.12/site-packages")
    try:
        import yaml
        import cantera
    except ImportError as exc:
        die("PyYAML and Cantera 3.2.0 are required (%s)" % exc)

    class MechanismLoader(yaml.SafeLoader):
        """a loader that leaves the species names NO and ON as strings."""

    MechanismLoader.yaml_implicit_resolvers = {
        key: [(tag, regexp) for tag, regexp in resolvers
               if tag != "tag:yaml.org,2002:bool"]
        for key, resolvers in yaml.SafeLoader.yaml_implicit_resolvers.items()
    }

    with open(yaml_path, "r", encoding="utf-8") as handle:
        document = yaml.load(handle, Loader=MechanismLoader)

    gas = cantera.Solution(yaml_path)
    phase = document["phases"][0]
    species_names = list(phase["species"])
    species_index = {name: i for i, name in enumerate(species_names)}
    if sorted(species_names) != sorted(sp.name for sp in gas.species()):
        die("yaml and cantera species lists disagree")
    weights = {sp.name: sp.molecular_weight / 1000.0 for sp in gas.species()}

    entries = [parse_reaction_entry(e) for e in document["reactions"]]
    if len(entries) != gas.n_reactions:
        die("yaml reaction count %d disagrees with cantera %d"
            % (len(entries), gas.n_reactions))

    reactions = []
    falloff_indices = []
    for index, entry in enumerate(entries):
        cantera_reaction = gas.reaction(index)
        reaction_type = cantera_reaction.reaction_type

        # cantera's own third-body description decides the collider's
        # contribution: a NAMED collider contributes only its own
        # concentration at an efficiency of one, and every other species
        # contributes nothing, while a generic 'M' collider takes the
        # efficiency table and a default efficiency of one.
        third_body = cantera_reaction.third_body
        if third_body is None:
            default_efficiency = 0.0
            efficiencies = {}
        else:
            default_efficiency = float(third_body.default_efficiency)
            efficiencies = {
                name: float(value) for name, value in third_body.efficiencies.items()
            }

        # cantera's stoichiometry leaves the collider out of both sides,
        # so the table follows cantera's dictionaries rather than the
        # equation text, which spells the collider out for a named one.
        record = {
            "index": index,
            "equation": entry["equation"],
            "reactants": [
                (species_index[name], float(coefficient))
                for name, coefficient in cantera_reaction.reactants.items()
            ],
            "products": [
                (species_index[name], float(coefficient))
                for name, coefficient in cantera_reaction.products.items()
            ],
            "reversible": entry["reversible"],
            "duplicate": entry["duplicate"],
            "efficiencies": efficiencies,
            "default_efficiency": default_efficiency,
            "kind": None,
            "falloff": None,
            "troe": None,
        }
        if reaction_type == "Arrhenius":
            record["kind"] = "Arrhenius"
            rate = cantera_reaction.rate
            record["forward"] = {
                "a": rate.pre_exponential_factor,
                "n": rate.temperature_exponent,
                "ea": rate.activation_energy,
            }
        elif reaction_type == "three-body-Arrhenius":
            record["kind"] = "ThreeBody"
            rate = cantera_reaction.rate
            record["forward"] = {
                "a": rate.pre_exponential_factor,
                "n": rate.temperature_exponent,
                "ea": rate.activation_energy,
            }
        elif reaction_type in ("falloff-Troe", "falloff-Lindemann"):
            record["kind"] = "Falloff"
            record["falloff"] = "Troe" if reaction_type == "falloff-Troe" else "Lindemann"
            low = cantera_reaction.rate.low_rate
            high = cantera_reaction.rate.high_rate
            record["low"] = {
                "a": low.pre_exponential_factor,
                "n": low.temperature_exponent,
                "ea": low.activation_energy,
            }
            record["high"] = {
                "a": high.pre_exponential_factor,
                "n": high.temperature_exponent,
                "ea": high.activation_energy,
            }
            troe = entry["troe"]
            record["troe"] = None if troe is None else {
                "a": troe["A"],
                "t3": troe["T3"],
                "t1": troe["T1"],
                "t2": troe.get("T2", 0.0),
            }
            falloff_indices.append(index)
        else:
            die("unexpected reaction type %r for %r"
                % (reaction_type, entry["equation"]))
        reactions.append(record)

    # ------------------------------------------------------------------
    # rust tables
    # ------------------------------------------------------------------
    out = io.StringIO()
    out.write("// generated by crates/nufor-core/tools/generate_gri30_tables.py,"
              " do not edit.\n")
    out.write("//\n")
    out.write("// source mechanism: gri30.yaml from cantera %s.\n" % cantera.__version__)
    out.write("// source sha256: %s\n" % sha256_of(yaml_path))
    out.write("//\n")
    out.write("// concentrations are kmol/m3, activation energies are j/kmol,\n")
    out.write("// rate constants are m3/(kmol s) for bimolecular laws and\n")
    out.write("// m6/(kmol2 s) for trimolecular ones, and the molar gas constant\n")
    out.write("// is 8314.46261815324 j/(kmol k).\n\n")

    out.write("/// the molar gas constant in j/(kmol k).\n")
    out.write("pub const GAS_CONSTANT: f64 = %s;\n\n" % ruslit(RU_KMOL))

    out.write("/// the sha256 of the mechanism yaml these tables came from.\n")
    out.write("pub const SOURCE_SHA256: &str = \"%s\";\n\n" % sha256_of(yaml_path))

    out.write("/// the number of species in the mechanism.\n")
    out.write("pub const SPECIES_COUNT: usize = %d;\n\n" % len(species_names))

    out.write("/// the number of reactions in the mechanism.\n")
    out.write("pub const REACTION_COUNT: usize = %d;\n\n" % len(reactions))

    out.write("/// one chemical species.\n")
    out.write("pub struct Species {\n")
    out.write("    /// the species name as cantera spells it.\n")
    out.write("    pub name: &'static str,\n")
    out.write("    /// the molecular weight in kg/mol.\n")
    out.write("    pub molecular_weight: f64,\n")
    out.write("}\n\n")

    out.write("/// the arrhenius rate parameters in si units.\n")
    out.write("#[derive(Clone, Copy)]\n")
    out.write("pub struct Arrhenius {\n")
    out.write("    /// the pre-exponential factor.\n")
    out.write("    pub a: f64,\n")
    out.write("    /// the temperature exponent.\n")
    out.write("    pub n: f64,\n")
    out.write("    /// the activation energy in j/kmol.\n")
    out.write("    pub activation_energy: f64,\n")
    out.write("}\n\n")

    out.write("/// the troe broadening parameters in cantera's convention.\n")
    out.write("#[derive(Clone, Copy)]\n")
    out.write("pub struct TroeParameters {\n")
    out.write("    /// the mixture weight between the two exponentials.\n")
    out.write("    pub a: f64,\n")
    out.write("    /// the low-temperature characteristic temperature in kelvin.\n")
    out.write("    pub t3: f64,\n")
    out.write("    /// the high-temperature characteristic temperature in kelvin.\n")
    out.write("    pub t1: f64,\n")
    out.write("    /// the temperature of the added correction term in kelvin.\n")
    out.write("    pub t2: f64,\n")
    out.write("}\n\n")

    out.write("/// one stoichiometric entry of a reaction.\n")
    out.write("pub struct Stoichiometry {\n")
    out.write("    /// the species index into the species table.\n")
    out.write("    pub species: usize,\n")
    out.write("    /// the stoichiometric coefficient.\n")
    out.write("    pub coefficient: f64,\n")
    out.write("}\n\n")

    out.write("/// one third-body efficiency.\n")
    out.write("pub struct Efficiency {\n")
    out.write("    /// the species index into the species table.\n")
    out.write("    pub species: usize,\n")
    out.write("    /// the collisional efficiency.\n")
    out.write("    pub efficiency: f64,\n")
    out.write("}\n\n")

    out.write("/// the kind of kinetic description a reaction carries.\n")
    out.write("pub enum ReactionKind {\n")
    out.write("    /// a plain arrhenius rate law.\n")
    out.write("    Arrhenius,\n")
    out.write("    /// an arrhenius rate law with an explicit third body.\n")
    out.write("    ThreeBody,\n")
    out.write("    /// a pressure dependent falloff rate law.\n")
    out.write("    Falloff,\n")
    out.write("}\n\n")

    out.write("/// the falloff broadening a reaction uses.\n")
    out.write("pub enum FalloffKind {\n")
    out.write("    /// the troe broadening.\n")
    out.write("    Troe,\n")
    out.write("    /// the plain lindemann blend with no broadening.\n")
    out.write("    Lindemann,\n")
    out.write("}\n\n")

    out.write("/// one reaction of the mechanism.\n")
    out.write("pub struct Reaction {\n")
    out.write("    /// the index of this reaction in the mechanism.\n")
    out.write("    pub index: usize,\n")
    out.write("    /// the reaction equation as the mechanism spells it.\n")
    out.write("    pub equation: &'static str,\n")
    out.write("    /// the kinetic description of this reaction.\n")
    out.write("    pub kind: ReactionKind,\n")
    out.write("    /// the reactants in the order the equation writes them.\n")
    out.write("    pub reactants: &'static [Stoichiometry],\n")
    out.write("    /// the products in the order the equation writes them.\n")
    out.write("    pub products: &'static [Stoichiometry],\n")
    out.write("    /// whether the reaction reverses.\n")
    out.write("    pub reversible: bool,\n")
    out.write("    /// whether the reaction duplicates an earlier one.\n")
    out.write("    pub duplicate: bool,\n")
    out.write("    /// whether the reaction consumes a third body.\n")
    out.write("    pub third_body: bool,\n")
    out.write("    /// the third-body efficiencies, empty when there is no third body.\n")
    out.write("    pub efficiencies: &'static [Efficiency],\n")
    out.write("    /// the efficiency a species outside the table contributes.\n")
    out.write("    ///\n")
    out.write("    /// a generic collider takes one, so every species joins the bath,\n")
    out.write("    /// while a named collider takes zero, so only the collider itself\n")
    out.write("    /// joins the bath.\n")
    out.write("    pub default_efficiency: f64,\n")
    out.write("    /// the forward arrhenius rate law, unset for falloff reactions.\n")
    out.write("    pub forward: Option<Arrhenius>,\n")
    out.write("    /// the low-pressure arrhenius rate law, set only for falloff reactions.\n")
    out.write("    pub low: Option<Arrhenius>,\n")
    out.write("    /// the high-pressure arrhenius rate law, set only for falloff reactions.\n")
    out.write("    pub high: Option<Arrhenius>,\n")
    out.write("    /// the falloff broadening, set only for falloff reactions.\n")
    out.write("    pub falloff: Option<FalloffKind>,\n")
    out.write("    /// the troe parameters, set only for troe falloff reactions.\n")
    out.write("    pub troe: Option<TroeParameters>,\n")
    out.write("}\n\n")

    out.write("/// the species of the mechanism in mechanism order.\n")
    out.write("pub static SPECIES: &[Species] = &[\n")
    for name in species_names:
        out.write("    Species { name: \"%s\", molecular_weight: %s },\n"
                  % (name, ruslit(weights[name])))
    out.write("];\n\n")

    out.write("/// the reactions of the mechanism in mechanism order.\n")
    out.write("pub static REACTIONS: &[Reaction] = &[\n")
    for record in reactions:
        out.write("    Reaction {\n")
        out.write("        index: %d,\n" % record["index"])
        out.write("        equation: \"%s\",\n"
                  % record["equation"].replace("\\", "\\\\").replace("\"", "\\\""))
        out.write("        kind: ReactionKind::%s,\n" % record["kind"])
        out.write("        reactants: &[\n")
        for species, coefficient in record["reactants"]:
            out.write("            Stoichiometry { species: %d, coefficient: %s },\n"
                      % (species, ruslit(coefficient)))
        out.write("        ],\n")
        out.write("        products: &[\n")
        for species, coefficient in record["products"]:
            out.write("            Stoichiometry { species: %d, coefficient: %s },\n"
                      % (species, ruslit(coefficient)))
        out.write("        ],\n")
        out.write("        reversible: %s,\n"
                  % ("true" if record["reversible"] else "false"))
        out.write("        duplicate: %s,\n"
                  % ("true" if record["duplicate"] else "false"))
        out.write("        third_body: %s,\n"
                  % ("true" if record["kind"] in ("ThreeBody", "Falloff") else "false"))
        out.write("        default_efficiency: %s,\n" % ruslit(record["default_efficiency"]))
        efficiencies = record["efficiencies"]
        if efficiencies:
            out.write("        efficiencies: &[\n")
            for name in sorted(efficiencies):
                out.write("            Efficiency { species: %d, efficiency: %s },\n"
                          % (species_index[name], ruslit(efficiencies[name])))
            out.write("        ],\n")
        else:
            out.write("        efficiencies: &[],\n")
        if "forward" in record:
            forward = record["forward"]
            out.write("        forward: Some(Arrhenius { a: %s, n: %s, activation_energy: %s }),\n"
                      % (ruslit(forward["a"]), ruslit(forward["n"]),
                         ruslit(forward["ea"])))
            out.write("        low: None,\n")
            out.write("        high: None,\n")
        else:
            low = record["low"]
            high = record["high"]
            out.write("        forward: None,\n")
            out.write("        low: Some(Arrhenius { a: %s, n: %s, activation_energy: %s }),\n"
                      % (ruslit(low["a"]), ruslit(low["n"]), ruslit(low["ea"])))
            out.write("        high: Some(Arrhenius { a: %s, n: %s, activation_energy: %s }),\n"
                      % (ruslit(high["a"]), ruslit(high["n"]), ruslit(high["ea"])))
        if record["kind"] == "Falloff":
            out.write("        falloff: Some(FalloffKind::%s),\n" % record["falloff"])
        else:
            out.write("        falloff: None,\n")
        if record["troe"] is not None:
            troe = record["troe"]
            out.write("        troe: Some(TroeParameters { a: %s, t3: %s, t1: %s, t2: %s }),\n"
                      % (ruslit(troe["a"]), ruslit(troe["t3"]),
                         ruslit(troe["t1"]), ruslit(troe["t2"])))
        else:
            out.write("        troe: None,\n")
        out.write("    },\n")
    out.write("];\n")
    write_if_changed(RS_OUT, format_rust(out.getvalue()))

    # ------------------------------------------------------------------
    # cantera oracle fixture
    #
    # the broadening factor comes from a clean cantera rate whose low-
    # and high-pressure branches are identical, so cantera's reduced
    # pressure is exactly the argument handed to falloff_function and
    # the reference factor carries no dependence on the transcribed
    # rate constants. the blended forward rate comes from the physical
    # reaction, assembled by cantera from its own converted branches.
    # ------------------------------------------------------------------
    csv_buffer = io.StringIO()
    writer = csv.writer(csv_buffer, lineterminator="\n")
    writer.writerow([
        "reaction_index",
        "equation",
        "temperature_K",
        "pressure_Pa",
        "mole_fraction_csv",
        "efficiency_weighted_M_mol_per_L",
        "cantera_third_body_M_mol_per_L",
        "reduced_pressure_pr",
        "broadening_factor_cantera_clean",
        "k_fwd_cantera_SI",
        "k_fwd_projected_SI",
    ])

    for index in falloff_indices:
        record = reactions[index]
        troe_parameters = record["troe"]
        clean = {
            "type": "falloff",
            "low-P-rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0},
            "high-P-rate-constant": {"A": 1.0, "b": 0.0, "Ea": 0.0},
        }
        if troe_parameters is not None:
            clean["Troe"] = {
                "A": troe_parameters["a"],
                "T3": troe_parameters["t3"],
                "T1": troe_parameters["t1"],
                "T2": troe_parameters["t2"],
            }
        clean_rate = cantera.ReactionRate.from_dict(clean)
        physical = gas.reaction(index).rate

        for temperature in TEMPERATURES_K:
            for pressure_atm in PRESSURES_ATM:
                pressure = pressure_atm * cantera.one_atm
                for mixture, _label in MIXTURES:
                    gas.TPX = temperature, pressure, mixture.replace(";", ",")
                    mole_fractions = list(gas.X)
                    concentrations = list(gas.concentrations)
                    bath = 0.0
                    for position, name in enumerate(gas.species_names):
                        efficiency = record["efficiencies"].get(
                            name, record["default_efficiency"]
                        )
                        bath += efficiency * concentrations[position]
                    cantera_bath = float(gas.third_body_concentrations[index])
                    if abs(bath - cantera_bath) > 1e-12 * cantera_bath:
                        die("third-body convention mismatch on reaction %d" % index)
                    k_low = branch_rate(physical.low_rate, temperature)
                    k_high = branch_rate(physical.high_rate, temperature)
                    reduced = bath * k_low / k_high
                    factor = float(clean_rate.falloff_function(temperature, reduced))
                    k_forward = float(gas.forward_rate_constants[index])
                    physical_factor = float(physical.falloff_function(temperature, bath))
                    k_projected = k_high * (reduced / (1.0 + reduced)) * physical_factor
                    writer.writerow([
                        index,
                        record["equation"],
                        "%.17g" % temperature,
                        "%.17g" % pressure,
                        ";".join(
                            "%s:%.17g" % (gas.species_names[p], mole_fractions[p])
                            for p in range(len(gas.species_names))
                            if mole_fractions[p] > 0.0
                        ),
                        "%.17g" % bath,
                        "%.17g" % cantera_bath,
                        "%.17g" % reduced,
                        "%.17g" % factor,
                        "%.17g" % k_forward,
                        "%.17g" % k_projected,
                    ])
    write_if_changed(CSV_OUT, csv_buffer.getvalue())

    print("generated %s" % RS_OUT)
    print("  sha256 %s" % sha256_of(RS_OUT))
    print("generated %s" % CSV_OUT)
    print("  sha256 %s" % sha256_of(CSV_OUT))
    print("species %d, reactions %d, falloff reactions %d"
          % (len(species_names), len(reactions), len(falloff_indices)))
    print("mechanism sha256 %s" % sha256_of(yaml_path))


if __name__ == "__main__":
    main(sys.argv)
