"""Command-line integration checks for parameter metadata generation."""

import json
import os
import re
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

@pytest.fixture(scope="module")
def template_codegen():
    executable = Path(
        os.environ.get(
            "CADET_TEMPLATE_CODEGEN", ROOT / "build/src/build-tools/templateCodeGen"
        )
    ).resolve()
    if not executable.is_file():
        pytest.fail("Build templateCodeGen or set CADET_TEMPLATE_CODEGEN to its path")
    return executable


def declared_block(source):
    """Return the single codegen block of a source file as parsed JSON."""
    blocks = re.findall(r"/\*<codegen>(.*?)</codegen>\*/", source.read_text(), re.DOTALL)
    assert len(blocks) == 1
    return json.loads(blocks[0])


def test_cpp_generation_preserves_surrounding_source(template_codegen, tmp_path):
    source = tmp_path / "model.cpp"
    template = tmp_path / "template.in"
    output = tmp_path / "generated.cpp"
    source.write_text('before\n/*<codegen>{"name":"Model"}</codegen>*/\nafter\n')
    template.write_text("ignored\n/* <codegentemplate> */class {{ name }} {};")

    subprocess.run([template_codegen, template, source, output], check=True)

    generated = output.read_text()
    assert generated.startswith(
        "/***********************************************************"
    )
    assert generated.endswith("before\nclass Model {};\nafter\n")


@pytest.mark.parametrize(
    ("block", "error"),
    [
        ({"inputs": [], "name": "Handler"}, "cannot be combined with legacy field name"),
        ({"inputs": [], "parameters": []}, "cannot be combined with legacy field parameters"),
        ({"inputs": [], "implementation": "Handler"}, "implementation must be an object"),
    ],
)
def test_inputs_reject_legacy_block_fields(template_codegen, tmp_path, block, error):
    source = tmp_path / "model.cpp"
    template = tmp_path / "classes.in"
    output = tmp_path / "classes.txt"
    source.write_text(f"/*<codegen>{json.dumps(block)}</codegen>*/")
    template.write_text("unused")
    output.write_text("previous output")

    result = subprocess.run(
        [template_codegen, template, source, output],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert error in result.stderr
    assert output.read_text() == "previous output"


def _input(conf_name, descriptor="ScalarComponentDependentParameter", **fields):
    entry = {
        "confName": conf_name,
        "dtype": "double",
        "dims": ["component"] if descriptor == "ScalarComponentDependentParameter" else [],
        "sensitivity": True,
        "externalFunction": True,
        "sectionDependent": False,
    }
    if descriptor is not None:
        entry["implementation"] = {"descriptor": descriptor, "varName": conf_name.lower()}
    entry.update(fields)
    return entry


def test_inputs_are_routed_by_external_function_flag(template_codegen, tmp_path):
    source = tmp_path / "model.cpp"
    template = tmp_path / "lists.in"
    output = tmp_path / "lists.txt"
    block = {
        "implementation": {"handler": "Handler", "externalHandler": "ExternalHandler"},
        "inputs": [
            _input("A"),
            _input("B", descriptor="ScalarParameter", externalFunction=False),
            _input("C", descriptor=None),
            _input("D", descriptor="ScalarParameter"),
        ],
    }
    source.write_text(f"/*<codegen>{json.dumps(block)}</codegen>*/")
    template.write_text(
        "{% for p in parameters %}{{ p/confName }}:{{ p/type }} {% endfor %}| "
        "{% for p in constantParameters %}{{ p/confName }}:{{ p/type }}{% endfor %}"
    )

    subprocess.run([template_codegen, template, source, output], check=True)

    # C has no descriptor: read by handwritten code, so no generated reader
    assert output.read_text().endswith(
        "A:ScalarComponentDependentParameter D:ScalarParameter | B:ScalarParameter"
    )


@pytest.mark.parametrize(
    ("entry", "error"),
    [
        (_input("A", sensitivity=False), "contradicts descriptor"),
        (_input("A", sectionDependent=True), "contradicts descriptor"),
        (_input("A", dims=[]), "contradicts descriptor"),
        (_input("A", dtype="int"), "contradicts descriptor"),
        (_input("A", descriptor="ScalarParameter", dims=["component"]), "contradicts descriptor"),
        (_input("A", descriptor="UnknownParameter"), "no contract for descriptor"),
        (
            {k: v for k, v in _input("A").items() if k != "sensitivity"},
            "needs boolean field sensitivity",
        ),
        (_input("A", externalFunction="yes"), "needs boolean field externalFunction"),
        (
            _input("A", descriptor=None, type="ScalarParameter"),
            "field type belongs under implementation",
        ),
        (_input("A", varName="a"), "field varName belongs under implementation"),
    ],
)
def test_input_contradicting_its_reader_fails(template_codegen, tmp_path, entry, error):
    source = tmp_path / "model.cpp"
    template = tmp_path / "classes.in"
    output = tmp_path / "classes.txt"
    source.write_text(f"/*<codegen>{json.dumps({'inputs': [entry]})}</codegen>*/")
    template.write_text("unused")
    output.write_text("previous output")

    result = subprocess.run(
        [template_codegen, template, source, output],
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert error in result.stderr
    assert output.read_text() == "previous output"


@pytest.mark.parametrize(
    ("source_name", "model", "expected_names"),
    [
        ("LangmuirBinding", "MULTI_COMPONENT_LANGMUIR", "MCL_KA MCL_KD MCL_QMAX"),
        (
            "StericMassActionBinding",
            "STERIC_MASS_ACTION",
            "SMA_KA SMA_KD SMA_NU SMA_SIGMA SMA_LAMBDA SMA_REFC0 SMA_REFQ",
        ),
    ],
)
def test_declarations_cover_each_model_specific_input_once(
    template_codegen, tmp_path, source_name, model, expected_names
):
    declared = declared_block(ROOT / f"src/libcadet/model/binding/{source_name}.cpp")

    assert declared["model"] == model
    entries = declared["inputs"]
    names = [parameter["confName"] for parameter in entries]
    assert len(names) == len(set(names))
    # Core-only fields appear only under "implementation".
    core_only = {"name", "externalName", "handler", "externalHandler", "type", "descriptor", "varName"}
    assert not core_only & set(declared)
    for parameter in entries:
        assert not core_only & set(parameter)
    assert set(names) == set(expected_names.split())

    # Documentation-only properties must not change generated readers.
    reader_fields = {
        "implementation",
        "type",
        "descriptor",
        "varName",
        "confName",
        "objName",
        "confPrefix",
        "skipConfig",
        "default",
        # Checked against the descriptor contract or used to route inputs
        "dtype",
        "dims",
        "sensitivity",
        "externalFunction",
        "sectionDependent",
    }
    readers_only = dict(declared)
    readers_only["inputs"] = [
        {key: value for key, value in entry.items() if key in reader_fields}
        for entry in entries
    ]
    temporary_source = tmp_path / "model.cpp"
    generated = tmp_path / "generated.cpp"
    cpp_template = ROOT / "src/libcadet/model/binding/ExternalFunctionTemplate.cpp"
    baseline = None
    for block in (
        declared,
        readers_only,
    ):
        temporary_source.write_text(f"/*<codegen>{json.dumps(block)}</codegen>*/")
        subprocess.run(
            [template_codegen, cpp_template, temporary_source, generated], check=True
        )
        if baseline is None:
            baseline = generated.read_bytes()
        else:
            assert generated.read_bytes() == baseline


def test_declarations_state_ranges_and_optional_inputs():
    declarations = {}
    for name in ("Langmuir", "StericMassAction"):
        declared = declared_block(ROOT / f"src/libcadet/model/binding/{name}Binding.cpp")
        declarations[name] = {entry["confName"]: entry for entry in declared["inputs"]}

    langmuir = declarations["Langmuir"]
    assert langmuir["MCL_KA"]["range"] == {"ge": 0}
    assert langmuir["MCL_QMAX"]["range"] == {"gt": 0}
    sma = declarations["StericMassAction"]
    for name, unit in (("SMA_REFC0", "mol m_MP^-3"), ("SMA_REFQ", "mol m_SP^-3")):
        assert sma[name]["required"] is False
        assert sma[name]["default"] == 1.0
        assert sma[name]["dims"] == []
        assert sma[name]["unit"] == unit


@pytest.mark.parametrize(
    ("source_name", "expected_names"),
    [
        ("model/binding/BindingModelBase.cpp", ["IS_KINETIC"]),
        ("SimulatorImpl.cpp", ["RELTOL", "RELTOL_SENS", "INIT_STEP_SIZE"]),
    ],
)
def test_handwritten_inputs_pass_input_checks(
    template_codegen, tmp_path, source_name, expected_names
):
    # These blocks are not compiled through the generator, so run them through the checks here
    declared = declared_block(ROOT / "src/libcadet" / source_name)
    block = tmp_path / "model.cpp"
    template = tmp_path / "lists.in"
    output = tmp_path / "lists.txt"
    block.write_text(f"/*<codegen>{json.dumps(declared)}</codegen>*/")
    template.write_text("{{ length(parameters) }}")

    subprocess.run([template_codegen, template, block, output], check=True)

    # Read by handwritten code: described, but no generated reader
    assert output.read_text().endswith("0")
    assert [entry["confName"] for entry in declared["inputs"]] == expected_names
    for entry in declared["inputs"]:
        assert "implementation" not in entry
