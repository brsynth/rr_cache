"""
Pytest rewrite of the original unittest-based rrCache tests.

This version uses native pytest parametrization so each cspace/case is
collected and reported as a separate test item instead of being grouped under
unittest subTest blocks.
"""

from __future__ import annotations

from json import load as json_load
from logging import Logger
from os import path as os_path, remove as os_rm

import pytest

from brs_utils import check_file_size, create_logger, extract_gz
from rr_cache import rrCache

HERE = os_path.dirname(os_path.abspath(__file__))
DATA_PATH = os_path.join(HERE, "data")
# CSPACES = ["mnx3.1", "mnx4.4", "rr2026"]
CSPACES = ["rr2026"]
DATABASES = ["metanetx", "rhea"]
DATASETS = ["compounds", "metrics", "retrorules", "reactions"]


def _load_json(filepath: str):
    with open(filepath, "r") as handle:
        return json_load(handle)


# Load reference data at collection time so pytest can parametrize individual
# cases and report them separately.
REFERENCE_DATA = {
    cspace: {
        dataset: _load_json(os_path.join(DATA_PATH, f"{dataset}_{cspace}.json"))
        for dataset in DATASETS
    }
    for cspace in CSPACES
}

ALL_ATTR_CASES = [
    pytest.param(cspace, attr, metric["length"], id=f"{cspace}-{attr}")
    for cspace in CSPACES
    for attr, metric in REFERENCE_DATA[cspace]["metrics"].items()
]

COMPOUND_CASES = [
    pytest.param(
        cspace,
        next(iter(REFERENCE_DATA[cspace]["compounds"])),
        id=f"{cspace}-compound",
    )
    for cspace in CSPACES
]

REACTION_CASES = [
    pytest.param(
        cspace,
        next(iter(REFERENCE_DATA[cspace]["reactions"])),
        id=f"{cspace}-reaction",
    )
    for cspace in CSPACES
]

RULE_CASES = [
    pytest.param(
        cspace,
        next(iter(REFERENCE_DATA[cspace]["retrorules"])),
        id=f"{cspace}-rule",
    )
    for cspace in CSPACES
]


@pytest.fixture(scope="session")
def logger() -> Logger:
    return create_logger(__name__, "ERROR")


@pytest.fixture(scope="session")
def reference_data() -> dict[str, dict[str, dict]]:
    return REFERENCE_DATA


@pytest.fixture(scope="session")
def caches(logger: Logger) -> dict[str, rrCache]:
    return {
        cspace: rrCache(
            cspace=cspace, databases=DATABASES, interactive=False, logger=logger
        )
        for cspace in CSPACES
    }


@pytest.mark.parametrize("cspace, attr, expected_length", ALL_ATTR_CASES)
def test_all_attr(caches, cspace: str, attr: str, expected_length: int):
    """Load all rrCache attributes and validate their expected lengths."""
    assert len(caches[cspace].get(attr)) == expected_length


@pytest.mark.skip(reason="Too long, not in standard tests")
@pytest.mark.parametrize("cspace", CSPACES, ids=CSPACES)
def test_generate_cache(cspace: str, logger: Logger, reference_data):
    """Generate rrCache files and validate the extracted file sizes."""
    outdir = f"cache-{cspace}"
    rrCache.generate_cache(outdir, interactive=False, logger=logger)
    metrics = reference_data[cspace]["metrics"]

    for name, meta in metrics.items():
        filepath = os_path.join(outdir, f"{name}.json.gz")
        outfile = extract_gz(filepath, outdir)
        try:
            assert check_file_size(outfile, meta["file_size"], logger)
        finally:
            if os_path.exists(outfile):
                os_rm(outfile)


@pytest.mark.parametrize("cspace, cmpd_id", COMPOUND_CASES)
def test_get_compound(caches, reference_data, cspace: str, cmpd_id: str):
    compounds = reference_data[cspace]["compounds"]
    compound = caches[cspace].get_compound(cmpd_id)
    compound["xref"] = caches[cspace].get_compound_xref(cmpd_id)
    print(compound)
    print()
    print(compounds[cmpd_id])
    assert compound == compounds[cmpd_id]


@pytest.mark.parametrize("cspace, cmpd_id", COMPOUND_CASES)
def test_get_list_of_compounds(caches, reference_data, cspace: str, cmpd_id: str):
    compound_ids = caches[cspace].get_list_of_compounds()
    expected_length = reference_data[cspace]["metrics"]["cid_strc"]["length"]

    assert cmpd_id in compound_ids
    assert len(compound_ids) == expected_length


@pytest.mark.parametrize("cspace, rxn_id", REACTION_CASES)
def test_get_reaction(caches, reference_data, cspace: str, rxn_id: str):
    reactions = reference_data[cspace]["reactions"]
    assert caches[cspace].get_reaction(rxn_id) == reactions[rxn_id]


@pytest.mark.parametrize("cspace, rxn_id", REACTION_CASES)
def test_get_list_of_reactions(caches, reference_data, cspace: str, rxn_id: str):
    reaction_ids = caches[cspace].get_list_of_reactions()
    expected_length = reference_data[cspace]["metrics"]["template_reactions"]["length"]

    assert rxn_id in reaction_ids
    assert len(reaction_ids) == expected_length


@pytest.mark.parametrize("cspace, rule_id", RULE_CASES)
def test_get_reaction_rule(caches, reference_data, cspace: str, rule_id: str):
    retrorules = reference_data[cspace]["retrorules"]
    assert caches[cspace].get_reaction_rule(rule_id) == retrorules[rule_id]


@pytest.mark.parametrize("cspace, rule_id", RULE_CASES)
def test_get_list_of_reaction_rules(caches, reference_data, cspace: str, rule_id: str):
    rule_ids = caches[cspace].get_list_of_reaction_rules()
    expected_length = reference_data[cspace]["metrics"]["rr_reactions"]["length"]

    assert rule_id in rule_ids
    assert len(rule_ids) == expected_length
