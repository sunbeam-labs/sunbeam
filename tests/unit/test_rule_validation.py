import pytest

from sunbeam.ai.rule_creator import _snakemake_rules_validator


def test_validator_accepts_valid_rules():
    rules = (
        "rule test:\n"
        "    output: 'out.txt'\n"
        "    shell:\n"
        "        'echo hi > {output}'\n"
    )
    _snakemake_rules_validator(rules)


def test_validator_rejects_invalid_rules():
    bad_rules = (
        "rule test\n"  # missing colon after rule name
        "    output: 'out.txt'\n"
        "    shell:\n"
        "        'echo hi > {output}'\n"
    )
    with pytest.raises(ValueError):
        _snakemake_rules_validator(bad_rules)
