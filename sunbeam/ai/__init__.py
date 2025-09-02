"""AI utilities for Sunbeam (optional).

This package is only available when the optional dependency group
`ai` is installed, e.g. `pip install sunbeamlib[ai]`.

Public API:
- `create_rules_from_prompt` – draft Snakemake rules from a natural-language prompt.
"""

from .prompts import RULE_GEN_SYSTEM_PROMPT, CONDA_ENV_GEN_SYSTEM_PROMPT
from .rule_creator import create_rules_from_prompt, RuleCreationResult

__all__ = [
    "create_rules_from_prompt",
    "RuleCreationResult",
    "RULE_GEN_SYSTEM_PROMPT",
    "CONDA_ENV_GEN_SYSTEM_PROMPT",
]
