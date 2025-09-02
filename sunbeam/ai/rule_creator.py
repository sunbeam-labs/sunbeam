from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from sunbeam.ai import CONDA_ENV_GEN_SYSTEM_PROMPT, RULE_GEN_SYSTEM_PROMPT
from typing import Optional, Sequence


@dataclass
class RuleCreationResult:
    """Result of an AI drafting request.

    Attributes:
        rules_text: The generated Snakemake rules in `.smk` syntax.
        env_texts: Mapping of conda environment names to their YAML contents.
        written_to: Path where rules were saved (if any).
        model_name: The underlying LLM identifier, if available.
    """

    rules_text: str
    env_texts: dict[str, str]
    written_to: Optional[Path] = None
    model_name: Optional[str] = None


def _require_ai_dependencies() -> None:
    """Eagerly validate optional AI deps and raise a clear error if missing."""

    try:
        import langchain  # noqa: F401
        import langgraph  # noqa: F401
        import langchain_openai  # noqa: F401
    except Exception as exc:  # pragma: no cover - import error path
        raise RuntimeError(
            "AI features require the optional 'ai' extras. Install with\n"
            "  pip install 'sunbeamlib[ai]'\n"
            "Missing dependency details: " + str(exc)
        )


def _default_model_name() -> str:
    # Keep generic to avoid hard coupling to a specific vendor.
    # Users can override with their preferred model string.
    return "gpt-4o-mini"


def create_rules_from_prompt(
    prompt: str,
    *,
    context_files: Optional[Sequence[Path]] = None,
    write_to: Optional[Path] = None,
    model: Optional[str] = None,
    api_key: Optional[str] = None,
) -> RuleCreationResult:
    """Draft Snakemake rules from a natural-language prompt using an LLM.

    This function requires the optional `ai` extras. It imports AI libraries
    lazily and raises a clear error if they are not installed.

    Args:
        prompt: User's high-level request describing desired workflow behavior.
        context_files: Optional paths whose contents provide additional context
            to the LLM (e.g., existing rules or configs). Read as text.
        write_to: If provided, writes the generated rules to this path.
        model: Optional model identifier; falls back to a generic default.
        api_key: Optional API key for OpenAI-compatible providers. If omitted,
            provider-specific defaults are used (e.g., env vars).

    Returns:
        RuleCreationResult with the generated rules and metadata.
    """

    _require_ai_dependencies()

    # Import after dependency check to keep base imports light.
    from langchain_openai import ChatOpenAI
    from langchain_core.messages import SystemMessage, HumanMessage

    system = RULE_GEN_SYSTEM_PROMPT

    # Assemble context payloads.
    context_blobs: list[str] = []
    if context_files:
        for p in context_files:
            try:
                blob = Path(p).read_text()
            except Exception:
                continue
            context_blobs.append(f"CONTEXT FILE {Path(p).name}:\n{blob}")

    context_text = ("\n\n".join(context_blobs)).strip()
    user_content = (
        prompt if not context_text else (f"Context:\n{context_text}\n\nTask:\n{prompt}")
    )

    model_name = model or _default_model_name()
    llm = ChatOpenAI(model=model_name, api_key=api_key)

    # Run rule generation
    resp = llm.invoke(
        [
            SystemMessage(content=system),
            HumanMessage(content=user_content),
        ]
    )

    rules_text = (
        (resp.content or "").strip() if hasattr(resp, "content") else str(resp).strip()
    )

    # Run envs generation
    envs = get_envs_from_rules(rules_text)
    env_texts = {}

    for env_name, rule in envs.items():
        context = f"Generate a conda environment YAML file for the following Snakemake rule:\n\n{rule}"
        resp = llm.invoke(
            [
                SystemMessage(content=CONDA_ENV_GEN_SYSTEM_PROMPT),
                HumanMessage(content=context),
            ]
        )

        env_texts[env_name] = resp.content if hasattr(resp, "content") else str(resp)

    written_path: Optional[Path] = None
    if write_to:
        # Make ext dir
        out_path = Path(write_to)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Write rules
        out_path.write_text(rules_text)
        # Write envs
        envs_path = out_path.parent / "envs"
        envs_path.mkdir(exist_ok=True)
        for env_name, env_text in env_texts.items():
            (envs_path / f"{env_name}.yaml").write_text(env_text)
        written_path = out_path

    return RuleCreationResult(
        rules_text=rules_text,
        env_texts=env_texts,
        written_to=written_path,
        model_name=model_name,
    )


def get_envs_from_rules(
    rules_text: str,
) -> dict[str, str]:
    """Return rules with conda environments that need generation."""

    rules = []
    current_rule = []

    for line in rules_text.splitlines(keepends=True):
        if line.strip().startswith("rule:"):
            # If we're already collecting a rule, save it
            if current_rule:
                rules.append("".join(current_rule).rstrip("\n"))
                current_rule = []
        current_rule.append(line)

    # Save the last rule (if any)
    if current_rule:
        rules.append("".join(current_rule).rstrip("\n"))

    rules = [r for r in rules if "conda:" in r]

    # Extract env names
    envs = {}
    for rule in rules:
        for line in rule.splitlines():
            if line.strip().startswith("conda:"):
                parts = line.split(":", 1)
                if len(parts) == 2:
                    env_name = parts[1].strip().strip('"').strip("'").strip(".yaml")
                    envs[env_name] = rule
                else:
                    print("Warning: Malformed conda line:", line)

    return envs
