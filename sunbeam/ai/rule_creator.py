from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence


@dataclass
class RuleCreationResult:
    """Result of an AI drafting request.

    Attributes:
        rules_text: The generated Snakemake rules in `.smk` syntax.
        written_to: Path where rules were saved (if any).
        model_name: The underlying LLM identifier, if available.
    """

    rules_text: str
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


def _build_system_prompt() -> str:
    return (
        "You are an expert Snakemake engineer.\n"
        "Generate valid Snakemake `.smk` rules only.\n"
        "Constraints:\n"
        "- Output ONLY rules and required Python blocks for Snakemake.\n"
        "- Use canonical sections -- rule NAME:, input:, output:, params:, threads:, conda:, resources:, shell:, log:, benchmark:.\n"
        "- Do not include prose, markdown, or triple backticks. However, each rule should include a docstring with the rule's purpose and any additional context.\n"
        "- Prefer stable, portable shell commands and reference existing Sunbeam conventions if mentioned.\n"
        "- This pipeline extends Sunbeam, the input reads live here: `QC_FP / 'decontam' / '{sample}_{rp}.fastq.gz'`. Default to paired end if there's ambiguity.\n"
        "Some common Sunbeam conventions:\n"
        "- Other extensions may use similar rules and rule names; avoid collisions by prefixing each rule name with the extension name (e.g., `myext_rule_name`).\n"
        "- Use `log: LOG_FP / 'rule_name_{sample}.log'` to capture standard out/err for each sample. Expand over wildcards as necessary to match the output. In the shell command, try to include everything in a subshell and redirect everything that doesn't go into outputs to the log file.\n"
        "- Use `benchmark: BENCHMARK_FP / 'rule_name_{sample}.tsv'` to capture resource usage for each sample.\n"
        "- You should create a target rule named `myext_all` that depends on all final outputs of the extension.\n"
        "Some important Sunbeam variables available in rules:\n"
        "`Cfg` is a configuration dictionary holding the content of `sunbeam_config.yml`. You will probably not use this nor make your own config. If there are obvious configurable parameters for a rule, define them in code at the top of the file.\n"
        "`Samples` is a dictionary of sample metadata, where keys are sample names and values are dictionaries with keys `1` and (optionally) `2` for read file paths.\n"
        "`Pairs` is a list. If the run is paired end, it is ['1', '2']. If single end, it is ['1'].\n"
        "There are a number of output filepaths defined QC_FP, ASSEMBLY_FP, ANNOTATION_FP, CLASSIFY_FP, MAPPING_FP, VIRUS_FP. All outputs should live in one of these directories. If none of these fit the theme of the new extension, define your own at the top of the file with `SOMETHING_FP = output_subdir(Cfg, 'something')`.\n"
    )


def _default_model_name() -> str:
    # Keep generic to avoid hard coupling to a specific vendor.
    # Users can override with their preferred model string.
    return "gpt-4o-mini"


def _simple_rules_validator(text: str) -> None:
    """Basic safety/structure checks for generated rules.

    Raises ValueError if any quick structural checks fail. This is not
    a full parser; it enforces minimal sanity before writing to disk.
    """

    if not text.strip():
        raise ValueError("No content generated for rules.")

    lowered = text.lower()
    if "```" in text:
        raise ValueError("Output contains markdown fences; expected raw `.smk` rules.")
    if "rule " not in lowered:
        raise ValueError("No `rule` blocks detected in generated content.")
    if any(k in lowered for k in ["<todo>", "[fill", "<fill", "tbd"]):
        raise ValueError("Placeholder text found; refusing to write incomplete rules.")


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

    system = _build_system_prompt()

    # Assemble context payloads.
    context_blobs: list[str] = []
    if context_files:
        for p in context_files:
            try:
                blob = Path(p).read_text()
            except Exception:
                continue
            context_blobs.append(f"FILE: {Path(p).name}\n{blob}")

    context_text = ("\n\n".join(context_blobs)).strip()
    user_content = (
        prompt if not context_text else (f"Context:\n{context_text}\n\nTask:\n{prompt}")
    )

    model_name = model or _default_model_name()
    llm = ChatOpenAI(model=model_name, api_key=api_key)

    resp = llm.invoke(
        [
            SystemMessage(content=system),
            HumanMessage(content=user_content),
        ]
    )

    rules_text = (
        (resp.content or "").strip() if hasattr(resp, "content") else str(resp).strip()
    )

    _simple_rules_validator(rules_text)

    written_path: Optional[Path] = None
    if write_to:
        out_path = Path(write_to)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(rules_text)
        written_path = out_path

    return RuleCreationResult(
        rules_text=rules_text,
        written_to=written_path,
        model_name=model_name,
    )
