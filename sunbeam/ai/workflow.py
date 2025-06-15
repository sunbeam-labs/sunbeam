"""Utilities to generate and run AI-driven workflows."""

from __future__ import annotations

from dataclasses import dataclass

from langchain.prompts import ChatPromptTemplate
from langchain.schema.runnable import Runnable
from langchain.schema.output_parser import StrOutputParser

try:
    from langchain_openai import ChatOpenAI
except ImportError:  # pragma: no cover - optional dependency
    ChatOpenAI = None


def default_llm():
    if ChatOpenAI is None:
        raise ImportError(
            "langchain_openai is required for the default planner. Install"
            " the 'ai' extra and set OPENAI_API_KEY"
        )
    return ChatOpenAI(temperature=0)


@dataclass
class WorkflowPlanner:
    """Plan bioinformatics workflows using an LLM."""

    llm: Runnable | None = None

    def __post_init__(self):
        if self.llm is None:
            self.llm = default_llm()
        self.chain = self._build_chain(self.llm)

    @staticmethod
    def _build_chain(llm: Runnable) -> Runnable:
        prompt = ChatPromptTemplate.from_messages(
            [
                (
                    "system",
                    "You are an expert bioinformatics workflow planner. Given a description"
                    " of the data, desired outputs, and allowed steps, produce a concise"
                    " plan with shell commands to run each step.",
                ),
                (
                    "user",
                    "Input: {input_desc}\nDesired output: {output_desc}\nAllowed steps: {steps_hint}",
                ),
            ]
        )
        return prompt | llm | StrOutputParser()

    def plan(self, input_desc: str, output_desc: str, steps_hint: str = "") -> str:
        """Return a text plan for accomplishing the task."""

        return self.chain.invoke(
            {
                "input_desc": input_desc,
                "output_desc": output_desc,
                "steps_hint": steps_hint,
            }
        )
