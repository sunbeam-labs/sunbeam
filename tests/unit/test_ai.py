from sunbeam.ai.workflow import WorkflowPlanner
from langchain_community.llms.fake import FakeListLLM


def test_planner_returns_plan():
    llm = FakeListLLM(responses=["plan"])
    planner = WorkflowPlanner(llm=llm)
    plan = planner.plan("reads", "assembly", "QC, assembly")
    assert plan == "plan"
