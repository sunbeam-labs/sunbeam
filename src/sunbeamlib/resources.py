import os

MIN_MEM_MB = os.getenv("SUNBEAM_MIN_MEM_MB", 1000)
MIN_RUNTIME = os.getenv("SUNBEAM_MIN_RUNTIME", 20)

# def adapter_removal_mem_mb(attempts):


def sample_intake_mem_mb(wildcards, attempt):
    return MIN_MEM_MB


def sample_intake_runtime(wildcards, attempt):
    return MIN_RUNTIME
