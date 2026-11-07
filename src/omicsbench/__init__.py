"""Discovery, validation and local analysis of compact omics benchmarks."""
__version__ = "11.0.0"

__all__ = ["Benchmark", "Collection", "__version__"]


def __getattr__(name):
    if name == "Benchmark":
        from .client import Benchmark
        return Benchmark
    if name == "Collection":
        from .client import Collection
        return Collection
    raise AttributeError(name)
