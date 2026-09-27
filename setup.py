from pathlib import Path
from shutil import copy2, copytree

from setuptools import setup
from setuptools.command.build_py import build_py as setuptools_build_py


class build_py(setuptools_build_py):
    def run(self):
        super().run()
        source = Path(__file__).parent / "datasets"
        collection = Path(self.build_lib) / "omicsbench" / "_collection"
        destination = collection / "datasets"
        copytree(source, destination, dirs_exist_ok=True)
        collection.mkdir(parents=True, exist_ok=True)
        for name in ("CITATION.cff", "LICENSE", "LICENSE-METADATA", "NOTICE"):
            copy2(Path(__file__).parent / name, collection / name)


setup(cmdclass={"build_py": build_py})
