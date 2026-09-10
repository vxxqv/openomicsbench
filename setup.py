from pathlib import Path
from shutil import copytree

from setuptools import setup
from setuptools.command.build_py import build_py as setuptools_build_py


class build_py(setuptools_build_py):
    def run(self):
        super().run()
        source = Path(__file__).parent / "datasets"
        destination = Path(self.build_lib) / "omicsbench" / "_collection" / "datasets"
        copytree(source, destination, dirs_exist_ok=True)


setup(cmdclass={"build_py": build_py})
