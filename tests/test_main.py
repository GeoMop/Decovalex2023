import os
import shutil
import subprocess
from pathlib import Path

import pytest

from decodfn import main

script_dir = os.path.dirname(os.path.realpath(__file__))


class workdir:
    """
    Context manager for creation and usage of a workspace dir.
    """
    def __init__(self, path: str="sandbox", clean=False):
        self.work_dir = Path(path)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self._clean = clean
        self._orig_dir = os.getcwd()

    def __enter__(self):
        os.chdir(self.work_dir)
        return self.work_dir

    def __exit__(self, type, value, traceback):
        os.chdir(self._orig_dir)
        if self._clean:
            shutil.rmtree(self.work_dir)

def eq_h5(a, b, reltol=0.01):
    cmd = ["h5diff", f"-p {reltol}", f"{a}", f"{b}"]
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if len(result.stdout) > 0:
        print("diff cmd: ", ["h5diff", f"-p {reltol}", f"{a}", f"{b}"])
        print(result.stdout)
        return False
    return True

def check_output_dir(dir: Path):
    all = True
    for f in ["isotropic_k", "porosity"]:
        f = Path(f"{f}.h5")
        all = all and eq_h5(dir / ".." / "ref" / f.name, f)
    return all


@pytest.mark.parametrize(
    "case",
    ["dfn_2", "dfn_5", "repo_cube", "surface_cube"]
)
def test_dfn(case):
    with workdir(Path(script_dir) / "test_data" / case / "output"):
        main.main(Path("..").absolute())
        assert(check_output_dir(Path(".")))
