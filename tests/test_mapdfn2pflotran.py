import os
import shutil
import subprocess
from pathlib import Path

import mapdfn2pflotran

script_dir = os.path.dirname(os.path.realpath(__file__))


class workdir:
    """
    Context manager for creation and usage of a workspace dir.

    name: the workspace directory
    inputs: list of files and directories to copy into the workspaceand
        TODO: fine a sort of robust ad portable reference
    clean: if true the workspace would be deleted at the end of the context manager.
    TODO: clean_before / clean_after
    TODO: File constructor taking current workdir environment, openning virtually copied files.
    TODO: Workdir would not perform change of working dir, but provides system interface for: subprocess, file openning
          only perform CD just before executing a subprocess also interacts with concept of an executable.
    portable reference and with lazy evaluation. Optional true copy possible.
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
    for f in dir.glob('*.h5'):
        all = all and eq_h5(dir / ".." / "ref" / f.name, f)
    return all


def test_dfn():
    """
    TODO:
    - comparison using HDF5
    """

    with workdir(Path(script_dir) / "test_data" / "dfn_2" / "output"):
        mapdfn2pflotran.main(Path(".."))
        assert(check_output_dir(Path(".")))