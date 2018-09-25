import pathlib
import json
import sys
import os
import logging
from subprocess import Popen, PIPE, run
logger = logging.getLogger(__name__)

VERSION = "0.3.3"


def as_str(path):
    if isinstance(path, str):
        return path
    else:
        return str(path.as_posix())


def check_julia(command):
    # Check to see if Julia is already in the path

    p = Popen([command, "-v"], stdout=PIPE, stderr=None)

    out, err = p.communicate()
    text = out.decode('utf8').strip()
    if text.startswith("julia version 0.6"):
        return True
    else:
        return False


class Config:
    _which = 'which'
    base = pathlib.Path.home().absolute() / '.BondGraphTools'
    file = base / 'config.json'

    def __init__(self, julia_executable=None, python_executable=None, **kwargs):

        self._julia = None
        self._de = None
        new_config = False
        rebuild = False

        if not julia_executable:
            self.julia_executable = self.find_julia()
            new_config = True
        else:
            self.julia_executable = pathlib.Path(julia_executable)

        if python_executable and python_executable is not sys.executable:
            rebuild=True
        elif python_executable:
            self.python_executable = pathlib.Path(python_executable)
        else:
            self.python_executable = pathlib.Path(sys.executable).resolve()

        if new_config:
            if not self.base.exists():
                self.base.mkdir()
                self.file.touch()
            self.install_dependencies(rebuild)
            self.save()

    def find_julia(self):
        if check_julia('julia'):
            p = Popen([self._which, 'julia'], stdout=PIPE, stderr=None)
            out, err = p.communicate()
            julia_dir = out.decode('utf8').strip()
            return pathlib.Path(julia_dir).resolve()
        else:
            raise NotImplementedError("Could not find Julia 0.6.4: please install and add to path")

    def find_conda(self):
        p = Popen([self._which, 'conda'], stdout=PIPE, stderr=None)
        out, err = p.communicate()
        conda_dir = out.decode('utf8').strip()
        if not conda_dir:
            return None
        else:
            return pathlib.Path(conda_dir).resolve()

    def install_dependencies(self, rebuild=False):
        # we assume julia and python are already in the path
        logger.debug('Installing Julia dependencies; this may take some time')
        env = os.environ
        env.update({"PYTHON": as_str(self.python_executable),
                    "JULIA": as_str(self.julia_executable)})
        conda = self.find_conda()
        if conda:
            env.update({"CONDA": as_str(conda)})

        julia_code = [
            "Pkg.init()\n",
            "Pkg.update()\n",
            """Pkg.add("PyCall")\n""",
            """Pkg.add("DifferentialEquations")\n"""]

        if rebuild:
            julia_code +=[
                """Pkg.build("PyCall")""",
                "using PyCall\n",
                "using DifferentialEquations\n"]

        temp = self.base / 'deps.jl'
        temp.touch()
        with open(as_str(temp), 'w') as fs:
            fs.writelines(julia_code)
        run([as_str(self.julia_executable), as_str(temp)], env=env)
        os.remove(temp)
        logger.debug("Complete")

    @property
    def julia(self):
        if not self._julia:
            self.start_julia()

        return self._julia

    def start_julia(self):

        import diffeqpy
        from diffeqpy import de

        self._julia = diffeqpy.setup()
        self._de = de

    @property
    def de(self):
        if not self._julia: self.start_julia()
        return self._de


    @staticmethod
    def load():
        try:
            with open(as_str(Config.file), 'r') as f:
                kwargs = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as ex:
            kwargs = {}

        if sys.platform.startswith('win'):
            return WinConfig(**kwargs)
        else:
            return Config(**kwargs)

    def save(self):

        config_dict = {
            'julia_executable': as_str(self.julia_executable),
            'python_executable': as_str(self.python_executable),
            'version': VERSION
        }
        with open(as_str(Config.file), 'w') as f:
            json.dump(config_dict, f)


class WinConfig(Config):
    _which = 'where'


config = Config.load()



