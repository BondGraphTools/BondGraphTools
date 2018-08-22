import pathlib
import json


def as_str(path):
    if isinstance(path, str):
        return path
    else:
        return path.as_posix()


class Config:

    base = pathlib.Path.home().absolute() / '.BondGraphTools'
    file = base / 'config.json'

    def __init__(self):

        self._julia = None
        if not self.base.exists():
            self.base.mkdir()
            self.file.touch()

    @property
    def julia(self):
        if not self._julia:
            self._julia = Julia_Config(self, 'julia', self.base / 'julia_packages')
        return self._julia

    def load(self):
        with open(str(self.file), 'r') as f:
            config_dict = json.load(f)

    def save(self):
        config_dict = {
            "julia": self.julia.to_dict()
        }


class Julia_Config:
    def __init__(self, parent, executable, packages):
        self.__config = parent
        self._executable = executable
        self._packages = as_str(packages)

    def to_dict(self):
        return {
            "executable": self._executable,
            "packages": self._packages
        }

    @property
    def executable(self):
        return self._executable

    @property
    def packages(self):
        return self._packages

    @packages.setter
    def packages(self, val):
        self._packages = as_str(val)
        self.__config.save()

    @executable.setter
    def executable(self, val):
        self._executable = as_str(val)
        self.__config.save()


config = Config()
if Config.file.exists():
    config.load()

