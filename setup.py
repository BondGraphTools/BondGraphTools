import os, sys, platform, shutil, pathlib, stat
import requests
from subprocess import Popen, PIPE, run, call




from BondGraphTools.config import config

_S3HOST = 'https://julialang-s3.julialang.org/bin'

_PYTHON = sys.executable if os.name == 'posix' else sys.executable.replace("""\\""", """/""")

julia_commands = [
    f"""ENV["PYTHON"]="{_PYTHON}"\n""",
    f"""ENV["JULIA_PKGDIR"]="{config.julia.packages}"\n""",
    "Pkg.init()\n",
    "Pkg.update()\n",
    """Pkg.add("PyCall")\n""",
    """Pkg.add("DifferentialEquations")\n""",
    "using PyCall\n"
    "using DifferentialEquations\n"
]


def setup_julia():

    setup_jl = config.base / "setup.jl"

    if not setup_jl.exists():
        setup_jl.touch()
    with open(str(setup_jl), 'w') as filestream:
        filestream.writelines(julia_commands)

    p = run([config.julia.executable, str(setup_jl.absolute())])
    call([sys.executable, '-m', 'pip', 'install', 'diffeqpy==0.4'])


def make_executable(file):
    if os.name == 'posix':
        # Make the resource executable
        mode = (os.stat(file).st_mode | 0o555) & 0o7777
        os.chmod(file, mode)


def web_copy(url, local_file):
    response = requests.get(url, stream=True)
    with open(local_file, 'wb') as filestream:
        response.raw.decode_content = True
        shutil.copyfileobj(response.raw, filestream)


def _install_osx():
    julia_path = config.base / "julia"
    julia_exec = julia_path / "bin" / "julia"

    temp_file = config.base / "file.dmg"
    julia_file = _S3HOST + '/mac/x64/0.6/julia-0.6.4-mac64.dmg'
    mount_point = pathlib.Path('/Volumes/Julia-0.6.4')

    app_base = mount_point / "Julia-0.6.app" / "Contents" \
               / "Resources" / "julia"

    if not temp_file.exists():
        web_copy(julia_file, temp_file)

    p = run(['hdiutil', 'mount', temp_file], stdout=None, stderr=None)
    shutil.copytree(app_base, julia_path)

    p = run(['hdiutil', 'unmount', mount_point], stdout=None, stderr=None)
    os.remove(temp_file)
    make_executable(julia_exec)
    config.julia.executable = str(julia_exec)


def _install_linux():
    julia_path = config.base / "julia"
    julia_exec = julia_path / "bin" / "julia"
    temp_file = config.base / "file.tar.gz"
    if platform.architecture()[0] is '64bit':
        julia_file = _S3HOST + '/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz'
    else:
        julia_file = _S3HOST + '/linux/x86/0.6/julia-0.6.4-linux-i686.tar.gz'

    extract_path = config.base / "julia-9d11f62bcb"

    if not temp_file.exists():
        web_copy(julia_file, temp_file)

    shutil.unpack_archive(temp_file, format='gztar')
    shutil.move(extract_path, julia_path)
    make_executable(julia_exec)
    config.julia.executable = str(julia_exec)


def _install_win32():
    julia_file = _S3HOST + '/winnt/x86/0.6/julia-0.6.4-win32.exe'
    temp_file = config.base /'julia-0.6.4-win32.exe'
    if not temp_file.exists():
        web_copy(julia_file, temp_file)
    p = run(str(temp_file))
    default_locn = pathlib.Path.home() / "AppData" / "Local" / "Julia-0.6.4" / "bin"


def _install_win64():
    julia_file = _S3HOST + '/winnt/x64/0.6/julia-0.6.4-win64.exe'
    temp_file = config.base / 'julia-0.6.4-win64.exe'

    if not temp_file.exists():
        web_copy(julia_file, temp_file)
    p = run(str(temp_file))


def check_julia(command):
    # Check to see if Julia is already in the path

    p = Popen([command, "-v"], stdout=PIPE, stderr=None)

    out, err = p.communicate()
    text = out.decode('utf8').strip()
    if text.startswith("julia version 0.6"):
        return True
    else:
        return False


if sys.platform == 'darwin':
    if not check_julia(config.julia.executable):
        _install_osx()


elif sys.platform == 'linux':
    if not check_julia(config.julia.executable):
        _install_linux()


elif sys.platform == 'win32':
    if not check_julia(config.julia.executable):
        if platform.architecture()[0] == '64bit':
            _install_win64()
        else:
            _install_win32()
else:
    raise NotImplementedError('Operating system not supported')

setup_julia()












