import os, sys, platform, shutil,pathlib, stat
import requests
from subprocess import Popen, PIPE, run, call

_S3HOST = 'https://julialang-s3.julialang.org/bin'

julia_path = pathlib.Path(__file__).parent.absolute() / \
             "BondGraphTools" / "julia"

julia_exec = julia_path / "bin" / "julia"

base = pathlib.Path(__file__).parent.absolute()
mask = stat.S_IXUSR | stat.S_IXGRP|stat.S_IXOTH

def web_copy(url, local_file):
    response = requests.get(url, stream=True)
    with open(local_file, 'wb') as filestream:
        response.raw.decode_content = True
        shutil.copyfileobj(response.raw, filestream)


def _install_osx():
    temp_file = base / "file.dmg"
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
    os.chmod(julia_exec, mask)


def _install_linux():
    temp_file = base / "file.tar.gz"
    if platform.architecture()[0] is '64bit':
        julia_file = _S3HOST + '/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz'
    else:
        julia_file = _S3HOST + '/linux/x86/0.6/julia-0.6.4-linux-i686.tar.gz'

    extract_path = base / "julia-9d11f62bcb"

    if not temp_file.exists():
        web_copy(julia_file, temp_file)

    shutil.unpack_archive(temp_file, format='gztar')
    shutil.move(extract_path, julia_path)
    os.chmod(julia_exec, mask)


def _install_win32():
    julia_file = _S3HOST + '/winnt/x86/0.6/julia-0.6.4-win32.exe'


def _install_win64():
    julia_file = _S3HOST + '/winnt/x64/0.6/julia-0.6.4-win64.exe'


def install_julia_deps():
    run([julia_exec, 'setup.jl'])


def install_python_deps():
    call([sys.executable, '-m', 'pip', 'install', 'diffeqpy==0.4'])


def check_julia():
    # Check to see if Julia is already in the path
    if not julia_path.exists():
        return False
    try:
        p = Popen([julia_exec, "-v"],
                  stdout=PIPE, stderr=None)
        out, err = p.communicate()
        text = out.decode('utf8').strip()
        if text.startswith("julia version 0.6"):
            return True
        else:
            return False
    except FileNotFoundError:
        return False


if not check_julia():
    if sys.platform == 'darwin':
        _install_osx()
    elif sys.platform == 'linux':
        _install_linux()
    elif sys.platform == 'windows':
        raise NotImplementedError(
            "Automatic Installation for Windows is not supported. "
            "Please manually install and add julia to the environment path."
        )

install_julia_deps()
install_python_deps()









