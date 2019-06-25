Installing BondGraphTools on Linux
===================================

A working installation of both Python and Julia is required to use all of
BondGraphTools features.
In this guide we walk step-by-step through the recommended method of
installation.

Commands in this section beginning with `$` are entered at the terminal prompt.
Commands beginning with `>>>` are enter at the Python interpreter prompt.

We assume Ubuntu 16.04 or greater. BondGraphTools _should_ work on all
distributions but it is only tested and supported on Ubuntu 16.04 and up.

Step 1: Install Python 3.6 or greater
-------------------------------------

Ubuntu 17.10 or Greater
+++++++++++++++++++++++
Skip this step as Python 3.6 is already installed.

Ubuntu 14.04 Only
+++++++++++++++++
A newer version of the C++ standard library is required. First, download it a
working directory with::

  $ wget -q -O libstdc++6 http://security.ubuntu.com/ubuntu/pool/main/g/gcc-5/libstdc++6_5.4.0-6ubuntu1~16.04.10_amd64.deb

then install using::

    $ sudo dpkg --force-all -i libstdc++6

Ubuntu 14.04 and 16.04
++++++++++++++++++++++
Add the 'deadsnakes' repository to the Ubuntu package manager::

    $ sudo add-apt-repository ppa:deadsnakes/ppa

Ubuntu 14.04 to 16.10
+++++++++++++++++++++
Update the package indices and install Python 3.6::

    $ sudo apt update
    $ sudo apt install python36 python36-pip

Step 2: Install Julia
---------------------
We must now install Julia 0.6.
Our examples follow the Julia platform instructions https://julialang.org/downloads/platform.html#generic-binaries .

1. Download Julia 0.6.4 generic binaries from https://julialang.org/downloads/oldreleases.html .
   For example using curl to save the archive to your downloads::

    $ curl -L https://julialang-s3.julialang.org/bin/linux/x64/0.6/julia-0.6.4-linux-x86_64.tar.gz > ~/Downloads/julia-0.6.4-linux-x86_64.tar.gz
2. Extract the archive `.tar.gz` to the desired directory. For example, here
   we first create the target directory, then extract the tarball::

    $ sudo mkdir /usr/local/opt/julia
    $ sudo tar xf ~/Downloads/julia-0.6.4-linux-x86_64.tar.gz --directory /usr/local/opt/julia --strip-components=1

3. Add julia to the PATH, for example by adding a symbolic link to the local bin
   folder::

    $ sudo ln -s /usr/local/opt/julia/bin/julia /usr/local/bin/julia

4. Verify that this has worked by typing `julia -v', which should result in
   `julia version 0.6.4`.

Step 3: Install BondGraphTools
------------------------------
Use python to install BondGraphTools with::

    $ pip3 install BondGraphTools

or alternatively::

    $ python3.6 -m pip install BondGraphTools

Step 4: Julia dependencies.
---------------------------
`BondGraphTools` will automatically download and install Julia dependencies the first time
a user attempts to run a simulation.

This can be manually triggered via the following commands in the python
interpreter::

    >>> from bondgraphtools.config import config
    >>> config.install_dependencies()

