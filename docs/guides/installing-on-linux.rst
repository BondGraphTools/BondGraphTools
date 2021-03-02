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

Step 2: Install Sundials
---------------------


Step 3: Install BondGraphTools
------------------------------
Use python to install BondGraphTools with::

    $ pip3 install BondGraphTools

or alternatively::

    $ python3.6 -m pip install BondGraphTools

It is optional but highly recommended to install the Jupyter packages via::

    $ pip3 install Jupyter
