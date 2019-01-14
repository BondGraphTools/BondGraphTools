Installing BondGraphTools on Mac OSX
====================================

A working installation of both Python and Julia is required to use all of BondGraphTools features.
In this guide we walk step-by-step through the reccomended method of installation.

Step 1: Instal Homebrew and Python
----------------------------------
First we must install python 3. We recommend doing this via Homebrew_.

1. Open a terminal instance.
2. Paste the following command into the terminal instance::

    /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

3. Once Homebrew has finished, install python by typing the following command into the terminal::

    brew install python

Python 3.7 is now installed.

Step 2: Install Julia
---------------------
We must now install Julia 0.6.

1. Download Julia 0.6.4 from https://julialang.org/downloads/oldreleases.html
2. Once the download is complete open the file. It will contain a link to th Applications folder, and the Julia-0.6
   package. Drag the julia package into the applications folder. Once this is done, julia has been installed and you
   can remove the downloaded file.
3. Add julia to your profile path by typing the following command into the terminal::

    echo 'export PATH=/Applications/Julia-0.6.app/Contents/Resources/julia/bin/:$PATH' >> ~/.bash_profile

Step 3: Install BondGraphTools
------------------------------

1. Install BondGraphTools from PyPI by typing the following command into the terminal::

    pip3 install BondGraphTools

2. Install the recommended additional packages (jupyter) by typing the following command into the terminal::

    pip3 install jupyter

`BondGraphTools` is now installed.

Step 4: Julia dependencies.
---------------------------
`BondGraphTools` will automatically download and install Julia dependencies the first time
a user attempts to run a simulation.

This can be manually triggered via the following command::

    from bondgraphtools.config import config
    config.install_dependencies()

.. _Homebrew: https://brew.sh/



