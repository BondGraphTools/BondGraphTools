Installing BondGraphTools on Mac OSX
====================================

A working installation of both Python and Julia is required to use all of BondGraphTools features.
In this guide we walk step-by-step through the recommended method of installation.

Step 1: Install Homebrew and Python
----------------------------------
First we must install python 3. We recommend doing this via Homebrew_.

1. Open a terminal instance.
2. Paste the following command into the terminal instance::

    $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

3. Once Homebrew has finished, install python by typing the following command into the terminal::

    $ brew install python

Python 3.7 is now installed.

Step 2: Install Sundials
------------------------
Next, install Sundials_ and other scikit.odes dependencies
(see https://scikits-odes.readthedocs.io/en/latest/installation.html#requirements-before-install for details).

Step 3: Install BondGraphTools
------------------------------

1. Install BondGraphTools from PyPI by typing the following command into the terminal::

    $ pip3 install BondGraphTools

2. Install the recommended additional packages (jupyter) by typing the following command into the terminal::

    $ pip3 install jupyter

`BondGraphTools` is now installed.


.. _Homebrew: https://brew.sh/
.. _Sundials: https://computing.llnl.gov/projects/sundials


