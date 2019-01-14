Troubleshooting
===============

Julia Errors
------------
If you encounter julia errors, you may need to rebuild pycall, or re-install your julia libraries.

1. To rebuild pycall and related libraries enter the following commands in your python interpreter::

        from bondgraphtools.config import config
        config.install_dependencies(rebuilt=True)

2. If that fails, you may need to reset your julia library. **Warning** this will delete any
   additional libraries you may have installed. This can be done by removing the directory '~/.julia/lib'
   in Linux/MacOSX, or 'C:\\users\\<username>\\.julia\\lib' in windows.



