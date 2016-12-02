## mdbexport python package

This package provides a simple function, mdbexport, that extracts all tables from a Microsoft Access database file and saves each table to a separate comma-delimited text file.

This function was written and tested in a Windows environment, and requires the appropriate database driver (Microsoft Access Driver (*.mdb)) in order to function.  

It may be possible to get this running on Linux/Mac with the unixODBC driver, but I have not tested that.  I prefer to use the mdbtools library on those platforms, since it already has the unixODBC capabilities bundled.

I've added a setup script to this folder, in case any users want to install the module locally.  To do so, from the command prompt, navigate to this folder ([your-install-location]\ecopath_matlab\mdbexport) and run

`python setup.py install`  

However, if you're only using this function via its internal call from mdb2ecopathmodel via Matlab, then the local install isn't necessary.  That function is designed to add this folder to the python system path within Matlab and then read the function directly.