## mdbexport python package

This package provides a simple function, mdbexport, that extracts all tables from a Microsoft Access database file and saves each table to a separate comma-delimited text file.

This function was written and tested in a Windows environment, and requires the appropriate database driver (Microsoft Access Driver (*.mdb)) in order to function.  

It may be possible to get this running on Linux/Mac with the unixODBC driver, but I have not tested that.  I prefer to use the mdbtools library on those platforms, since it already has the unixODBC capabilities bundled.

### Installation instructions

From the command prompt, navigate to this folder ([your-install-location]\ecopath_matlab\mdbexport) and run

`python setup.py install`  