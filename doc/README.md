DOCUMENTATION
=============
The documentation for LSL consists primarily of Python docstrings in the LSL modules themselves.  However, a sphinx-based documentation system has been included.  sphinx extracts the docstrings from the installed LSL package and builds a variety of manual formats, including HTML and PDF.

To build the documentation sphnix needs to be installed on your system.  To list all avaliable documentation styles, run:

    make

To build the HTML documentation, run:

    make html

This command create a complete set of HTML files in the build/html directory.