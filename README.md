pyBioSig is tool for the identification of putative biosignatures from transcriptomical data. It employs Genetic Algorithms as gene selector and optimizes intra-group cohesion and inter-group dispersion. Python is the main programming language and R is the statistical framework. The project portal is under development.

= Installing =
*pyBioSig* was developed using Python and R. So, the first thing you must do in order to have it working is ensure that all dependencies are installed in your computer.

===The interface is based on GTK===
  # *GTK* (http://www.gtk.org/)

===All programming was done with Python and its modules===
  # *Python* >= 2.6 ( [http://www.python.org] )
  # *PyGTK* >= 2.17 (http://www.pygtk.org/)
  # *Rpy* (http://rpy.sourceforge.net/rpy.html)
  # *pyEvolve* >= 0.5 (http://pyevolve.sourceforge.net/)

===All stats need R and pvclust===
  # *R* (http://www.r-project.org/)
  # *pvclust* (http://www.is.titech.ac.jp/~shimo/prog/pvclust/)


Once you have everything installed, download the pyBioSig package and uncompress it in any folder you desire. 

= Running pyBioSig =
Access the folder you uncompressed pyBioSig and type:

{{{python BioSig_interface.py}}}, in order to use it via the graphical interface

or

{{{python script_name.py}}} in order to use it via python scripts

An example of how to write a python script for pyBioSig can be found at _example_biosig_script.py_



_Please inform any error you detect in this or other page of pyBioSig portal_


