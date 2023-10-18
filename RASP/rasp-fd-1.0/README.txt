BUILD
-----

For building the package, run in your CLI:

$ make

If you want to optimize the final programme run:

$ make CFLAG=<optimization flag>

where <optimization flag> could be: "-O", "-O2", "-O3". For example:
"CFLAG=-O3".

For re-building, do:

$ make clean
$ make

All executables will be created in a local subfolder named "bin".

Before executing the software for the first time, please read the section below.

EXECUTION
---------

Before executing for the first time, you should set the environment variable "RASP":

    export RASP=/path/to/the/rasp-fd/folder   (in bash or ksh)
    setenv RASP /path/to/the/rasp-fd/folder   (in csh or tcsh)

Set this variable in your shell configuration file, commonly a hidden file with the
end "rc" located in your home folder (e.g: ".tcshrc"). After that do not forget to
either launch a new terminal or to run:

source ~/.tcshrc  (in case your file is ".tcshrc").

If this configuration file does not exist, you must create it. If you do not know how
to do this, please ask your system administrator.

Execute the software without any parameters and their mode of usage will be displayed
in the screen.


EXAMPLE
-------

In order to test the building of the package, enter the folder "example" and run the
following instructions:

$ ../bin/rasp_fd -e all -p 1a9n.pdb

You should get on screen the results saved in the file "output_all.txt".

$ ../bin/rasp_fd -e bbr -p 1a9n.pdb

You should get on screen the results saved in the file "output_bbr.txt". Note that the
three last values may be different since there is a randomization process.

$ ../bin/rasp_profile_fd -e all -r -p 1a9n.pdb

You should get on screen the results saved in the file "profile_all.txt". A file called
"profile.scr" will also be created.


Run "rasmol" as follows so as to visualize the structure and the energies for each residue.

$ rasmol -script profile.scr 1a9n.pdb



For questions, please write to:

Tomas Norambuena A. <tanoramb@puc.cl>
Molecular Bioinformatics Laboratory
Pontificia Universidad Catolica de Chile
