@echo off

set arg1=%1 
::sciezka az do folderu z R-em
set arg2=%2
::nazwa R-a

python Main.py -pl %arg1%\pdb\pdblist -model RNA3DCNN_MD.hdf5 -local 0 >%arg2%_rna3dcnn.txt