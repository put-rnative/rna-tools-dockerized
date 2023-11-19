@echo off
title Konfiguracja folderow
set arg1=%1
set arg2=%2
echo %arg1%
docker cp casp-15-pdbs.tar.gz sweet_proskuriakova:/home/3dRNAscore/
docker exec --workdir /home/3dRNAscore sweet_proskuriakova tar xvf casp-15-pdbs.tar.gz
docker exec --workdir /home/3dRNAscore sweet_proskuriakova ls
pause
::R1107/pdb <-- folder z tymi tymi