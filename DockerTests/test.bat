@echo off
title To jest podejscie testowe
set arg1=%1
set arg2=%2
echo %arg1%
::docker cp test3.txt eager_benz:/home/dfire/dane
::docker exec --workdir /home/dfire eager_benz tar xvf casp.tar.gz
::docker exec --workdir /home/dfire eager_benz ls
::R1107/pdb <-- folder z tymi tymi
docker exec --workdir /home/dfire/%arg2% eager_benz ../../bin/DFIRE_RNA %arg1% >dfire.txt
::pause