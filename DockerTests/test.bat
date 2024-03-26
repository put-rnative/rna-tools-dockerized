@echo off
title dfire
set arg1=%1
set arg2=%2
echo %arg1%
::docker cp test3.txt eager_benz:/home/dfire/dane
::docker exec --workdir /home/dfire eager_benz tar xvf casp.tar.gz
::docker exec --workdir /home/dfire eager_benz ls
::R1107/pdb <-- folder z tymi tymi
::docker cp randstr hungry_lalande:/home/dfire/randstr
docker exec --workdir /opt/dfire/%arg2% funny_chaplygin ../../../bin/DFIRE_RNA %arg1% >dfire.txt
::pause