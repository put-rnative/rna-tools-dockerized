@echo off
title To jest podejscie testowe
set arg1=%1
set arg2=%2
::echo %arg1%
::docker cp casp-15-pdbs.tar.gz sweet_proskuriakova:/home/3dRNAscore/
::docker exec --workdir /home/dfire eager_benz tar xvf casp.tar.gz
::docker exec --workdir /home/dfire eager_benz ls
::R1107/pdb <-- folder z tymi tymi

::docker exec --workdir /opt/dfire/%arg2% funny_chaplygin ../../../bin/DFIRE_RNA %arg1% >dfire.txt
docker exec --workdir /home/3dRNAscore/%arg2% fervent_shamir /home/3dRNAscore/bin/3dRNAscore -s %arg1% >3dRNA.txt

::docker exec --workdir /home/3dRNAscore/example sweet_proskuriakova /home/3dRNAscore/bin/3dRNAscore -s 1A9L.pdb >3dRNA.txt
::pause