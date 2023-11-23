@echo off
::docker container ps -l
title To jest podejscie testowe
set arg1=%1
set arg2=%2
::echo %arg1%
docker cp casp-15-pdbs.tar.gz nostalgic_dijkstra:/home/RNA-BRiQ/
docker exec --workdir /home/RNA-BRiQ nostalgic_dijkstra tar xvf casp-15-pdbs.tar.gz
docker exec --workdir /home/RNA-BRiQ nostalgic_dijkstra ls
::R1107/pdb <-- folder z tymi tymi
::docker exec --workdir /home/3dRNAscore/%arg2% nostalgic_dijkstra /home/3dRNAscore/bin/3dRNAscore -s %arg1% >3dRNA.txt

::docker exec --workdir /home/3dRNAscore/example sweet_proskuriakova /home/3dRNAscore/bin/3dRNAscore -s 1A9L.pdb >3dRNA.txt
pause