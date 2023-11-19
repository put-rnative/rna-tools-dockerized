@echo off
title To jest podejscie testowe
set arg1=%1
set arg2=%2
::echo %arg1%
::docker cp casp-15-pdbs.tar.gz laughing_lamport:/home/cgRNASP/cgRNASP/
::docker exec --workdir /home/cgRNASP/cgRNASP/ laughing_lamport tar xvf casp-15-pdbs.tar.gz
::docker exec --workdir /home/cgRNASP/cgRNASP/ laughing_lamport ls
::R1107/pdb <-- folder z tymi tymi

docker exec -i --workdir /home/cgRNASP/cgRNASP/%arg2%  laughing_lamport /home/cgRNASP/cgRNASP/example/cgRNASP  ./ %arg1% energy_list.txt >cgrnasdp.txt
docker exec -i --workdir /home/cgRNASP/cgRNASP/%arg2%  laughing_lamport ls
docker exec -i --workdir /home/cgRNASP/cgRNASP/%arg2%  laughing_lamport cat energy_list.txt >cgrnasp.txt
::docker cp laughing_lamport:/home/cgRNASP/cgRNASP/%arg2%/energy_list.txt ./cgrnasp.txt 
::docker exec --workdir /home/3dRNAscore/example laughing_lamport /home/3dRNAscore/bin/3dRNAscore -s 1A9L.pdb >3dRNA.txt
pause