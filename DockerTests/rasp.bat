@echo off
title Rasp testing
set arg1=%1
set arg2=%2

::jakby zwracal ze costam costam /bash bad interpreter:
::sed -i -e 's/\r$//' cgRNASP

::echo %arg1%
::docker cp casp-15-pdbs  admiring_swanson:/home/rasp-fd-1.0/rasp-fd-1.0/
::docker exec --workdir /home/cgRNASP/cgRNASP/ laughing_lamport tar xvf casp-15-pdbs.tar.gz
::docker exec --workdir /home/cgRNASP/cgRNASP/ laughing_lamport ls
::R1107/pdb <-- folder z tymi tymi

::docker cp test.pdb laughing_lamport:/home/cgRNASP/cgRNASP/
docker exec -i --workdir /home/  admiring_swanson rasp-fd-1.0/rasp-fd-1.0/bin/rasp_fd -e all -p rasp-fd-1.0/rasp-fd-1.0/casp-15-pdbs/%arg2%/%arg1% >rasp.txt
::docker exec -i --workdir /home/cgRNASP/cgRNASP/%arg2%  laughing_lamport ls
::docker exec -i --workdir /home/cgRNASP/cgRNASP/%arg2%  laughing_lamport cat energy_list.txt >cgrnasp.txt


::docker cp laughing_lamport:/home/cgRNASP/cgRNASP/%arg2%/energy_list.txt ./cgrnasp.txt 
::docker exec --workdir /home/3dRNAscore/example laughing_lamport /home/3dRNAscore/bin/3dRNAscore -s 1A9L.pdb >3dRNA.txt
::pause