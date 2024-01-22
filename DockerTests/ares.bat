

@echo off
title kopiowanie na aresa
set arg1=%1
::set arg2=%2

::jakby zwracal ze costam costam /bash bad interpreter:
::sed -i -e 's/\r$//' cgRNASP
::docker run -d --ipc="host" -it adamczykb/ares_qa
::python3 -m ares.predict data/sample/remoteTest data/weights.ckpt wynik.csv -f pdb --nolabels --num_workers=8
::docker exec --workdir /home/3dRNAscore/%arg2% priceless_margulis /home/3dRNAscore/bin/3dRNAscore -s %arg1% >3dRNA.txt
::echo %arg1%
::docker cp ares.sh priceless_margulis:/app/release/ares.sh 
docker cp reduce/reduce priceless_margulis:/app/ares_release/data/sample/reduce
docker exec --workdir /app/ares_release/ priceless_margulis python3 -m ares.predict data/sample/remoteTest data/weights.ckpt wynik.csv -f pdb --nolabels --num_workers=8
docker exec --workdir /app/ares_release/ priceless_margulis cat wynik.csv >ares.txt
docker exec --workdir /app/ares_release/data/sample/remoteTest/ rm %arg1%
