#! /bin/bash
for target in 3dRNAscore BRiQ cgRNASP cgRNASP-CN DFIRE-RNA RASP RNA3DCNN rsRNASP; do
	(
		cd ${target}
		docker build -t ${target,,*} .
	)
done

docker pull adamczykb/ares_qa
