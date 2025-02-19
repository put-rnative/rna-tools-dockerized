ARES:
we have an image https://hub.docker.com/r/adamczykb/ares_qa, however it doesn't work with my Docker on Windows.
Webserwer (DOESN'T WORK), go to page and upload a PDB file to get results.
 https://drorlab.stanford.edu/ares.html
 

RNA-BRiQ-main:
DOESN'T WORK


lociPARSE:
WORKS, calculated, with minor changes to code (added ",map_location='cuda:0'" to torch.load() call)

cgRNASP:
WORKS, calculated 

cgRNASP-CN:
WORKS, calculated

cgRNASP-C:
WORKS, calculated

cgRNASP-PC:
WORKS, calculated

rsRNASP:
WORKS, calculated

3dRNAscore:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb 3drnascore:latest bash -c "perl /opt/3dRNAscore/lib/format.pl /tmp/input.pdb > /tmp/formatted.pdb; 3dRNAscore -s /tmp/formatted.pdb"
```
The command mounts your PDB file into the container, formats it using the provided script, and runs the scoring algorithm.

rasp-fd-1.0:
WORKS, calculated,  but with minor changes to code (env variable with path on machine changed to just a string added by hand)

dfire_rna-master:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb dfire-rna:latest bash -c "DFIRE_RNA /tmp/input.pdb"
```
The command mounts your PDB file into the container and runs the DFIRE-RNA scoring algorithm.

RNA3DCC:
WORKS, calculated
