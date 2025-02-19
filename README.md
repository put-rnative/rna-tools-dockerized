ARES:
we have an image https://hub.docker.com/r/adamczykb/ares_qa, however it doesn't work with my Docker on Windows.
Webserwer (DOESN'T WORK), go to page and upload a PDB file to get results.
 https://drorlab.stanford.edu/ares.html
 

RNA-BRiQ-main:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb briq:latest bash -c "BRiQ_AssignSS /tmp/input.pdb /tmp/input.briq; sed -i '1 i pdb /tmp/input.pdb' /tmp/input.briq; BRiQ_Energy /tmp/input.briq >/dev/null 2>/tmp/output; cat /tmp/output"
```
The command:
1. Mounts your PDB file into the container
2. Assigns secondary structure using BRiQ_AssignSS
3. Prepares the input file for energy calculation
4. Runs BRiQ_Energy and captures the output


lociPARSE:
WORKS, calculated, with minor changes to code (added ",map_location='cuda:0'" to torch.load() call)

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb lociparse:latest lociPARSE /tmp/input.pdb
```
The command mounts your PDB file into the container and runs the lociPARSE scoring algorithm.

cgRNASP:
WORKS, calculated 

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp:latest cgRNASP /tmp
```
The command mounts your PDB file into the container and runs the cgRNASP scoring algorithm.

cgRNASP-CN:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp-cn:latest cgRNASP-CN /tmp
```
The command mounts your PDB file into the container and runs the cgRNASP-CN scoring algorithm.

cgRNASP-C:
WORKS, calculated

Docker usage:
Uses the same Docker image as cgRNASP but with a different command:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp:latest cgRNASP-C /tmp
```

cgRNASP-PC:
WORKS, calculated

Docker usage:
Uses the same Docker image as cgRNASP but with a different command:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp:latest cgRNASP-PC /tmp
```

rsRNASP:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb rsrnasp:latest rsRNASP /tmp/input.pdb
```
The command mounts your PDB file into the container and runs the rsRNASP scoring algorithm.

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

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb rasp:latest bash -c "rasp_fd -p /tmp/input.pdb"
```
The command mounts your PDB file into the container and runs the RASP scoring algorithm.

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

RNA3DCNN:
WORKS, calculated

Docker usage:
```bash
# Assuming $pdb is the path to your PDB file:
# Using MD model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MD.hdf5 -local 0

# Or using MDMC model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MDMC.hdf5 -local 0
```
The command mounts your PDB file into the container and runs RNA3DCNN with either the MD or MDMC model.
