# Building

Clone the repository and build all images using `./build.sh`

# ARES

Prerequisites:

1. Install reduce from https://github.com/rlabduke/reduce
2. Process your PDB files with reduce before scoring

Docker usage:

```bash
# First process your PDB with reduce:
# Assuming $pdb is the path to your input PDB file and $pdbdir is where to store processed PDBs
reduce ${pdb} > ${pdbdir}/output.pdb

# Then run ARES in Docker:
# Assuming $pdbdir points to directory with processed PDBs
docker run -it --rm -v ${pdbdir}:/tmp/pdb adamczykb/ares_qa:latest

# Run this in the container's shell!
python3 -m ares.predict /tmp/pdb data/weights.ckpt output.csv -f pdb --nolabels --num_workers=$(nproc)
cat output.csv
```

# RNA-BRiQ

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb briq:latest bash -c "BRiQ_AssignSS /tmp/input.pdb /tmp/input.briq; sed -i '1 i pdb /tmp/input.pdb' /tmp/input.briq; BRiQ_Energy /tmp/input.briq >/dev/null 2>/tmp/output; cat /tmp/output"
```

# lociPARSE

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb lociparse:latest lociPARSE /tmp/input.pdb
```

# cgRNASP, cgRNASP-C, cgRNASP-PC

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp:latest cgRNASP /tmp
```

You can substitute `cgRNASP` with `cgRNASP-C` or `cgRNASP-PC` to run a different variant.

# cgRNASP-CN

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp-cn:latest cgRNASP-CN /tmp
```

Note! This uses a different Docker image than cgRNASP, cgRNASP-C and cgRNASP-CN

# rsRNASP

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb rsrnasp:latest rsRNASP /tmp/input.pdb
```

# 3dRNAscore

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb 3drnascore:latest bash -c "perl /opt/3dRNAscore/lib/format.pl /tmp/input.pdb > /tmp/formatted.pdb; 3dRNAscore -s /tmp/formatted.pdb"
```

# RASP

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb rasp:latest bash -c "rasp_fd -p /tmp/input.pdb"
```

# DFIRE

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb dfire-rna:latest bash -c "DFIRE_RNA /tmp/input.pdb"
```

# RNA3DCNN

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
# Using MD model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MD.hdf5 -local 0

# Or using MDMC model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MDMC.hdf5 -local 0
```
