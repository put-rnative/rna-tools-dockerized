# RNA Structure Scoring Tools

This repository provides Docker images for various RNA structure scoring methods. The recommended way to use these tools is through our unified wrapper image that incorporates all scoring methods except ARES (which requires special preprocessing and must be used separately):

```bash
# Assuming $pdb_dir contains your PDB files:
docker run --rm -v ${pdb_dir}:/data rna-tools bash -c 'scoring-wrapper.py /data/*.pdb'
```

The wrapper provides several features:
- Parallel processing of files and methods
- Automatic checkpointing to resume interrupted runs
- Unified CSV output format
- Progress tracking with completion estimates

## Wrapper Usage

```bash
scoring-wrapper.py [-h] [--scoring-method {3dRNAscore,DFIRE,RASP,RNA-BRiQ,RNA3DCNN_MD,RNA3DCNN_MDMC,cgRNASP,cgRNASP-C,cgRNASP-CN,cgRNASP-PC,lociPARSE,rsRNASP}] 
                  [--output FILE.csv]
                  pdb_files [pdb_files ...]
```

The wrapper will:
1. Process multiple PDB files in parallel
2. Save progress automatically to a checkpoint file (unique per set of files/methods)
3. Resume from checkpoint if interrupted
4. Output results in a CSV format with one row per PDB file and columns for each scoring method

Example output:

| PDB File | 3dRNAscore | DFIRE | RASP | RNA3DCNN_MD | RNA3DCNN_MDMC | cgRNASP | cgRNASP-C | cgRNASP-CN | cgRNASP-PC | lociPARSE | rsRNASP | RNA-BRiQ |
|----------|------------|-------|------|-------------|---------------|----------|-----------|------------|------------|-----------|----------|----------|
| 1a9nR_M10.pdb | 11.605 | -10349.440 | -7606.76 | 15.971 | 17.630 | -169.709 | -332.194 | -81.094 | -100.894 | 0.67 | -2465.760 | 167.642 |
| 1a9nR_M11.pdb | 13.646 | -9617.315 | -6666.93 | 16.434 | 18.388 | -163.088 | -356.332 | -82.149 | -100.366 | 0.69 | -2186.164 | 238.550 |
| 1a9nR_M12.pdb | 10.691 | -10420.018 | -8273.36 | 15.640 | 17.739 | -174.046 | -334.797 | -87.357 | -96.707 | 0.69 | -2454.089 | 110.188 |

## Individual Method Images

If you prefer to use individual scoring methods separately, the following Docker images are available:

# ARES

ARES (Atomic Rotationally Equivariant Scorer) is a geometric deep learning method for RNA structure assessment described in:

Townshend, R. J. L., Eismann, S., Watkins, A. M., Rangan, R., Karelina, M., Das, R., & Dror, R. O. (2021). Geometric deep learning of RNA structure. Science, 373(6558), 1047–1051. https://doi.org/10.1126/science.abe5650

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

RNA-BRiQ (RNA Base-pair and Residue Interaction Quality) is a high-resolution statistical potential for RNA structure assessment described in:

Xiong, P., Wu, R., Zhan, J., & Zhou, Y. (2021). Pairing a high-resolution statistical potential with a nucleobase-centric sampling algorithm for improving RNA model refinement. Nature Communications, 12(1), 2777. https://doi.org/10.1038/s41467-021-23100-4

Source code: https://github.com/Jian-Zhan/RNA-BRiQ

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb briq:latest bash -c "BRiQ_AssignSS /tmp/input.pdb /tmp/input.briq; sed -i '1 i pdb /tmp/input.pdb' /tmp/input.briq; BRiQ_Energy /tmp/input.briq >/dev/null 2>/tmp/output; cat /tmp/output"
```

# lociPARSE

lociPARSE (Locality-aware Invariant Point Attention Model for Scoring RNA 3D Structures) is described in:

Tarafder, S., & Bhattacharya, D. (2024). lociPARSE: A Locality-aware Invariant Point Attention Model for Scoring RNA 3D Structures. Journal of Chemical Information and Modeling, 64(22), 8655–8664. https://doi.org/10.1021/acs.jcim.4c01621

Source code: https://github.com/Bhattacharya-Lab/lociPARSE

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb lociparse:latest lociPARSE /tmp/input.pdb
```

# cgRNASP, cgRNASP-C, cgRNASP-PC

cgRNASP (Coarse-Grained RNA Statistical Potentials) is a family of statistical potentials for RNA structure evaluation described in:

Tan, Y.-L., Wang, X., Yu, S., Zhang, B., & Tan, Z.-J. (2023). cgRNASP: Coarse-grained statistical potentials with residue separation for RNA structure evaluation. NAR Genomics and Bioinformatics, 5(1), lqad016. https://doi.org/10.1093/nargab/lqad016

Source code: https://github.com/Tan-group/cgRNASP

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp:latest cgRNASP /tmp
```

You can substitute `cgRNASP` with `cgRNASP-C` or `cgRNASP-PC` to run a different variant.

# cgRNASP-CN

cgRNASP-CN is a minimal coarse-grained representation-based statistical potential described in:

Song, L., Yu, S., Wang, X., Tan, Y.-L., & Tan, Z.-J. (2022). cgRNASP-CN: A minimal coarse-grained representation-based statistical potential for RNA 3D structure evaluation. Communications in Theoretical Physics, 74(7), 075602. https://doi.org/10.1088/1572-9494/ac7042

Source code: https://github.com/Tan-group/cgRNASP-CN

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb cgrnasp-cn:latest cgRNASP-CN /tmp
```

Note! This uses a different Docker image than cgRNASP, cgRNASP-C and cgRNASP-CN

# rsRNASP

rsRNASP (Residue-Separation-based RNA Statistical Potential) is described in:

Tan, Y.-L., Wang, X., Shi, Y.-Z., Zhang, W., & Tan, Z.-J. (2022). rsRNASP: A residue-separation-based statistical potential for RNA 3D structure evaluation. Biophysical Journal, 121(1), 142–156. https://doi.org/10.1016/j.bpj.2021.11.016

Source code: https://github.com/Tan-group/rsRNASP

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb rsrnasp:latest rsRNASP /tmp/input.pdb
```

# 3dRNAscore

3dRNAscore is a distance and torsion angle dependent evaluation function described in:

Wang, J., Zhao, Y., Zhu, C., & Xiao, Y. (2015). 3dRNAscore: A distance and torsion angle dependent evaluation function of 3D RNA structures. Nucleic Acids Research, 43(10), e63–e63. https://doi.org/10.1093/nar/gkv141

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

DFIRE-RNA is an all-atom knowledge-based potential for RNA structure discrimination described in:

Zhang, T., Hu, G., Yang, Y., Wang, J., & Zhou, Y. (2020). All-Atom Knowledge-Based Potential for RNA Structure Discrimination Based on the Distance-Scaled Finite Ideal-Gas Reference State. Journal of Computational Biology, 27(6), 856–867. https://doi.org/10.1089/cmb.2019.0251

Source code: https://github.com/tcgriffith/dfire_rna

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
docker run --rm -v ${pdb}:/tmp/input.pdb dfire-rna:latest bash -c "DFIRE_RNA /tmp/input.pdb"
```

# RNA3DCNN

RNA3DCNN is a deep learning-based method for RNA 3D structure quality assessment described in:

Li, J., Zhu, W., Wang, J., Li, W., Gong, S., Zhang, J., & Wang, W. (2018). RNA3DCNN: Local and global quality assessments of RNA 3D structures using 3D deep convolutional neural networks. PLOS Computational Biology, 14(11), e1006514. https://doi.org/10.1371/journal.pcbi.1006514

Source code: https://github.com/lijunRNA/RNA3DCNN

Docker usage:

```bash
# Assuming $pdb is the path to your PDB file:
# Using MD model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MD.hdf5 -local 0

# Or using MDMC model:
docker run --rm -v ${pdb}:/tmp/input.pdb rna3dcnn:latest python /opt/RNA3DCNN/Main.py -pn /tmp/input.pdb -model /opt/RNA3DCNN/RNA3DCNN_MDMC.hdf5 -local 0
```
