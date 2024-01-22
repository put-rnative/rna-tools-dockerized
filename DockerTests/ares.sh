#!/bin/bash
python3 -m ares.predict data/sample/remoteTest data/weights.ckpt wynik.csv -f pdb --nolabels --num_workers=8