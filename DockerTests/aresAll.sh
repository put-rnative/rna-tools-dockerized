#!/bin/bash


for f in $1/*;
do
        echo "$f"
        #nazwa pliku
        filename=$(basename "$f")
        #przekopiuj do data/sample/remoteTest
        cp $1/$filename ./data/sample/remoteTest
        ls ./data/sample/remoteTest
        #policz wynik i zapisz w remotetest.csv
        python3 -m ares.predict data/sample/remoteTest data/weights.ckpt remotetest.csv -f pdb --nolabels --num_workers=8
        #wyczysc linijke z naglowkami
        sed -i '1d' remotetest.csv
        #otworz remotetest.csv, zapisz do remoteresult.csv wyczysc remotetest.csv
        cat remotetest.csv >> remoteresult.csv
        # "" >remotetest.csv
        #usun
        rm ./data/sample/remoteTest/$filename
done