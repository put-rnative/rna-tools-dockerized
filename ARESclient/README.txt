This client app sends pdb files to ARES tool. 
Results are sent to e-mail adres in config.txt
You can edit this file manually or use
python arescli.py chm [newemail]


To send single file use 
python arescli.py [filename]


To send multiple files, place them in one folder and use
python arescli.py dirsend [foldername]

Note: pdb files in this project are taken from cgRNASP-CN example folder 
https://github.com/Tan-group/cgRNASP-CN/tree/main/example