import docker
import subprocess
import os

#Dfire zwraca "atomtype error U_OP3 not found" dla niektórych plików, ale liczy wynik i tak

print("reste")

def csv_dfire(numerPliku):
    i =0
    nazwa_csv = "dfireR1"+numerPliku+".csv"
    folder_w_maszynie = "R1"+numerPliku+"/pdb"
    for f in os.listdir("casp-15-pdbs/"+folder_w_maszynie):
        # if os.path.isfile(f):
        i+=1
        komenda = "test.bat" + " "+f+" "+folder_w_maszynie
        os.system(komenda)
        val = ""
        with open("dfire.txt", "r") as myfile:
            val = myfile.read()
        with open(nazwa_csv, "a") as myfile:
            myfile.write(val)
        with open("dfire.txt", "w") as myfile:
            myfile.write("")
    print("number of files: "+ str(i))

def csv_dfireFolder(folder):
    
    for d in os.listdir(folder):
        print("folder " + folder + "d: "+ d)
        i =0
        nazwa_csv = "dfire_"+d+".csv"
        folder_w_maszynie = folder+"/"+d
        for f in os.listdir(folder_w_maszynie):
            # if os.path.isfile(f):
            i+=1
            komenda = "test.bat" + " "+f+" "+folder_w_maszynie
            os.system(komenda)
            val = ""
            with open("dfire.txt", "r") as myfile:
                val = myfile.read()
            with open(nazwa_csv, "a") as myfile:
                myfile.write(val)
            with open("dfire.txt", "w") as myfile:
                myfile.write("")
        print("number of files: "+ str(i))

def csv_3drna(numerPliku):
    i =0
    nazwa_csv = "3drna"+numerPliku+".csv"
    folder_w_maszynie = "R1"+numerPliku+"/pdb"
    for f in os.listdir("casp-15-pdbs/"+folder_w_maszynie):
        # if os.path.isfile(f):
        i+=1
        komenda = "3dRNA.bat" + " "+f+" "+folder_w_maszynie
        os.system(komenda)
        val = ""
        with open("3dRNA.txt", "r") as myfile:
            val = myfile.read()
        with open(nazwa_csv, "a") as myfile:
            myfile.write(val)
        with open("3dRNA.txt", "w") as myfile:
            myfile.write("")
    print("number of files: "+ str(i))

def csv_3drnaFolder(folder):
    
    for d in os.listdir(folder):
        print("folder " + folder + "d: "+ d)
        i =0
        nazwa_csv = "3dRNA"+d+".csv"
        folder_w_maszynie = folder+"/"+d
        for f in os.listdir(folder_w_maszynie):
            # if os.path.isfile(f):

            #check if file is pdb

            i+=1
            komenda = "3dRNA.bat" + " "+f+" "+folder_w_maszynie
            os.system(komenda)
            val = ""
            with open("3dRNA.txt", "r") as myfile:
                val = myfile.read()
            with open(nazwa_csv, "a") as myfile:
                myfile.write(val)
            with open("3dRNA.txt", "w") as myfile:
                myfile.write("")
        print("number of files: "+ str(i))

# def csv_cgrnasp(numerPliku):
#     i =0
#     nazwa_csv = "cgrnasp"+numerPliku+".csv"
#     folder_w_maszynie = "R1"+numerPliku+"/pdb"
#     for f in os.listdir("casp-15-pdbs/"+folder_w_maszynie):
#         # if os.path.isfile(f):
#         i+=1
#         komenda = "cgrnasp.bat" + " "+f+" "+folder_w_maszynie
#         os.system(komenda)
#         val = ""
#         with open("cgrnasp.txt", "r") as myfile:
#             val = myfile.read()
#         with open(nazwa_csv, "a") as myfile:
#             myfile.write(val)
#         with open("3dRNA.txt", "w") as myfile:
#             myfile.write("")
#     print("number of files: "+ str(i))

def csv_cgrnasp():
    i=0
    for f in os.listdir("casp-15-pdbs/"):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnasp.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))

def csv_rsrnasp(numerPliku):
    i =0
    nazwa_csv = "rsrnasp"+numerPliku+".csv"
    folder_w_maszynie = "R1"+numerPliku
    for f in os.listdir("casp-15-pdbs/"+folder_w_maszynie+"/pdb"):
        # if os.path.isfile(f):
        i+=1
        komenda = "rsrnasp.bat" + " "+folder_w_maszynie+" "+f
        os.system(komenda)
        val = ""
        with open("rsrnasp.txt", "r") as myfile:
            val = myfile.read()
        with open(nazwa_csv, "a") as myfile:
            myfile.write(val)
        with open("rsrnasp.txt", "w") as myfile:
            myfile.write("")
    print("number of files: "+ str(i))  

def csv_rsrnaspFolder(folder):
    print("number of files: "+ str(i))  
    for d in os.listdir(folder):
        print("folder " + folder + "d: "+ d)
        i =0
        nazwa_csv = "rsrnasp_"+d+".csv"
        folder_w_maszynie = folder+"/"+d
        for f in os.listdir(folder_w_maszynie):
            # if os.path.isfile(f):

            #check if file is pdb

            i+=1
            komenda = "rsrnasp.bat" + " "+folder_w_maszynie+" "+f
            os.system(komenda)
            val = ""
            with open("rsrnasp.txt", "r") as myfile:
                val = myfile.read()
            with open(nazwa_csv, "a") as myfile:
                myfile.write(val)
            with open("rsrnasp.txt", "w") as myfile:
                myfile.write("")
        print("number of files: "+ str(i))




def csv_cgrnasp_Folder(folder):
    i=0
    for f in os.listdir(folder):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnasp.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))
def csv_cgrnasp_cn():
    i=0
    for f in os.listdir("casp-15-pdbs/"):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnaspcn.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))

def csv_cgrnasp_cnFolder(folder):
    i=0
    for f in os.listdir(folder):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnaspcn.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))
def csv_cgrnasp_c():
    i=0
    for f in os.listdir("casp-15-pdbs/"):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnaspc.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))
def csv_cgrnasp_pc():
    i=0
    for f in os.listdir("casp-15-pdbs/"):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnasppc.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))

def csv_cgrnasppc_Folder(folder):
    i=0
    for f in os.listdir(folder):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnasppc.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))

def csv_cgrnaspc_Folder(folder):
    i=0
    for f in os.listdir(folder):
        # if os.path.isfile(f):
        print(f)
        i+=1
        komenda = "cgrnaspc.bat" + " "+f
        os.system(komenda)
    print("number of files: "+ str(i))


def rasp(numerPliku):
    i =0
    nazwa_csv = "rasp"+numerPliku+".csv"
    folder_w_maszynie = "R1"+numerPliku+"/pdb"
    for f in os.listdir("casp-15-pdbs/"+folder_w_maszynie):
        # if os.path.isfile(f):
        i+=1
        komenda = "rasp.bat" + " "+f+" "+folder_w_maszynie
        os.system(komenda)
        val = ""
        with open("rasp.txt", "r") as myfile:
            val = myfile.read()
        with open(nazwa_csv, "a") as myfile:
            myfile.write(val)
        with open("rasp.txt", "w") as myfile:
            myfile.write("")
    print("number of files: "+ str(i))


def ares(numerPliku):
    i =0
    nazwa_csv = "ares/ares.csv"
    for f in os.listdir("reduce/reduce/"):
        # if os.path.isfile(f):
        i+=1
        komenda = "ares.bat" + " "+f
        os.system(komenda)
        val = ""
        with open("ares/ares.txt", "r") as myfile:
            val = myfile.read()
        with open(nazwa_csv, "a") as myfile:
            myfile.write(val)
        with open("ares/ares.txt", "w") as myfile:
            myfile.write("")
    print("number of files: "+ str(i))

# rasp("107")
# #csv_cgrnasp_pc()

# rasp("108")
# rasp("116")
# rasp("117")
# rasp("126")
# rasp("128")
# rasp("136")
# rasp("138")
# rasp("149")
# rasp("156")
# rasp("189")
# rasp("190")

csv_3drnaFolder("randstr/decoys")

# komenda = "cgrnasp.bat" + " t R1107/pdb "
# os.system(komenda)
#csv_dfire("136")
# csv_dfire("138")
# csv_dfire("149")
# csv_dfire("156")
# csv_dfire("189")
# csv_dfire("190")
# print('test')
#result =subprocess.call("docker run -it dfire2 bash; cd home; ls", shell=True)#, stdout=output, stderr=output
# result = subprocess.run(["docker run -it dfire2 bash; cd home; ls"], shell=True, capture_output=True, text=True)
# s = """
# echo hello piwo
# echo piwo piwo
# echo sesja dziekan"""
# result = subprocess.run([s],shell=True,)

# print(repr(result.stdout))


# cmd=['cd',r'C:\Program Files (x86)\Notepad++','&&','notepad','LICENSE','&&',r'D:\Program\Tools\Putty.exe','-v']
# d=subprocess.Popen(cmd, shell=True)
#p = subprocess.run('docker exec --workdir /home eager_benz ls', shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
# subprocess.run([r"test.bat"])
#print(p)
# p = subprocess.Popen(["python", "--help","&&","echo","piwo"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
# # p = subprocess.Popen(["python, --help"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
# output, errors = p.communicate()

# print(output)
# with open("output.log", "a") as output:
#     #result =subprocess.call("docker run -it dfire2 bash; cd home; ls", shell=True, stdout=output, stderr=output)
#     result = subprocess.run(["docker run -it dfire2 bash"], shell=False, capture_output=True, text=True)
#     #subprocess.run('docker run -it dfire2 ', shell=True, executable="/bin/bash", stdout=output)
#     #subprocess.run("docker run -it dfire2 bash", shell=True, stdout=output, stderr=output)
#     print(result)
