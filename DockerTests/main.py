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

csv_dfire("136")
csv_dfire("138")
csv_dfire("149")
csv_dfire("156")
csv_dfire("189")
csv_dfire("190")
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
