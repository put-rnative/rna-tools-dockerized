import os


def rna3dcnn_to_csv(filename):
    fnamesplit = filename.split(".")
    newfilename = fnamesplit[0]+".csv"
    csv = open(newfilename, "a")
    with open(filename) as file:
        for line in file:
            if(line.find("Total score for")!=-1):
                line2 = line.partition("Total score for ")
                line3 = line2[2].partition(" is  ")
                csv.write(line3[0]+"\t"+line3[2])

def rna3dcnn_folder(folder):
    for f in os.listdir(folder):
        print(f)
        rna3dcnn_to_csv(folder+"/"+f)

#rna3dcnn_to_csv("R1108_rna3dcnnMC.txt")
rna3dcnn_folder("rna3DCNN")