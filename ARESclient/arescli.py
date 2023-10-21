import requests
import sys
import os


# chm [email] changes email 
if sys.argv[1] == 'chm':
    f = open("config.txt","w")
    f.write(sys.argv[2])
    f.close()
    exit()



f = open("config.txt", "r")
email = f.readline()
f.close()

url = 'http://167.99.175.117/run/score?email='+email+'&method=ARES'

#if we want to send all files contained in a directory, then this is used
if sys.argv[1] == 'dirsend':
    l = os.listdir(sys.argv[2])
    for fn in l:
        files = {'file': open(sys.argv[2]+'/'+fn, 'rb')}
        r = requests.post(url, files=files)
        r.text
    exit()

#otherwise we send the file given as the argument
files = {'file': open(sys.argv[1], 'rb')}

r = requests.post(url, files=files)
r.text