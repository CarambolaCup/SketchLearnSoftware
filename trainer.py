import csv
import os
import sys
import subprocess
import time

from matplotlib.pyplot import close

result = open('D://workspace//SketchLearnSoftware//result.csv','w')
with result:
    fnames = ['eta','e1','e2','e3','e4']
    writer = csv.DictWriter(result,fieldnames=fnames)
    writer.writeheader()

    for i in range(500,1000,50):
        eta = i / 1000
        cmd = 'cmd /c D://workspace//SketchLearnSoftware//train18000.exe "{0}"'.format(eta)
        print(eta)
        p = subprocess.Popen(cmd,shell=True)
        return_code=p.wait()  #等待子进程结束，并返回状态码；
        time.sleep(49)
        testout = open('D://workspace//SketchLearnSoftware//out.txt','r')
        with testout:
            line = testout.readline()
            line = line.strip()
            a = line.split(',')
            writer.writerow({'eta':a[0],'e1':a[1],'e2':a[2],'e3':a[3],'e4':a[4]})
            print(i)




