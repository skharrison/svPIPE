import os 
import sys 

outDir = sys.argv[1] 
pDir = outDir + "/pindel/"
statDir = outDir + "/stats"
bamDir = outDir + "/bams/"
for file in os.listdir(statDir):
    bam = "."
    if file.endswith(".stats"):
        mySample = file.split(".")[0]
        theBam = file.split(".")
        theBam = theBam[:-1]
        theBam.append("bam")
        bam = bam.join(theBam)
        bam = bamDir + bam
        myStat = statDir + "/" + file 
        file1 = open(myStat, "r")
        for line in file1.readlines():
            if line.startswith("SN"):
                if line.split("\t")[1] == "insert size average:":
                    insertSize = line.rstrip().split("\t")[2]
        myText = bam + "\t" + insertSize + "\t" + mySample
        newName = pDir + mySample + "_pindel.cfg"
        pFile = open(newName, "w")
        pFile.write(myText)
        pFile.close()
        file1.close() 