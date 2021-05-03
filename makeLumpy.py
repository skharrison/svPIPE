import os 
import sys 

outDir = sys.argv[1] 



myFile = "#!/bin/bash" + "\n"
outDir = sys.argv[1] 
lDir = outDir + "/lumpy/"
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
                if line.split("\t")[1] == "insert size standard deviation:":
                    sd = line.rstrip().split("\t")[2]
                if line.split("\t")[1] == "average length:":
                    readLen = line.rstrip().split("\t")[2]
        myText = "lumpy -mw 4 -tt 0 -pe id:" + mySample + ",bam_file:" + bamDir + mySample + ".discordants.bam," + "histo_file:" + bamDirc + mySample + ".lib1.histo,mean:" + insertSize + ",stdev:" + sd + ",read_length:" + readLen + ",min_non_overlap:" + readLen + ",discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20" + " -sr id:" + mySample + ",bam_file:" + bamDir + mySample + ".splitters.bam" + ",back_distance:10,weight:1,min_mapping_threshold:20" + " > " + lDir + mySample + "_lumpy.vcf" + "\n"
        myFile += myText

lFile = open(lDir + "runLumpy.sh", "w")
lFile.write(myFile)
lFile.close()
















for file in os.listdir("/home/sharrison/data/bams/"):
    bam = "."
    if file.endswith(".stats"):
        mySample = file.split(".")[0]
        theBam = file.split(".")
        theBam = theBam[:-1]
        bam = bam.join(theBam)
        myBam = "/home/sharrison/data/bams/" + file 
        file1 = open(myBam, "r")
        for line in file1.readlines():
            if line.startswith("SN"):
                if line.split("\t")[1] == "insert size average:":
                    insertSize = line.rstrip().split("\t")[2]
                if line.split("\t")[1] == "insert size standard deviation:":
                    sd = line.rstrip().split("\t")[2]
                if line.split("\t")[1] == "average length:":
                    readLen = line.rstrip().split("\t")[2]
        myText = "lumpy -mw 4 -tt 0 -pe id:" + mySample + ",bam_file:" + mySample + ".discordants.bam," + "histo_file:" + mySample + ".lib1.histo,mean:" + insertSize + ",stdev:" + sd + ",read_length:" + readLen + ",min_non_overlap:" + readLen + ",discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20" + " -sr id:" + mySample + ",bam_file:" + mySample + ".splitters.bam" + ",back_distance:10,weight:1,min_mapping_threshold:20" + " > " + mySample + "_lumpy.vcf" + "\n"
        myFile += myText


pFile = open("runLumpy.sh", "w")
pFile.write(myFile)
pFile.close()