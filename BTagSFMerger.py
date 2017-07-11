from ROOT import TChain
import ROOT
import glob
import xml.etree.cElementTree as ET
import os
from datetime import datetime


# Define time variable
#date = "test"
date = "170711"

scalesys = "nominal"
#scalesys = "JERup"

sys = ["central","up","down"]
systype = "mujets"


#Define path where btag histos are stored
pathNonMerged = "BTagHistos/"+date+"/"
pathMerged = "BTagHistos/"+date+"/"+"Merged/"

if not os.path.exists(pathMerged):
    os.makedirs(pathMerged)

# get filenames from the xml!!!
tree = ET.ElementTree(file='config/topWidth_MC_loc.xml')


# get the list of dataset
root = tree.getroot()
datasets = root.find('datasets')
print "found  "  + str(len(datasets)) + " datasets"
datasetNames = []

print ""
# loop over the datasets to be added and fill the "topTrees" vector
for d in datasets:
    if d.attrib['add'] == '1':
        print "found dataset to be added..." + str(d.attrib['name'])

        # select a subset of the existing root file
        datasetNames.append(str(d.attrib['name']))
        print str(d.attrib['name'])


listOfZombies= []

# loop over data set to search root files
for n in datasetNames:
    for s in sys:
        filenames = glob.glob(pathNonMerged + "/*" + n + "*" + s + ".root")
        hadd = "hadd -f " + pathMerged + "/BTagSFs_"+ n + "_" + systype + "_" + s + ".root"

        if (len(filenames) == 0):
            print "no root files found in directory" , pathNonMerged ,  " for dataset " , n , " and sys " , s , " !!"
        else :
            # loop over root files
            for f in filenames:
                file=ROOT.TFile(f,"read")
                # check if the file is a zombie
                if (file.IsZombie()):
                    print "File" , f, "is a Zombie.... Skipping"
                    listOfZombies.append(f)
                else:
                    print f
                    hadd = hadd + " " + f
            print "Merging btag SF histos for " + n + " " + s
            os.system(hadd)

print "\n\n"

# print the list of zombies
print "The total number of zombie files is ", len(listOfZombies)
if (len(listOfZombies) > 0):
    outfile = open (pathMerged+"/Zombie_"+sys+".txt", 'a')
    print "And the list of the zombie is put in "+pathMerged+"/Zombie_"+sys+".txt "
    for zombie in listOfZombies:
        print >> outfile, zombie

print "\n\n"


# Further merging of data if there are multiple data collections
mergeData=False   # exclude TT_nominal_backup (already in TT_nominal) and TT_widthxX

if (mergeData):
    for s in sys:
        # combining all the dataset histos in one
        dataList=glob.glob(pathMerged+"*"+s+".root")

        cmd = "hadd " + pathMerged + "/PlotsForBTagSFs_" + scalesys + "_" + s + ".root"
        for data in dataList:
            cmd = cmd + " " + data
        print "Merging btag SF histos for sys " + s
        os.system(cmd)
