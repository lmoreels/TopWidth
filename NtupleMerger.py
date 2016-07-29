from ROOT import TChain
import ROOT
import glob
import xml.etree.cElementTree as ET
import os
from datetime import datetime


# Define time variable
now = datetime.now()
dd = str(now.day).zfill(2)
mm = str(now.month).zfill(2)
yyyy = str(now.year)
yy = str(now.year-2000)
# pick one of the two above
#date = "160602"
#date = "17_1_2016"
#date = yy+mm+dd
date = "160727"

#channels = ["mu","el"] 
channels = ["mu"]

for chan in channels:
    
    #Define path where ntuples are stored
    pathNonMerged = "NtupleOutput/"+chan+"/"+date+"/"
    pathMerged = "NtupleOutput/MergedTuples/"+chan+"/"+date+"/"
    
    if not os.path.exists(pathMerged):
        os.makedirs(pathMerged)
    
    # get filenames from the xml!!!    
    if "mu" in chan:
        tree = ET.ElementTree(file='config/topWidth_localgrid.xml')
    elif "el" in chan:
        tree = ET.ElementTree(file='config/topWidth_el_localgrid.xml')
    else:
        print "No tree has been loaded!!! Make sure the correct xml file are in the right directories!!!"
        sys.exit()

    #tree = ET.ElementTree(file='config/FullMcBkgdSamplesV9.xml')
    #tree = ET.ElementTree(file='config/DisplacedTopsSignal.xml')
    #tree = ET.ElementTree(file='config/DataSamples.xml')

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
            if not "over" in str(d.attrib['name']) :
                datasetNames.append(str(d.attrib['name']))
                print str(d.attrib['name'])
    
    
    listOfZombies= []
    
    # loop over data set to search root files
    for n in datasetNames:
        filenames = glob.glob(pathNonMerged + "/*" + n + "*.root")
        hadd = "hadd -f " + pathMerged + "/Ntuples_"+ n + ".root"

        if (len(filenames) == 0):
            print "no root files found in directory" , pathNonMerged ,  " for dataset " , n , " !!"
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
            print "Merging ntuples for " + n
            os.system(hadd)
        
    print "\n\n"
    
    # print the list of zombies
    print "The total number of zombie files is ", len(listOfZombies)
    if (len(listOfZombies) > 0):
        outfile = open (pathMerged+"/Zombie_"+chan+".txt", 'a')
        print "And the list of the zombie is put in "+pathMerged+"/Zombie_"+chan+".txt "
        for zombie in listOfZombies:
            print >> outfile, zombie
    
    
    # Further merging of data if there are multiple data collections
    mergeData=False
    
    if (mergeData):
    # combining all the Data in one
        dataList=glob.glob(pathMerged+"*data*.root")
    
	cmd = "hadd " + pathMerged + "/Ntuples_"+ "data.root"
        for data in dataList:
            cmd = cmd + " " + data
        os.system(cmd)
