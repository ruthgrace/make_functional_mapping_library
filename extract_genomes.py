#############################################################################
# These functions allow one to retrieve NCBI taxon ID using the GI number,
#  and then download the corresponding full bacterial genome.
# The input is a csv BLAST output file, downloaded from the NCBI BLAST webtool.
# Ruth Grace Wong
# ruthgracewong@gmail.com
# October 1, 2015
#############################################################################

import urllib2
import re
import random
import os
from bs4 import BeautifulSoup
import ftplib
from ftplib import FTP

# Extracts taxon ID using GI number, using information found on NCBI website.
def getTaxon(giNumber):
    page = urllib2.urlopen('http://www.ncbi.nlm.nih.gov/nuccore/' + giNumber).read()
    soup = BeautifulSoup(page,"html.parser")
    divText = soup.find(submit_url=re.compile("ORGANISM"))
    search = re.search('ORGANISM=([0-9]*)',str(divText))
    if (search):
        return search.group(1)
    return ""

# Extracts species name using GI number, using information found on NCBI website.
#  Return value is an array of length 3 where the first 2 items are the binomial species name and the 3rd is the strain.
def getSpeciesName(giNumber):
    page = urllib2.urlopen('http://www.ncbi.nlm.nih.gov/nuccore/' + str(giNumber)).read()
    soup = BeautifulSoup(page,"html.parser")
    titleWords = soup.find("title").getText().split()
    return titleWords[:3]

# Extracts GI numbers from BLAST output.
# In this case I wanted the all the matches >= 98% identity or the best match
#  for an OTU if none were > 98% identity
def getGI(gi,filename):
    otu = set()
    blastResults = open(filename,'r')
    for line in blastResults:
        line = line.strip()
        items = line.split(",")
        if len(items) == 12:            
            # if you just want to get all the GI numbers, you would replace the
            #  rest of this for loop with:
            # newGI = items[1].split("|")[1]
            # if newGI not in gi:
            #     gi.add(newGI)
            newOTU = items[0]
            pident = float(items[2])
            newGI = items[1].split("|")[1]
            if newOTU not in otu:
                otu.add(newOTU)
                if newGI not in gi:
                    gi.add(newGI)
            elif pident >= 98:
                if newGI not in gi:
                    gi.add(newGI)
    blastResults.close()
    return gi

# Extracts all the binomial species names for all the GI numbers
def getAllSpeciesNames(gi, filename):
    # stores how many strains we've found for each species
    speciesNames = {}
    strainOutput = open(filename,'w')
    for i in gi:
        success = False
        while not success:
            try:
                names = getSpeciesName(i)
                print("gi " + i + " species name " + str(names))
                newSpeciesName = names[0] + " " + names[1]
                newStrainName = newSpeciesName + " " + names[2]
                if newSpeciesName not in speciesNames:
                    tempList = {}
                    tempList[newStrainName] = ""
                    speciesNames[newSpeciesName] = tempList
                    strainOutput.write(newStrainName+"\n")
                else:
                    if newStrainName not in speciesNames[newSpeciesName]:
                        speciesNames[newSpeciesName][newStrainName] = ""
                        strainOutput.write(newStrainName+"\n")
                success=True
            except urllib2.HTTPError:
                print("Trouble accessing website for gi " + i)
    strainOutput.close()
    return speciesNames

# Extracts all the taxon ID for all the GI numbers
def getAllTaxon(gi,filename):
    taxon = set()
    taxonOutput = open(filename,'w')
    for i in gi:
        success = False
        while not success:
            try:
                newTaxon = getTaxon(i)
                print("gi "+i+" taxon "+newTaxon)
                if newTaxon not in taxon and newTaxon!="":
                    taxon.add(newTaxon)
                    taxonOutput.write(newTaxon+"\n")
                success=True
            except urllib2.HTTPError:
                print("Trouble accessing website for gi " + i)
    taxonOutput.close()
    return taxon

# Downloads all genomes associated with the taxon into the specified folder
def getGenomesForTaxon(taxon, foldername):
    # This file was downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
    genbank = open('data/assembly_summary_genbank.txt','r')
    ftp = FTP('ftp.ncbi.nlm.nih.gov')     # connect to host, default port
    ftp.login()
    for line in genbank:
        if line[0]!='#':
            items = line.split("\t")
            if len(items) == 20:
                currentTaxon  = items[5]
                if currentTaxon in taxon:
                    # Get the FTP link to the genome for the taxa
                    genomeFolder = items[19].strip()
                    if genomeFolder!="":
                        # removing ftp://ftp.ncbi.nlm.nih.gov from beginning of genome folder path
                        genomeFolder = genomeFolder[26:]
                        filename = genomeFolder.split("/")[-1].strip()
                        if filename!="":
                            print("attempting to retrieve " + filename + "_genomic.fna.gz" + " from " + genomeFolder + " for taxon " + currentTaxon)
                            ftp.cwd(genomeFolder)
                            localGenome = open(foldername + filename + "_genomic.fna.gz", 'wb')
                            ftp.retrbinary('RETR %s' % filename + "_genomic.fna.gz", localGenome.write)
                            localGenome.close()
    genbank.close()

def getAllOtherStrainsForTaxon(taxon, taxonFilename, otherStrainsFilename):
    # This file was downloaded from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt
    ftpReference = 'data/assembly_summary_genbank.txt'
    genbank = open(ftpReference,'r')
    taxonOutput = open(taxonFilename, 'w')
    species = {}
    for line in genbank:
        if line[0]!='#':
            items = line.split("\t")
            if len(items) == 20:
                currentTaxon  = items[5]
                if currentTaxon in taxon:
                    names = items[7].split()
                    newSpeciesName = names[0] + " " + names[1]
                    newStrainName = items[7].strip()
                    genomeFolder = items[19].strip()
                    # removing ftp://ftp.ncbi.nlm.nih.gov from beginning of genome folder path
                    genomeFolder = genomeFolder[26:]
                    if len(genomeFolder) > 10:
                        if newSpeciesName in species:
                            species[newSpeciesName][newStrainName] = genomeFolder
                        else:
                            tempDict = {}
                            tempDict[newStrainName] = genomeFolder
                            species[newSpeciesName] = tempDict
                        taxonOutput.write(newStrainName + " " + genomeFolder + "\n")
    taxonOutput.close()
    otherStrains = {}
    otherStrainsOutput = open(otherStrainsFilename, 'w')
    genbank.close()
    genbank = open(ftpReference,'r')
    for line in genbank:
        if line[0]!='#':
            items = line.split("\t")
            if len(items) == 20:
                names = items[7].split()
                if (len(names) >= 2):
                    newSpeciesName = names[0] + " " + names[1]
                    newStrainName = items[7].strip()
                    genomeFolder = items[19].strip()
                    # removing ftp://ftp.ncbi.nlm.nih.gov from beginning of genome folder path
                    genomeFolder = genomeFolder[26:]
                    if len(genomeFolder) > 10:
                        if newSpeciesName in species:
                            if newStrainName not in species[newSpeciesName]:
                                if newSpeciesName in otherStrains:
                                    if newStrainName not in otherStrains[newSpeciesName]:
                                        otherStrains[newSpeciesName][newStrainName] = genomeFolder
                                        otherStrainsOutput.write(newStrainName + " " + genomeFolder+"\n")
                                else:
                                    tempDict = {}
                                    tempDict[newStrainName] = genomeFolder
                                    otherStrains[newSpeciesName] = tempDict
                                    otherStrainsOutput.write(newStrainName + " " + genomeFolder+"\n")
    otherStrainsOutput.close()
    genbank.close()
    returnList = []
    returnList.append(species)
    returnList.append(otherStrains)
    return returnList


# Downloads ten random strains for each species name
def getAllGenomesForTaxon(taxon, foldername):
    if foldername[-1]!='/':
        foldername += "/";
    data = getAllOtherStrainsForTaxon(taxon, "data/species_with_FTP_link.txt", "data/other_strains_with_FTP_link.txt")
    species = data[0]
    otherStrains = data[1]
    ftp = FTP('ftp.ncbi.nlm.nih.gov')     # connect to host, default port
    ftp.login()
    selectedStrainsOutput = open("data/selected_strains_with_FTP_link.txt",'w')
    #select 10 strains from each species (or as many as exist if there are less than 10)
    for name in species:
        selectedStrains = {}
        if len(species[name]) > 10:
            randomStrains = random.sample(species[name].keys(), 10)
            selectedStrains = {k: species[name][k] for k in randomStrains}
        elif len(species[name]) == 10:
            selectedStrains = species[name]
        else:
            selectedStrains = species[name].copy()
            if name in otherStrains:
                if len(otherStrains[name]) <= (10 - len(species[name])):
                    selectedStrains.update(otherStrains[name])
                else:
                    randomStrains = random.sample(otherStrains[name].keys(), 10 - len(species[name]))
                    selectedStrains.update({k: otherStrains[name][k] for k in randomStrains})
        for strain in selectedStrains:
            selectedStrainsOutput.write(strain+" "+selectedStrains[strain]+"\n")
            filenameBase = selectedStrains[strain].split("/")[-1].strip()
            print("attempting to retrieve " + filenameBase + "_genomic.fna.gz from folder " + selectedStrains[strain] + " for species " + strain)
            ftp.cwd(selectedStrains[strain])
            fnaFilename = selectedStrains[strain] + "/" + filenameBase + "_genomic.fna.gz"
            localGenome = open(foldername + filenameBase + "_genomic.fna.gz", 'wb')
            ftp.retrbinary('RETR %s' % fnaFilename, localGenome.write)
            localGenome.close()
            try:
                featureFilename = selectedStrains[strain] + "/" + filenameBase + "_feature_table.txt.gz"
                localGenome = open(foldername + filenameBase + "_feature_table.txt.gz", 'wb')
                ftp.retrbinary('RETR %s' % featureFilename, localGenome.write)
                localGenome.close()
            except ftplib.error_perm:
                localGenome.close()
                os.remove(foldername + filenameBase + "_feature_table.txt.gz")
                gffFilename = selectedStrains[strain] + "/" + filenameBase + "_genomic.gff.gz"
                localGenome = open(foldername + filenameBase + "_genomic.gff.gz", 'wb')
                ftp.retrbinary('RETR %s' % gffFilename, localGenome.write)
                localGenome.close()
    selectedStrainsOutput.close()

### Calling the functions
gi = set()
# A bug in NCBI BLAST made it necessary for me to search the wgs database for
#  bacterial draft genomes separately from the complete genomes
getGI(gi,'./output/draft-genomes-Alignment-HitTable.csv')
print("Extracted " + str(len(gi)) + " gi numbers from wgs output")
getGI(gi,'./output/complete-genomes-Alignment-HitTable.csv')
print("Extracted " + str(len(gi)) + " gi numbers from wgs and complete genome output")
taxon = getAllTaxon(gi, "./data/taxon.txt")
print("Output all taxon names")
getAllGenomesForTaxon(taxon,'data/genomes')
print("Complete")
