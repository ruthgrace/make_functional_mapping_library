import os
import re

# for test purposes
# featureTableFileName = "~/Downloads/test/GCA_000432675.1_MGS50_feature_table.txt"
# fnaFilename = "~/Downloads/test/GCA_000432675.1_MGS50_genomic.fna"

def outputCodingSequencesFromFeatureTable(featureTableFilename, fnaFilename, outputFilename):
    fna = open(fnaFilename, 'r')
    sequences = {}
    isId = True
    id = ""
    for line in fna:
        if (isId):
            id = line.split()[0]
            # get rid of the ">" at the beginning of each line
            id = id[1:]
            isId = False
        else:
            sequences[id] = line.strip()
            id = ""
            isId = True
    fna.close()
    featureTable = open(featureTableFilename,'r')
    output = open(outputFilename,'w')
    for line in featureTable:
        if line[:3] == "CDS":
            items = line.split("\t")
            seqid = items[6].strip()
            protid = items[10].strip()
            start = items[7].strip()
            end = items[8].strip()
            sequence = sequences[seqid][(start-1):end].strip()
            output.write(">genomic_accession|"+seqid+"|product_accession|"+protid+"\n")
            output.write(sequence+"\n")
    featureTable.close()
    output.close()

def outputCodingSequencesFromGFF(gffFilename, fnaFilename, cdsFilename):
    # FILL THIS OUT !!

def getGenomeNamesPerGenus(strains, strainFilename):
    strainFile = open(strainFilename,'r')
    for line in strainFile:
        line = line.strip()
        strains[line.split()[-1][13:]] = re.sub(" [^ ]*$","",line)
    return strains
    strainFile.close()

strains = {}
strains = getGenomeNamesPerGenus(strains,"data/selected_strains_with_FTP_link.txt")

genus = {}
for filename in os.listdir("data/genomes"):
    if filename.endswith(".txt") or filename.endswith(".gff"):
        currentGenus = strains[filename].split()[0]
        if currentGenus not in genus:
            os.mkdir("data/genomes/" + currentGenus)
            genus[currentGenus] = {}
        filepath = "data/genomes/" + filename
        if filename.endswith(".txt"):
            fnaFilepath = "data/genomes/" + re.sub("_feature_table.txt$","_genomic.fna",filename)
            cdsFilepath = "data/genomes/" + currentGenus + "/" + re.sub("_feature_table.txt$","_coding_sequences.fa",filename)
            outputCodingSequencesFromFeatureTable(filepath, fnaFilepath, cdsFilepath)
            genus[currentGenus][strains[filename]] = cdsFilepath
        elif file.endswith(".gff"):
            fnaFilepath = "data/genomes/" + re.sub(".gff$",".fna",filename)
            cdsFilepath = "data/genomes/" + currentGenus + "/" + re.sub("_feature_table.txt$","_coding_sequences.fa",filename)
            outputCodingSequencesFromGFF(filepath, fnaFilepath, cdsFilepath)
            genus[currentGenus][strains[filename]] = cdsFilepath

genusFile = open("data/genus.txt",'w')
for g in genus.keys():
    genusFile.write(g+"\n")
    filenames = genus[g].values()
    with open("data/genomes/" + g + "/" + g + ".fa", 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
genusFile.close()