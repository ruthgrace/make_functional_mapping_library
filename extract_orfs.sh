import os
import re

# for test purposes
# featureTableFileName = "~/Downloads/test/GCA_000432675.1_MGS50_feature_table.txt"
# fnaFilename = "~/Downloads/test/GCA_000432675.1_MGS50_genomic.fna"

def outputCodingSequencesFromFeatureTable(featureTableFilename, fnaFilename, outputFilename):
    fna = open(fnaFilename, 'r')
    sequences = {}
    id = ""
    for line in fna:
        if (line[0]=='>'):
            if id!="":
                sequences[id] = "".join(sequences[id])
            id = line.split()[0]
            # get rid of the ">" at the beginning of each line
            id = id[1:]
            sequences[id] = []
        else:
            if id!="":
                sequences[id].append(line.strip())
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
            sequence = sequences[seqid][(int(start)-1):int(end)].strip()
            output.write(">genomic_accession|"+seqid+"|product_accession|"+protid+"\n")
            output.write(sequence+"\n")
    featureTable.close()
    output.close()

def outputCodingSequencesFromGFF(gffFilename, fnaFilename, cdsFilename):
    fna = open(fnaFilename, 'r')
    sequences = {}
    id = ""
    for line in fna:
        if (line[0]=='>'):
            if id!="":
                sequences[id] = "".join(sequences[id])
            id = line.split()[0]
            # get rid of the ">" at the beginning of each line
            id = id[1:]
            sequences[id] = []
        else:
            if id!="":
                sequences[id].append(line.strip())
    fna.close()
    gff = open(gffFilename,'r')
    output = open(outputFilename,'w')
    for line in gff:
        if line[0]!='#':
            items = line.split("\t")
            if len(items > 5):
                seqid = items[6].strip()
                protid = items[10].strip()
                start = items[7].strip()
                end = items[8].strip()
                sequence = sequences[seqid][(int(start)-1):int(end)].strip()
                output.write(">genomic_accession|"+seqid+"|product_accession|"+protid+"\n")
                output.write(sequence+"\n")
    gff.close()
    output.close()

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
        rootName = re.sub("_genomic.gff","",re.sub("_feature_table.txt","",filename))
        currentGenus = strains[rootName].split()[0]
        if currentGenus not in genus:
            os.mkdir("data/genomes/" + currentGenus)
            genus[currentGenus] = {}
        filepath = "data/genomes/" + filename
        if filename.endswith(".txt"):
            fnaFilepath = "data/genomes/" + re.sub("_feature_table.txt$","_genomic.fna",filename)
            cdsFilepath = "data/genomes/" + currentGenus + "/" + re.sub("_feature_table.txt$","_coding_sequences.fa",filename)
            outputCodingSequencesFromFeatureTable(filepath, fnaFilepath, cdsFilepath)
            genus[currentGenus][strains[rootName]] = cdsFilepath
        elif file.endswith(".gff"):
            fnaFilepath = "data/genomes/" + re.sub(".gff$",".fna",filename)
            cdsFilepath = "data/genomes/" + currentGenus + "/" + re.sub("_feature_table.txt$","_coding_sequences.fa",filename)
            outputCodingSequencesFromGFF(filepath, fnaFilepath, cdsFilepath)
            genus[currentGenus][strains[rootName]] = cdsFilepath

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