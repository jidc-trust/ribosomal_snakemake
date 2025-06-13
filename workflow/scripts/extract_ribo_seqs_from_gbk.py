#!/usr/bin/python

from Bio import SeqIO
from Bio.Seq import Seq
import re
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import sys
import datetime

ribosomal_proteins = ['L14','L16','L18','L2','L22','L24','L3','L4','L5','L6','S10','S17','S19','S3','S8','L15']
ribo_dic = {'L14': 'rplN', 'L16': 'rplP', 'L18': 'rplR', 'L2': 'rplB', 'L22': 'rplV', 'L24': 'rplX', 'L3': 'rplC', 'L4': 'rplD', 'L5': 'rplE', 'L6': 'rplF', 'S10': 'rpsJ', 'S17': 'rpsQ', 'S19': 'rpsS', 'S3': 'rpsC', 'S8': 'rpsH', 'L15': 'rplO'}

if sys.argv[1].endswith('.gz'):
    filehandle = sys.argv[1][:-3]
    print(filehandle)
else:
    filehandle = sys.argv[1]
handle = open(filehandle, 'r')
params = sys.argv[2]
print("parameters for protein or dna are:", params)
outname = datetime.date.today().strftime("%d-%m-%y")+"extracted.fasta"
records = list(SeqIO.parse(handle, "genbank"))
seen_before = []
sequences = []

#Loop through annotation to look for an annotation match
#Save correct matches to file
ribo_count = 0

for rec in records:
    for feature in rec.features:
        if feature.type == "CDS":
          try:
            product = ''.join(feature.qualifiers['product'])
            if re.search(r'\bribosomal( subunit)? protein\b', product, re.IGNORECASE):
                print(feature.qualifiers['locus_tag'], feature.qualifiers['product'])
                print("ribo_count is", ribo_count, "seen_before", seen_before)
                if ribo_count > 16:
                    with open('strains_missing_ribos.txt', 'a') as f:
                            print("breaking out")
                            f.write("{}\n".format(filehandle))
                    break
                else:
                    elements = ''.join(feature.qualifiers['product']).split()
                    for rib in ribosomal_proteins:
                        if elements[-1] == rib:
                                ribo_count+=1
                                if rib in seen_before:
                                    with open('strains_missing_ribos.txt', 'a') as f:
                                           print("breaking out as seen before")
                                           f.write("{}\n".format(filehandle))
                                           break
                                seen_before.append(rib)
#                          print("annotation match")
                                if params.lower() == 'dna':
                            #DNA parameters so extract the DNA sequence and save to a file with the appropriate ID
                                    newrec = SeqRecord(feature.location.extract(rec.seq), id="{}_{}|{}".format(elements[-1], ribo_dic[elements[-1]], filehandle), description="")
                                    print("length of the translate", len(feature.location.extract(rec.seq).translate(to_stop=True)), "seq length", len(feature.location.extract(rec.seq))//3)
                                    if len(feature.location.extract(rec.seq).translate(to_stop=True)) < (len(feature.location.extract(rec.seq))//3)-1:
                                         pass
                                    else:
                                         sequences.append(newrec)
                                elif params.lower() == 'protein':
                            #protein parameters so extract the protein sequence, (by translating the DNA sequence) and save to file with appropriate ID
                                     newrec = SeqRecord(feature.location.extract(rec.seq).translate(to_stop=True), id="{}_{}|{}".format(elements[-1], ribo_dic[elements[-1]], filehandle), description="")
                                     print('>{}\n{}'.format(newrec.id, newrec.seq))
                                     sequences.append(newrec)
                                else:
                                     print("unknown parameters for DNA or protein, please choose either dna or protein")
          except:
               print("problem {}".format(feature))
outfile = open(outname, 'a+')
#save chosen sequences to file with date stamp

SeqIO.write(sequences, outfile, "fasta")
