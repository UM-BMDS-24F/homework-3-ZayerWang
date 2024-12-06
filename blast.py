'''
Chose blastp as we are querying protein sequences, thus we should be using blastp (protein) instead of blastn (nucleotide)
The set evalue and substitution matrix (1*10^-12 and PAM30) were chosen through testing.
My goal was to output (around) 100 or less of the closest homologous sequences (if only 1 or 2 were very close, then only output those 1 or 2). 
Additionally, I set a requirement that there should be at least one result for every query.

Through testing, I found that some sequences in the query have a lot of similar matches while other has far fewer. 
Depending on the parameters, I could get up to a range of 300 for number of sequences that met the evalue threshold.
This was true for all the BLOSUM matrices I tested, and for any evalue with magnitude value greater than 10^-12.
PAM30 (or PAM70) with a evalue in the magnitude of 10^-13 gave me the desired result of 100 max homologous sequences 

Additionally from a holistic view, such a low evalue and using PAM makes sense. While mice and humans are very different,
genetically they should be relatively similar. Both are mammals that split recently (evolutionarily), 
PAM is usually for closely related (and global) alignments and the lower the evalue the more closely the sequences must match.
Lots of the human and mouse protein sequences should align (globally) and so a substitution matrix that is used for high global alignment makes sense in context.
Since these sequences should be fairly similar, we want to lower the evalue so only the highest of alignments will be presented.
'''

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

#Running on local blast, please modify to so that this program can see whatever local version of blast you are using
blastpLocation = "./ncbi-blast-2.16.0+/bin/blastp"
queryFileName = "human.fa"
#If you only have the mouse.fa file please the the following command to get the database processed
#./ncbi-blast-2.16.0+/bin/makeblastdb -in mouse.fa -dbtype "prot" -out mouse_db
databaseName = "mouse_db"

query_dict = []
for seq_record in SeqIO.parse(queryFileName, 'fasta'):
    query_dict.append(seq_record.id)

blastp_cline = NcbiblastpCommandline(
    cmd = blastpLocation, 
    query = queryFileName, 
    db = databaseName, 
    evalue = 5 * pow(10,-13), 
    matrix = "PAM30", 
    outfmt = 5,
    out = "blast_result.xml"
    )

blastp_cline()

result_handle = open("blast_result.xml")
blast_records = NCBIXML.parse(result_handle)

outputFile = open("results.txt", 'w')

index = 0
for blast_record in blast_records:
    outputFile.write("Alignments for Human Sequence ID: " + str(query_dict[index]) + "\n")
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            split = alignment.title.split("|") 
            outputFile.write("\tMouse Sequence ID: " + str(split[3]) + "\n")
            outputFile.write("\tQuery sequence:          " + str(hsp.query) + "\n")
            outputFile.write("\tCorresponding alignment: " + str(hsp.match) + "\n")
            outputFile.write("\tE-value: " + str(hsp.expect) + "\n")
            outputFile.write("\tBit Score: " + str(hsp.score) + "\n\n")
    
    index += 1


outputFile.close()
