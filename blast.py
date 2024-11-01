from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

human_dict = []
for seq_record in SeqIO.parse("human.fa", 'fasta'):
    human_dict.append(seq_record.id)

blastpLocation = input("Path to blast program: ")
queryFileName = input("Query file name: ")
databaseName = input("Database Name: ")

blastp_cline = NcbiblastpCommandline(
    cmd = blastpLocation, 
    query = queryFileName, 
    db = databaseName, 
    evalue = 5 * pow(10,-13), 
    matrix = "PAM30", 
    outfmt = 5,
    out = "blast_result.xml"
    )

'''
Chose blastp as we are querying protein sequences, thus we should be using blastp (protein) instead of blastn (nucleotide)
The set evalue and substitution matrix (1*10^-12 and PAM30) were chosen through testing.
My goal was to output around 100 or less of the closest homologous sequences. 
Through testing, it was found that some sequences in the query have a lot of similar matches, 
depending on the parameters I could get up to a range of 300 for number of sequences that met the parameter requirements.
This was true for all the BLOSUM matrices I tested, and for any evalue with magnitude value greater than 10^-12
PAM30 (or PAM70) with a evalue in the magnitude of 10^-13 gave me the desired results

Additionally, such a low evalue and using PAM makes sense, as while mice and humans are very different,
genetivally they should be fairly similar. Both are mammals that split recently (evolutionarily), 
PAM is usually used for higher matching sequences and the lower the evalue the more closely the sequences match.
Thus for such similar sequences, using PAM and a very low evalue makes sense.
'''

blastp_cline()

result_handle = open("blast_result.xml")
blast_records = NCBIXML.parse(result_handle)

outputFile = open("results.txt", 'w')

index = 0
for blast_record in blast_records:
    outputFile.write("Alignments for Human Sequence ID: " + str(human_dict[index]) + "\n")
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