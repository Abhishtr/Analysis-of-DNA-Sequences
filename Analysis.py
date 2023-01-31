#importing required modules
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#defining a function to analyse a single sequence of DNA
def Analysis(DNA_record):
	
	#Printing general information aboout the sequence such as ID, length and DNA Sequence.
	print("ID: %s\n" % (DNA_record.id))
	print("length: %d\n" % (len(DNA_record.seq)))
	print("DNA: %s\n" % (DNA_record.seq))
	
	#Printing the GC Content using the gc_fraction module
	print("GC Content: %0.3f\n" %(gc_fraction(DNA_record.seq)))
	
	#printing the complement and reverse complement of the DNA sequence
	print("Complement: %s\n" % (DNA_record.seq.complement()))
	print("Reverse Complement: %s\n" % (DNA_record.seq.reverse_complement()))
	
	#Transcribing the DNA to RNA
	print("RNA: %s\n" % (DNA_record.seq.transcribe()))
	
	#Checking whether the DNA Sequence is unambiguous
	if str.count(str(DNA_record.seq),'N') == 0:
		#Translating the DNA to protein
		protein_sequence = DNA_record.seq.translate(to_stop=True)
		print("Protein: %s\n" % (protein_sequence))
		
		#Doing protein Analysis and getting molecular weight, aromaticity, instability index and isoelectric point of protein
		ps = ProteinAnalysis(protein_sequence)
		print("Molecular Weight of Protein: %0.3f\n" % (ps.molecular_weight()))
		print("Aromaticity of Protein: %0.3f\n" % (ps.aromaticity()))
		print("Instability Index of Protein: %0.3f\n" % (ps.instability_index()))
		print("Isoelectric Point of Protein: %0.3f\n" % (ps.isoelectric_point()))
		
		#Printing the count and percentage of amino acids in the protein
		print("Amino Acids, their count and percentage:")
		for i in ps.get_amino_acids_percent():
			print("%s : %d, %0.3f " % (i,ps.count_amino_acids()[i],ps.get_amino_acids_percent()[i]))

if __name__ == "__main__":
	
	#Reading the fasta file
	DNA_Records = list(SeqIO.parse("Oikopleura_transcripts_reference_v1.0.fa.txt", "fasta"))
	
	#taking the input for choice
	choice = int(input("Type \'1\' for Analysis for one DNA Sequence and any other number for analysis of all sequences :"))
	print("");
	
	if choice == 1 :
		Analysis(DNA_Records[0])
	else:
		for DNA_record in DNA_Records:
			Analysis(DNA_record)
