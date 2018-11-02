#!/usr/bin/python

import argparse
import sys

# Declare tuple of UMIs
UMIList = ("AACGCCAT","AAGGTACG","AATTCCGG","ACACAGAG","ACACTCAG","ACACTGTG",
           "ACAGGACA","ACCTGTAG","ACGAAGGT","ACGACTTG","ACGTCAAC","ACGTCATG",
           "ACTGTCAG","ACTGTGAC","AGACACTC","AGAGGAGA","AGCATCGT","AGCATGGA",
           "AGCTACCA","AGCTCTAG","AGGACAAC","AGGACATG","AGGTTGCT","AGTCGAGA",
           "AGTGCTGT","ATAAGCGG","ATCCATGG","ATCGAACC","ATCGCGTA","ATCGTTGG",
           "CAACGATC","CAACGTTG","CAACTGGT","CAAGTCGT","CACACACA","CAGTACTG",
           "CATCAGCA","CATCGTTC","CCAAGGTT","CCTAGCTT","CGATTACG","CGCCTATT",
           "CGTTCCAT","CGTTGGAT","CTACGTTC","CTACTCGT","CTAGAGGA","CTAGGAAG",
           "CTAGGTAC","CTCAGTCT","CTGACTGA","CTGAGTGT","CTGATGTG","CTGTTCAC",
           "CTTCGTTG","GAACAGGT","GAAGACCA","GAAGTGCA","GACATGAG","GAGAAGAG",
           "GAGAAGTC","GATCCTAG","GATGTCGT","GCCGATAT","GCCGATTA","GCGGTATT",
           "GGAATTGG","GGATAACG","GGCCTAAT","GGCGTATT","GTCTTGTC","GTGATGAG",
           "GTGATGTC","GTGTACTG","GTGTAGTC","GTTCACCT","GTTCTGCT","GTTGTCGA",
           "TACGAACC","TAGCAAGG","TAGCTAGC","TAGGTTCG","TATAGCGC","TCAGGACT",
           "TCCACATC","TCGACTTC","TCGTAGGT","TCGTCATC","TGAGACTC","TGAGAGTG",
           "TGAGTGAG","TGCTTGGA","TGGAGTAG","TGTGTGTG","TTCGCCTA","TTCGTTCG")

# List of chromosomes encountered 
chrs = []

def get_arguments():
    '''Takes one input file'''
    parser = argparse.ArgumentParser(description = "Deduper Inputs")
    parser.add_argument("-f", "--file", help= "Designates the Input File", required=True, type=str)
    parser.add_argument("-p", "--paired", help= "(Unsupported) Designates if paired end", required=False, type=str)
    parser.add_argument("-u", "--umi", help= "(Unsupported) Designates the UMI list", required=False, type=str)
    return parser.parse_args()

def TestUMI(UMI):
    '''Tests if UMI contains 8nt, no N's, and is in list of 96 UMIs'''
    if "N" not in UMI and len(UMI) == 8 and UMI in UMIList:
        return True
    else:
        return False
    
def ForOrRev(FLAG):
    '''Determines from FLAG if the read is on the forward strand (True) or reverse strand (False)'''
    FLAG = int(FLAG)
    Forward = True
    if (FLAG & 16 == 16):
        Forward = False
    return Forward

def AdjustPosFor(CIGAR, POS):
    '''Adjust the position of the alignment based on the cigar for the forward strand'''
    
    # NOTE: No I, D, or N will be reached before the first M
    # This keeps track of when the M is reached when reading through the cigar string
    M_Reached = False  
    POS = int(POS)
    for i in range(len(CIGAR)):
        if M_Reached == True:
            break
        elif CIGAR[i] == "M":
            M_Reached = True
        elif CIGAR[i] == "S":
            POS -= int(CIGAR[:i])
        elif CIGAR[i] == "H":
            POS -= int(CIGAR[:i])
            
    return POS

def AddSpaces(CIGAR):
    '''Add spaces after characters to allow for splitting'''
    # Determine if each character is present
    CIGAR = CIGAR.replace("M", "M ")
    CIGAR = CIGAR.replace("I", "I ")
    CIGAR = CIGAR.replace("N", "N ")
    CIGAR = CIGAR.replace("D", "D ")
    CIGAR = CIGAR.replace("S", "S ")
    CIGAR = CIGAR.replace("H", "H ")
    return CIGAR

def AdjustPosRev(CIGAR, POS):
    '''Adjust the position of the alignment based on the cigar for the reverse strand'''
    
    # This keeps track of when M is reached when reading through the cigar string
    M_Reached = False
    
    POS = int(POS)
    ADJ = 0
    
    CIGAR = AddSpaces(CIGAR)
    CIGAR = CIGAR.split()
    for pos in range(len(CIGAR)):
        if "M" in CIGAR[pos]:
            ADJ += int(CIGAR[pos].replace("M", ""))
            M_Reached = True
        elif "N" in CIGAR[pos]:
            ADJ += int(CIGAR[pos].replace("N", ""))
        elif "D" in CIGAR[pos]: 
            ADJ += int(CIGAR[pos].replace("D", ""))
        elif "S" in CIGAR[pos] and M_Reached == True:
            ADJ += int(CIGAR[pos].replace("S", ""))
        elif "H" in CIGAR[pos] and M_Reached == True:
            ADJ += int(CIGAR[pos].replace("H", ""))
    
    POS += ADJ
    return POS    

        
def main():
    '''Removes PCR duplicates from SAM file, assuming SAMTOOLS sort has been used first'''
    
    # Initialize forward and reverse dictionaries for duplicates
    chr_f = {}
    chr_r = {}
    
    # Get file arguments
    args = get_arguments()
    File = args.file
    
    if args.paired is not None or args.umi is not None:
        print("Hey friend! That's not supported. You can literally only use the -f flag which is required.")
        sys.exit()
   
    In_File = open(File, "r")
    
    # Create File Names
    deduped = File.replace(".sam", "") + "_Deduped.sam"
    
    # Open Files
    Deduped = open(deduped, "w")
    
    # Read all the lines in file
    while True:
        ToWrite = In_File.readline()
        Line = ToWrite.strip("\n")
        if not Line:
            break
        if Line.startswith("@") != True:
            Line = Line.split()
            
            # Get the variables for each line and test UMI
            UMI = Line[0].split(":")[7]
            if TestUMI(UMI) == True:
                CHR = Line[2]
                POS = Line[3]
                CIGAR = Line[5]

                # Purge both dictionaries if new chromosome
                if CHR not in chrs:
                    chrs.append(CHR)
                    chr_f = {}
                    chr_r = {}

                # Determine if the read is forward or reverse
                Forward = ForOrRev(Line[1])

                # Adjust the position, write to appropriate file, and keep record in dictionary
                if Forward == True:
                    POS = AdjustPosFor(CIGAR, POS)
                    key = str(CHR) + "_" + str(UMI) + "_" + str(POS)
                    if key not in chr_f:
                        chr_f[key] = 1
                        Deduped.write(ToWrite)             

                elif Forward == False:
                    POS = AdjustPosRev(CIGAR, POS)
                    key = str(CHR) + "_" + str(UMI) + "_" + str(POS)
                    if key not in chr_r:
                        chr_r[key] = 1
                        Deduped.write(ToWrite)    

    Deduped.close()
            
main()