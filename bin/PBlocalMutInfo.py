# -*- coding: utf-8 -*-
""" Mutual information matrix based on Protein Blocks (PBs)

Calculator and graphical representation on Mutual Information (MI) matrix based
on PBxplore PBs assignment and GSA-tools open-source functions for matrix
calculations.
The first and last 2 aminoacids (aa), thus undetermined ("Z") in PB structural 
alphabet are discarded for analysis. Please, consider it to adjust between PB
and aa numbers.


@author: Benedicte Noblet
started on Wed Sep  8 16:55:01 2021

Program produced alone for Master 2 BioInformatics (M2BI) programming training 
course in more or less 1 week. Consider your results carefully.
"""


# librairies

import sys  # to exit when issues
import os   # to check for folder and files existence
import argparse  # arguments retrieval

import math          # logarithm of different numbers
import numpy as np   # arrays
import pandas as pd  # data frame

import time   # evaluate time, not used so far



# a constant
PB_ALPHABET = "abcdefghijklmnop"


# functions
def get_frame_num(fastaheadline):
     """
     Parse header line of a fasta sequence(s) from PBassign
     output to get sequence rank.
     
     Parameters
     ----------
     fastaheadline : string
         Line starting with '>' corresponding to descriptor line
         in the fasta file output of PBassign. It handles both 
         output from a trajectory or a list of pdb files.

     Returns
     -------
     [number] : list(int)
         Lis of length 1 containing frame number only.

     """
     frame = fastaheadline.strip().split(sep = " | ")[1]
     number = int(frame.split(sep = " ")[-1])
     
     return [number]


def mutual_information(freq_dataframe):
    """
    Mutual Information calculation as implemented by GSA tools authors and 
    remembered by Pandini et al., 2012 (FASEB Journal).    
   
    Parameter
    ----------
    frequency_dataframe : a pandas dataframe
        Squared dataframe corresponding to the 16 PBs frequency between two 
        sequence position.

    Returns
    -------
    mut_info_val : float
        Value of Mutual Information for corresponding frequency dataframe

    """
    mut_info_val = 0.0

    # val1 en ligne
    # val2 en colonne
    # voir comment faire avec pandas
    for val1 in freq_dataframe.index:
        for val2 in freq_dataframe.columns:
            if freq_dataframe.loc[val1][val2] != 0:
                mut_info_val += (freq_dataframe.loc[val1, val2]
                * math.log( freq_dataframe.loc[val1, val2] /
                    (freq_dataframe.loc[val1].sum(axis=1) * freq_dataframe[val2].sum()),
                    2 )
                )
    
    if mut_info_val > 0:
        return mut_info_val
    else:
        return 0


def get_reduced_position_string(sequences_dataframe, seqposition, step = 10):
    """
    Retrieve a sequence for defined position with all or selected frames.  
    
    Parameters
    ----------
    sequences_dataframe : pandas dataframe
        Dataframe containing all Protein Blocks sequences, without undefined \
        "Z" values. Rows stand for one pdb file or molecular dynamic frame, \
        columns matches position in structural alphabet protein sequence. The \
        first two undertermined "Z" PBs motifs are considered as 0 and 1 \
        position respectively but do not exist in dataframe.
    seqposition : integer
        Value to select an aminoacid position. Must be superior or equal to 2 \
        considering dataframe column numerotation.
    step : positive not null integer, optional
        Step between two frames and pdb snapshots. The default is set to 10 \
        but this option is aimed to be used with -s or --step value defined \
        by user (args.step).

    Returns
    -------
    thesequence : string.
        The sequence corresponding to PBs transitions for seqposition and \
        some frames depending on step value.    
    """
    # indexes selection
    reduced = range(1, sequences.dataframe.shape[1] + 1, step)
    # column selection
    partialcol = sequences_dataframe.loc[reduced, seqposition]
    # create a string
    thesequence = list(partialcol).join()
    
    return thesequence
    

def column_mutual_information(sequences_dataframe, pos1, pos2, step,
                              alphabet = list(PB_ALPHABET)):
    """
    Calculate Mutual Information for two columns in fasta alignment
    
    Parameters
    ----------
    sequences_dataframe : pandas dataframe
        Dataframe containing all Protein Blocks sequences, without undefined \
        "Z" values. Rows stand for one pdb file or molecular dynamic frame, \
        columns matches position in structural alphabet protein sequence. The \
        first two undertermined "Z" PBs motifs are considered as 0 and 1 \
        position respectively but do not exist in dataframe.
    pos1, pos2 : positive integers, >= 2
        Two protein blocks number positions to compute mutual information \
        between them. pos1 and pos2 can be equal.
    step : positive not null integer
        Step between two frames and pdb snapshots. The default is set to 10 \
        but this option is aimed to be used with -s or --step value defined \
        by user (args.step).
    alphabet : list
        All available letters in a structural alphabet. By default, parameter \
        set to Protein Blocks (PBs) alphabet.
    
    Returns
    -------
    mut_info_val : float
        Value of Mutual Information for corresponding set of positions using \
        mutual_information() function defined above    
    """
    mut_info_val = 0.0

    # retrieve sequences
    seq1 = get_reduced_position_string(sequences_dataframe, pos1, step)
    seq2 = get_reduced_position_string(sequences_dataframe, pos2, step)

    # initialise co-occurence frequency dataframe (to have names in margin ^^)
    alpha_matrix = np.zeros( (len(alphabet), len(alphabet)), dtype = int )
    counts_dataframe = pd.DataFrame(data = alpha_matrix,
                                   columns= alphabet,
                                   index = alphabet)
    
    for tps in range(len(seq1)):
        # frame / row PB motif for sequence in position1
        for l1 in counts_dataframe.index:
            if l1 == seq1[tps]:
                break
            
        # frame / row PB motif for sequence in position2  
        for l2 in counts_dataframe.columns:
            if l2 == seq2[tps]:
                break
        
        # assert sentences removed as undefined PBs ("Z") were discarded
        
        # update co-occurrence in counts dataframe
        counts_dataframe.loc[l1, l2] += 1
    
    # generate co-occurrence frequencies dataframe then compute MI
    
    
    mut_info_val = mutual_information(freq_dataframe)  
    
    return mut_info_val

}

# END Raw extracts from g_sa_analyze.c (GSA-tools script on GitHub)




if __name__ == "__main__":
    
    print("Program starting...")
    
    # Communication between command line and Python program
    parser = argparse.ArgumentParser(
        description='Read PDB structures, assign protein blocs (PBs) \
            using PBxplore then generate alignement and Mutual \
            Information (MI) matrix in plain text files along with MI \
            graphical representation.')
    
    ## PBxplore arguments
    parser.add_argument("-p", "--pdbfolder",
                        help="path of a directory containing pdb files.")
    parser.add_argument("-x", "--trajectory",
                        help="name of the trajectory file, \
                            if used topology file required.")
    parser.add_argument("-t", "--topology",
                        help="name of the topology file, \
                            if used trajectory file is required.")
    
    ## Complementary arguments                    
    parser.add_argument("-o", "--outputfolder",
                        help="path for the results folder writing. Defaut, \
                            same level as pdb inputfolder")
    parser.add_argument("-on", "--outputname",
                        help="Prefix for results files", required = True)

    ## Mutual Information arguments
    parser.add_argument("-s", "--step",
                        help="step to select frames for allosteric study, \
                            one file will be kept every s files. Value must be \
                            a positive integer.")
    parser.add_argument("-v", "--verbose",
                        help="make the script verbose",
                        action="store_true")
    
    # Arguments retrieval
    args = parser.parse_args()
    
    
    ## Some checks before going on
    if args.pdbfolder is None and args.trajectory is None:
        sys.exit("Program needs at least path to folder or file(s)")
    
    if args.pdbfolder:
        if not os.path.exists(args.pdbfolder):
            sys.exit("Folder not found, please check it or give absolute path")
        else:
            if not os.listdir(args.pdbfolder):
                sys.exit("No file found inside input folder")
            
    if args.trajectory and not args.topology:
        sys.exit("Program needs both trajectory and topology files")
    
    if not args.outputname:
        sys.exit("Please specify a prefix for your results using -on option")          
    
    if args.outputfolder and os.path.exists(args.outputfolder):
        print(f"Results folder is {args.outputfolder:s}")
    if args.outputfolder is None:
        if args.pdbfolder:
            args.outputfolder = os.path.split(args.pdbfolder)[0]
        else:
            args.outputfolder = os.path.split(args.trajectory)[0]
        print(f"Results folder will be written in {args.outputfolder:s}")
    
    if not args.step:
        args.step = 1
    else if float(args.step) <= 0 or type(args.step) != int:
        sys.exit("Step value, is invalid must be a positive integer")
        
    # Let's use PBassign
    if args.verbose:
        print("Calling for PBassign function")
    
    if args.pdbfolder:
        assigned_o = os.path.join(args.outputfolder, args.outputname)
        os.system(f'PBassign -p {args.pdbfolder:s} -o {assigned_o:s}')
    if args.trajectory and args.topology:        
        assigned_o = os.path.join(args.outputfolder, args.outputname)
        os.system(f'PBassign -x {args.trajectory:s} -g {args.topology} \
                  -o {assigned_o:s}')           
    
    # for devopment ---
    outputfolder = "/Users/bened/Documents/M2BI/1_Programmation_et_Gestion_Projets/Project"
    outputname = "WT_allFrames"
    assigned_o = os.path.join(outputfolder, outputname)
    # --- for development end
    
    fastafile = f"{assigned_o:s}.PB.fasta"
    
    if not os.path.isfile(fastafile):
        sys.exit("PBassign output file not found: did it display any error?")
    
    
    # File parsing
    pb_sequences = []
    pb_frames = []
    gotone = False
    current = []
    
    with open(fastafile, 'r') as PBoutput:
        
        for line in PBoutput:
            if line.startswith(">"):
                if gotone:
                    # remove non-PBs letters, save time and calculation
                    current[1] = current[1].replace("Z","")  
                    # store in big lists
                    pb_sequences.append(list(current[1]))
                    pb_frames.append(current[0])
                    # and empty for next sequence
                    current = get_frame_num(line) + [""]
                else:
                    current = get_frame_num(line) + [""]
                    gotone = True
            else:
                # get PBs sequence itself
                current[1] += line.strip()


# ajouter un paramètre pour conserver aminoacid number ?            

    # Insert into a dataframe
    pb_sequences_df = pd.DataFrame(data = pb_sequences,
                                   columns = range(2, len(pb_sequences[0])+2),
                                   index = pb_frames)
    pb_sequences_df.sort_index(inplace=True)
    
    if args.verbose:
        print(f"{pb_sequences_df.shape[0]:d} sequences found.")
        print(f"{pb_sequences_df.shape[1]:d} PBs per sequence.")






# Raw extracts from g_sa_analyze.c (GSA-tools script on GitHub)

    /** allocate and initialize Mutual Information matrix */
    MIMat = alloc_float_matrix(MIMat, sequenceLength, sequenceLength);
    initialise_float_matrix(MIMat, sequenceLength, sequenceLength, 0.0);
    eeMIMat = alloc_float_matrix(eeMIMat, sequenceLength, sequenceLength);
    initialise_float_matrix(eeMIMat, sequenceLength, sequenceLength, 0.0);
    JentropyMat = alloc_float_matrix(JentropyMat, sequenceLength, sequenceLength);
    initialise_float_matrix(JentropyMat, sequenceLength, sequenceLength, 0.0);

    /*________________________________________________________________________*/
    /** allocate and initialize probability matrix */
    initialize_probability_matrix(&probMat, &iCodeSet, &jCodeSet);

    /*________________________________________________________________________*/
    /** calculate Mutual Information */
		int zz;
		int completion_i = 0;
		double completion;

		MPI_Barrier(MPI_COMM_WORLD);

		if (my_rank == 0) fprintf(stdout, "\nMutual Information and Joint Entropy\n");
		fflush(stdout);
    for(i = 0, zz = 0, completion_i = 0; i < sequenceLength; ++i) {
        for(j = i; j < sequenceLength; ++j) {
			   /* print progress */
				++ zz; 
				completion = (long double)zz / ((long double)(sequenceLength*(sequenceLength-1)) / 2) * 100;
				if ((int)completion > completion_i) {
					completion_i = (int)completion;
					if (my_rank == 0) {
						fprintf(stdout, "\t%3d%%\r", completion_i);
						fflush(stdout);
					}
				}

            MIMat[i][j] = column_mutual_information(&inputSequenceSet, &probMat, i, j, &eeMI);
            eeMIMat[i][j] = eeMI;
            if (i != j) {
                MIMat[j][i] = MIMat[i][j];
                eeMIMat[j][i] = eeMIMat[i][j];
            }
        }
    }


        /*________________________________________________________________________*/
    /** output Mutual Information matrices */
    /** output MI */
    if (my_rank == 0) {
        MIFile = safe_open(opt2fn("-MImat", asize(fnm), fnm), "w");
        for(i = 0; i < sequenceLength; ++i) {
            for(j = 0; j < sequenceLength; ++j) {
                fprintf(MIFile, "%8.3f", MIMat[i][j]);
            }
            fprintf(MIFile, "\n");
        }
        fclose(MIFile);
    }

# END Raw extracts from g_sa_analyze.c (GSA-tools script on GitHub)