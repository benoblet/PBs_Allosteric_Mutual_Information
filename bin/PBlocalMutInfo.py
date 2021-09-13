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

import time   # evaluate time for Mutual Information matrix
import matplotlib.pyplot as plt


# a constant
PB_ALPHABET = "abcdefghijklmnop"


# functions
def arguments_retrieval():
    """
    A function to group arguments parsing and make main code easier to read.
    
    Returns
    -------
    args : list
        If no error was raised, valid argument list is returned.
    """    
    
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
                        help="path to the trajectory file, \
                            if used topology file required.")
    parser.add_argument("-t", "--topology",
                        help="path to the topology file, \
                            if used trajectory file is required.")
    parser.add_argument("-pb", "--pbfile",
                        help="path to PBassign output file if already done.")
    
    ## Additional Pxwplore argument with adaptation             
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
    parser.add_argument("-faa", "--firstaminoacid",
                        help="number of the first aminoacid in the complete \
                            sequence protein. Value must be a positive integer.")
    parser.add_argument("-v", "--verbose",
                        help="make the script verbose",
                        action="store_true")
    
    # Arguments list
    args = parser.parse_args()
    
    
    # Some checks before going on
    ## An existing and valid input
    if args.pdbfolder is None and args.trajectory is None \
        and args.pbfile is None:
        sys.exit("Program needs at least path to folder or file(s)")
    
    if args.pdbfolder:
        if not os.path.exists(args.pdbfolder):
            sys.exit("Folder not found, please check it or give absolute path")
        else:
            if not os.listdir(args.pdbfolder):
                sys.exit("No file found inside input folder")
        
    if args.trajectory and not args.topology:
        sys.exit("Program needs both trajectory and topology files")

    if args.pbfile:
        if not os.path.isfile(args.pbfile):
            sys.exit("PB.fasta file not found")

    ## A valid output        
    if not args.outputname:
        sys.exit("Please specify a prefix for your results using -on option")          
    
    if args.outputfolder and os.path.exists(args.outputfolder):
        print(f"Results folder is {args.outputfolder:s}")
    if args.outputfolder is None:
        if args.pdbfolder:
            args.outputfolder = os.path.split(args.pdbfolder)[0]
        elif args.trajectory:
            args.outputfolder = os.path.split(args.trajectory)[0]
        else:
            args.outputfolder = os.path.split(args.pbfile)[0]
        print(f"Results folder will be written in {args.outputfolder:s}")
    
    ## Existing additionnal arguments 
    if not args.step:
        args.step = 1
    else:
        try:
            args.step = int(args.step)
        except:
            sys.exit("Step value is invalid, must be an integer")
    if args.step <= 0:
        sys.exit("Step value is invalid, must be a positive integer")
        
            
    if not args.firstaminoacid:
        args.firstaminoacid = 1
    else:
        try:
            args.firstaminoacid = int(args.firstaminoacid)
        except:
            sys.exit("Aminoacid number is invalid, must be an integer")
    if args.firstaminoacid <= 0:
        sys.exit("Aminoacid number is invalid, must be a positive integer")
            
    # Finished, ready to go on
    return args


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


def get_reduced_frames_df(sequences_dataframe, step):
    """
    Generate a dataframe for only studied time or pdb frames. 
    
    Parameters
    ----------
    sequences_dataframe : pandas dataframe
        Dataframe containing all Protein Blocks sequences, without undefined \
        "Z" values. Rows stand for one pdb file or molecular dynamic frame, \
        columns matches position in structural alphabet protein sequence. The \
        first two undertermined "Z" PBs motifs are considered as 0 and 1 \
        position respectively but do not exist in dataframe.
    step : positive not null integer
        Step between two frames and pdb snapshots. No default value is set as \
        this option is meant to be used with -s or --step value defined by \
        user (args.step).

    Returns
    -------
    partial_df : pandas datarame
        The original dataframe limited to regularly spaced selected frames or \
        pdb files depending on step value.    
    """
    # indexes selection
    reduced = range(1, sequences_dataframe.shape[0] + 1, step)
    # reduced dataframe selection
    partial_df = sequences_dataframe.loc[reduced, ]
    
    return partial_df


def get_combinationseq_in_position(sequences_dataframe, seqposition):
    """
    Retrieve a sequence for defined position in the input dataframe.  
    
    Parameters
    ----------
    sequences_dataframe : pandas dataframe
        Dataframe containing all Protein Blocks sequences, without undefined \
        "Z" values. Rows stand for one pdb file or molecular dynamic frame, \
        columns matches position in structural alphabet protein sequence. The \
        first two undertermined "Z" PBs motifs are considered as 0 and 1 \
        position respectively but do not exist in dataframe.
    seqposition : integer >= 2
        Value to select an aminoacid position. Must be superior or equal to 2 \
        considering dataframe column numerotation.

    Returns
    -------
    thesequence : string.
        The sequence corresponding to PBs transitions for seqposition.    
    """
    # column selection
    wantedcol = sequences_dataframe[seqposition]
    # create a string
    thesequence = "".join(list(wantedcol))
    
    return thesequence


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

    for val1 in freq_dataframe.index:
        for val2 in freq_dataframe.columns:
            if freq_dataframe.loc[val1, val2] != 0:
                mut_info_val += (freq_dataframe.loc[val1, val2]
                * math.log( freq_dataframe.loc[val1, val2] /
                    (freq_dataframe.loc[val1].sum() * freq_dataframe[val2].sum()),
                    2 )
                )
    
    if mut_info_val > 0:
        return mut_info_val
    else:
        return 0
    

def column_mutual_information(sequences_dataframe, pos1, pos2,
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
    alphabet : list
        All available letters in a structural alphabet. By default, parameter \
        set to Protein Blocks (PBs) alphabet.
    
    Returns
    -------
    mut_info_val : float
        Value of Mutual Information for corresponding set of positions using \
        mutual_information() function defined above    
    """

    # retrieve sequences
    seq1 = get_combinationseq_in_position(sequences_dataframe, pos1)
    seq2 = get_combinationseq_in_position(sequences_dataframe, pos2)

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
    freq_dataframe = counts_dataframe / len(seq1)
    mut_info_val = mutual_information(freq_dataframe)  
    
    return mut_info_val


# Function from Matplolib documentation with some slight changes
# https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html?highlight=heatmap
def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=+90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=+90, ha="right",
             rotation_mode="default")

    # Turn spines off and create white grid.
    #ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    #ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar
# End of functions from Matplolib documentation




if __name__ == "__main__":
    
    args = arguments_retrieval()     
    print("Program starting...")  
    
    # Let's check if we already have a PBassign output file
    if args.pbfile:
        fastafile = args.pbfile
        assigned_o = fastafile.replace(".PB.fasta", "")
    else:
        # Let's use PBassign
        if args.verbose:
            print("Calling for PBassign function")
        
        if args.pdbfolder:
            assigned_o = os.path.join(args.outputfolder, args.outputname)
            print(f'PBassign -p {args.pdbfolder:s} -o {assigned_o:s}')
      
        if args.trajectory and args.topology:        
            assigned_o = os.path.join(args.outputfolder, args.outputname)
            os.system(f'PBassign -x {args.trajectory:s} -g {args.topology} \
                      -o {assigned_o:s}')      
        
        fastafile = f"{assigned_o:s}.PB.fasta"             
        
    if not os.path.isfile(fastafile):
        sys.exit("PBassign output file not found: did it display any error or \
                 file extension '.PB.fasta' changed?")
    
    
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

    print(f"{len(pb_sequences):d} sequences lues")
    print(f"{len(pb_frames):d} numéros de frames lus")

    ## Insert into a dataframe
    pb_sequences_df = pd.DataFrame(data = pb_sequences,
                                   columns = range(2, len(pb_sequences[0])+2),
                                   index = pb_frames)    
    pb_sequences_df.sort_index(inplace = True)
    
    if args.verbose:
        print(f"{pb_sequences_df.shape[0]:d} sequences found.")
        print(f"{pb_sequences_df.shape[1]:d} PBs per sequence.")
    
    ## Print sequences into text file
    if args.verbose:
        print("Writting sequences in text file.")
    with open(assigned_o + ".allPBsequences.txt", 'w') as alignoutput:
        for seq in pb_sequences:
            alignoutput.write("".join(seq) + "\n")
    
    # Reduce dataframe size to limited frames
    pb_sequences_df = get_reduced_frames_df(pb_sequences_df,
                                            step = args.step)
    if args.verbose:
        if args.step == 1:
            print("All frames will be considered for Mutual information \
            calculation.")
        else:
            print(f"Only {pb_sequences_df.shape[0]:d} frames will be studied \
            as --step option as been set.")


    # Initialise Mutual Information matrix 
    pos_matrix = np.zeros((pb_sequences_df.shape[1], pb_sequences_df.shape[1]),
                          dtype = float )
    mutinfo_df = pd.DataFrame(data = pos_matrix,
                              columns= pb_sequences_df.columns,
                              index = pb_sequences_df.columns)
    
    
    # Compute Mutual Information dataframe
    starting_time = time.time()
    
    for pb1 in mutinfo_df.index:
        for pb2 in mutinfo_df.columns:
            mutinfo_df.loc[pb1, pb2] = column_mutual_information(pb_sequences_df,
                                                                 pb1, pb2)
            if pb1 != pb2:
                mutinfo_df.loc[pb2, pb1] = mutinfo_df.loc[pb1, pb2]

    ending_time = time.time()
    
    # Display human readable time
    duration_sec = ending_time - starting_time
    duration_min = duration_sec // 60
    duration_remaining_sec = duration_sec % 60
    if args.verbose:
        print(f"It took {duration_min:.0f} minute(s) and \
        {duration_remaining_sec:.2f} seconds to compute MI matrix.")
        
    # Save Mutual Information Matrix (MI) to file
    tsvfilename = assigned_o + f"_MIMatrix_{args.step:d}spacedFrames.tsv"
    if args.verbose:
        print("Writting Mutual information values in tsv file.")
    mutinfo_df.to_csv(os.path.join(args.outputfolder, tsvfilename),
                      sep = "\t", float_format = "%.3f",
                      header = False, index = False)



    # Make graphical MImatrix representation
    
    ## Set space size and axes
    figure, axes = plt.subplots(figsize = [20, 12],
                                dpi = 500)
    
    mutinfo_array = np.array(mutinfo_df)
    image, _ = heatmap(mutinfo_array,
                       list(mutinfo_df.columns + args.firstaminoacid),
                       list(mutinfo_df.index + args.firstaminoacid),
                       cmap = "GnBu", vmin = 0,
                       vmax = math.ceil(mutinfo_array.max()),
                       cbarlabel = "mutual information values")
    
    ## Add a main title and fit nicely graph
    plt.suptitle(t = "Mutual Information Matrix for studied protein sequence \
    in PBs", ha = "center", va = "top", size = 15)
    plt.tight_layout()
    
    ## Save graphical reprensetation (figure) in a file
    #plt.show()
    if args.verbose:
        print("Writting Mutual information Matrix plot in png file.")
    figfilename = assigned_o + f"_MIMatrix_{args.step:d}spacedFrames.png"
    plt.savefig(fname = figfilename, format = "png")

    print("... Program ends")