# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 16:55:01 2021

@author: Benedicte
"""


# librairies
import sys  # to exit when issues
import os   # to check for folder and files existence
import argparse  # arguments retrieval

import numpy as np   # arrays
import pandas as pd  # data frame

import time   # evaluate time, not used so far










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
                        help="path of a directory containing pdb files")
    parser.add_argument("-x", "--trajectory",
                        help="name of the trajectory file, \
                            if used topology file required")
    parser.add_argument("-t", "--topology",
                        help="name of the topology file, \
                            if used trajectory file is required")
    
    ## Complementary arguments                    
    parser.add_argument("-o", "--outputfolder",
                        help="path for the results folder writing. Defaut, \
                            same level as pdb inputfolder")
    parser.add_argument("-on", "--outputname",
                        help="Prefix for results files", required = True)

    ## Mutual Information arguments
    parser.add_argument("-s", "--step",
                        help="step to select frames for allosteric study, \
                            one file will be kept every s files")
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
        print(f"Results folder will be written in {args.outputfolder:s}.")
            
        
    # Let's use PBassign
    if args.verbose:
        print("Calling for PBassign function")
    
    if args.pdbfolder:
        assigned_o = os.path.join(args.outputfolder, args.outputname, )
        os.system(f'PBassign -p {args.pdbfolder:s} -o {assigned_o:s}')
    if args.trajectory and args.topology:        
        assigned_o = os.path.join(args.outputfolder, args.outputname, )
        os.system(f'PBassign -x {args.trajectory:s} -g {args.topology} \
                  -o {assigned_o:s}')
                  
    
    # File parsing
    projectfolder = "/Users/bened/Documents/M2BI/1_Programmation_et_Gestion_Projets/Project/"
    fastafile = "/Users/bened/Documents/M2BI/1_Programmation_et_Gestion_Projets/Project/WT_allFrames.PB.fasta"
    dossiervide="/Users/bened/Documents/M2BI/PBs_Allosteric_Mutual_Information/data"
    outname="monpremieressai"
    m2bifolder="/Users/bened/Documents/M2BI"
    module1 ="1_Programmation_et_Gestion_Projets"
    #os.path.split(fastafile)
    
    pb_sequences = []
    with open(fastafile, 'r') as PBoutput:
        
        # starting with one avoid comparison in first if
        current = ["model", "number", "sequence"]
        
        for line in PBoutput:
            
            if line.startswith(">"):
                
                # store in big list
                pb_sequences.append(current)
                
                # empty information in little list
                current = []
                
                # get sequence info
                model = line.strp().split(sep = " | ")[1]
                current.append(model)
                number = int(model.split(sep = " ")[-1])
                current.append(number)
                current.append("")
            
            else:
                # get PBs sequenc itself
                current[2] += line.strip()
                
    # Insert into a dataframe
    pb_sequences_df = pd.DataFrame(data = pb_sequences,
                                   columns = pb_sequences[0],
                                   index = [row[1] for row in pb_sequences])
    
            