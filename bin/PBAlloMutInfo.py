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
    
    # Communication between command line and Python program
    parser = argparse.ArgumentParser()
    
    ## Positional arguments
    parser.add_argument("pdbfolder",
                        help="path to the folder containing trajectory \
                        pdb files.")
    
    ## Additionnal arguments                    
    parser.add_argument("-s", "--step",
                        help="step to select frames for allosteric study, \
                            one file will be kept every s files.")
    parser.add_argument("-o", "--outputfolder",
                        help="path for the results folder writing. Defaut, \
                            same level as pdb inputfolder.")
    parser.add_argument("-v", "--verbose",
                        help="make the script verbose",
                        action="store_true")
    
    # Arguments retrieval
    args = parser.parse_args()
    
    ## Some checks before going on
    if args.pdbfolder is None:
        sys.exit("Program needs at least path to folder contaning pdb files")
    if not os.path.exists(args.pdbfolder):
        sys.exit("Folder not found, please check it or give absolute path")
    if len(os.listdir(args.pdbfolder)) == 0:
        sys.exit("No file found inside input folder")


    if args.outputfolder is None:
        sys.exit("Program needs at least path to folder containgn pdb files")
    elif os.path.exists(args.outputfolder):
        print(f"Results folder is {args.outputfolder:s}")
    else:
        sys.exit("Destination folder not found")



    # Let's use PBassign
        ## TODO
        
    
    # File parsing
    fastafile = "/Users/bened/Documents/M2BI/1_Programmation_et_Gestion_Projets/Project/WT_allFrames.PB.fasta"
    
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
    
                