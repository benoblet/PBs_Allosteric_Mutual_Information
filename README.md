# PBs_Allosteric_Mutual_Information
Master 2 BioInformatcis short project repository for Mutual Information determination with Protein Blocks (PBs) structural alphabet.
The aim of this project is to study allosteric communication in proteins.

## Content
The script allows to read pdb files in a folder or trajectory and topology files, thanks to [PBxplore](https://github.com/pierrepo/PBxplore) program and **PBassign** in particular. Alternatively, **PBassign** output can be used as input if already done.
After aminoacids conversion into Protein Blocks structural alphabet sequences for all sequences found in a folder or a trajectory, they are processed to obtain the mutual information matrix.
There are 3 to 4 outputs:
- PBassign output in PB.fasta format,
- parsed PBs sequences, without undetermined PBs, in plain text format,
- Mutual Information matrix (MI) as a tsv file (dot as decimal separator),
- Heatmap of the obtained MI with protein position if *-faa* option is used.

## Installation

```bash
$ conda env create -f binder/environnement_menu.yml
```
Alternatively, replace `conda` by if `mamba` is installed on your computer.
This environnment is ready to use `PBlocalMutInfo.py`.

Please consider that under Windows Powershell, an error seems to come from PBxplore: **PBassign** reads sequences in a folder twice. This is not observed under Ubuntu 20.04 system (even installed thanks to Windows Subsystem for Linux, WSL).

## Basic usages

To read from a folder containing pdb files:
```bash
$ python PBlocalMutInfo.py -p path/to/pdbfolder -on "prefix_for_files"
```

To read from **PBassign** output file, select only 1 frame every 10 frames, set first aminoacid number (protein rank of the first aminoacid used in MD) and display some information (including MI computational time) and follow program progression:
```bash
$ python PBlocalMutInfo.py --pbfile path/to/PBassignFile.PB.fasta \
                           --outputname "prefix_for_files" \
                           --step 10 --firstaminoacid 1824 \
                           --verbose
```

For more details, get help details with:
```bash
$ python bin/PBlocalMutInfo.py --help
```

## Context details 
This code and documentation were developped in more or less one week, some issues may have not been tested considering lack of time. Deadline for code and manuscript submissions is Tuesday, September 14th at 12:00 am.  

As Molecular Dynamics(MD) and Mathematics are not my specialities (coming from Molecular Biology and Genomics), please consider results with care.  
Nonetheless, I used in my script 2 functions for mutual information computation from C-coded scripts coming from [GSA tools](https://github.com/AllosterIt/GSAtools), with Java programming knowledge and internet.

Graphical representation is based on *heatmap()* function found in [matplolib tutorial for heatmap](https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html?highlight=heatmap).


## License
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation.