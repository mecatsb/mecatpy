# mecatpy

## requirements

* Python 2.7
* KEGG data (raw text files of KEGG REACTION, COMPOUND, RCLASS and ENZYME entries) 
* graph-tool
* numpy
* matplotlib


## KEGG data
the KEGG raw data (not included) can be parsed with the script KEGG\KEGGreader.py

## models
The models used in the publication are contained in Models (without KEGG datasets).
Further models can be built by adding the KEGG organism abbreviation to the file organisms.txt and running the script MECAT\buildHostModels.py (parsed KEGG data required)
