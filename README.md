# cfml-maskrc
Masks recombinant regions in an alignment based on ClonalFrameML output

## Author
Jason Kwong (@kwongjc)  
GitHub: [kwongj](https://github.com/kwongj)  

## Dependencies
* [Python 2.x](https://www.python.org/downloads/)
* [BioPython](http://biopython.org/wiki/Main_Page)
* [ete2](http://etetoolkit.org/)

## Usage
`$ cfml-maskrc.py -h`  
```
usage: 
  cfml-maskrc.py --aln FASTA --out OUTFILE --symbol ? <CFML_PREFIX>

Script to mask recombination from ClonalFrameML output

positional arguments:
  CFML_PREFIX    prefix used for CFML output files

optional arguments:
  -h, --help     show this help message and exit
  --aln FASTA    multiFASTA alignment used as input for CFML
  --out OUTFILE  output file for masked alignment (default="maskrc.aln")
  --symbol ?     symbol to use for masking (default="?")
  --version      show program's version number and exit
```

**Requires:**
* Output from ClonalFrameML (specify prefix used to name output files)
* MultiFASTA alignment used as input for ClonalFrameML

**Options:**
* Specify output file using `--out`
* Specify symbol to use in the alignment for masking `--symbol`

##Bugs
Please submit via the GitHub issues page: [https://github.com/kwongj/cfml-maskrc/issues](https://github.com/kwongj/cfml-maskrc/issues)  

##Software Licence
GPLv2: [https://github.com/kwongj/cfml-maskrc/blob/master/LICENSE](https://github.com/kwongj/cfml-maskrc/blob/master/LICENSE)

## Other links
* [ClonalFrameML](https://github.com/xavierdidelot/clonalframeml)
