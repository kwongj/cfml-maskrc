# cfml-maskrc
Masks recombinant regions in an alignment based on ClonalFrameML output  
Option to draw SVG of recombinant regions

## Author
Jason Kwong (@kwongjc)  
GitHub: [kwongj](https://github.com/kwongj)  

## Dependencies
* [Python 3.x](https://www.python.org/downloads/)
* [BioPython](http://biopython.org/wiki/Main_Page)
* [ete3](http://etetoolkit.org/)
* [svgwrite](https://pypi.python.org/pypi/svgwrite/)

## Usage
`$ cfml-maskrc.py -h`  
```
usage: 
  cfml-maskrc.py --aln FASTA --out OUTFILE <CFML_PREFIX>

Script to mask recombination from ClonalFrameML output and draw SVG of recombinant regions

positional arguments:
  CFML_PREFIX          prefix used for CFML output files (required)

optional arguments:
  -h, --help           show this help message and exit
  --aln FASTA          multiFASTA alignment used as input for CFML (required)
  --out OUTFILE        output file for masked alignment (default="maskrc.aln")
  --symbol ?           symbol to use for masking (default="?")
  --regions FILE       output recombinant regions to file
  --svg FILE           draw SVG output of recombinant regions and save as specified file
  --svgsize WIDExHIGH  specify width and height of SVG in pixels (default="800x600")
  --svgorder FILE      specify file containing list of taxa (1 per line) in desired order
  --svgcolour COLOUR   specify colour of recombination regions in HEX format (default=black)
  --consensus          add consensus row of recombination hotspots
  --version            show program's version number and exit
```

**Requires:**
* Output from ClonalFrameML (specify prefix used to name output files)
* MultiFASTA alignment used as input for ClonalFrameML

**Options:**
* Specify output file using `--out OUTFILE`
* Specify symbol to use in the alignment for masking `--symbol ?`
* Save tab-separated file of recombinant region coordinates `--regions FILE`
* Draw SVG of recombinant regions `--svg FILE`
* Specify size of SVG in pixels eg. 800x600 `--svgsize WIDTHxHEIGHT`
* Specify desired order of taxa in SVG `--svgorder FILE`
* Specify colour to show extant recombination `--svgcolour COLOUR`
* Add consensus row of recombination hotspots to SVG `--consensus`

## Bugs
Please submit via the GitHub issues page: [https://github.com/kwongj/cfml-maskrc/issues](https://github.com/kwongj/cfml-maskrc/issues)  

## Software Licence
GPLv2: [https://github.com/kwongj/cfml-maskrc/blob/master/LICENSE](https://github.com/kwongj/cfml-maskrc/blob/master/LICENSE)

## Other links
* [ClonalFrameML](https://github.com/xavierdidelot/clonalframeml)
