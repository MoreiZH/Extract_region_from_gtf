# Installation
This tool is written with python3.8, the packages used are listed in requirements.txt, install them by command line: *pip install -r requirements.txt*
# Usage
Use the tool by command line: *python extract_region.py -f path/to/gtf/file* <br>
usage: extract_region.py [-h] [-f FILE]<br>

This tool helps to extract region from gtf file<br>

optional arguments:<br>
  -h, --help            show this help message and exit<br>
  -f FILE, --file FILE  gtf file to be extracted
# Output
The output file will be created in current work directory, which include all gene, isoform, coding exons, 5'UTR and 3'UTR coordinates in BED format. Bed file header includes ['seq_id', 'start', 'end', 'gene_id', 'attribute'/'transcript_id']. <br>For 5' and 3' UTR regions, there are two versions 1) Gene level and 2) Isoform level. For gene level include all UTRs of all the isoforms of a gene -merge any overlapping isoform UTRs.



