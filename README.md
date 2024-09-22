# Generate-Kraken2-report
This script allows you to generate Kraken2-style reports.
Your initial file should be in a tsv format and should contain 2 columns: NCBI taxid and Read amount for a given taxid. 
You can specify column names in a command line.

Also, in a code add path to directory with <em>names.dmp</em> and <em>nodes.dmp</em> files from NCBI taxdmp. 

Example of input file
| taxid | scaffold_count |
| ------ | ------ |
| 1 | 15555 |
| 9606 | 1254|
| 131567 | 7 |


```sh
python generate_kraken_report.py -i input.tsv -o kraken_report.tsv
```
