This is a NOMAD parser for [Wannier90](http://www.wannier.org/). It will
read Wannier90 input and output files and provide all information in NOMAD's unified
Metainfo based Archive format.

For Wannier90 please provide at least the files from this table if applicable to your
calculations (remember that you can provide more files if you want):

| Input Filename | Description |
| --- | --- |
| `*.wout` | output text file w/ arbitrary name |
| `*.win` | input text file |
| `*_hr.dat` | hopping matrices (written if write_hr *.win is true) |