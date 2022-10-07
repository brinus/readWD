<a href="https://brinus.github.io/readWD/">
  <img align="center" src="https://github.com/brinus/readWD/actions/workflows/docs.yml/badge.svg" />
</a>

# readWD
Cpp script to read raw data from DRS Evaluation Board and WaveDream Board. To run the script you must first compile it:
```
 $ make
 $ ./readWD
```
You can type for help:
```
$ ./readWD -h

***********************************************************************
*                                                                     *
*           Welcome to the WaveDREAM data analysis tool               *
*                  Version: 0.9.4 Beta (22.06.2020)                   *
*                                                                     *
***********************************************************************

 Use to read data from binary WaveDREAM output files to root files

 ./readWD [args] [-o {outputPath}] file1.dat [file2.dat ... fileN.dat]

 Reads the files file1.dat ... fileN.dat and creates a root file for each

   -o {outputPath}        set the path where the rootfiles are created
                          default is ./

 Possible [args] are : 
   -d, --DEBUG            print debug information
   -f, --ForceOverwrite   enforce the overwriting of existing rootfiles
   -h, --help             print this help text and quit
   -p, --pos_wf           use for positive waveforms
   -c, --config           enter configuration interface, requires user input
   -s, --subtractNoise    subtract sine noise from waveforms
```

