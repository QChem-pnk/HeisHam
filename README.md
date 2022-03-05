##HeisHam

This is a Fortran program to obtain the energies of a system with a Heisenberg Dirac van Vleck Hamiltonian. This program also assigns the S number to each energy.

# Compiling and running the code

Compiling the code is really easy. A makefile is included
in the main directory to compile the different source files
included in the src folder. This makefile contains several
rules, with the following commands:

- `make`: default compiling.
- `make clean`: delete all compiled files.
- `make print`: print info.
- `make debug`: execution of rules clean and print and
compilation with debug mode. This mode incorporates
several flags to the compiler, used for debugging.

For the purpose of using the program, the command make
should be used. This will create an executable in the bin
folder.

For the execution of the program, a bash script is also
incorporated. The command to use this script is:

`./run.sh input [output]`

The input name can be provided with or without extension.
For example, to run file s1.inp, the command can be
`./run.sh s1` and the output files will be s1.out and
s1.resume.out. Another location apart from inp can also
be provided.

If no output name is given, the program will create by
default two output files in the ./out folder, one with the
full output, and one with a resume of the results, with
extensions .out and .resume.out respectively. If a name
is provided, the full output will have that name and the
resume filename will be based on that.

Also, the program can be run directly by executing the
resulting compiled code with the following command:

`./bin/magnetic.exe INPUT [[-o] [OUTPUT]]`

The option `-o` can be used alone and let the program decide
the name of the output file, or an output name can be
provided. Without either of these options, only the resume
output will be created and the full results will be printed in
the standard input.

In resume, the commands for running it easily are:

`make`
`./run.sh s1`

#Format of the input

The format of the input is the following:
```
F # T for verbose mode , for debugging 
4 # Number of sites
3 # Number of couplings
1 2 −10.d0 # Site 1 is coupled with site 2 by −10
2 3 −5.d0 # Site 2 is coupled with site 3 by −5
3 4 −1.d0 # Site 3 is coupled with site 4 by −1
 ( whatever the unit : cm−1, meV, . . . )
```

Input files for all the systems in the homework are named
s1 to s10. Input files for all the systems in the Winter
School presentation are also included, a1 to a10.
