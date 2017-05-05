Performing your own simulations (this is essentially the same text as
given in the htm-file 'Introduction / readme'

PC-DOS programs are supplied that allow repetition of the simulations.
In these simulations, parameters such as the production, decay and
diffusion rates can be changed. The results can be visualized in
different display modes. Guided tours are available that allow to
perform complex simulations only by pressing the <RETURN> key. These
programs have been written in Power-Basic for the DOS mode. This compiler
can be obtained from http://www.powerbasic.com/products/pbdos/. The
slightly modified programs SPMS.BAS and XYMS.BAS can be modified and
executed under Microsoft Professional Basic 7.0 (QBX) or QB 4.5.

The programs should run on any PC. MAC computers require a program that
allows a PC emulation. The programs can started from the CD, e.g., by
the explorer and a click to the program names given below. The programs
are in the subdirectory DOS-SIMU. Go to that subdirectory, e.g., with
the Explorer and start the one of programs listed below (.exe files)

SP      (this program that allows the most frequently used space-time
plots)

XY       (this program allows simulations in a two-dimensional field of
cells)

OLIVA     (this is a simplified program for Oliva porphyria like Fig.
6.3, with many comments)


For those in hurry: a quick start to the computer program

Most figure captions in this book contain commands to reproduce the
corresponding simulations. Such commands are prefixed with the letter S
(Simulation) followed by the figure number. For instance, S61 will
produce the simulation shown in Figure 6.1. Equations not explicitly
given in the book may be found by reading in the parameters of the
corresponding simulation (for instance, r61) and using the command PE
(PrintEquation). In most cases the equations are also provided together
with the animated simulations on the CD-ROM.

Q <RETURN> causes program termination. In the hopefully never occurring
case of a hang up, try to finish the program by <Ctrl><break>. A list of
the most frequently used commands is given on the initial screen. Other
commands and parameter changes can be selected from a menu by F1. The
initiation screen is available by pressing the <ESC> or <RETURN> key
once or twice.

Each parameter has a name. The names are the same as those used in the
equations. Entering its name may change the value of a parameter. For
example, type DA (or da) to change the value of the diffusion constant
of the activator. The computer will first display the current value. You
can either type a new value or <ESC> to leave the parameter unchanged.
Pressing <RETURN> will set the parameter to zero. Press F1 for help and
menu to see the list of parameters and their significance.

The first 14 parameters (KT, KP, ...) are integer numbers and are used to
control the program flow. They include number of iterations, the number
of iterations between updating the display, type of display, the type of
equation or initial conditions (see F1 and 'Integer parameters').

No parameters are given in the book. For those how wish to have a
printed list of the parameters, a file PARLIST.DOC exists. A printout of
this file would generate the corresponding list in the same format as
they appear in the program. It contains, in addition, the program code
for the interactions actually used. This code is given only once and
precedes the set of parameters in which the interaction is first
employed. This list as well as the program codes can be also reached
from the index page of the animated simulations.

Simulations that require more complex inputs from the keyboard have been
simplified using "GUIDED TOURS". In this case, the input is read from a
file step by step by pressing the <RETURN> key. For instance, the
command GT32 leads to the simulations shown in Figure 3.2. A list of the
GUIDED TOURS can be obtained by the command GT.

An additional minimal program OLIVA is included that contains many
comments. It reproduces roughly Fig. 6.3 of the book. This may be
convenient to study the general structure of the programs. In contrast
to the main programs SP and XY, it runs also with the QBASIC-
interpreter, but this can be very slow. The program codes can be
inspected with the Animated Simulations by the corresponding links
from the index page (chapter 10)


Possible problems in the DOS mode:

On some systems, DOS programs may no longer work. PerfectSync provides a
program called DOSBOX that allows you to run DOS programs on a Windows 32
environment. A shareware version is supplied on this disk. A registered
version can be obtained from http://perfectsync.com/DevelopmentTools.htm
(it starts and runs somewhat faster). Instead of starting the program by
commands SP or XY, go o the directory DOS-SIMU and use one of the
following commands (e.g., via   START /  EXECUTE)

dosbox sp
dosbox xy
dosbox oliva

No installation required:

The animated simulations and the programs maybe started from the CD.
For more frequent use it is convenient to copy the CD onto the hard
disk. Make a new directory of your choice (for instance 'SP'), go into
that directory and copy the contents of the CD to this directory with
all (!) subdirectories.


A word of caution...

This software was written with care but the customer knows that software
cannot be produced without errors. Springer-Verlag and the author
disclaim any liability for software errors or their consequences.
Exclusively the author holds all user rights to the software.
Springer-Verlag does not originate data and programs but makes them
available. Any commercial use or distribution independently of the book
is not permitted. [(C) 1995-2002 Hans Meinhardt].



