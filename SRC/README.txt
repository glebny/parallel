To compile & run the program on Windows,
do the following:

COMPILATION
1) Find make.cmd file among those in the folder;

2) Click twice on that file.
This will run the compilation script
in the command line;

3) Press Enter when you see 'pause' in the line.
This means compilation finished successfully;


RUNNING
4) Open Windows PowerShell.
It can be found by search in your Start menu;

5) Print 'cd <path to SRC folder>'
This will lead you to the folder your code is located;

6) Print 'ls' and make sure there is a prog.exe file
in the folder now;

7) Print 'set OMP_NUM_THREADS=2' or '...=4' or else
to set number of cores you want to use for calculations;

8) Print './prog.exe' to run the program.


HOW PROGRAM WORKS

When it starts, print '10' to continue.
If you print anything else and repeat this 10 times,
the game will be over and you will have to start ./prog.exe
again to try getting your data.

When you print '10' successfully, the window closes and
the prog runs for about 3 seconds. After that you receive
your files, which are located in the 'results' folder.
Time elapsed will be given you in the new window.
You can watch it and press enter to close the window
and end the prog.