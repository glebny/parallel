#include "cursedpart.h"

void printHello()
{
    const int width = 50;
    const int height = 20;

    if (!initscr())
        {
            fprintf(stderr, "Error initialising ncurses.\n");
            exit(1);
        }

    keypad (stdscr, TRUE);
    curs_set(0);
    refresh();

    int offsetx = (COLS - width) / 2;
    int offsety = (LINES - height) / 2;

    WINDOW *win = newwin(LINES, COLS, offsetx, offsety);

    char hello[] = "It all reminds that episode of The Simpsons about Ned Flanders' house ahahehe...\n";

    mvaddstr(LINES/2, (COLS-strlen(hello))/2, hello);

    wrefresh(win);
    wgetch(win);

    delwin(win);
    endwin();
}

void printWorkTime(double elapsed)
{
    const int width = 50;
    const int height = 20;

    if (!initscr())
        {
            fprintf(stderr, "Error initialising ncurses.\n");
            exit(1);
        }

    keypad (stdscr, TRUE);
    curs_set(0);
    refresh();

    int offsetx = (COLS - width) / 2;
    int offsety = (LINES - height) / 2;

    WINDOW *win = newwin(LINES, COLS, offsetx, offsety);

    char timeStr[] = "Time measured (seconds): ";
    //Time measured: %.3f seconds.\n", elapsed

    move (LINES/2, (COLS-strlen(timeStr))/2 - 10);
    for (int i = 0; i < 25; ++i) //25 is the length of timeStr. Yes, I was that lazy.
        wprintw (win, "%c", timeStr[i]);
    wprintw(win, "%d\n", elapsed);

    wrefresh(win);
    wgetch(win);

    delwin(win);
    endwin();
}

void askMatrixSize()
{

}
