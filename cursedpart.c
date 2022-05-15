#include "cursedpart.h"

const int width = 50;
const int height = 20;

int welcome()
{
    if (!initscr())
        {
            fprintf(stderr, "Error initialising ncurses.\n");
            exit(1);
        }

    char hello[] = "Welcome! Press 10 to postavit' otl Glebu>>>\n";
    mvaddstr(LINES/2, (COLS-strlen(hello))/2, hello);
    refresh();

    int Nres = 0;
    int ch = 0;
    int i = 1;
    int flag = 0;

    while(1)
    {
        noecho();
        ch = getch();

        if (ch == 49)
            flag = 1;

        else
        {
            ++i;
            if (i > 10)
            {
                char gameOver[] = "GAME OVER.\n";
                mvaddstr(LINES/2 + i, (COLS-strlen(gameOver))/2, gameOver);
                refresh();

                getch();
                exit(-1);
            }

            char tryAgain[] = "Ne, tak ne poydyot. Try again.\n";
            mvaddstr(LINES/2 + i, (COLS-strlen(tryAgain))/2, tryAgain);
            refresh();
        }

        if (flag == 1)
        {
            ch = getch();
            if (ch == 48)
            {
                char goFurther[] = "Nice, press any key to continue>>>\n";
                mvaddstr(LINES/2 + i, (COLS-strlen(goFurther))/2, goFurther);
                refresh();
                getch();

                char askN1[] = "Enter N1: \n";
                mvaddstr(LINES/2 + i + 1, (COLS-strlen(goFurther))/2, askN1);
                refresh();

                int j = 1;
                while(flag != 2)
                {
                    ch = getch();
                    if (ch >= 48 && ch <= 57)
                        Nres = Nres*10 + (ch - 48);

                    else flag = 2;
                }

                char Nstr[] = "Entered Nres is: ";
                mvaddstr(LINES/2 + i + 2, (COLS-strlen(goFurther))/2, Nstr);
                printw("%d\n", Nres);
                refresh();
                getch();

                break;
            }

            else
            {
                char smthWrong[] = "Something went wrong. Try again.\n";
                mvaddstr(LINES/2 + i, (COLS-strlen(smthWrong))/2, smthWrong);
                refresh();

                flag = 0;
            }
        }
    }
    endwin();
    return Nres;
}
void printWorkTime(double elapsed)
{
    if (!initscr())
        {
            fprintf(stderr, "Error initialising ncurses.\n");
            exit(1);
        }
    refresh();

    char timeStr[] = "Time measured (seconds): ";

    mvaddstr(LINES/2, (COLS-strlen(timeStr))/2, timeStr);
    printw("%.3f\n", elapsed);

    refresh();
    getch();

    endwin();
}
