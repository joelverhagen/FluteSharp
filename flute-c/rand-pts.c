#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if __GLIBC__ 
#include <unistd.h>
#else
#include <process.h>
#endif

int main(int ac, char *av[])
{
    int d=10, tmp, i;
    int PNUM = 0;

    for (i=1; i<ac; i++) {
        if (strcmp(av[i], "-r")==0)  // random
#if __GLIBC__ 
            srandom((int) getpid());
#else
            srand((int) _getpid());
#endif
        else if (strncmp(av[i], "-s", 2)==0)  // set random seed
#if __GLIBC__ 
            srandom(atoi(av[i] + 2));
#else
            srand(atoi(av[i] + 2));
#endif
        else if (strcmp(av[i], "-n")==0)  // print # of points first
            PNUM=1;
        else if (sscanf(av[i], "%d", &tmp))  // set # of points
            d = tmp;
        else {
            printf("Usage: %s [-r] [-s<S>] [-n] [<D>]\n", av[0]);
            printf("  Output <D> random points ");
            printf("as <D> lines of coordinate pairs.\n");
            printf("  Default <D> is 10.\n");
            printf("  -r\t Randomize. Use getpid() as seed.\n");
            printf("  -s<S>\t Set random seed to <S>.\n");
            printf("  -n\t Write <D> first before the random points.\n");
            exit(-1);
        }
    }
    
    if (PNUM)
        printf("%d\n", d);
    for (i=1; i<=d; i++)
#if __GLIBC__ 
        printf("%4d %4d\n", (int) random()%10000, (int) random()%10000);
#else
        printf("%4d %4d\n", (int) rand()%10000, (int) rand()%10000);
#endif
}
