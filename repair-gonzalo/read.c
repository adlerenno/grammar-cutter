
#include <stdio.h>

main()

  { FILE *f = fopen("ver.R","r");
    int a,b;
    int alph;
    char c;
    fread(&alph,4,1,f);
    while (alph--)
       { fread(&c,1,1,f);
         printf ("%c ",c);
       }
    printf("\n");
    while (!feof(f))
       { fread(&a,4,1,f);
         fread(&b,4,1,f);
         printf ("%i %i\n",a,b);
       }
  }
