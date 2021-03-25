#include"body.h"


int findxiny(int x, int *y);
int outpart (struct part ans,double n,int N);
int updateG(struct graph *G, int x, int y);
void extrafG(struct part ans,struct graph G,int *s);
void gcopy(struct graph G,struct graph *Ga);
struct part RG(struct graph G, int ke,unsigned int *seed);
