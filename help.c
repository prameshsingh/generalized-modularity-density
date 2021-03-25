#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"body.h"
extern double cp;


int randint(int K,unsigned int *seed)
{
	return ((int)(rand_r(seed)/(RAND_MAX*1.0+1)*K));
}

double rand1(unsigned int* seed)
{
	return (rand_r(seed)/(RAND_MAX*1.0+1));
}

double abso(double x)
{	if (x<0)  return (-x); return (x);	}

double compQ(struct graph G)
{
    double ans=0;
    int i;
	double den;
    for (i=0;i<G.com;i++)
    {
	if (G.ml[i][0]==0)
		continue;
       if (G.ml[i][0]==1)
		den=0;
	else 
		den=2*G.wl[i][0]/G.ml[i][0]/(G.ml[i][0]-1);
	ans+=(2*G.wl[i][0]-G.dl[i]*1.0*G.dl[i]/(2*G.n))*pow(den,cp);
    }
    return (ans);
}


void lcopy(int *p, int *q, int N)
{
	int i;
	for (i=0;i<N;i++)
		p[i]=q[i];
}

void erout(FILE *fo,char *message)
{
	fprintf(fo,"%s", message);
	exit(0);
}
