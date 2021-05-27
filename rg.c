#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"body.h"
#include"help.h"
extern double cp;


int findxiny(int x, int *y)
{
    int j;
    for (j=1;j<=y[0];j++)
    
    if (y[j]==x) return (j);
    return (-1);
}


int outpart(struct part ans,double n, int N, char *fname)
{
    int i;
    printf("%d\t%lf\n",ans.com,ans.Q/(2*n));
	FILE *fo;
	fo=fopen(fname,"w");
    	for (i=0;i<N;i++)
    		fprintf(fo,"%d\n",ans.pa[i]);
	fclose(fo);
    return 0;
}



int updateG(struct graph *G, int x, int y)
{
    int i,k;
    int *newlist;
    double *newl2;

    newlist=(int*)malloc(sizeof(int)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    newl2  =(double*)malloc(sizeof(double)*(G[0].nl[x][0]+G[0].nl[y][0]+1));
    k=0;
    newl2[0]=G->wl[x][0]+G->wl[y][0];
    for (i=1;i<=G[0].nl[x][0];i++)
    {
        if (G[0].nl[x][i]==y){	k=-1;	newl2[0]+=G[0].wl[x][i];}
        else{
            newlist[i+k]=G[0].nl[x][i];
            newl2[i+k]=G[0].wl[x][i];
        }
    }
    newlist[0]=G[0].nl[x][0]-1;

    for (i=1;i<=G[0].nl[y][0];i++)
    if (G[0].nl[y][i]!=x)
    {
        int sg,k1,k2;
        sg=G[0].nl[y][i];
        k1=findxiny(x,G[0].nl[sg]);
        k2=findxiny(y,G[0].nl[sg]);
        if (k1>0)
        {
            G[0].wl[sg][k1]+=G[0].wl[sg][k2];
            int test;
            test=findxiny(sg,newlist);
            if (test>0)
                newl2[test]=G[0].wl[sg][k1];
            else
            {
                printf("Error2!\n");
                exit(0);
            }
            
            if (k2!=G[0].nl[sg][0])
            {
                G[0].nl[sg][k2]=G[0].nl[sg][G[0].nl[sg][0]];
                G[0].wl[sg][k2]=G[0].wl[sg][G[0].nl[sg][0]];
            }
            G[0].nl[sg][0]--;
        }
        else
        {
            if (k2>0)
            G[0].nl[sg][k2]=x;
            else
            {
                printf("Error3!\n");   
		printf("sg=%d,k1=%d,k2=%d,x=%d,y=%d,i=%d,y[i]=%d\n",sg,k1,k2,x,y,i,G->nl[y][i]);
		for (k=0;k<=G->nl[y][0];k++)
		printf("%d ",G->nl[y][k]);
		printf("\n");
		for (k=0;k<=G->nl[sg][0];k++)
		printf("%d ",G->nl[sg][k]);
		printf("\n");
		for (k=0;k<G->com;k++)
		printf("%d ",G->nl[k][0]);
		printf("\n");
		exit(0);
            }
            newlist[0]++;
            newlist[newlist[0]]=sg;
            newl2[newlist[0]]=G[0].wl[sg][k2];
        }
    }
    free(G[0].nl[x]);	free(G[0].wl[x]);
    G->nl[x]=newlist;
    G->wl[x]=newl2;
    
    
    free(G[0].nl[y]);	free(G[0].wl[y]);
    G->nl[y]=NULL;	G->wl[y]=NULL;
    if (y!=G[0].com-1)
    {
        for (i=1;i<=G[0].nl[G[0].com-1][0];i++)
        {
            int k1,k2;
            k1=G[0].nl[G[0].com-1][i];
            k2=findxiny(G[0].com-1,G[0].nl[k1]);
            G[0].nl[k1][k2]=y;
        }
        G->nl[y]=G->nl[G->com-1];
        G->wl[y]=G->wl[G->com-1];
        G->nl[G->com-1]=NULL;
        G->wl[G->com-1]=NULL;
    }
    
    G[0].dl[x]+=G[0].dl[y];
    if (y!=G[0].com-1)
    {
        G[0].dl[y]=G[0].dl[G[0].com-1];
    }
    G->ml[x][0]+=G->ml[y][0];
    if (y!=G->com-1)
    {
	G->ml[y][0]=G->ml[G->com-1][0];
    }
    return 0;
}

void extrafG(struct part ans,struct graph G,int *s)
{
    int i,j;
    for (i=0;i<G.com;i++)
    {
	if (G.ml[i]!=NULL)
		{
        	for (j=1;j<=G.ml[i][0];j++)
		if (G.ml[i][j]<=G.N&&G.ml[i][j]>=1)
        		ans.pa[G.ml[i][j]-1]=s[i]+1;
		else
			{printf("extrafG wrong 1\n");
			exit(1);}
		}
		else
		{
		printf("extrafG wrong 2\n");
		exit(1);
		}
    }
}

void gcopy(struct graph G,struct graph *Ga)
{
    Ga->N=G.N;	Ga->com=G.com;	Ga->n=G.n;
    int i,j;
    Ga->ml=(int**)malloc(sizeof(int*)*G.com);
    Ga->nl=(int**)malloc(sizeof(int*)*G.com);
    Ga->wl=(double**)malloc(sizeof(double*)*G.com);
    Ga->dl=(double*)malloc(sizeof(double)*G.com);
    for (i=0;i<G.com;i++)
    {
        Ga->ml[i]=(int*)malloc(sizeof(int));
        Ga[0].ml[i][0]=G.ml[i][0];
        
        Ga->nl[i]=(int*)malloc(sizeof(int)*(G.nl[i][0]+1));
        for (j=0;j<=G.nl[i][0];j++)
        Ga->nl[i][j]=G.nl[i][j];
        
        Ga[0].wl[i]=(double*)malloc(sizeof(double)*(G.nl[i][0]+1));
        for (j=0;j<=G.nl[i][0];j++)
        Ga[0].wl[i][j]=G.wl[i][j];
        Ga[0].dl[i]=G.dl[i];
    }
    
    
}

void freeG(struct graph G)
{
    int i;
    for (i=0;i<G.com;i++)
    {
        free(G.ml[i]);
        free(G.nl[i]);
        free(G.wl[i]);
    }
    free(G.dl);
    free(G.ml);
    free(G.nl);
    free(G.wl);
}

double movedQ(struct graph G, int child, int home, int school, int *s,double *deg,int *nodes, double *links)
{
	int i,j,mem;
	double ans;
	double L1t=0,L2t=0;
	for (i=1;i<=G.nl[child][0];i++)
	{
		mem=s[G.nl[child][i]];
		if (mem==home)
			L1t+=G.wl[child][i];
		else if (mem==school)
			L2t+=G.wl[child][i];
	}
	double q1,q2,q1p,q2p,d1,d2,d1p,d2p;
	q1=2*links[home]-deg[home]*deg[home]/(2*G.n);
	q2=2*links[school]-deg[school]*deg[school]/(2*G.n);
	if (nodes[home]<=1)
	d1=0;
	else 
	d1=2*links[home]/nodes[home]/(nodes[home]-1);
	if (nodes[school]<=1)
	d2=0;
	else
	d2=2*links[school]/nodes[school]/(nodes[school]-1);
	q1p=2*(links[home]-G.wl[child][0]-L1t)-(deg[home]-G.dl[child])*(deg[home]-G.dl[child])/(2*G.n);
	q2p=2*(links[school]+G.wl[child][0]+L2t)-(deg[school]+G.dl[child])*(deg[school]+G.dl[child])/(2*G.n);
	if (nodes[home]-G.ml[child][0]<=1)
	d1p=0;
	else
	d1p=2*(links[home]-G.wl[child][0]-L1t)/(nodes[home]-G.ml[child][0])/(nodes[home]-G.ml[child][0]-1);
	if (nodes[school]+G.ml[child][0]<=1)
	d2p=0;
	else
	d2p=2*(links[school]+G.wl[child][0]+L2t)/(nodes[school]+G.ml[child][0])/(nodes[school]+G.ml[child][0]-1);
	ans=q1p*pow(d1p,cp)+q2p*pow(d2p,cp)-q1*pow(d1,cp)-q2*pow(d2,cp);

	return  (ans);
}




struct part RG(struct graph G, int ke, unsigned int *seed)
{
    int i;
    struct part ans;
    struct graph Ga;
	
    gcopy(G,&Ga);
    ans.pa=(int*)malloc(sizeof(int)*G.N);;

    int *s,*ns,kee;
    double *deg;
	s=(int*)malloc(sizeof(int)*G.com);
	ns=(int*)malloc(sizeof(int)*G.com);
	deg=(double*)malloc(sizeof(double)*G.com);
    double Q;
    Q=compQ(Ga);

    double dQmax,dQ;

    int j,k,g1,g2,lr,randg;
    double Q1,Q2,den1,den2,den12;
    int temp;

    for (i=0;i<G.com;i++)
	{
		s[i]=i; 
		ns[i]=i;
		deg[i]=G.dl[i];   
    	}

    ans.com=G.com;
    ans.Q=Q;
    int maxiter=G.com-1;
    
    for (i=1;i<=maxiter;i++)
    {

        dQmax=-2*G.n-1; lr=0;
	if (Ga.com<ke)
		kee=Ga.com;
	else
		kee=ke;


        for (j=0;j<kee;j++)
        {
            randg=randint(Ga.com,seed);
            for (k=1;k<=Ga.nl[randg][0];k++)
            {
		temp=Ga.nl[randg][k];
		Q1=2*Ga.wl[randg][0]-Ga.dl[randg]*Ga.dl[randg]/(2*Ga.n);
		Q2=2*Ga.wl[temp][0]-Ga.dl[temp]*Ga.dl[temp]/(2*Ga.n);
		if (Ga.ml[randg][0]==1)
			den1=0;
		else
			den1=2*Ga.wl[randg][0]/Ga.ml[randg][0]/(Ga.ml[randg][0]-1);
		if (Ga.ml[temp][0]==1)
			den2=0;
		else
			den2=2*Ga.wl[temp][0]/Ga.ml[temp][0]/(Ga.ml[temp][0]-1);
		den12=2*(Ga.wl[randg][0]+Ga.wl[temp][0]+Ga.wl[randg][k])/(Ga.ml[randg][0]+Ga.ml[temp][0])/(Ga.ml[randg][0]+Ga.ml[temp][0]-1);

                dQ=2*(Ga.wl[randg][k]-Ga.dl[randg]*1.0*Ga.dl[Ga.nl[randg][k]]/(2*Ga.n))*pow(den12,cp)+Q1*(pow(den12,cp)-pow(den1,cp))+Q2*(pow(den12,cp)-pow(den2,cp));
                if (abso(dQ-dQmax)<Inf_Sma)
                {
                    lr++;
                    if (rand1(seed)*lr<1)
                    {
                        g1=randg;	g2=Ga.nl[randg][k];
                    }
                }
                else if (dQ>dQmax)
                {
                    lr=1;
                    dQmax=dQ;
                    g1=randg;	g2=Ga.nl[randg][k];
                }
            }
        }
        if (lr==0) break;
        else	{
		if (Ga.nl[g1][0]<Ga.nl[g2][0])
        	{
			g1=g1+g2;	g2=g1-g2;	g1=g1-g2;
		}    	

		updateG(&Ga,g1,g2);
		
            Ga.com--;
		for (j=0;j<G.com;j++)
			if (s[j]==g2)
				s[j]=g1;
		if (g2!=Ga.com)
		for (j=0;j<G.com;j++)
			if (s[j]==Ga.com)
				s[j]=g2;
        }
        Q+=dQmax;
        if (Q>ans.Q)
        {
            ans.com=Ga.com;
            ans.Q=Q;
		for (j=0;j<Ga.com;j++)
			deg[j]=Ga.dl[j];
            for (j=0;j<G.com;j++)
			ns[j]=s[j];
        }

    }

	int *nodes;
	double *links;
	nodes=(int*)malloc(sizeof(int)*ans.com);
	links=(double*)malloc(sizeof(double)*ans.com);
	for (i=0;i<ans.com;i++)
	{
		nodes[i]=0;	links[i]=0;
	}
	for (i=0;i<G.com;i++)
	{
		nodes[ns[i]]+=G.ml[i][0];
		for (j=1;j<=G.nl[i][0];j++)
		if (ns[i]==ns[G.nl[i][j]])
		links[ns[i]]+=G.wl[i][j];
	}
	for (i=0;i<ans.com;i++)
		links[i]=links[i]/2;
	for (i=0;i<G.com;i++)
		links[ns[i]]+=G.wl[i][0];

//refine section
	int change=1;
	while (change)
	{
		change=0;
		for (i=0;i<G.com;i++)
		{
			g1=ns[i];
			dQmax=-1;
			for (j=0;j<ans.com;j++)
				if (j!=g1)
				{
					dQ=movedQ(G,i,g1,j,ns,deg,nodes,links);
					if (dQ>dQmax)
					{
						dQmax=dQ;
						g2=j;
					}
	
			}
			if (dQmax>Inf_Sma)
			{
				ans.Q+=dQmax;

				deg[g1]-=G.dl[i];
				deg[g2]+=G.dl[i];

				nodes[g1]-=G.ml[i][0];
				nodes[g2]+=G.ml[i][0];
				for (j=1;j<=G.nl[i][0];j++)
				if (ns[G.nl[i][j]]==g1)
				links[g1]-=G.wl[i][j];
				else if (ns[G.nl[i][j]]==g2)
				links[g2]+=G.wl[i][j];
				links[g1]-=G.wl[i][0];
				links[g2]+=G.wl[i][0];
				ns[i]=g2;
				change=1;
			}
		}
	}

	free(links);
	free(nodes);
	extrafG(ans,G,ns);
//below is to reorder ans.partition
	j = 0;
	int *lab;
	lab=(int*)malloc(sizeof(int)*G.N);
	for (i=0;i<G.N;i++)
		lab[i]=0;
	for (i = 0;i < G.N;i++)
	{
		if (!lab[i])
		{
			j++;
			lab[i]=j;
		
		for (k = i + 1;k < G.N;k++)
			if (ans.pa[k] == ans.pa[i])
				lab[k] = lab[i];
		}
	}
	for (i=0;i<G.N;i++)
		ans.pa[i]=lab[i];
	ans.com=j;
	free(lab);
	free(s);
	free(ns);
	free(deg);
    freeG(Ga);
    
    return (ans);
}
