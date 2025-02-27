/************************************/
/*PARTICLE SWARM OPTIMIZATION */
/************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*CONSTANTS */
#define NOPS 100		/*NUMBER OF PARTICLES */
#define LIMITL 256		/*MAXIMUM SIZE OF THE ARRAY FOR PARTICLE DATA */
#define ILIMIT 150		/*MAXIMUM NUMBER OF ITERATION */
#define SEED 32767		/*INITIALIZATION OF RANDOM NUMBER */
#define W 0.7			/*INERTIA CONSTANT */
#define C1 1.4			/*CONSTANT OF ATTRACTION FROM PERSONAL BEST */
#define C2 1.4			/*CONSTANT OF ATTRACTION FROM GROUP BEST */
#define DMAX 7			/*THE NUMBER OF THE NODES */


/*STRUCTURE FOR ARBITRARY COORDINATE*/
struct point
{
  double x[DMAX];
};

/*STRUCTURE FOR PARTILE*/
struct particle
{
  struct point pos;
    /*POSITION*/ double value;	/*EVALUATED VALUE */
  struct point v;
    /*VELOCITY*/ struct point bestpos;	/*THE BEST POSITION */
  double bestval;		/*THE BEST EVALUATED VALUE */
};


/*PROTOTYPES OF FUNCTIONS */
void initps (struct particle ps[]);
/*INITIALIZATION OF PARTICLE SWARM*/
void printps (struct particle ps[]);
/*PRINT OUT OF PARTICLE SWARM DATA*/
double frand ();		/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE */
double calcval (double *x);	/*THE OBJECTIVE FUNCTION */
double calcresidue (double *x);	/*THE RESIDUE OF THE OPTIMIZATION */

void optimize (struct particle ps[]);
/*OPTIMIZATION*/ void setgbest (struct particle ps[]);	/*FIND THE GROUP BEST POSITION */

/*GLOBAL VARIABLES*/
struct point gbestpos;		/*THE GROUP BEST POSITION */
double gbestval;		/*THE BEST VALUE OF THE GROUP */
double YTARGET[DMAX * DMAX];	/*LOGARITHM OF THE FLOWS */
double ETARGET[DMAX * DMAX];	/*THE RAW DATA OF THE FLOWS ON THE EDGES */
int EDGE[DMAX * DMAX];		/*THE ADJACENCY BETWEEN THE NODES */
double SOL[DMAX];		/*WORK ARRAY */


/*DIV OPERATOR*/
double
div (int i, double *Y)
{
  double sum = 0.;
  for (int k = i + 1; k < DMAX; k++)
    {
      sum += Y[DMAX * i + k];
    }
  for (int k = 0; k < i; k++)
    {
      sum += -Y[DMAX * k + i];
    }
  return sum;
}


/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE */
double
doublerand (void)
{
  double result;
  while ((result = (double) rand () / RAND_MAX) >= 1);
  return result;
}

/*READ INPUT DATA OF FLOWS AND EDGES*/
void
readTARGET ()
{
  FILE *f1;
  f1 = fopen ("edge.txt", "r");
  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
      fscanf (f1, "%d%d%d%d%d%d%d",
	      &EDGE[k],
	      &EDGE[k + 1],
	      &EDGE[k + 2],
	      &EDGE[k + 3], &EDGE[k + 4], &EDGE[k + 5], &EDGE[k + 6]);
    }

  f1 = fopen ("exc.txt", "r");
  char d1[10], d2[10], d3[10];
  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
      fscanf (f1, "%s%s%s%lf%lf%lf%lf%lf%lf%lf",
	      d1, d2, d3,
	      &ETARGET[k],
	      &ETARGET[k + 1],
	      &ETARGET[k + 2],
	      &ETARGET[k + 3],
	      &ETARGET[k + 4], &ETARGET[k + 5], &ETARGET[k + 6]);
    }
  int k = 0;
  fscanf (f1, "%s%s%lf%lf%lf%lf%lf%lf%lf",
	  d1, d2,
	  &SOL[k],
	  &SOL[k + 1],
	  &SOL[k + 2], &SOL[k + 3], &SOL[k + 4], &SOL[k + 5], &SOL[k + 6]);

  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
      printf ("%s %s %s %lf %lf %lf %lf %lf %lf %lf\n",
	      d1, d2, d3,
	      ETARGET[k],
	      ETARGET[k + 1],
	      ETARGET[k + 2],
	      ETARGET[k + 3], ETARGET[k + 4], ETARGET[k + 5], ETARGET[k + 6]);
    }
  k = 0;
  printf ("%s %s %lf %lf %lf %lf %lf %lf %lf\n",
	  d1, d2,
	  SOL[k],
	  SOL[k + 1],
	  SOL[k + 2], SOL[k + 3], SOL[k + 4], SOL[k + 5], SOL[k + 6]);


  for (int i = 0; i < DMAX * DMAX; i++)
    {
      YTARGET[i] = -log (ETARGET[i]);
    }

  for (int i = 0; i < DMAX; i++)
    {
      k = DMAX * i;
      printf ("%lf %lf %lf %lf %lf %lf %lf\n",
	      YTARGET[k],
	      YTARGET[k + 1],
	      YTARGET[k + 2],
	      YTARGET[k + 3], YTARGET[k + 4], YTARGET[k + 5], YTARGET[k + 6]);
    }

  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  YTARGET[DMAX * i + j] += 2 * doublerand () - 1;
	  YTARGET[DMAX * j + i] = -YTARGET[DMAX * i + j];
	}
    }

}

int
main ()
{

  readTARGET ();

  struct particle ps[LIMITL];	/*THE SWARM OF PARTICLES */
  int i;			/*ITERATOR OF THE REPETITION */

/*INITIALIZE THE RANDOM NUMBER*/
  srand (SEED);

/*INITIALIZATION OF PARTICLE SWARM*/
  initps (ps);

  printps (ps);

   /*OPTIMIZATION*/ for (i = 0; i < ILIMIT; ++i)
    {
      optimize (ps);
      printf ("\n%d\n", i);
      printps (ps);
    }


  calcresidue (gbestpos.x);


  return 0;

}


void
optimize (struct particle ps[])
{
  int i, j;
  double r1, r2;		/*RANDOM NUMBER */

  for (i = 0; i < NOPS; ++i)
    {

/*SET RANDOM NUMBER*/
      r1 = frand ();
      r2 = frand ();

/*VELOCITY UPDATE*/
      for (j = 0; j < DMAX; j++)
	{
	  ps[i].v.x[j] =
	    W * ps[i].v.x[j] + C1 * r1 * (ps[i].bestpos.x[j] -
					  ps[i].pos.x[j]) +
	    C2 * r2 * (gbestpos.x[j] - ps[i].pos.x[j]);
	}

/*POSITION UPDATE*/
      for (j = 0; j < DMAX; j++)
	{
	  ps[i].pos.x[j] += ps[i].v.x[j];
	}
/*UPDATE THE BEST VALUE FOR EACH PARTICLE*/
      ps[i].value = calcval (ps[i].pos.x);

      if (ps[i].value < ps[i].bestval)
	{
	  ps[i].bestval = ps[i].value;
	  ps[i].bestpos = ps[i].pos;
	}
    }

/*UPDATE THE BEST VALUE IN THE GROUP*/
  setgbest (ps);

}



/*RESIDUE OF THE OPTIMIZATION*/
double
calcresidue (double *XIN)
{
  double XX[DMAX * DMAX];
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = 0; j < DMAX; j++)
	{
	  XX[i * DMAX + j] = XIN[i] - XIN[j];
	}
    }

  double ZZ[DMAX * DMAX];
  double norm = 0.0;
  for (int i = 0; i < DMAX * DMAX; i++)
    {
      ZZ[i] = YTARGET[i] - XX[i];
      norm += (XX[i] - YTARGET[i]) * (XX[i] - YTARGET[i]);
    }

  printf ("\n KER(DIV)?");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i; j < DMAX; j++)
	{
	  printf ("%lf ", div (i, ZZ));
	}
    }
  printf ("\n");
  printf ("2-NORM=%lf", norm);

  FILE *f1;
  f1 = fopen ("KERDIV.txt", "w");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  int k = i * DMAX + j;
	  fprintf (f1, "%lf ", ZZ[k]);
	}
    }
  fclose (f1);
  printf ("THE RESIDUE OF GRAD-S OPTIMIZATION IS WRITTEN IN KERDIV.txt.\n");
  printf ("PROCESS IT BY CURLSTAR OPTIMIZATION TO GET THE HARMONIC COMPONENT.\n");
  return norm;
}

/*THE OBJECTIVE FUNCTION*/
double
calcval (double *XIN)
{

  double XX[DMAX * DMAX];
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = 0; j < DMAX; j++)
	{
	  XX[i * DMAX + j] = XIN[i] - XIN[j];
	}
    }

  double ZZ[DMAX * DMAX];
  double norm = 0.0;
  for (int i = 0; i < DMAX * DMAX; i++)
    {
      ZZ[i] = YTARGET[i] - XX[i];
      norm += (XX[i] - YTARGET[i]) * (XX[i] - YTARGET[i]);
    }

  return norm;
}


/*FIND THE BEST IN THE GROUP*/
void
setgbest (struct particle ps[])
{
  int i, j;
  double besti;
  double x[DMAX];
  besti = ps[0].value;

  for (j = 0; j < DMAX; j++)
    {
      x[j] = ps[0].pos.x[j];
    }

  for (i = 0; i < NOPS; ++i)
/*FIND THE CURRENT BEST VALUE*/
    if (ps[i].value < besti)
      {
	besti = ps[i].value;
	for (j = 0; j < DMAX; j++)
	  {
	    x[j] = ps[i].pos.x[j];
	  }
      }

/*UPDATE IF POSSIBLE*/
  if (besti < gbestval)
    {
      gbestval = besti;
      for (j = 0; j < DMAX; j++)
	{
	  gbestpos.x[j] = x[j];
	}
    }

}


/*INITIALIZATION OF PARTICLE SWARM */
void
initps (struct particle ps[])
{

  int i, j;
  double x[DMAX];


  for (i = 0; i < NOPS; ++i)
    {
       /*POSITION*/ for (j = 0; j < DMAX; j++)
	{
	  x[j] = ps[i].pos.x[j] = 4 * (frand () * 2 - 1.0);
	}
/*EVALUATED VALUE*/
      ps[i].value = calcval (x);
       /*VELOCITY*/ for (j = 0; j < DMAX; j++)
	{
	  ps[i].v.x[j] = frand () * 2 - 1.0;
	}
/*THE BEST POSITION*/
      for (j = 0; j < DMAX; j++)
	{
	  ps[i].bestpos.x[j] = ps[i].pos.x[j];
	}
/*THE BEST EVALUATED VALUE*/
      ps[i].bestval = ps[i].value;


    }
/*THE BEST POSITION IN THE GROUP*/
  gbestval = ps[0].value;
  for (j = 0; j < DMAX; j++)
    {
      gbestpos.x[j] = ps[0].pos.x[j];
    }
  setgbest (ps);
}


/*PRINT OUT OF PARTICLE SWARM DATA */

void
printps (struct particle ps[])
{
  int debug = 0;
  if (debug == 1)
    {
      for (int i = 0; i < NOPS; ++i)
	{
	  printf ("%d ", i);
	  for (int j = 0; j < DMAX; j++)
	    {
	      printf ("%lf ", ps[i].pos.x[j]);
	    }
	  for (int j = 0; j < DMAX; j++)
	    {
	      printf ("%lf ", ps[i].v.x[j]);
	    }
	  printf ("\n");
	}
    }
  printf ("\n BEST [given by exp(-s)]: ");
  for (int j = 0; j < DMAX; j++)
    {
      printf ("%lf ", exp (-gbestpos.x[j] + gbestpos.x[0]));
    }
  printf ("\n BEST [s, as it is]: ");
  for (int j = 0; j < DMAX; j++)
    {
      printf ("%lf ", exp (-gbestpos.x[j] + gbestpos.x[0]));
    }
//for (int j=0;j<DMAX;j++){printf("%lf ",SOL[j]/SOL[0]);} 
  printf (" -> %lf\n", gbestval);
}


/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
double
frand (void)
{
  double result;
  while ((result = (double) rand () / RAND_MAX) >= 1);
  return result;
}
