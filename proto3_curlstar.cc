/************************************/
/*PARTICLE SWARM OPTIMIZATION */
/* 
pso.c 
*/
/* 
*/
/************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

/*CONSTANTS */
#define NOPS 250		/*NUMBER OF PARTICLES */
#define LIMITL 256		/*MAXIMUM SIZE OF THE ARRAY FOR PARTICLE DATA */
#define ILIMIT 250		/*MAXIMUM NUMBER OF ITERATION */
#define SEED 32767		/*INITIALIZATION OF RANDOM NUMBER */
#define W 0.7			/*INERTIA CONSTANT */
#define C1 1.4			/*CONSTANT OF ATTRACTION FROM PERSONAL BEST */
#define C2 1.4			/*CONSTANT OF ATTRACTION FROM GROUP BEST */


#define DMAX 7			/*THE NUMBER OF NODES */
#define TMAX 7*5		/*THE NUMBER OF THE REPRESENTATIVES OF 3-CLIQUES */

/*STRUCTURE FOR ARBITRARY COORDINATE*/
struct point
{
  double x[TMAX];
};


/*STRUCTURE FOR PARTICLE*/
struct particle
{
  struct point pos;
    /*POSITION*/ double value;	/*EVALUATED VALUE */
  struct point v;
    /*VELOCITY*/ struct point bestpos;	/*THE BEST POSITION */
  double bestval;		/*THE BEST EVALUATED VALUE */
};


/*PROTOTYPES OF FUNCTIONS */
void initps (struct particle ps[]);	/*INITIALIZATION OF PARTICLE SWARM */
void printps (struct particle ps[]);	/*PRINT OUT OF PARTICLE SWARM DATA */
double frand ();		/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE */
double calcval (double *x);	/*THE OBJECTIVE FUNCTION */
double calcresidue (double *x);	/*THE RESIDUE OF THE OPTIMIZATION */

void optimize (struct particle ps[]);
/*OPTIMIZATION*/ void setgbest (struct particle ps[]);	/*FIND THE GROUP BEST */

/*GLOBAL VARIABLES*/
struct point gbestpos;		/*THE BEST POSITION IN THE GROUP */
double gbestval;		/*THE BEST VALUE IN THE GROUP */

double YTARGET[DMAX * DMAX];	/*LOGARITHM OF THE EDGE FLOW */
double ETARGET[DMAX * DMAX];	/*THE RAW DATA OF THE EDGE FLOW */
double SOL[TMAX];
int EDGE[DMAX * DMAX];
/*ADJACENCY*/ int TCLIQUE[TMAX];
/*TCLIQUE[i] : DOES i-th 3-CLIQUE REALLY EXIST?
IF IT EXISTS, RETURN 1; IF NOT 0 
*/

double forphi[TMAX];		/*WORK ARRAY */

int
max (int a, int b)
{
  if (a > b)
    {
      return a;
    }
  return b;
}

int
min (int a, int b)
{
  if (a < b)
    {
      return a;
    }
  return b;
}

/*OPERATOR CURL*/
double
curl (int i, int j, int k, double *Y)
{
  int tmin, tmax;
  double sum = 0.;
  printf ("\n %d %d %d\n", i, j, k);
  printf ("\n%d %d %d ", EDGE[DMAX * i + j], EDGE[DMAX * j + k],
	  EDGE[DMAX * k + i]);
  printf ("%lf %lf %lf \n", Y[DMAX * i + j], Y[DMAX * j + k],
	  Y[DMAX * k + i]);

  sum += Y[DMAX * i + j] * EDGE[DMAX * i + j];
  sum += Y[DMAX * j + k] * EDGE[DMAX * j + k];
  sum += Y[DMAX * k + i] * EDGE[DMAX * k + i];
  return sum;
}

/*READ THE DATA OF THE EDGES AND THE FLOW*/
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

  for (int i = 0; i < DMAX; i++)
    {
      for (int j = 0; j < DMAX; j++)
	{
	  if (i == j)
	    {
	      EDGE[DMAX * i + j] = 0;
	    }
	  int k1 = DMAX * i + j;
	  int k2 = DMAX * j + i;
	  EDGE[k2] = EDGE[k1];
	}
    }

  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
      printf ("%d%d%d%d%d%d%d\n",
	      EDGE[k],
	      EDGE[k + 1],
	      EDGE[k + 2],
	      EDGE[k + 3], EDGE[k + 4], EDGE[k + 5], EDGE[k + 6]);
    }



  f1 = fopen ("KERDIV2.txt", "r");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = 0; j < DMAX; j++)
	{
	  fscanf (f1, "%lf", &YTARGET[DMAX * i + j]);
	  YTARGET[DMAX * j + i] = -YTARGET[DMAX * i + j];
	}
    }
  printf ("\nYTARGET\n");
  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
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
	  YTARGET[DMAX * i + j] *= EDGE[DMAX * i + j];
	  YTARGET[DMAX * j + i] = -YTARGET[DMAX * i + j];
	}
    }


  for (int i = 0; i < DMAX; i++)
    {
      int k = DMAX * i;
      printf ("%lf %lf %lf %lf %lf %lf %lf\n",
	      YTARGET[k],
	      YTARGET[k + 1],
	      YTARGET[k + 2],
	      YTARGET[k + 3], YTARGET[k + 4], YTARGET[k + 5], YTARGET[k + 6]);
    }

  int ic = 0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  for (int k = j + 1; k < DMAX; k++)
	    {
	      int l1 = DMAX * i + j;
	      int l2 = DMAX * j + k;
	      int l3 = DMAX * k + i;
	      TCLIQUE[ic] = EDGE[l1] * EDGE[l2] * EDGE[l3];
	      ic++;
	    }
	}
    }
}


/*THE SEQUENTIAL NUMBER OF THE CLIQUE (i0,j0,k0)*/
int
Tindex (int i0, int j0, int k0)
{
  int i, j, k, ic;
  ic = 0;
  for (i = 0; i < DMAX; i++)
    {
      for (j = i + 1; j < DMAX; j++)
	{
	  for (k = j + 1; k < DMAX; k++)
	    {
	      if (i == i0 && j == j0 && k == k0)
		{
		  return ic;
		}
	      ic++;
	    }
	}
    }
  return -1;
}

int
BubSort (int x[], int n)
{
  int i, j, temp;

  for (i = 0; i < n - 1; i++)
    {
      for (j = n - 1; j > i; j--)
	{
	  if (x[j - 1] > x[j])
	    {
	      temp = x[j];
	      x[j] = x[j - 1];
	      x[j - 1] = temp;
	    }
	}
    }
}


double
getphi (int i, int j, int k, double *PHI)
{
/*
GET THE VALUE OF 3-CLIQUE(i,j,k) WITH THE CORRECT SIGN
FROM THE REPRESENTATIVE PHI[]
*/
  int ijk = (i - j) * (j - k) * (k - i);
  int ip;
  if (ijk > 0)
    {
      ip = 1;
    };
  if (ijk == 0)
    {
      ip = 0;
    };
  if (ijk < 0)
    {
      ip = -1;
    };
  int x[3];
  x[0] = i;
  x[1] = j;
  x[2] = k;
  BubSort (x, 3);
  return ip * PHI[Tindex (x[0], x[1], x[2])] *
    TCLIQUE[Tindex (x[0], x[1], x[2])];
}


double
grad (int i, int j, double *X)
{
  int k = DMAX * i + j;
  return (X[j] - X[i]) * EDGE[k];
}

double
div (int i, double *Y)
{
  double sum = 0.;
  for (int k = i + 1; k < DMAX; k++)
    {
      int l = DMAX * i + k;
      sum += Y[DMAX * i + k] * EDGE[l];
    }
  for (int k = 0; k < i; k++)
    {
      int l = DMAX * k + i;
      sum += -Y[DMAX * k + i] * EDGE[l];
    }
  return sum;
}

double
curlstar (int i, int j, double *PHI)
{
  double sum = 0.;
  for (int k = 0; k < DMAX; k++)
    {
      sum += getphi (i, j, k, PHI);
    }
  return sum;
}


int
deltaone (double *Y)
{
  double s[DMAX];
  double Z[DMAX * DMAX];
  for (int i = 0; i < DMAX * DMAX; i++)
    {
      Z[i] = 0.;
    }

  printf ("CURL\n");
  int ic = 0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  for (int k = i + 1; k < DMAX; k++)
	    {
	      forphi[ic] = curl (i, j, k, YTARGET);
	      printf ("%lf ", forphi[ic]);
	      ic++;

	    }
	}
    }


  double A[DMAX * DMAX];
  for (int i = 0; i < DMAX * DMAX; i++)
    {
      A[i] = 0.;
    }
  printf ("\nCURLSTAR*CURL\n");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  Z[DMAX * i + j] = curlstar (i, j, forphi);
	  A[DMAX * i + j] = curlstar (i, j, forphi);
//    A[DMAX*j+i]=-curlstar(i,j,forphi);
	  printf ("%lf ", curlstar (i, j, forphi));
	}
    }

  printf ("\nDIV*CURLSTAR -- MUST BE ZERO UP TO THE MACHINE PRECISION!\n");

  for (int i = 0; i < DMAX; i++)
    {
      s[i] = div (i, A);
      printf ("%lf ", s[i]);

    }



  printf ("\ndiv\n");
  for (int i = 0; i < DMAX; i++)
    {
      s[i] = div (i, Y);
      printf ("%lf ", s[i]);

    }


  for (int i = 0; i < DMAX * DMAX; i++)
    {
      A[i] = 0.;
    }

  printf ("\ngrad*div\n");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  Z[DMAX * i + j] -= grad (i, j, s);
	  A[DMAX * i + j] = grad (i, j, s);
	  A[DMAX * j + i] = -grad (i, j, s);

	  printf ("%d %d %lf ", i, j, grad (i, j, s));
	}
    }

  printf ("\nCURL*GRAD -- MUST BE ZERO UP TO THE MACHINE PRECISION!\n");
  ic = 0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  for (int k = i + 1; k < DMAX; k++)
	    {
	      forphi[ic] = curl (i, j, k, A);
	      printf ("%lf ", forphi[ic]);
	      ic++;

	    }
	}
    }


  printf ("\n");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  printf ("%lf ", Z[DMAX * i + j]);
	}
    }
  return 0;


}

FILE *monitor;
int
main ()
{

  monitor = fopen ("monitor.txt", "w");
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

      for (j = 0; j < TMAX; j++)
	{
	  ps[i].v.x[j] =
	    W * ps[i].v.x[j] + C1 * r1 * (ps[i].bestpos.x[j] -
					  ps[i].pos.x[j]) +
	    C2 * r2 * (gbestpos.x[j] - ps[i].pos.x[j]);
	}

/*POSITION UPDATE*/
      for (j = 0; j < TMAX; j++)
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


/*OBJECTIVE FUNCTION*/
double
calcval (double *XIN)
{

  double XX[DMAX * DMAX];
  for (int i = 0; i < DMAX; i++)
    {
      XX[i * DMAX + i] = 0.0;
      for (int j = i + 1; j < DMAX; j++)
	{
	  XX[i * DMAX + j] = curlstar (i, j, XIN);
	}
    }

  double norm = 0.0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  int k = i * DMAX + j;
	  norm += (XX[k] - YTARGET[k]) * (XX[k] - YTARGET[k]);
	}
    }
  return norm;

}

/*THE RESIDUE OF THE OPTIMIZATION*/
double
calcresidue (double *XIN)
{
  printf ("CHECK OF RESIDUE\n");
  double XX[DMAX * DMAX];
  printf ("CURLSTAR \n");
  for (int i = 0; i < DMAX; i++)
    {
      XX[i * DMAX + i] = 0.0;
      for (int j = i + 1; j < DMAX; j++)
	{
	  XX[i * DMAX + j] = curlstar (i, j, XIN);
	  printf ("%lf ", XX[i * DMAX + j]);
	  XX[j * DMAX + i] = -XX[i * DMAX + j];
	}
    }
  printf ("\n");


  double norm = 0.0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  int k = i * DMAX + j;
	  norm += (XX[k] - YTARGET[k]) * (XX[k] - YTARGET[k]);
	}
    }

  double normY;
  printf ("DOES YTARGET LIE IN KER(CURL)?");
  double curlsaveY[TMAX];
  int ic = 0;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  for (int k = j + 1; k < DMAX; k++)
	    {
	      curlsaveY[ic] = curl (i, j, k, YTARGET);
	      normY += curlsaveY[ic] * curlsaveY[ic];
	      printf (" %lf", curlsaveY[ic]);
	      ic++;
	    }
	}
    }
  printf ("\n");
  printf ("DOES RESIDUE LIE IN KER(CURL)?");
  ic = 0;
  double normR;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  for (int k = j + 1; k < DMAX; k++)
	    {
	      curlsaveY[ic] -= curl (i, j, k, XX);
	      normR += curlsaveY[ic] * curlsaveY[ic];
	      printf (" %lf", curlsaveY[ic]);
	      ic++;
	    }
	}
    }
  printf ("\n");
  printf ("2-norm of Ytarget = %10.5lf \n", normY);
  printf ("2-norm of Residue = %10.5lf \n", normR);

  return norm;

}



/*FIND THE GROUP BEST*/
void
setgbest (struct particle ps[])
{
  int i, j;
  double besti;
  double x[TMAX];
  besti = ps[0].value;

  for (j = 0; j < TMAX; j++)
    {
      x[j] = ps[0].pos.x[j];
    }

  for (i = 0; i < NOPS; ++i)
/*FIND THE CURRENT BEST*/
    if (ps[i].value < besti)
      {
	besti = ps[i].value;
	for (j = 0; j < TMAX; j++)
	  {
	    x[j] = ps[i].pos.x[j];
	  }
      }

/*UPDATE, IF POSSIBLE*/
  if (besti < gbestval)
    {
      gbestval = besti;
      for (j = 0; j < TMAX; j++)
	{
	  gbestpos.x[j] = x[j];
	}
    }

}


/*INITIALIZATION OF PARTICLE SWARM*/
void
initps (struct particle ps[])
{

  int i, j;
  double x[TMAX];


  for (i = 0; i < NOPS; ++i)
    {
       /*POSITION*/ for (j = 0; j < TMAX; j++)
	{
	  x[j] = ps[i].pos.x[j] = 4 * (frand () * 2 - 1.0);
	}
/*EVALUATED VALUE*/
      ps[i].value = calcval (x);

       /*VELOCITY*/ for (j = 0; j < TMAX; j++)
	{
	  ps[i].v.x[j] = frand () * 2 - 1.0;
	}
/*THE BEST POSITION*/
      for (j = 0; j < TMAX; j++)
	{
	  ps[i].bestpos.x[j] = ps[i].pos.x[j];
	}
/*THE BEST EVALUATED VALUE*/
      ps[i].bestval = ps[i].value;


    }
/*THE BEST IN THE GROUP*/
  gbestval = ps[0].value;
  for (j = 0; j < TMAX; j++)
    {
      gbestpos.x[j] = ps[0].pos.x[j];
    }
  setgbest (ps);
}


/*PRINT OUT OF PARTICLE SWARM DATA*/
void
printps (struct particle ps[])
{
  int debug = 0;
  if (debug == 1)
    {
      int i;
      for (i = 0; i < NOPS; ++i)
	{
	  printf ("%d ", i);
	  for (int j = 0; j < TMAX; j++)
	    {
	      printf ("%lf ", ps[i].pos.x[j]);
	    }
	}
    }
  printf ("\n BEST: ");
  for (int j = 0; j < TMAX; j++)
    {
      printf ("%lf ", gbestpos.x[j]);
    }
  printf (" -> %lf\n", gbestval);

  printf ("curlstar\n");
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  printf ("%10.5lf ", curlstar (i, j, gbestpos.x));
	}
    }
  printf ("\n");

  double norm = 0.;
  for (int i = 0; i < DMAX; i++)
    {
      for (int j = i + 1; j < DMAX; j++)
	{
	  norm += curlstar (i, j, gbestpos.x) * curlstar (i, j, gbestpos.x);
	}
    }
  fprintf (monitor, "%lf\n", norm);

}


/*UNIFORM REAL RANDOM NUMBER FROM ZERO TO ONE*/
double
frand (void)
{
  double result;
  while ((result = (double) rand () / RAND_MAX) >= 1);
  return result;
}
