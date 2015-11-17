#include "function.h"
#include "u.c"

int 
ReadFouInp (char *space_group, T_CellParam * CellParam, int *grid_dim,
	    int dev, int *limit,
	    double *iso_value, double *escale,
	    double *Scale, double *B,
	    int *perm_limit, int *atom_number, double *f0,
	    double *weight, double *PAWidth, char *code_name)
{


  FILE *in;
  TOKEN testtok;
  int kint, n;
  char *ptr, *comand;
  int i, j;
  char buff[8000];
  int valid_lines, count;
  int *h_ind, *k_ind, *l_ind;

  if ((in = fopen ("sp.inp", "r")) == NULL)
    {
      fprintf (stdout, "Cannot open sp.inp\n");
      exit (1);
    }


  n = 0;
  valid_lines = 0;

  if ((in = fopen ("sp.inp", "r")) == NULL)
    exit (1);
  i = 0;
  while (EOF != fgetline (in, buff))
    if (buff[0] != '#')
      {
	n = 0;

	testtok.actual_token = testtok.next_token = buff;
	while (NULL != (ptr = str_token (&testtok)))
	  {
	    /* fprintf(stdout,"%s %d %s\n", buff, n, ptr); */
	    if (n == 0)
	      comand = ptr;
	    if ((strcmp (comand, "SpaceGroup") == 0) && (n == 1))
	      for (kint = 0; kint <= testtok.tok_length; kint++)
		{
		  space_group[kint] = ptr[kint];
		}
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 1))
	      sscanf (ptr, " %lf", &CellParam->a);
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 2))
	      sscanf (ptr, " %lf", &CellParam->b);
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 3))
	      sscanf (ptr, " %lf", &CellParam->c);
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 4))
	      sscanf (ptr, " %lf", &CellParam->alpha);
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 5))
	      sscanf (ptr, " %lf", &CellParam->beta);
	    if ((strcmp (comand, "UnitCell") == 0) && (n == 6))
	      sscanf (ptr, " %lf", &CellParam->gamma);
	    if ((strcmp (comand, "MaxRvalue") == 0) && (n == 1))
	      sscanf (ptr, " %lf", &MaxRvalue[n]);
	    if ((strcmp (comand, "MaxRvalue") == 0) && (n == 2))
	      sscanf (ptr, " %lf", &MaxRvalue[n]);
	    if ((strcmp (comand, "MaxRvalue") == 0) && (n == 3))
	      sscanf (ptr, " %lf", &MaxRvalue[n]);
	    if ((strcmp (comand, "EhklScale") == 0) && (n == 1))
	      sscanf (ptr, " %lf", escale);
	    if ((strcmp (comand, "Fscale") == 0) && (n == 1))
	      sscanf (ptr, " %lf", Scale);
	    if ((strcmp (comand, "AtomNumber") == 0) && (n == 1))
	      sscanf (ptr, " %d", atom_number);
	    if ((strcmp (comand, "ScattFac") == 0) && (n == 1))
	      sscanf (ptr, " %lf", f0);
	    if ((strcmp (comand, "TempFac") == 0) && (n == 1))
	      sscanf (ptr, " %lf", B);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 1))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 2))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 3))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 4))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 5))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "GridDimension") == 0) && (n == 6))
	      sscanf (ptr, " %d", &grid_dim[n]);
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 1))
	      for (kint = 0; kint <= testtok.tok_length; kint++)
		{
		  Comp->atom[0][kint] = ptr[kint];
		}
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 2))
	      sscanf (ptr, " %d", &Comp->atom_ratio[0]);
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 3))
	      for (kint = 0; kint <= testtok.tok_length; kint++)
		{
		  Comp->atom[1][kint] = ptr[kint];
		}
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 4))
	      sscanf (ptr, " %d", &Comp->atom_ratio[1]);
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 5))
	      for (kint = 0; kint <= testtok.tok_length; kint++)
		{
		  Comp->atom[2][kint] = ptr[kint];
		}
	    if ((strcmp (comand, "ChemicalComp") == 0) && (n == 6))
	      sscanf (ptr, " %d", &Comp->atom_ratio[2]);
	    if ((strcmp (comand, "PseudoAtom") == 0) && (n == 1))
	      sscanf (ptr, " %lf", PAWidth);
	    if ((strcmp (comand, "PermCode") == 0) && (n == 1))
	      for (kint = 0; kint <= testtok.tok_length; kint++)
		{
		  code_name[kint] = ptr[kint];
		}
	    if ((strcmp (comand, "End") == 0) && (n == 0))
	      {
		fprintf (stdout, "@Q@@@@@@@\n");
		valid_lines++;
	      }
	    if (valid_lines > 0)
	      valid_lines++;
	    n++;
	  }
      }

  fclose (in);
  valid_lines /= 6;


  if (BuildSgInfo (&SgInfo, space_group) != 0)
    fprintf (stderr, "%s\n", SgError);
  else
    {
      ListSgInfo (&SgInfo, 1, 0, stdout);
      if (SgError)
	fprintf (stderr, "%s\n", SgError);
    }
  fprintf (stdout, " acentric? %3d \n", SgInfo.Centric);


  if ((in = fopen ("sp.inp", "r")) == NULL)
    exit (1);

  h_ind = malloc ((valid_lines + 1) * sizeof (int));
  k_ind = malloc ((valid_lines + 1) * sizeof (int));
  l_ind = malloc ((valid_lines + 1) * sizeof (int));
  stat[0] = 0;

  i = 1;
  count = -1;
  while (EOF != fgetline (in, buff))
    {

      if (buff[0] != '#')
	{
	  n = 0;

	  testtok.actual_token = testtok.next_token = buff;
	  while (NULL != (ptr = str_token (&testtok)))
	    {
	      if (n == 0)
		comand = ptr;
	      if ((strcmp (comand, "End") == 0) && (n == 0))
		count++;
	      if (count > -1)
		{


		  switch (count)
		    {
		    case 1:
		      sscanf (ptr, " %d", &stat[i]);
		      break;
		    case 2:
		      sscanf (ptr, " %d", &h_ind[i]);
		      break;
		    case 3:
		      sscanf (ptr, " %d", &k_ind[i]);
		      break;
		    case 4:
		      sscanf (ptr, " %d", &l_ind[i]);
		      break;
		    case 5:
		      sscanf (ptr, " %lf", &ehkl[i]);
		      break;
		    case 6:
		      sscanf (ptr, " %lf", &phase0[i]);
		      break;
		    default:
		      fprintf (stdout, "Check reflection list\n");
		    }
		  /* fprintf(stdout," %d %d %s %d %s\n", i, count, buff, n, ptr); */
		  if (count == 6)
		    {
		      count = 0;
		      i++;
		    }
		  count++;
		}
	      n++;
	    }

	}
    }

  for (j = 1; j <= valid_lines; j++)
    {
      BuildEq_hkl (&SgInfo, 1, &Eq_hkl[j], h_ind[j],
		   k_ind[j], l_ind[j]);
      fprintf (stdout, "%4d %4d %4d %4d %4d  %8.3f %8.3f\n", valid_lines,
	       stat[j], h_ind[j], k_ind[j],
	       l_ind[j], ehkl[j],
	       phase0[j]);
    }
  fclose (in);
  return (1);
}


 /* 
    Returns the number of symmetry dependent reflections and 
    reads the lattice parameter, hkl, E(h), and phases from the
    hkl.inp file and build up the structure Eq_hkl
  */


int 
CountHklInpLine (void)
{

  TOKEN testtok;
  FILE *in;
  int i, n;
  char buff[500];
  char *ptr, *comand;
  int count;

  if ((in = fopen ("sp.inp", "r")) == NULL)
    exit (1);
  count = -1;
  i = 0;
  while (EOF != fgetline (in, buff))
    {

      if (buff[0] != '#')
	{
	  n = 0;

	  testtok.actual_token = testtok.next_token = buff;
	  while (NULL != (ptr = str_token (&testtok)))
	    {
	      if (n == 0)
		comand = ptr;
	      if ((strcmp (comand, "End") == 0) && (n == 0))
		count++;
	      if (count > -1)
		{
		  if (count == 6)
		    {
		      count = 0;
		      i++;
		    }
		  count++;
		}
	      n++;
	    }
	}

    }
  fclose (in);
  return (i);
}


/*
   Reads the possible phases from phase.inp and returns their number
 */


int 
ReadPhaseInp (double *perm_phase)
{

  FILE *ph;
  int i;
  char buff[500];
  int valid_lines;


  i = 0;

  valid_lines = CountPhaseInpLine ();

  if ((ph = fopen ("phase.inp", "r")) == NULL)
    {
      fprintf (stdout, "Cannot open phase.inp\n");
      exit (1);
    }

  while (EOF != fgetline (ph, buff))
    if (buff[0] != '#')
      {
	i++;
	sscanf (buff, "%lf", perm_phase + i);
	/* fprintf(stdout, "phase[%4d] %4.3f\n", i, perm_phase[i] );
	 */

      }

  fclose (ph);

  return (valid_lines);

}

/*
   Reads the possible phases from phase_comb.inp 
 */


void 
ReadPhaseCombInp (double *phase_comb)
{

  FILE *ph;
  int i;

  if ((ph = fopen ("phase_comb.inp", "r")) == NULL)
    {
      fprintf (stdout, "Cannot open phase_comb.inp\n");
      exit (1);
    }

  for (i = 0; i < 112; i++)
    fscanf (ph, "%lf", phase_comb + i);

  fclose (ph);

  return;

}

int 
CountPhaseInpLine (void)
{

  FILE *ph;
  char buff[500];
  int valid_lines = 0;

  if ((ph = fopen ("phase.inp", "r")) == NULL)
    {
      fprintf (stdout, "Cannot open phase.inp");
      exit (1);
    }


  while (EOF != fgetline (ph, buff))
    {
      if (buff[0] != '#')
	valid_lines++;
    }

  fclose (ph);

  return (valid_lines);

}



void 
FillHklAlloc (void)
{
  int count, i1, i2;
  count = 0;

  for (i1 = 1; i1 <= NumIndepRefl; i1++)
    {
      for (i2 = 0; i2 < Eq_hkl[i1].N; i2++)
	{
	  count++;
	  /* fprintf( stdout, "count = %d\n", count ); */
	  h[count] = Eq_hkl[i1].h[i2];
	  k[count] = Eq_hkl[i1].k[i2];
	  l[count] = Eq_hkl[i1].l[i2];
	}
    }
  return;

}

static int 
str_ibegin (const char *s1, const char *s2)	/* string ignore-case */
{				/* begin              */

  char u1, u2;


  while (*s1 && *s2)
    {
      u1 = toupper (*s1++);
      u2 = toupper (*s2++);
      if (u1 < u2)
	return -1;
      else if (u1 > u2)
	return 1;

    }

  if (*s2)
    return -1;
  return 0;

}

int 
BuildSgInfo (T_SgInfo * SgInfo, const char *SgName)

{

  int VolLetter;
  const T_TabSgName *tsgn;

  /* look for "VolA", "VolI", or "Hall"

   */

  while (*SgName && isspace (*SgName))
    SgName++;
  VolLetter = -1;

  if (isdigit (*SgName))
    VolLetter = 'A';
  else if (str_ibegin (SgName, "VolA") == 0)
    {
      VolLetter = 'A';
      SgName += 4;
    }
  else if (str_ibegin (SgName, "VolI") == 0
	   || str_ibegin (SgName, "Vol1") == 0)
    {

      VolLetter = 'I';
      SgName += 4;

    }

  else if (str_ibegin (SgName, "Hall") == 0)

    {
      VolLetter = 0;
      SgName += 4;
    }

  while (*SgName && isspace (*SgName))
    SgName++;

  /* default is "VolA"

   */

  if (VolLetter == -1)
    VolLetter = 'A';

  /* if we don't have a Hall symbol do a table look-up

   */

  tsgn = NULL;
  if (VolLetter)

    {

      tsgn = FindTabSgNameEntry (0, SgName, VolLetter);
      if (tsgn == NULL)
	return -1;		/* no matching table entry */
      SgName = tsgn->HallSymbol;

    }


  /* Allocate memory for the list of Seitz matrices and
     a supporting list which holds the characteristics of
     the rotation parts of the Seitz matrices
   */

  SgInfo->MaxList = 192;	/* absolute maximum number of symops */
  SgInfo->ListSeitzMx
    = malloc (SgInfo->MaxList * sizeof (*SgInfo->ListSeitzMx));

  if (SgInfo->ListSeitzMx == NULL)
    {
      SetSgError ("Not enough core");
      return -1;

    }

  SgInfo->ListRotMxInfo
    = malloc (SgInfo->MaxList * sizeof (*SgInfo->ListRotMxInfo));

  if (SgInfo->ListRotMxInfo == NULL)
    {
      SetSgError ("Not enough core");
      return -1;

    }



  /* Initialize the SgInfo structure

   */

  ResetSgInfo (SgInfo);
  SgInfo->TabSgName = tsgn;	/* in case we know the table entry */

  /* Translate the Hall symbol and generate the whole group

   */


  ParseHallSymbol (SgName, SgInfo);
  if (SgError != NULL)
    return -1;

  /* Do some book-keeping and derive crystal system,  group,
     and - if not already set - find the entry in the internal
     table of space group symbols

   */
  return (0);

}


/*
   Slow fourier transformation, not very clever, sorry
   It can be fixed, which reflection are invoved for the calculation
 */

void 
fourier (int *grid_dim, double *roxyz, int begin, int end)
{

  int UNIT_PART = 1;

  int iD;
  int i, j, m, n, t;
  double realpart, phase;
  double x, y, z, STEPX, STEPY, STEPZ;
  int hkl_line, first_hkl;

  if (end == 0)
    hkl_line = CountHklInpLine ();
  else
    hkl_line = end;

  if (begin == 0)
    first_hkl = 1;
  else
    first_hkl = begin;


  x = 0.0;
  y = 0.0;
  z = 0.0;

  STEPX = 1 / (((double) grid_dim[1] - 1) * UNIT_PART);
  STEPY = 1 / (((double) grid_dim[2] - 1) * UNIT_PART);
  STEPZ = 1 / (((double) grid_dim[3] - 1) * UNIT_PART);

  /* fprintf(stdout,"reflection number : %4d  %4d\n", first_hkl, hkl_line); */



  for (i = 0; i < grid_dim[1] * grid_dim[4]; i++)
    {
      if (i > 0)
	x = x + STEPX;
      y = 0.0;


      for (j = 0; j < grid_dim[2] * grid_dim[5]; j++)
	{
	  if (j > 0)
	    y = y + STEPY;
	  z = 0.0;

	  for (m = 0; m < grid_dim[3] * grid_dim[6]; m++)
	    {

	      if (m > 0)
		z = z + STEPZ;
	      realpart = 0.0;


	      for (n = first_hkl; n <= hkl_line; n++)

		for (t = 0; t < Eq_hkl[n].N; t++)

		  {


		    phase = 2 * PI * (Eq_hkl[n].h[t] * x + Eq_hkl[n].k[t] * y
				      + Eq_hkl[n].l[t] * z - (phase0[n]
					+ (double) Eq_hkl[n].TH[t] / STBF));
		    realpart = realpart + ehkl[n] * cos (phase);

		  }

	      iD = (grid_dim[3] * grid_dim[6] * j + m) + (grid_dim[2] * grid_dim[5]
					   * grid_dim[3] * grid_dim[6] * i);
	      roxyz[iD] = 2 * realpart;


	    }
	}
    }
  /* fprintf(stdout,"\n\n\n"); 
     for(t=0; t < Eq_hkl[hkl_line].N ; t++)    
     fprintf(stdout," phases in fourier : %6.2f\n", (phase0[4]
     + (double)Eq_hkl[hkl_line].TH[t]/STBF)); */
  return;
}


/****************Peakiness calculation*****************************/

double 
Peakiness (int *grid_dim, double *roxyz)
{

  int UNIT_PART = 1;

  int iD;
  int i, j, m, n, t;
  double realpart, phase;
  double x, y, z, STEPX, STEPY, STEPZ;
  double peakiness;

  peakiness = 0.0;
  x = 0.0;
  y = 0.0;
  z = 0.0;

  STEPX = 1 / (((double) grid_dim[1] - 1) * UNIT_PART);
  STEPY = 1 / (((double) grid_dim[2] - 1) * UNIT_PART);
  STEPZ = 1 / (((double) grid_dim[3] - 1) * UNIT_PART);

  /* fprintf(stdout,"reflection number : %4d  %4d\n", first_hkl, hkl_line); */



  for (i = 0; i < grid_dim[1]; i++)
    {
      if (i > 0)
	x = x + STEPX;
      y = 0.0;


      for (j = 0; j < grid_dim[2]; j++)
	{
	  if (j > 0)
	    y = y + STEPY;
	  z = 0.0;

	  for (m = 0; m < grid_dim[3]; m++)
	    {

	      if (m > 0)
		z = z + STEPZ;
	      realpart = 0.0;


	      for (n = 2; n <= NumIndepRefl; n++)

		for (t = 0; t < Eq_hkl[n].N; t++)

		  {


		    phase = 2 * PI * (Eq_hkl[n].h[t] * x + Eq_hkl[n].k[t] * y
				      + Eq_hkl[n].l[t] * z - (phase0[n]
					+ (double) Eq_hkl[n].TH[t] / STBF));
		    realpart = realpart + ehkl[n] * cos (phase);

		  }

	      iD = (grid_dim[3] * j + m) + (grid_dim[2] * grid_dim[3] * i);
	      roxyz[iD] = 2 * realpart;


	      peakiness = peakiness + roxyz[iD] * roxyz[iD] * roxyz[iD];


	    }
	}
    }

  return peakiness;
}





/*
   Puts pout a fourier map readable by Cerius module "Isosurface"
 */

void 
OutputFourier (double *roxyz, T_CellParam * CellParam, int *grid_dim, int map_num)
{

  FILE *dens;
  int i;
  char simon[20];

  simon[0] = (65 + map_num);
  simon[1] = '.';
  simon[2] = 'g';
  simon[3] = 'r';
  simon[4] = 'd';
  if ((dens = fopen (simon, "w")) == NULL)
    exit (1);
  fprintf (dens, "\n\n");
  fprintf (dens, "%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", grid_dim[4] * CellParam->a,
	   grid_dim[5] * CellParam->b, grid_dim[6] * CellParam->c,
	   CellParam->alpha, CellParam->beta, CellParam->gamma);
  fprintf (dens, "%4d %4d %4d \n", grid_dim[1] * grid_dim[4] - 1,
	   grid_dim[2] * grid_dim[5] - 1, grid_dim[3] * grid_dim[6] - 1);
  fprintf (dens, "\n");

  for (i = 0; i < grid_dim[1] * grid_dim[4] * grid_dim[2] * grid_dim[5]
       * grid_dim[3] * grid_dim[6]; i++)

    fprintf (dens, " %6.3f\n", roxyz[i]);
  fclose (dens);
  return;
}




int 
IsTriplATripl (int h1, int k1, int l1,
	       int h2, int k2, int l2,
	       int h3, int k3, int l3)
{



  if (((h1 + h2 + h3 == 0) && (k1 + k2 + k3 == 0) && (l1 + l2 + l3 == 0)) ||
      ((h1 + h2 - h3 == 0) && (k1 + k2 - k3 == 0) && (l1 + l2 - l3 == 0)) ||
      ((h1 - h2 + h3 == 0) && (k1 - k2 + k3 == 0) && (l1 - l2 + l3 == 0)) ||
      ((h1 - h2 - h3 == 0) && (k1 - k2 - k3 == 0) && (l1 - l2 - l3 == 0)) ||
   ((-h1 + h2 + h3 == 0) && (-k1 + k2 + k3 == 0) && (-l1 + l2 + l3 == 0)) ||
   ((-h1 + h2 - h3 == 0) && (-k1 + k2 - k3 == 0) && (-l1 + l2 - l3 == 0)) ||
   ((-h1 - h2 + h3 == 0) && (-k1 - k2 + k3 == 0) && (-l1 - l2 + l3 == 0)) ||
   ((-h1 - h2 - h3 == 0) && (-k1 - k2 - k3 == 0) && (-l1 - l2 - l3 == 0)) ||
      ((-h1 - h2 - h3 == 0) && (-k1 - k2 - k3 == 0) && (-l1 - l2 - l3 == 0)))
    return (1);

  else
    return (0);

}

int 
IsQuartAQuart (int h1, int k1, int l1,
	       int h2, int k2, int l2,
	       int h3, int k3, int l3,
	       int h4, int k4, int l4)
{



  if (((h1 + h2 + h3 + h4 == 0) && (k1 + k2 + k3 + k4 == 0) && (l1 + l2 + l3 + l4 == 0)) ||
      ((h1 + h2 + h3 - h4 == 0) && (k1 + k2 + k3 - k4 == 0) && (l1 + l2 + l3 - l4 == 0)) ||
      ((h1 + h2 - h3 + h4 == 0) && (k1 + k2 - k3 + k4 == 0) && (l1 + l2 - l3 + l4 == 0)) ||
      ((h1 + h2 - h3 - h4 == 0) && (k1 + k2 - k3 - k4 == 0) && (l1 + l2 - l3 - l4 == 0)) ||
      ((h1 - h2 + h3 + h4 == 0) && (k1 - k2 + k3 + k4 == 0) && (l1 - l2 + l3 + l4 == 0)) ||
      ((h1 - h2 + h3 - h4 == 0) && (k1 - k2 + k3 - k4 == 0) && (l1 - l2 + l3 - l4 == 0)) ||
      ((h1 - h2 - h3 + h4 == 0) && (k1 - k2 - k3 + k4 == 0) && (l1 - l2 - l3 + l4 == 0)) ||
      ((h1 - h2 - h3 - h4 == 0) && (k1 - k2 - k3 - k4 == 0) && (l1 - l2 - l3 - l4 == 0)))
    /* ((-h1+h2+h3+h4 == 0)&& (-k1+k2+k3+k4==0)  && (-l1+l2+l3+l4==0)) ||
       ((-h1+h2+h3-h4 == 0) &&  (-k1+k2+k3-k4==0) && (-l1+l2+l3-l4==0)) ||
       ((-h1+h2-h3+h4 == 0) &&  (-k1+k2-k3+k4==0) && (-l1+l2-l3+l4==0)) ||
       ((-h1+h2-h3-h4 == 0) &&  (-k1+k2-k3-k4==0)  && (-l1+l2-l3-l4==0)) ||   
       ((-h1-h2+h3+h4 == 0) &&  (-k1-k2+k3+k4==0) &&  (-l1-l2+l3+l4==0)) ||
       ((-h1-h2+h3-h4 == 0) &&  (-k1-k2+k3-k4==0)  &&  (-l1-l2+l3-l4==0)) ||
       ((-h1-h2-h3+h4 == 0) &&  (-k1-k2-k3+k4==0)  &&  (-l1-l2-l3+l4==0)) ||
       ((-h1-h2-h3-h4 == 0) &&  (-k1-k2-k3-k4==0)  &&  (-l1-l2-l3-l4==0))) */



    return (1);

  else
    return (0);

}

void 
MakeNormalization (double *SayreE, double *ObsE, int ClustMem)
{
  int i;


  double SumSayre;
  double SumObs;



  SumSayre = 0;
  SumObs = 0;



  for (i = 1; i <= ClustMem; i++)
    {

      SumSayre = SumSayre + SayreE[i];
      SumObs = SumObs + ObsE[i];
    }



  for (i = 1; i <= ClustMem; i++)
    {

      ObsE[i] = (SumObs / SumSayre) * SayreE[i] * (SayreE[i] / fabs (SayreE[i]));


    }

  return;
}


/* returns the  reziprocal d-value */

double 
CalculateH (int h, int k, int l)
{

  double a, b, c, alpha, beta, gamma;
  double a2, b2, c2;
  double q, q1, q2;
  double sin2alpha, sin2beta, sin2gamma;
  double cosalpha, cosbeta, cosgamma;
  double cos2alpha, cos2beta, cos2gamma;
  int h2, k2, l2;

  a = CellParam->a;
  b = CellParam->b;
  c = CellParam->c;
  alpha = (CellParam->alpha / 360) * TWO_PI;
  beta = (CellParam->beta / 360) * TWO_PI;
  gamma = (CellParam->gamma / 360) * TWO_PI;

  a2 = a * a;
  b2 = b * b;
  c2 = c * c;

  h2 = h * h;
  k2 = k * k;
  l2 = l * l;

  sin2alpha = sin (alpha) * sin (alpha);
  cosalpha = cos (alpha);
  sin2beta = sin (beta) * sin (beta);
  cosbeta = cos (beta);
  sin2gamma = sin (gamma) * sin (gamma);
  cosgamma = cos (gamma);
  cos2alpha = cosalpha * cosalpha;
  cos2beta = cosbeta * cosbeta;
  cos2gamma = cosgamma * cosgamma;

  /* fprintf(stdout,"%6.3f %6.3f %6.3f\n%6.3f %6.3f %6.3f", 
     alpha, beta, gamma, cosalpha, cosbeta, cosgamma); */

  q1 = b2 * c2 * sin2alpha * h2 + c2 * a2 * sin2beta * k2 + a2 * b2 * sin2gamma * l2 +
    2 * a * b * c2 * (cosalpha * cosbeta - cosgamma) * h * k +
    2 * a * b2 * c * (cosalpha * cosgamma - cosbeta) * h * l +
    2 * a2 * b * c * (cosbeta * cosgamma - cosalpha) * k * l;

  q2 = a2 * b2 * c2 * (1 - cos2alpha - cos2beta - cos2gamma + 2 * cosalpha * cosbeta * cosgamma);

  if (q2 < 1.0e-4)
    {
      fprintf (stdout, "Check your lattice parameter\n");
      exit (1);
    }
  q = sqrt (q1 / q2);
  return q;

}




void 
ToAtomsOneDim (double *SfacSquare, double *Sfac)
{
  int i;
  double Radius, Radius2;


  for (i = 1; i <= NumIndepRefl; i++)
    {
      Radius = CalculateH (Eq_hkl[i].h[0],
			   Eq_hkl[i].k[0], Eq_hkl[i].l[0]);

      Radius2 = Radius * Radius;
      Sfac[i] = ehkl[i] * 0.353553 * exp (-0.5 * PI * Radius2);
      SfacSquare[i] = ehkl[i] * 0.125 * exp (-0.25 * PI * Radius2);

    }
}


void 
PointToAtoms (double *SfacSquare, double *Sfac, double PAWidth)
{
  int i;
  double Radius, Radius2;
  double b;

  b = SQRT_2 / (PI * PAWidth);

  for (i = 1; i <= NumIndepRefl; i++)
    {
      Radius = CalculateH (Eq_hkl[i].h[0],
			   Eq_hkl[i].k[0], Eq_hkl[i].l[0]);

      Radius2 = Radius * Radius;

      Sfac[i] = ehkl[i] * (2 / (sqrt (b) * b)) * exp (-(1 / b) * PI * Radius2);
      SfacSquare[i] = ehkl[i] * (1 / (SQRT_2 * b * sqrt (b))) * exp (-(1 / (2 * b)) * PI * Radius2);

    }
  return;
}


/* structur factor calculation for Sayre Woolfson Approach */

void 
PointToAtomsWoolfson (double *SfacP, double *SfacPSquare, double *SfacPCub, double weightP,
       double *SfacQ, double *SfacQSquare, double *SfacQCub, double weightQ)
{
  int i;
  double Radius, Radius2;
  double expfac1, expfac2, expfac3, weightP1, weightQ1;


  for (i = 1; i <= NumIndepRefl; i++)
    {
      Radius = CalculateH (Eq_hkl[i].h[0],
			   Eq_hkl[i].k[0], Eq_hkl[i].l[0]);

      Radius2 = Radius * Radius;
      /*scattering factors of the atom, the squared and the cubes Atom */

      expfac1 = exp (-10.47197551 * Radius2);
      expfac1 = expfac1 * expfac1 * expfac1;
      expfac2 = exp (-5.235987758 * Radius2);
      expfac2 = expfac2 * expfac2 * expfac2;
      expfac3 = exp (-3.490658504 * Radius2);
      expfac3 = expfac3 * expfac3 * expfac3;
      weightP1 = weightP * weightP * weightP;
      weightQ1 = weightQ * weightQ * weightQ;


      SfacP[i] = 6.085806191 * weightP1 * expfac1;

      SfacPSquare[i] = weightP1 * weightP1;
      SfacPSquare[i] *= 2.151657416 * expfac2;

      SfacPCub[i] = weightP1 * weightP1 * weightP1;
      SfacPCub[i] *= 1.171213947 * expfac3;



      SfacQ[i] = 6.085806191 * weightQ1 * expfac1;

      SfacQSquare[i] = weightQ1 * weightQ1;
      SfacQSquare[i] *= 2.151657416 * expfac2;

      SfacQCub[i] = weightQ1 * weightQ1 * weightQ1;
      SfacQCub[i] *= 1.171213947 * expfac3;

    }


}

void 
AtomsToPoints (double *SfacSquare, double *Sfac)
{
  int i;
  double Radius, Radius2;
  double expfac1, expfac2;
  double *edecon;

  edecon = malloc ((NumIndepRefl + 1) * sizeof (double));

  for (i = 1; i <= NumIndepRefl; i++)
    if (stat[i] > 4)
      {
	Radius = CalculateH (Eq_hkl[i].h[0],
			     Eq_hkl[i].k[0], Eq_hkl[i].l[0]);

	Radius2 = Radius * Radius;
	expfac1 = exp (-10.4719 * Radius2);
	expfac2 = exp (-5.235 * Radius2);
	edecon[i] = SfacSquare[i] / (2.151 * expfac2 * expfac2 * expfac2);
	Sfac[i] = edecon[i] * 6.08 * expfac1 * expfac1 * expfac1;

	/* fprintf (stdout, " in AtomsToPoints  %5.3f  %6.4f   %6.4f   %6.4f   %6.4f \n", Sfac[i], SfacSquare[i]); */


      }
  return;
}

void 
RemoveTempScattPart (double B, double Scale, int ENum, int AtKind)
{
  double *SfacBuff;
  double H, H2 = 0.0;
  int i = 0;

  SfacBuff = malloc ((NumIndepRefl + 1) * sizeof (double));
  for (i = 1; i <= NumIndepRefl; i++)
    {
      H = CalculateH (Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0]);
      H2 = H * H;

      if (i > 1)
	SfacBuff[i] = (ehkl[i] * Scale) / (exp (-B * 0.25 * H2));
      else
	SfacBuff[i] = ehkl[i] / (exp (-B * 0.25 * H2));
    }

  for (i = 1; i <= NumIndepRefl; i++)
    ehkl[i] = SfacBuff[i];

  free (SfacBuff);

  return;

}


int 
CountSymmDepRefl (void)
{

  int count, i, j;

  count = 0;
  for (i = 1; i <= NumIndepRefl; i++)
    for (j = 0; j < Eq_hkl[i].N; j++)
      count++;

  return (count);
}


void 
WriteNewHklInp (T_CellParam * CellParam, double *Sfac)
{

  int i;
  double phase_ratio;

  num_plus = num_minus = 0;

  fprintf (stdout, "# Title\n");
  fprintf (stdout, "# Title\n");
  fprintf (stdout, "%9.3f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
	   CellParam->a, CellParam->b, CellParam->c,
	   CellParam->alpha, CellParam->beta, CellParam->gamma);
  fprintf (stdout, "#       h   k   l   |F(hkl)|    phase\n");

  for (i = 1; i <= NumIndepRefl; i++)
    {
      if (phase0[i] == 0.0)
	num_plus++;
      if (phase0[i] == 0.5)
	num_minus++;
      fprintf (stdout, "%4d %6d %4d %4d %9.3f %6.2f\n",
	       stat[i], Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
	       ehkl[i], phase0[i]);

    }
  phase_ratio = (double) num_plus / num_minus;	/*avoids the minus or plus catastroph */
  fprintf (stdout, "PhaseSet");
  for (i = 1; i <= NumIndepRefl; i++){
    if (stat[i] == 1){
      fprintf (stdout, "%6.2f", phase0[i]);}}


  fprintf (stdout, "\n");


  return;
}





void 
PositivityProbab (int atom_number)
{

  int i1, i2, i3;
  double positivity, E_sum;
  int count_trip, count;
  double *intern_ehkl;
  int *stat_intern;

  intern_ehkl = malloc ((NumDepRefl + 1) * sizeof (double));
  stat_intern = malloc ((NumDepRefl + 1) * sizeof (int));

  for (i1 = 2; i1 <= NumIndepRefl; i1++)
    for (i2 = 0; i2 < Eq_hkl[i1].N; i2++)
      {
	count++;
	intern_ehkl[count] = norm_ehkl[i1] *
	  cos ((phase0[i1] + Eq_hkl[i1].TH[i2] / STBF) * TWO_PI);
	stat_intern[count] = stat[i1];

      }

  for (i1 = 2; i1 <= NumIndepRefl; i1++)
    {

      E_sum = 0;
      count_trip = 0;

      for (i2 = 2; i2 <= NumDepRefl; i2++)
	for (i3 = 2; i3 <= NumDepRefl; i3++)
	  {

	    /* fprintf(stdout," ##%7.3f %7.3f \n", intern_ehkl[i2], intern_ehkl[i3]); */

	    if (IsTriplATripl (Eq_hkl[i1].h[0], Eq_hkl[i1].k[0],
			       Eq_hkl[i1].l[0],
			       h[i2], k[i2], l[i2],
			       h[i3], k[i3], l[i3]) == 1)
	      {
		if ((stat[i1] != 0) && (stat_intern[i2] == 0) && (stat_intern[i3] == 0))
		  {
		    E_sum = E_sum + intern_ehkl[i2] * intern_ehkl[i3];
		    count_trip++;
		  }
	      }
	  }

      positivity = 0.5 + (0.5 * tanh ((1 / sqrt (atom_number)) *
				      fabs (ehkl[i1]) * E_sum));


      if (count_trip > 0)

	{
	  if (positivity <= 0.5)
	    phase0[i1] = 0.5;
	  else
	    phase0[i1] = 0.0;
	}


      fprintf (stdout, "prob : %4d %4d %4d %4.2f %4d %7.3f %7.3f %4d\n",
	       Eq_hkl[i1].h[0], Eq_hkl[i1].k[0],
	       Eq_hkl[i1].l[0],
	       positivity, count_trip,
	       E_sum, phase0[i1], stat[i1]);



    }


  free (intern_ehkl);
  free (stat_intern);
  return;
}



/*****************************************************************************
 test of the program using the data of the woolfson paper (1957)
 The same notation as in the paper is used
******************************************************************************/



void 
MakeWoolfson (double escale, double volume,
	      double *fp, double *fp2, double *fp3,
	      double *fq, double *fq2, double *fq3,
	      int Reflexstat1, int Reflexstat2, int Reflexstat3,
	      int Reflexstat4)

{



  int count, count_trip, count_quart;
  int i, j, m, n;
  double imm_sumA, imm_sumB;
  double *Ahkl, *Bhkl;
  double E_SumASquare, E_SumBSquare, E_SumACub, E_SumBCub;
  double tmp, tmp1, tmp2;
  double *PhaseSayre;
  int *stat_intern;
  double Hh, Gh, Bs, As, Fh, R, Robs;
  double dval;

  count = 0;
  count_trip = count_quart = 0;
  E_SumASquare = E_SumACub = 0.0;
  E_SumBSquare = E_SumBCub = 0.0;
  R = Robs = 0.000;
  PhaseSayre = malloc ((NumDepRefl + 1) * sizeof (double));
  Ahkl = malloc ((NumDepRefl + 1) * sizeof (double));
  Bhkl = malloc ((NumDepRefl + 1) * sizeof (double));
  stat_intern = malloc ((NumDepRefl + 1) * sizeof (int));

  for (i = 1; i <= NumIndepRefl; i++)
    for (j = 0; j < Eq_hkl[i].N; j++)
      {
	count++;
	Ahkl[count] = ehkl[i] * fq[i] * cos ((phase0[i] + Eq_hkl[i].TH[j] / STBF) * TWO_PI);
	Bhkl[count] = ehkl[i] * fq[i] * sin ((phase0[i] + Eq_hkl[i].TH[j] / STBF) * TWO_PI);
	stat_intern[count] = stat[i];
      }

  fprintf (stdout, "                  numt nubq     As       -Bs       Gh       Hh         fp       ehkl       As*Gh\n");

  count = 0;
  for (i = 1; i <= NumIndepRefl; i++)
    {

      dval = CalculateH (Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0]);
      dval = dval * dval;
      /* fprintf(stdout," fp %9.3f fp2 %9.3f fp3 %9.3f fq  %9.3f fq2 %9.3f  fq3 %9.3f\n", fp[i], fp2[i],
         fp3[i], fq[i], fq2[i], fq3[i]); */

      if ((stat[i] != Reflexstat1) && (stat[i] != Reflexstat2) &&
	  (stat[i] != Reflexstat3) && (stat[i] != Reflexstat4))
	{

	  for (j = 1; j <= NumDepRefl; j++)

	    if ((stat_intern[j] != Reflexstat1) && (stat_intern[j] != Reflexstat2) &&
		(stat_intern[j] != Reflexstat3) && (stat[i] != Reflexstat4))

	      for (n = 1; n <= NumDepRefl; n++)

		if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
		(stat_intern[n] != Reflexstat3) && (stat[i] != Reflexstat4))


		  if (IsTriplATripl (Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
				     h[j], k[j], l[j],
				     h[n], k[n], l[n]) == 1)


		    {

		      /* fprintf(stdout," %4d%4d%4d  %4d%4d%4d   %4d%4d%4d\n",
		         Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
		         h[j], k[j], l[j],
		         h[n], k[n], l[n]);  */

		      count_trip++;



		      imm_sumA = Ahkl[j] * Ahkl[n] - Bhkl[j] * Bhkl[n];
		      imm_sumB = Ahkl[j] * Bhkl[n] + Bhkl[j] * Ahkl[n];
		      /*        fprintf(stdout,"#%4d%4d %9.3f %9.3f \n", j, n, imm_sumA); */

		      E_SumASquare = E_SumASquare + imm_sumA;
		      E_SumBSquare = E_SumBSquare + imm_sumB;
		      /* fprintf(stdout,"#%9.3f %4d \n", E_SumASquare, i); */
		    }


	  if (fq[i] != fp[i])
	    {

	      for (j = 1; j <= NumDepRefl; j++)

		if ((stat_intern[j] != Reflexstat1) && (stat_intern[j] != Reflexstat2) &&
		(stat_intern[j] != Reflexstat3) && (stat[i] != Reflexstat4))

		  for (n = 1; n <= NumDepRefl; n++)
		    if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
			(stat_intern[n] != Reflexstat3) && (stat[i] != Reflexstat4))


		      for (m = 1; m <= NumDepRefl; m++)
			if ((stat_intern[n] != Reflexstat1) && (stat_intern[n] != Reflexstat2) &&
			    (stat_intern[n] != Reflexstat3) && (stat[i] != Reflexstat4))

			  if (IsQuartAQuart (Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
					     h[j], k[j], l[j],
					     h[n], k[n], l[n],
					     h[m], k[m], l[m]) == 1)
			    {
			      count_quart++;
			      imm_sumA = Ahkl[j] * Ahkl[n] * Ahkl[m] - Ahkl[j] * Bhkl[n] * Bhkl[m]
				- Bhkl[j] * Ahkl[n] * Bhkl[m] - Bhkl[j] * Bhkl[n] * Ahkl[m];
			      imm_sumB = Bhkl[j] * Bhkl[n] * Bhkl[m] + Ahkl[j] * Ahkl[n] * Bhkl[m]
				+ Ahkl[j] * Bhkl[n] * Ahkl[m] + Bhkl[j] * Ahkl[n] * Ahkl[m];

			      E_SumACub = E_SumACub + imm_sumA;
			      E_SumBCub = E_SumBCub + imm_sumB;
			    }
	    }


	  if ((Eq_hkl[i].h[0] == 0) && (Eq_hkl[i].k[0] == 0) && (Eq_hkl[i].l[0] == 0))
	    {
	      Gh = (E_SumASquare * 2 - (Ahkl[1] * Ahkl[1])) / volume;
	      Hh = (E_SumACub * 2 - Ahkl[1] * Ahkl[1] * Ahkl[1]) / (volume * volume);
	    }

	  else
	    {
	      Gh = E_SumASquare / (volume);
	      Hh = E_SumACub / (volume * volume);
	    }



	  E_SumASquare = E_SumACub = E_SumBCub = E_SumBSquare = 0.0;
	  tmp = fp2[i] / fq2[i];
	  tmp1 = fq[i] * tmp;
	  tmp2 = fq3[i] * tmp;
	  /* fprintf(stdout,"(fp[i] - tmp1) %9.3f  (fp3[i] - tmp2) %9.3f \n",fp[i] -
	     tmp1,fp3[i] - tmp2 ); */

	  if (tmp != 1.0000)
	    Bs = (fp[i] - tmp1) / (fp3[i] - tmp2);
	  else
	    Bs = 0.00;

	  As = (fq[i] - Bs * fq3[i]) / fq2[i];

	  /* fprintf(stdout," As%9.3f  Bs%9.3f Gh%9.3f Hh %9.3f\n",As, -Bs, Gh, Hh); */
	  Fh = As * Gh + Bs * Hh;
	  if (i == 1)
	    {
	      escale = Fh / ((fq[i] + fp[i]) / 2.0 * ehkl[i]);
	      /* fprintf(stdout, "escale %9.3f    ehkl %9.3f   Fh %9.3f \n", escale,
	         ehkl[i],Fh); */ 
	    }
	  Fh /= escale;

	  /*   if(i==1) 
	     escale = Fh/ehkl[i];


	     Fh = Fh/escale;
	   */

	  if (Fh < 0)
	    PhaseSayre[i] = 0.5;
	  else
	    PhaseSayre[i] = 0.0;

	  fprintf (stdout, " E(%2d %4d %4d) %4d %4d  %9.3f %9.4f %9.3f %9.3f %9.3f %9.3f %9.3f\n",
		   Eq_hkl[i].h[0], Eq_hkl[i].k[0], Eq_hkl[i].l[0],
		   count_trip, count_quart,
		   As, -Bs, Gh, Hh, ehkl[i] * (fq[i] + fp[i]) / 2.0, Fh, PhaseSayre[i]);


	  R = R + fabs (fabs (Fh) - fabs (ehkl[i] * (fq[i] + fp[i]) / 2.0));
	  Robs = Robs + fabs (ehkl[i] * (fq[i] + fp[i]) / 2.0);

	  count_trip = 0;
	  count_quart = 0;

	}
    }
  printf (" Relative deviation of E(values)  : %6.2f \n ", R / Robs);

  R = Robs = 0.0;



  free (Ahkl);
  free (Bhkl);

  return;

}

double 
CalcVolume (T_CellParam * CellParam)
{
  double a, b, c, alpha, beta, gamma;
  double salpha, calpha, sbeta, cbeta, sgamma, cgamma;
  double V;
  double alphar;
  double salphar;

  a = CellParam->a;
  b = CellParam->b;
  c = CellParam->c;
  alpha = CellParam->alpha * RAD;
  beta = CellParam->beta * RAD;
  gamma = CellParam->gamma * RAD;

  salpha = sin (alpha);
  calpha = cos (alpha);
  sbeta = sin (beta);
  cbeta = cos (beta);
  sgamma = sin (gamma);
  cgamma = cos (gamma);

  alphar = acos ((cbeta * cgamma - calpha) / (sbeta * sgamma));
  salphar = sin (alphar);

  V = a * b * c * salphar * sbeta * sgamma;
  return V;
}
