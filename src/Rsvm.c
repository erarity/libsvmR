#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "svm.h"


static int DEBUG_SVMLIB = 0;

void set_svm_debug(int d){
  DEBUG_SVMLIB = d;
}

static void free_svm_nodes(SEXP s){
  if(TYPEOF(s) != EXTPTRSXP){
    error("argument not external pointer");
  }
  struct svm_node *nodes = (struct svm_node *) R_ExternalPtrAddr(s);
  Free(nodes);
  if(DEBUG_SVMLIB){
    printf("Freeing problem");
  }
}

void svmtrain (double *x, int *r, int *c,
	       double *y,
	       int    *rowindex, int *colindex,
	       int    *svm_type,
	       int    *kernel_type,
	       int    *degree,
	       double *gamma,
	       double *coef0,
	       double *cost,
	       double *nu,
	       int    *weightlabels,
	       double *weights,
	       int    *nweights,
	       double *cache,
	       double *tolerance,
	       double *epsilon,
	       int    *shrinking,
	       int    *cross,
	       int    *sparse,
	       int    *probability,
	       int    *seed,
	       
	       int    *nclasses,
	       int    *nr,
	       int    *index,
	       int    *labels,
	       int    *nSV,
	       double *rho,
	       double *coefs,
	       double *sigma,
	       double *probA,
	       double *probB,

	       double *cresults,
	       double *ctotal1,
	       double *ctotal2,
	       char   **error)
{
  struct svm_parameter par;
  struct svm_problem   prob;
  struct svm_model    *model = NULL;
  int i, ii;
  const char* s;
    
  /* set parameters */
  par.svm_type    = *svm_type;
  par.kernel_type = *kernel_type;
  par.degree      = *degree;
  par.gamma       = *gamma;
  par.coef0       = *coef0;
  par.cache_size  = *cache;
  par.eps         = *tolerance;
  par.C           = *cost;
  par.nu          = *nu;
  par.nr_weight   = *nweights;
  if (par.nr_weight > 0) {
    par.weight      = (double *) malloc (sizeof(double) * par.nr_weight);
    memcpy(par.weight, weights, par.nr_weight * sizeof(double));
    par.weight_label = (int *) malloc (sizeof(int) * par.nr_weight);
    memcpy(par.weight_label, weightlabels, par.nr_weight * sizeof(int));
  }
  par.p           = *epsilon;
  par.shrinking   = *shrinking;
  par.probability = *probability;

  /* set problem */
  prob.l = *r;
  prob.y = y;
    
  if (*sparse > 0)
    prob.x = transsparse(x, *r, rowindex, colindex);
  else
    prob.x = sparsify(x, *r, *c);
    


  /* check parameters & copy error message */
  s = svm_check_parameter(&prob, &par);
  if (s) {
    strcpy(*error, s);
  } else {
    /* set seed */
    srand(*seed);

    /* call svm_train */
    model = svm_train(&prob, &par);
    
    /* set up return values */
    for (ii = 0; ii < model->l; ii++)
      for (i = 0; i < *r;	i++)
	if (prob.x[i] == model->SV[ii]) index[ii] = i+1;
	
    *nr  = model->l;
    *nclasses = model->nr_class;
    memcpy (rho, model->rho, *nclasses * (*nclasses - 1)/2 * sizeof(double));

    if (*probability && par.svm_type != ONE_CLASS) {
      if (par.svm_type == EPSILON_SVR || par.svm_type == NU_SVR)
	*sigma = svm_get_svr_probability(model);
      else {
	memcpy(probA, model->probA,
	       *nclasses * (*nclasses - 1)/2 * sizeof(double));
	memcpy(probB, model->probB,
	       *nclasses * (*nclasses - 1)/2 * sizeof(double));
      }
    }

    for (i = 0; i < *nclasses-1; i++)
      memcpy (coefs + i * *nr, model->sv_coef[i],  *nr * sizeof (double));
	
    if (*svm_type < 2) {
      memcpy (labels, model->label, *nclasses * sizeof(int));
      memcpy (nSV, model->nSV, *nclasses * sizeof(int));
    }
	
    /* Perform cross-validation, if requested */
    if (*cross > 0)
      do_cross_validation (&prob, &par, *cross, cresults,
			   ctotal1, ctotal2);

    /* clean up memory */
    svm_free_and_destroy_model(&model);
  }
    
  /* clean up memory */
  if (par.nr_weight > 0) {
    free(par.weight);
    free(par.weight_label);
  }
    
  for (i = 0; i < *r; i++) free (prob.x[i]);
  free (prob.x);
}
	     
void svmpredict  (int    *decisionvalues,
		  int    *probability,

		  double *v, int *r, int *c,
		  int    *rowindex,
		  int    *colindex,
		  double *coefs,
		  double *rho,
		  int    *compprob,
		  double *probA,
		  double *probB,
		  int    *nclasses,
		  int    *totnSV,
		  int    *labels,
		  int    *nSV,
		  int    *sparsemodel,

		  int    *svm_type,
		  int    *kernel_type,
		  int    *degree,
		  double *gamma,
		  double *coef0,

		  double *x, int *xr,
		  int    *xrowindex,
		  int    *xcolindex,
		  int    *sparsex,
		  
		  double *ret,
		  double *dec,
		  double *prob)
{
  struct svm_model m;
  struct svm_node ** train;
  int i;
    
  /* set up model */
  m.l        = *totnSV;
  m.nr_class = *nclasses;
  m.sv_coef  = (double **) malloc (m.nr_class * sizeof(double));
  for (i = 0; i < m.nr_class - 1; i++) {
    m.sv_coef[i] = (double *) malloc (m.l * sizeof (double));
    memcpy (m.sv_coef[i], coefs + i*m.l, m.l * sizeof (double));
  }
    
  if (*sparsemodel > 0)
    m.SV   = transsparse(v, *r, rowindex, colindex);
  else
    m.SV   = sparsify(v, *r, *c);
    
  m.rho      = rho;
  m.probA    = probA;
  m.probB    = probB;
  m.label    = labels;
  m.nSV      = nSV;

  /* set up parameter */
  m.param.svm_type    = *svm_type;
  m.param.kernel_type = *kernel_type;
  m.param.degree      = *degree;
  m.param.gamma       = *gamma;
  m.param.coef0       = *coef0;
  m.param.probability = *compprob;

  m.free_sv           = 1;

  /* create sparse training matrix */
  if (*sparsex > 0)
    train = transsparse(x, *xr, xrowindex, xcolindex);
  else
    train = sparsify(x, *xr, *c);

  /* call svm-predict-function for each x-row, possibly using probability
     estimator, if requested */
  if (*probability && svm_check_probability_model(&m)) {
    for (i = 0; i < *xr; i++)
      ret[i] = svm_predict_probability(&m, train[i], prob + i * *nclasses);
  } else {
    for (i = 0; i < *xr; i++)
      ret[i] = svm_predict(&m, train[i]);
  }

  /* optionally, compute decision values */
  if (*decisionvalues)
    for (i = 0; i < *xr; i++)
      svm_predict_values(&m, train[i], dec + i * *nclasses * (*nclasses - 1) / 2);

  /* clean up memory */
  for (i = 0; i < *xr; i++)
    free (train[i]);
  free (train);

  for (i = 0; i < *r; i++)
    free (m.SV[i]);
  free (m.SV);
    
  for (i = 0; i < m.nr_class - 1; i++)
    free(m.sv_coef[i]);
  free(m.sv_coef);
}	     
		
void svmwrite (double *v, int *r, int *c,
	       int    *rowindex,
	       int    *colindex,
	       double *coefs,
	       double *rho,
	       double *probA,
	       double *probB,
	       int    *nclasses,
	       int    *totnSV,
	       int    *labels,
	       int    *nSV,
	       int    *sparsemodel,

	       int    *svm_type,
	       int    *kernel_type,
	       int    *degree,
	       double *gamma,
	       double *coef0,

	       char **filename)

{
  struct svm_model m;
  int i;
  char *fname = *filename;

  /* set up model */
  m.l        = *totnSV;
  m.nr_class = *nclasses;
  m.sv_coef  = (double **) malloc (m.nr_class * sizeof(double));
  for (i = 0; i < m.nr_class - 1; i++) {
    m.sv_coef[i] = (double *) malloc (m.l * sizeof (double));
    memcpy (m.sv_coef[i], coefs + i*m.l, m.l * sizeof (double));
  }
    
  if (*sparsemodel > 0)
    m.SV   = transsparse(v, *r, rowindex, colindex);
  else
    m.SV   = sparsify(v, *r, *c);
    
  m.rho      = rho;
  m.label    = labels;
  m.nSV      = nSV;
  m.probA    = probA;
  m.probB    = probB;

  /* set up parameter */
  m.param.svm_type    = *svm_type;
  m.param.kernel_type = *kernel_type;
  m.param.degree      = *degree;
  m.param.gamma       = *gamma;
  m.param.coef0       = *coef0;

  m.free_sv           = 1;

  /* write svm model */
  svm_save_model(fname, &m);

  for (i = 0; i < m.nr_class - 1; i++)
    free(m.sv_coef[i]);
  free(m.sv_coef);


}


void exit_input_error(int line_num)
{
	fprintf(stderr,"Wrong input format at line %d\n", line_num);
	exit(1);
}


SEXP svm_read_problem(const char **fnameptr)
{
  struct svm_problem prob;
  struct svm_node *x_space;
  int elements, max_index, inst_max_index, i, j;
  const char *filename = fnameptr[0];
  FILE *fp = fopen(filename,"r");
  char *endptr;
  char *idx, *val, *label;

  if(fp == NULL)
    {
      error("can't open input file %s\n",filename);
      
    }

  prob.l = 0;
  elements = 0;

  printf("Starting processing\n");

  int max_line_len = 1024;
  char *line = Calloc(max_line_len,char);
  while(readline(fp)!=NULL)
    {
      char *p = strtok(line," \t"); // label

      // features
      while(1)
	{
	  p = strtok(NULL," \t");
	  if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
	    break;
	  ++elements;
	}
      ++elements;
      ++prob.l;
    }
  rewind(fp);

  SEXP y;			/* Easier than coercing raw double* to SEXP */
  PROTECT(y = NEW_NUMERIC(prob.l));
  prob.y = NUMERIC_POINTER(y);  
  prob.x = Calloc(prob.l,struct svm_node *);
  x_space = Calloc(elements, struct svm_node);

  max_index = 0;
  j=0;
  for(i=0;i<prob.l;i++) {
    inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
    readline(fp);
    prob.x[i] = &x_space[j];
    label = strtok(line," \t\n");
    if(label == NULL) // empty line
      exit_input_error(i+1);

    prob.y[i] = strtod(label,&endptr);
    if(endptr == label || *endptr != '\0')
      exit_input_error(i+1);

    int errno;			/* Not sure this is right */
    while(1)
      {
	idx = strtok(NULL,":");
	val = strtok(NULL," \t");

	if(val == NULL)
	  break;

	errno = 0;
	x_space[j].index = (int) strtol(idx,&endptr,10);
	if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
	  exit_input_error(i+1);
	else
	  inst_max_index = x_space[j].index;

	errno = 0;
	x_space[j].value = strtod(val,&endptr);
	if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
	  exit_input_error(i+1);

	++j;
      }

    if(inst_max_index > max_index)
      max_index = inst_max_index;
    x_space[j++].index = -1;
  }

  fclose(fp);
  
  /* R Code */
  SEXP features, nodes;
  PROTECT(features = R_MakeExternalPtr(prob.x, R_NilValue, R_NilValue));  
  PROTECT(nodes = R_MakeExternalPtr(x_space, R_NilValue, R_NilValue));  
  /* Register finalizer */
  R_RegisterCFinalizer(nodes, free_svm_nodes);
  
  SEXP l;
  PROTECT(l = NEW_INTEGER(1));
  (INTEGER_POINTER(l))[0] = prob.l;

  SEXP result;
  PROTECT(result = allocVector(VECSXP,4));
  SET_VECTOR_ELT(result, 0, l);
  SET_VECTOR_ELT(result, 1, y);
  SET_VECTOR_ELT(result, 2, features);
  SET_VECTOR_ELT(result, 3, nodes);

  UNPROTECT(5);
  return result;
}
