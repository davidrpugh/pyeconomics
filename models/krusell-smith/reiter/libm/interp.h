void hunt_double(const double *xx, int n, double x, int *jlo);
int index2n(const int *index, const int *ngrid, int dimen);
void piksr2(int n, double arr[], int brr[]);
int InterpWeightsSimplex(int *NumbGrid,
			 double *Probs,
			 double *State,
			 double **nodes,
			 int *nnodes,
			 int dim);

int InterpWeightsTensor(int *NumbGrid,
			double *Probs,
			double *State,
			double **nodes,
			int *nPerDim,
			int dim,
			int iLinear);
void myassert(int cond,const char *msg);

