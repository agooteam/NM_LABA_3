#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
typedef double TYPE;

#define FRM_STR_OUT "%.15lf\n"
#define FRM_STR_EPS "%e\n"

class matrix{
	protected:
		int N;
		int maxiter;
		TYPE eps;
		TYPE alfa;
		TYPE betta;
		int *ig;
		int *jg;
		TYPE *ggl;
		TYPE *ggu;
		TYPE *di;
		TYPE *cggl;
		TYPE *cggu;
		TYPE *cdi;
		TYPE *pr;
		TYPE *x;
		TYPE *r;
		TYPE *z;
		TYPE *p;
		TYPE *temp;
		TYPE *temp2;
	public:
		void read_kuslau();
		void read_matrix_data();
		void mul_matrix_vector(TYPE *v);
		void dec_LU_sq();
		void allocation_memory();
		void clear_memory();
		void calc_start_values(int solve);
		TYPE calc_otn_nevazka();
		void LOS_clean();
		void LOS_diag();
		int get_maxiter();
		TYPE get_eps();
		void write_result(int total,TYPE nevyazka);
};
