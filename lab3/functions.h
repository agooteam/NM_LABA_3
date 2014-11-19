#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;
typedef double TYPE;


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
		TYPE *pr;
		TYPE *x;
		TYPE *r;
		TYPE *z;
		TYPE *p;
	public:
		void read_kuslau();
		void read_matrix_data();
		TYPE mul_matrix_vector(TYPE *v);
		void dec_LU_sq();
		void allocation_memory();
		void clear_memory();
};
