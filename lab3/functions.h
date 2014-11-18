#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;
typedef double TYPE;

class vector{
	protected:
		int N;
		TYPE *elem;
	public:
		TYPE square_norm();
		TYPE square_norm_other(vector *v);
		TYPE sum_vector(vector *v);
		TYPE mul_scalar(TYPE *s);
};

class matrix{
	protected:
		int N;
		int maxiter;
		TYPE eps;
		TYPE *ig;
		TYPE *jg;
		TYPE *ggl;
		TYPE *ggu;
		TYPE *di;
		TYPE *pr;
	public:
		void read_kuslau();
		void read_matrix_data(TYPE *mas);
		TYPE mul_matrix_vector(vector *v);
		void dec_LU_sq();
		void allocation_memory();
		void clear_memory();
};
