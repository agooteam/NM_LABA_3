#include "functions.h"

TYPE square_norma(TYPE *v, TYPE *m , int size){
	TYPE summ = 0;
	for(int i = 0; i < size ; i++ ) summ += v[i]*m[i];
	return summ;
};

TYPE norma(TYPE *v, TYPE *m , int size){
	TYPE summ = 0;
	for(int i = 0; i < size ; i++ ) summ += v[i]*m[i];
	summ = sqrt(summ);
	return summ;
};

TYPE summ_vector(TYPE *v, TYPE *m , TYPE *res, int size){
	TYPE summ = 0;
	for(int i = 0; i < size ; i++ ) res[i] = v[i]+m[i];
	return *res;
};

TYPE mul_scalar_vector(TYPE scalar, TYPE *v , int size){
	TYPE summ = 0;
	for(int i = 0; i < size ; i++ ) v[i] = scalar * v[i];
	return *v;
};

void matrix::read_kuslau(){
	 ifstream input("kuslau.txt");
	 input >> N;
	 input >> maxiter;
	 input >> eps;
	 input.close();
};

void matrix::allocation_memory(){
	ifstream input("ig.txt");
	ig = new int[N+1];
	for(int i = 0; i <= N; i++){ input >> ig[i] ; ig[i]--;}
	int count_elem = ig[N] - ig[0];
	jg = new int[count_elem];
	ggl = new TYPE[count_elem];
	ggu = new TYPE[count_elem];
	pr = new TYPE[N];
	di = new TYPE[N];
	x = new TYPE[N];
	r = new TYPE[N];
	z = new TYPE[N];
	p = new TYPE[N];
	temp = new TYPE[N]();
};

void matrix::clear_memory(){
	delete [] ig;
	delete [] jg;
	delete [] ggu;
	delete [] ggl;
	delete [] pr;
	delete [] di;
	delete [] x;
	delete [] r;
	delete [] z;
	delete [] p;
	delete [] temp;
};

void matrix::read_matrix_data(){
	ifstream input;
	int count_elem = ig[N] - ig[0];
	input.open("di.txt");
	for(int i = 0; i < N ;i++) input >> di[i];
	input.close();
	input.open("pr.txt");
	for(int i = 0; i < N ;i++) input >> pr[i];
	input.close();
	input.open("jg.txt");
	for(int i = 0; i < count_elem; i++){ input >> jg[i] ; jg[i]--; }
	input.close();
	input.open("ggl.txt");
	for(int i = 0; i < count_elem; i++) input >> ggl[i];
	input.close();
	input.open("ggu.txt");
	for(int i = 0; i < count_elem; i++) input >> ggu[i];
	input.close();
};

void matrix::mul_matrix_vector(TYPE *v){
	for(int i = 0 ; i < N ; i++){
		temp[i] += di[i]*v[i];
		int count_elem = ig[i+1] - ig[i];
		for(int p = 0 ; p < count_elem; p++){
			int m = ig[i]+p;
			int column = jg[m];
			temp[i] += ggl[m]*v[column];
			temp[column] += ggu[m] * v[i];
		}
	}
};
		