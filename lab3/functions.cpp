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
	for(int i = 0 ; i < N ;i++) temp[i] = 0;
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

void matrix::calc_start_values(){
	for(int i = 0; i < N ;i++){
		x[i] = 0;
		r[i] = pr[i];
		z[i] = r[i];
	}
	mul_matrix_vector(r);
	for(int i = 0; i < N ;i++) p[i] = temp[i];
};

TYPE matrix::calc_otn_nevazka(){
	TYPE q1,q2;
	q1 = norma(r,r,N);
	q2 = norma(pr,pr,N);
	q1 /= q2; 
	return q1;
};

void matrix::LOS(){
	alfa = square_norma(p,r,N)/square_norma(p,p,N);
	for(int i = 0; i < N ; i++){
		x[i] += alfa*z[i];
		r[i] -= alfa*p[i];
	}
	mul_matrix_vector(r);
	betta = -square_norma(p,temp,N)/square_norma(p,p,N);
	for(int i = 0; i < N ; i++){
		z[i] = r[i] + betta*z[i];
		p[i] = temp[i] + betta*p[i];
	}
};
int matrix::get_maxiter(){
	return maxiter;
};
TYPE matrix::get_eps(){
	return eps;
};

void matrix::write_result(int total,TYPE nevyazka){
	FILE * fp = fopen("result.txt","w");
	fprintf(fp,"%d\n",total);
	for(int i = 0; i < N ; i++) fprintf(fp,FRM_STR_OUT,x[i]);
	for(int i = 1; i <= N ; i++) fprintf(fp,FRM_STR_EPS,x[i-1] - i);
	fprintf(fp,FRM_STR_EPS,nevyazka);
	fclose(fp);
};