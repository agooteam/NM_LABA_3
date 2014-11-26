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
	temp = new TYPE[N];
	temp2 = new TYPE[N];
	cggl = new TYPE[count_elem];
	cggu = new TYPE[count_elem];
	cdi = new TYPE[N];
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
	delete [] temp2;
	delete [] cggu;
	delete [] cggl;
	delete [] cdi;
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


TYPE matrix::calc_otn_nevazka(){
	TYPE q1,q2;
	q1 = norma(r,r,N);
	q2 = norma(pr,pr,N);
	q1 /= q2; 
	return q1;
};

void matrix::LOS_clean(){
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

void matrix::LOS_diag(){
	alfa = square_norma(p,r,N)/square_norma(p,p,N);
	for(int i = 0; i < N ; i++){
		x[i] += alfa*z[i];
		r[i] -= alfa*p[i];
	}
	for(int i = 0; i < N;i++) temp2[i] = r[i]/sqrt(di[i]);
	mul_matrix_vector(temp2);
	for(int i = 0; i < N;i++) temp[i] = temp[i]/sqrt(di[i]);
	betta = -square_norma(p,temp,N)/square_norma(p,p,N);
	for(int i = 0; i < N ; i++){
		z[i] = r[i]/sqrt(di[i]) + betta*z[i];
		p[i] = temp[i] + betta*p[i];
	}
};

void matrix::LOS_LU_sq(){
	alfa = square_norma(p,r,N)/square_norma(p,p,N);
	for(int i = 0; i < N ; i++){
		x[i] += alfa*z[i];
		r[i] -= alfa*p[i];
	}
	reverse(r,temp2);
	mul_matrix_vector(temp2);
	direct(temp,temp2);
	betta = -square_norma(p,temp2,N)/square_norma(p,p,N);
	reverse(r,temp);
	for(int i = 0; i < N ; i++){
		z[i] = temp[i] + betta*z[i];//
		p[i] = temp2[i] + betta*p[i];
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

void matrix::calc_start_values(int solve){
	if(solve == 1){
		for(int i = 0; i < N ;i++){
			x[i] = 0;
			r[i] = pr[i];
			z[i] = r[i];
		}
		mul_matrix_vector(r);
		for(int i = 0; i < N ;i++) p[i] = temp[i];
	}
	else if(solve == 2){
		for(int i = 0; i < N ;i++){
			x[i] = 0;
			r[i] = pr[i]/sqrt(di[i]);
			z[i] = r[i]/sqrt(di[i]);
		}
		mul_matrix_vector(z);
		for(int i = 0; i < N ;i++) p[i] = temp[i]/sqrt(di[i]);
	}
	else if(solve == 3){
		for(int i = 0; i < N ;i++) x[i] = 0;
		direct(pr,r);
		reverse(r,z);
		mul_matrix_vector(z);
		direct(temp,p);
	}
};

void matrix::dec_LU_sq(){
	int i,j,kol,kl=0,ku=0;
	for(i = 0 ; i < N; i++){
		kol = ig[i+1] - ig[i];
		for(j = 0; j < kol; j++){
			cggl[kl] = (ggl[kl] - dec_calc_elem(i,j,kl)) / cdi[jg[kl]];
			kl++;
		}
		for(j = 0;j < kol; j++){
			cggu[ku] = (ggu[ku] - dec_calc_elem(i,j,ku)) / cdi[jg[ku]];
			ku++;
		}
		cdi[i] = sqrt(di[i] - dec_calc_diag(j,kl));
	}

};

TYPE matrix::dec_calc_elem(int i, int j, int current_elem){
	TYPE s = 0;
	int k,J=jg[current_elem],p;
	for( k = j ; k > 0; k--){
		for(p = ig[J];p < ig[J+1]; p++) if(jg[p] == jg[current_elem - k]) s+= cggl[current_elem - k]*cggu[p];
	}
	return s;
};

TYPE matrix::dec_calc_diag(int j, int current_elem){
	TYPE s=0;
	int k;
	for( k = current_elem-j; k < current_elem; k++) s += cggl[k]*cggu[k];
	return s;
};

void matrix::direct(TYPE *in,TYPE *out){
	int k = 0,column,count;
	for(int i = 0; i < N; i++) out[i] = in[i];
	for(int i = 0 ; i < N; i++){	
		TYPE s=0;
		count = ig[i+1] - ig[i];
		for(int j = 0; j < count; j++, k++){ 	
			column = jg[k];
			s+= cggl[k]*out[column];
		}
		out[i] = (out[i]-s) / cdi[i];
	}
};

void matrix::reverse(TYPE *in,TYPE *out){
	int k = ig[N] - ig[0]-1,count,column;
	for(int i = 0; i < N; i++) out[i] = in[i];
	for(int i = N-1; i>=0; i--){	
		out[i] = out[i]/cdi[i];
		count = ig[i+1] - ig[i];
		for(int j = 0; j < count; j++, k--){	
			column = jg[k];
			out[column] -= cggu[k]*out[i];
		}	
	}

};