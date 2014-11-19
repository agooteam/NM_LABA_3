#include "functions.h"

matrix A;

void solve(int solve_type){
	//1 - clear
	//2 - diag
	//3- lu(sq)
	int maxiter = A.get_maxiter();
	TYPE eps = A.get_eps();
	int k = 0;
	for(; k < maxiter &&  A.calc_otn_nevazka() > eps;k++) A.LOS();
	TYPE mm =  A.calc_otn_nevazka();
	A.write_result(k-1,mm);
};

void main(){
	setlocale(LC_CTYPE,"rus");
	A.read_kuslau();
	A.allocation_memory();
	A.read_matrix_data();
	A.calc_start_values();
	solve(1);
	system("pause");
};