#include "functions.h"

matrix A;

void solve(int solve_type){
	int maxiter = A.get_maxiter();
	TYPE eps = A.get_eps();
	int k = 0;
	if(solve_type == 3)  A.dec_LU_sq();
	A.calc_start_values(solve_type);
	if(solve_type == 1) for(; k < maxiter &&  A.calc_otn_nevazka() > eps;k++) A.LOS_clean();//1 - clear
	if(solve_type == 2) for(; k < maxiter && A.calc_otn_nevazka() > eps; k++)A.LOS_diag();//2 - diag
	if(solve_type == 3) for(; k < maxiter && A.calc_otn_nevazka() > eps; k++)A.LOS_LU_sq();//3- lu(sq)
	TYPE mm =  A.calc_otn_nevazka();
	A.write_result(k-1,mm);
};

void main(){
	setlocale(LC_CTYPE,"rus");
	A.read_kuslau();
	A.allocation_memory();
	A.read_matrix_data();
	solve(3);
	system("pause");
};