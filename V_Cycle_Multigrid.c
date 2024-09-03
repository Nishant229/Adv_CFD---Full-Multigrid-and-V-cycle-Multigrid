#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double V_cycle_multigrid(int level_start, int level_max, int m, int k, int neu_1, int neu_2, double w, double C, double sigma, double f[level_max+1][m+1], double v[level_max+1][m+1], double r[level_max+1][m+1], double e[level_max+1][m+1], int relaxation_method );
double Weighted_GS_Multigrid_relaxation( int level_max, int m, int l, int neu, double sigma, int n, double w, double v[level_max+1][m+1], double f[level_max+1][m+1] );
double Weighted_Jacobi_Multigrid_relaxation( int level_max, int m, int l, int neu, double sigma, int n, double w, double v[level_max+1][m+1], double f[level_max+1][m+1] );	
double restriction_methods( int level_max, int m, int l, int n, double r[level_max+1][m+1], double f[level_max+1][m+1]);
double prolongation_methods(int level_max, int m, int l, int n, double e[level_max+1][m+1]);
double print_function(int level_max, int m, int l, int n, double v[level_max+1][m+1]);

int Weighted_Gauss_Seidel ( int m, double sigma, double w, double v_GS[m+1], double f[m+1], double epsilon, double residual_GS[1000] );
int Weighted_Jacobi ( int m, double sigma, double w, double v_WJ[m+1], double f[m+1], double epsilon, double residual_WJ[1000]  );

int main() {
	int N; // Grid size
	double w; // weighing factor
	double C, sigma;
	int neu, k; // Number of Jacobi/Gauss-Seidel relaxations

	FILE *fp1;
	fp1=fopen("input.dat","r");
	if (fp1 == NULL) {
		printf("------- File is empty of filename is incorrect !!!! -------\n ");
		exit(1);
	}
	// Taking input from file
	fscanf(fp1,"***** MultiGrid parameters *****\n\n");
	fscanf(fp1,"// Grid Size\n");
	fscanf(fp1,"\tn = %d\n\n",&N);
	fscanf(fp1,"// Jacobi Weight\n");
	fscanf(fp1,"\t w = %lf\n\n",&w);
	fscanf(fp1,"// Other constants\n");
	fscanf(fp1,"\t neu = %d\n",&neu);
	fscanf(fp1,"\t sigma = %lf\n",&sigma);
	fscanf(fp1,"\t k = %d\n",&k);
	fscanf(fp1,"\t C = %lf ( C = (pi*k)^2 + sigma )\n",&C);	
	printf("n = %d , w = %lf , neu = %d , sigma = %lf , k = %d , C = %lf \n",N,w,neu,sigma,k,C);
	
	// V-Cyle Algorithm
	int m=N;
	double u[N+1], v_initial[N+1], f_initial[N+1];
	// Homogeneous Boundary conditions
	u[0]=0; v_initial[0]=0;
	u[N]=0; v_initial[N]=0;
	f_initial[0] = 0, f_initial[N] = 0;
	for(int j=1;j<N;j++) {
		v_initial[j] = 0; // Initial guess
		u[j] = sin(k*j*M_PI / N); // Exact Solution
		f_initial[j] = C*sin(k*M_PI*j / N);
	}
	
	// Maximum number of levels
	int num_of_levels; // level of multigrid
	num_of_levels = log(N) / log(2);
	printf("Total number of levels: %d \n",num_of_levels);
	int level_max = num_of_levels;
	
	double f[level_max+1][N+1];
	double v[level_max+1][N+1];
	double r[level_max+1][N+1];
	double e[level_max+1][N+1];
	
	double v_old[N+1], v_new[N+1];
		
	// Initialization
	for(int i=0;i<=level_max;i++) {
		for(int j=0;j<=N;j++){
			f[i][j] = 0;
			v[i][j] = 0;
			r[i][j] = 0;
			e[i][j] = 0;
		}
	}
	
	// Initializing level 1
	for(int i=0;i<=N;i++) {
		f[1][i] = f_initial[i];
		v[1][i] = v_initial[i];
		v_old[i] = v_initial[i];
	}	
	
	int neu_1 = neu, neu_2 = neu;
	printf("neu_1 = %d  ,  neu_2 = %d \n",neu_1,neu_2);
	
	double epsilon = pow(10,-6);
	double error=1;
	double residual = 1;
	int V_cycle_iteration_count=0;
	double residual_V_Cycle[100];
	
	int relaxation_method;
	printf("Which Method You need to use for relaxation->    1 = Weighted Jacobi   OR  2 = Gauss Seidel Jacobi: ");
	scanf("%d",&relaxation_method);
	
	// Calling v_cycle
	while  ( residual > epsilon) { 
		error = 0;
		residual = 0;
		V_cycle_multigrid(1, level_max, m, k, neu_1, neu_2, w, C, sigma, f, v, r, e, relaxation_method );
		
		V_cycle_iteration_count++;
	
		// Error Calculation
		double h = 1.0/N;
		for(int i=1;i<m;i++) {
			v_new[i] = v[1][i];
			error = error + pow(v_new[i]-v_old[i],2); // Calculating error
			residual = residual + pow ( f[1][i] - (1/pow(h,2))*( - v[1][i-1] + (2.0+sigma*pow(h,2))*v[1][i] - v[1][i+1]  ) , 2);
			v_old[i] = v_new[i]; // updating the old V cycle iteration value
		}
		//error = pow(error/(m-1),0.5); // root Mean Square
		//residual = pow(residual/(m-1),0.5); // root Mean Square
		error = pow(error,0.5);
		residual = pow(residual,0.5);
		residual_V_Cycle[V_cycle_iteration_count] = residual;
		//printf("Iteration No. = %d \t Error = %0.10lf \n",V_cycle_iteration_count,error);
		printf("Iteration No. = %d \t Residual = %0.10lf \n",V_cycle_iteration_count,residual);
		
		// Testing
		//printf("After %d --> \n",V_cycle_iteration_count);
		//print_function(level_max, m, 1, N, v);
		
		
	}
		
	// Displaying result
	printf("The number of V Cycle Iterations: %d \n",V_cycle_iteration_count);
	
	
	// Displaying all 512  (N) grid points actual and numerical data
	/*printf("Actual Value \t Numerical Value \t Error \n\n");
	for(int i=0;i<=m;i++) {
		printf("i=%d \t u=%lf \t v=%lf \t error=%0.10lf \n",i,u[i], v_new[i], fabs(u[i] - v_new[i]));
	}*/
	
// Normal Iterative Methods  --->>> For Comparison of Normal Iterative and Multigrid Methods 	
	
	/*printf("-----------------------------------------------------------------------------------------------------------------------\n");
	//printf("-----------------------------------------------------------------------------------------------------------------------\n");
	//printf("-----------------------------------------------------------------------------------------------------------------------\n");
	
// Normal Iterative Methods  
	
	printf("\n------------Calculating Iterations and Residuals for Weighted Jacobi and Gauss Siedel Iterations----------\n");
	
	double v_GS[N+1],v_WJ[N+1],f_1[N+1],residual_GS[1000],residual_WJ[1000];
	for(int i=0;i<=N;i++) {
		f_1[i] = f_initial[i];
		v_WJ[i] = v_initial[i];
		v_GS[i] = v_initial[i];
		// printf("f = %lf , v_WJ = %lf , v_GS = %lf \n",f_1[i],v_WJ[i],v_GS[i]); // Testing
	}
	
	
	// Normal Gauss_Seidel (GS) Method (Weight = 1)
	int GS_iterations = Weighted_Gauss_Seidel ( m, sigma, 1.0, v_GS, f_1,  epsilon, residual_GS );
	
	// Weighted Jacobi (WJ) Method (Weight = 2/3)
	int WJ_iterations = Weighted_Jacobi ( m, sigma, w, v_WJ, f_1,  epsilon, residual_WJ );
	
	FILE *fp2;
	fp2 = fopen("Results.dat","w");
	fprintf(fp2,"Iterations \t V_Cycle \t Gauss_Seidel \t Weighted_Jacobi \n\n");
	int iterations = GS_iterations > WJ_iterations ? GS_iterations : WJ_iterations;
	int q=1,p=1,s=1;
	for(int i=1;i<=iterations;i++) {
		if(i<=10 || i%1000==0) {
			fprintf(fp2,"%d \t",i);
			if(i<=V_cycle_iteration_count) {
				fprintf(fp2,"%lf \t",residual_V_Cycle[q]);
				q++;
			}
			else {
				fprintf(fp2," 0.0 \t\t");
			}
			if(i<=GS_iterations) {
				fprintf(fp2,"%lf \t",residual_GS[p]);
				p++;
			}
			else {
				fprintf(fp2," 0.0 \t\t");
			}
			if(i<=WJ_iterations) {
				fprintf(fp2,"%lf \t",residual_WJ[s]);
				s++;
			}
			else {
				fprintf(fp2," 0.0 \t\t");
			}
			fprintf(fp2,"\n");
		}
	}
	fclose(fp2);*/
	
	fclose(fp1);
	
	return 1;
}



double V_cycle_multigrid(int level_start, int level_max, int m, int k, int neu_1, int neu_2, double w, double C, double sigma, double f[level_max+1][m+1], double v[level_max+1][m+1], double r[level_max+1][m+1], double e[level_max+1][m+1], int relaxation_method  ) {
		
		// Going down the V Cycle
		for(int level=level_start;level<level_max;level++) {
			int n=pow(2,level_max-(level-1));
			double h=1.0/n;
			if (level!=level_start) {
				for(int i=0;i<=m;i++) {
					v[level][i] = 0;
				}
			}
			// Relaxation Method
			if (relaxation_method == 2) {
				Weighted_GS_Multigrid_relaxation(level_max,m,level, neu_1, sigma, n,  1.0, v, f );
			}
			if (relaxation_method == 1) {
				Weighted_Jacobi_Multigrid_relaxation( level_max, m, level, neu_1, sigma, n, w, v, f );
			}
			r[level][0] = 0;
			r[level][n] = 0;
			for(int i=1;i<n;i++) {
				r[level][i] = f[level][i] - (1/pow(h,2))*( -v[level][i+1] + (2 + sigma*pow(h,2))*v[level][i] - v[level][i-1]);
			}	
			restriction_methods( level_max,m,level, n, r, f);
		}
		
		// Coursest level
		int n=pow(2,1);
		double h=1.0/n;
		e[level_max][0] = f[level_max][0];
		e[level_max][1] = pow(h,2)*f[level_max][1] / (2+sigma*pow(h,2));
		e[level_max][2] = f[level_max][2];
		// printf("Coursest Grid Solution = %lf \n",e[level_max][1]); //Testing
		
		 
		// Going up the V Cycle
		for(int level=level_max;level>=level_start; level-- ) {
			n=pow(2,level_max-(level-1));
			double h=1.0/n;
			if (level!=level_max){
				for(int i=0;i<=n;i++) {
					v[level][i] = v[level][i] + e[level][i];
				}
				if (relaxation_method == 2) {
					Weighted_GS_Multigrid_relaxation(level_max,m,level, neu_2, sigma, n,  1.0, v, f );
				}
				if (relaxation_method == 1) {
					Weighted_Jacobi_Multigrid_relaxation( level_max, m, level, neu_2, sigma, n, w, v, f );
				}
				for(int i=0;i<n;i++) {
					e[level][i] = v[level][i];
				}
				
			}
			if(level!=level_start) {
				prolongation_methods(level_max,m, level, n, e);
			}
		}
	
	return 1;
}

double Weighted_GS_Multigrid_relaxation( int level_max, int m, int l, int neu, double sigma, int n, double w, double v[level_max+1][m+1], double f[level_max+1][m+1] ) {
	int iteration_count = 0;
	double h = 1.0/n;
	double temp;
	while (iteration_count<neu) {
		for (int i=1;i<n;i++) {
			temp = (1/(2 + sigma*pow(h,2))) * ( v[l][i-1] + v[l][i+1] + pow(h,2)*f[l][i] );
			//printf("temp=%0.12lf \t",temp); // Testing
			v[l][i] = (1-w)*v[l][i] + w*temp;
			//printf("%0.12lf \t",v[l][i]); // Testing
		}
		iteration_count++;
	}
	// Testing
	// print_function(level_max, m, l, n, v);
	return 1;
}

double Weighted_Jacobi_Multigrid_relaxation( int level_max, int m, int l, int neu, double sigma, int n, double w, double v[level_max+1][m+1], double f[level_max+1][m+1] ) {
	int p = m / pow(2,l-1);
	double vel[n+1];
	for(int i=0;i<=n;i++) {
		vel[i] = v[l][i];
	} 
	int iteration_count = 0;
	double h = 1.0/n;
	double temp;
	while (iteration_count<neu) {
		for (int i=1;i<n;i++) {
			temp = (1/(2 + sigma*pow(h,2))) * ( vel[i-1] + vel[i+1] + pow(h,2)*f[l][i] );
			//printf("temp=%0.12lf \t",temp); // Testing
			v[l][i] = (1-w)*vel[i] + w*temp;
			//v[l][i] = temp;
			//printf("%0.12lf \t",v[l][i]); // Testing
		}
		iteration_count++;
		for (int i=1;i<n;i++) {
			vel[i] = v[l][i];	
		}
	}
	// Testing
	// print_function(level_max, m, l, n, v);
	return 1;
}

double restriction_methods( int level_max, int m, int l, int n, double r[level_max+1][m+1], double f[level_max+1][m+1]) {
	int p = (int)(n/2);
	r[l+1][0] = 0;
	r[l+1][p] = 0;
	for(int j=1;j<p;j++) {
		r[l+1][j] = 0.25*( r[l][2*j-1] + 2*r[l][2*j] + r[l][2*j+1] );
	}
	for(int j=0;j<=p;j++) {
		f[l+1][j] = r[l+1][j];
	}
	return 1;
} 

double prolongation_methods(int level_max, int m, int l, int n, double e[level_max+1][m+1]) {
	e[l-1][0] = 0;
	e[l-1][2*n] = 0;
	for(int j=0;j<n;j++) {
		e[l-1][2*j] = e[l][j];
		e[l-1][2*j+1] = 0.5*(e[l][j] +e[l][j+1]);
	}
	return 1;
}

int Weighted_Gauss_Seidel ( int m, double sigma, double w, double v_GS[m+1], double f[m+1], double epsilon, double residual_GS[1000]  ) {
	int iteration_count = 0;
	double h = 1.0/m;
	double temp;
	double residual=1.0;
	int p=0;
	while (residual>epsilon) {
		residual = 0.0;
		for (int i=1;i<m;i++) {
			temp = (1/(2 + sigma*pow(h,2))) * ( v_GS[i-1] + v_GS[i+1] + pow(h,2)*f[i] );
			temp = (1-w)*v_GS[i] + w*temp;
			v_GS[i] = temp;
		}
		for(int i=1;i<m;i++) {
			residual = residual + pow( f[i] - (1.0/pow(h,2))*( - v_GS[i-1] + (2.0+sigma*pow(h,2))*v_GS[i] - v_GS[i+1] ) ,2);
			//printf("%lf \n", residual);
		}
		residual = pow(residual,0.5); // L_2 Norm Residual
		iteration_count++;
		//printf("\nIteration = %d ==> residual = %lf \n",iteration_count,residual);
		if(iteration_count<=10 || iteration_count%1000==0) {
			p=p+1;
			residual_GS[p] = residual;
		}
	}
	printf("Weighted Gauss Seidel Iterations = %d \n",iteration_count);
	return iteration_count;
}

int Weighted_Jacobi ( int m, double sigma, double w, double v_WJ[m+1], double f[m+1], double epsilon, double residual_WJ[1000]  ) {
	int iteration_count = 0;
	double h = 1.0/m;
	double temp;
	double residual=1.0;
	double v_WJ_old[m+1];
	for(int i=0;i<=m;i++) {
		v_WJ_old[i] = v_WJ[i];
	}
	int p=0;
	while (residual>epsilon) {
		residual = 0.0;
		for (int i=1;i<m;i++) {
			temp = (1/(2 + sigma*pow(h,2))) * ( v_WJ_old[i-1] + v_WJ_old[i+1] + pow(h,2)*f[i] );
			temp = (1-w)*v_WJ_old[i] + w*temp;
			v_WJ[i] = temp;	
		}
		for(int i=1;i<m;i++) {
			residual = residual + pow( f[i] - (1.0/pow(h,2))*( - v_WJ[i-1] + (2.0+sigma*pow(h,2))*v_WJ[i] - v_WJ[i+1] ) ,2);
			//printf("%lf \n", residual);
		}
		residual = pow(residual,0.5); // L_2 Norm Residual
		iteration_count++;
		if(iteration_count<=10 || iteration_count%1000==0) {
			p=p+1;
			residual_WJ[p] = residual;
		}	
		//printf("\nIteration = %d ==> residual = %lf \n",iteration_count,residual);
		// Updating Calcluated values into old jacobi Matrix array
		for(int i=0;i<=m;i++) {
			v_WJ_old[i] = v_WJ[i];
		}
		
	}
	printf("Weighted Jacobi Iterations = %d \n",iteration_count);
	return iteration_count;
}

// For Testing Purpose
double print_function(int level_max, int m, int l, int n, double v[level_max+1][m+1]) {
	printf("-------------------------------------------------------\n");
	printf("The Matrix at level %d is as follows: \n",l);
	for(int i=0;i<=n;i++) {
		printf("%0.12lf \t",v[l][i]);
	}
	printf("\n\n");
}
