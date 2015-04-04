#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void read_input_data(FILE *input_file,
					 int *rows,
					 int *cols,
					 double *J,
					 double *t_min,
					 double *t_max,
					 int *n_taus,
					 int *n_steps,
					 int *n_measurments,
					 char *bc,
					 int lattice[40][40],
					 double magn_field[40][40],
					 int *corr_spin,
					 double random_numbers[100000][2]
					 ){
	char val;
	int n,rc,i,j;
	
	int spin;double mg;
	
	for(n=1;(rc=fscanf(input_file,"%c",&val))!=EOF;n++){
		
		//PASS HEADERS
		if((n<5) || (n>5 && n<15) || (n>15 && n<33) || (n>33 && n<48)) continue;
		//READ "DATA"
		if(n==5){
			fscanf(input_file,"%d	%d	%lf	%lf	%lf	%d	%c	%d	%d",rows,cols,J,t_min,t_max,n_taus,bc,n_steps,n_measurments);
		}
		//READ "ELEMENTS"
		if(n==15){
			for(i=0;i<(*rows);i++){
				for(j=0;j<(*cols);j++){
					fscanf(input_file,"%d	%lf",&lattice[i][j],&magn_field[i][j]);
				}
			}
		}
		//READ "CORRELATION SPIN"
		if(n==33){
			fscanf(input_file,"%d",corr_spin);
		}
		//READ "RANDOM NUMBERS"
		if(n==48){
			for(i=0;i<(*n_steps);i++){
					fscanf(input_file,"%lf	%lf",&random_numbers[i][0],&random_numbers[i][1]);
				}
		}
		
		
		}
}

//tested
double magn_field_interaction(int rows,
							 int cols ,
							 int lattice[40][40],
							 double magn_field[40][40]){
	int _r,_c;
	double S=0;
	for(_r=0;_r<rows;_r++){
			for(_c=0;_c<cols;_c++){
				S+=lattice[_r][_c]*magn_field[_r][_c];
			}
	}
	return S;
}
//tested
double spin_interaction(char boundry_conditions,
					   int rows,
					   int cols,
					   int lattice[40][40]){
		
		int _r,_c;
		double S=0;
		int steps=0;
		
		for(_r=0;_r<rows;_r++){
			for(_c=0;_c<cols-1;_c++){
				S+=lattice[_r][_c]*lattice[_r][_c+1];
				steps++;
			}
		}
		for(_c=0;_c<cols;_c++){
			for(_r=0;_r<rows-1;_r++){
				S+=lattice[_r][_c]*lattice[_r+1][_c];
				steps++;
			}
		}
		if(boundry_conditions=='c'){
				for(_c=0;_c<cols;_c++){
					S+=lattice[0][_c]*lattice[rows-1][_c];
					steps++;
				}
				for(_r=0;_r<rows;_r++){
					S+=lattice[_r][0]*lattice[_r][cols-1];
					steps++;
				}
		}
		return S;
}		
//tested
double E(char boundry_conditions,
		int rows,
		int cols,
		int lattice[40][40],
		double magn_field[40][40],
		double J){
			
	return -1*J*spin_interaction(boundry_conditions,rows,cols,lattice) - 
			magn_field_interaction(rows,cols,lattice,magn_field);
}
//tested
void flip_spin(int spin_row,int spin_col,int lattice[40][40]){
	lattice[spin_row][spin_col]*=-1;
}
//tested
//converts index in 2 dim array to its position using row/col
//index should be in range [0,rows*cols-1]
//0*cols ............. 1*cols-1
//1*cols ..............2*cols-1
//2*cols ..............3*cols-1
//                .
//                .
//                .
//(rows-1)*cols........rows*cols-1
void cell_index_to_row_col(int index, int rows, int cols ,int* row ,int* col){
	int r,c;
	r = index / cols; 
	c = index - cols*(r);
	
	
	*row=r;
	*col=c;
}
//tested
//picks two random numbers from the rows and cols
//higher random num should give a higher numbers
void choose_spin(double random_num,int rows, int cols ,int* row ,int* col){
	int max_index=rows*cols -1; //i count from zero
	int site_number=(int)(random_num*max_index);
	
	cell_index_to_row_col(site_number,rows,cols,row,col);
}

//tested
//copies source array into target arr
void copy_arr(int source[40][40],int target[40][40]){
	int _r,_c;
	for(_r=0;_r<40;_r++){
		for(_c=0;_c<40;_c++){
			target[_r][_c]=source[_r][_c];
		}
	}
}

//tested
//given an array with numbers and the size of it and some number "exponent"
//calculates SUM(array_item^exponent)/array_size
double thermal_average(int how_many,double *array,double exponent){
	double sum=0;
	int i;
	for(i=0;i<how_many;i++){
		sum+=pow(*array,exponent);
		array++;
	}
	return sum / (double)how_many;
}
//tested
//returns the POPULATION STDEV of the given array of size "how_many"
double thermal_stdev(int how_many,double *array){
	double th_avg_en_sq=thermal_average(how_many,array,2);
	double th_avg_sq=pow(thermal_average(how_many,array,1),2);
	double variance = th_avg_en_sq - th_avg_sq;
	double stdev=pow(variance ,0.5);
	if(stdev!=stdev){ //stdev is NaN,probably problem of double size precision
		return 0;
	}
	return stdev;
}

//tested
double next_tau(double tau,double t_min,double t_max,int n_taus){
	return tau + (t_max-t_min)/(double)(n_taus-1);
}
//tested(by eye)
void reset_system(int original_lattice[40][40],int active_lattice[40][40]){
	copy_arr(original_lattice,active_lattice);
	
}
//tested (by eye)
double M(int lattice[40][40],int rows,int cols){
	int _r,_c;
	double m=0;
	for(_r=0;_r<rows;_r++){
		for(_c=0;_c<cols;_c++){
			m+=lattice[_r][_c];
		}
	}
	return m/(double)(rows*cols);
}
//tested
void print_lattice(FILE *fout,int lattice[40][40],int rows,int cols){
	int _r,_c;
	for(_r=0;_r<rows;_r++){
		for(_c=0;_c<cols;_c++){
			fprintf(fout,"	%d",lattice[_r][_c]);
		}
		fprintf(fout, "\n");
	}
}
//tested
void print_magn_field(FILE *fout, double h[40][40],int rows,int cols){
	int _r,_c;
	for(_r=0;_r<rows;_r++){
		for(_c=0;_c<cols;_c++){
			fprintf(fout, "	%lf",h[_r][_c]);
		}
		fprintf(fout, "\n");
	}
}
double calculate_transition_probability(int spin_row,int spin_col,int rows,int cols,int lattice[40][40],double magn_field[40][40], double J,double tau,char boundry_conditions){
	double E_before,E_after;
	
	E_before=E(boundry_conditions,rows,cols,lattice,magn_field,J);
	flip_spin(spin_row,spin_col,lattice);
	E_after=E(boundry_conditions,rows,cols,lattice,magn_field,J);
	flip_spin(spin_row,spin_col,lattice);
	
	if(E_after > E_before){
		return  pow(M_E,-1*(E_after-E_before)/tau);
	}
	else{
		return 1;
	
	}
}
int need_to_print_measurment(int current_step,int total,int total_meas){
	return current_step % (total/(total_meas-1)) == 0;
}
int step_index(int current_step,int total,int total_meas){
	return current_step / (total/total_meas);
}
//make a step in the simulation (try to flip a spin)
void do_step(double rnd1,double rnd2, int rows, int cols,int lattice[40][40],double magn_field[40][40],double J,double tau,char bc ){
			int _r,_c;
			double probabilty_to_flip;

			choose_spin(rnd1,rows,cols,&_r,&_c);
			probabilty_to_flip = calculate_transition_probability(_r,_c,rows,cols,lattice,magn_field,J,tau,bc);
			if(rnd2<=probabilty_to_flip){
				flip_spin(_r,_c,lattice);
			}
}
//runs the bussiness logic
int start(FILE *input,FILE *meas, FILE *direc, FILE *corr){
	//given variables
	int rows,cols,n_taus,n_steps,n_measurments,corr_spin;
	double J,t_min,t_max;
	char boundry_conditions;
	int original_lattice[40][40];
	double magn_field[40][40];
	double random_numbers[100000][2];
	read_input_data(input,&rows,&cols,&J,&t_min,&t_max,&n_taus,&n_steps,&n_measurments,&boundry_conditions,original_lattice,magn_field,&corr_spin,random_numbers);
	
	//dynamic variables
	int step_counter; //in which step are we?
	int lattice[40][40];copy_arr(original_lattice,lattice); //duplicate the original lattice
	double tau;
	double energies[n_steps],magnetizations[n_steps];
	
	//meas.magn direc.spin
	fprintf(meas,"%d	%lf	%d	%d	%d	%d\n",n_taus,J,1,rows,cols,n_measurments);
	fprintf(direc,"%d	%lf	%d	%d	%d	%d\n",n_taus,J,1,rows,cols,n_measurments);
	fprintf(meas,"H=\n");
	fprintf(direc,"H=\n");
	print_magn_field(meas,magn_field,rows,cols);
	print_magn_field(direc,magn_field,rows,cols);
	
	for(tau=t_min;tau <= t_max;tau = next_tau(tau,t_min,t_max,n_taus)){
		reset_system(original_lattice,lattice);
		//meas.magn  direc.spin
		fprintf(meas, "%lf\n",tau);
		fprintf(direc, "%lf\n",tau);
		
		for(step_counter=0;step_counter<n_steps;step_counter++){
			
			do_step(random_numbers[step_counter][0],random_numbers[step_counter][1],rows,cols,lattice,magn_field,J,tau,boundry_conditions);
		
			energies[step_counter] = E(boundry_conditions,rows,cols,lattice,magn_field,J);
			magnetizations[step_counter] = M(lattice,rows,cols);
			
			if(need_to_print_measurment(step_counter,n_steps,n_measurments)){
				//meas.magn
				fprintf(meas,"%d	%lf	%lf	%lf	%lf\n",
									  step_index(step_counter,n_steps,n_measurments),
									  energies[step_counter],
									  thermal_stdev(step_counter+1,energies),
									  magnetizations[step_counter],
									  thermal_stdev(step_counter+1,magnetizations));
				//direc.spin
				fprintf(direc,"%d\n", step_index(step_counter,n_steps,n_measurments));
				print_lattice(direc,lattice,rows,cols);
			}
		}
		//meas.magn
		fprintf(meas, "result:\n");
		fprintf(meas, "%lf	%lf	%lf	%lf \n",
									thermal_average(n_steps,energies,1),
									thermal_stdev(n_steps,energies),
								    thermal_average(n_steps,magnetizations,1),
									thermal_stdev(n_steps,magnetizations));
	}
}
int main(int argc, char **argv)
{
	//prepare I/O files
	FILE *input_file = stdin; //fopen("/home/yardenst/workspace/ferro/input/4/input.dat","r");
	FILE *meas_output_file = fopen("./meas.magn","w"); //stdout;
	FILE *direc_output_file = fopen("./direc.spin","w"); //stdout;
	FILE *corr_output_file = fopen("./corr.corr","w"); //stdout;
	
	start(input_file, meas_output_file, direc_output_file, corr_output_file);
	
	
}
