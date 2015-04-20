#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//logs a message to the given file in a log info style
void logg(char msg[255],FILE *log_file){
	fprintf(log_file,"INFO: %s",msg);
	fprintf(log_file,"\n");	
}

//tested
//returns 1 if number is not an integer
//returns 0 if number is an integer
int validate_int(double number){
	return !(number==(int)number);
}
//returns 1 if number is not positive
//returns 0 if number is positive
int validate_positive(double number){
	return !(number>0);
}
//returns 1 if number is not a valid spin value (-1 or 1)
//returns 0 if number is a valid spin value (-1 or 1)
int validate_spin_value(double number){
	return !((number==1) || (number==-1));
}
//returns 1 if bc_value is not the char o or c, 0 otherwise
int validate_bc_value(char bc_value){
	return !(bc_value=='o' || bc_value=='c');
}

int validate_spin_range(int rows,int cols, int index){
		return !(index >=0 && index <=(rows*cols-1));
}

int validate_rnd_number(double number){
	return !(number>=0 && number<=1);
}

void log_section_error(FILE *log,int validator_count,char section_name[100]){
		if(validator_count>0){
				fprintf(log,"WARNING: %d errors on %s section\n", validator_count, section_name);
			}
}

//returns number of errors found on input FILE
int read_input_data(FILE *input_file, FILE *errors_file,
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
					 double random_numbers[300000][2]
					 ){
	char val;
	int n,rc,i,j;
	int validator_count=0,errors_before=0;
	
	//double values for int values
	double rows_d,cols_d,n_taus_d,n_steps_d,n_measurments_d,lattice_item_d,corr_spin_d;
	
	for(n=1;(rc=fscanf(input_file,"%c",&val))!=EOF;n++){
		
		//PASS HEADERS NAMES
		if((n<5) || (n>5 && n<15) || (n>15 && n<33) || (n>33 && n<48)) continue;
		//READ "DATA"
		if(n==5){
			fscanf(input_file,"%lf	%lf	%lf	%lf	%lf	%lf	%c	%lf	%lf",&rows_d,&cols_d,J,t_min,t_max,&n_taus_d,bc,&n_steps_d,&n_measurments_d);
			
			*rows=(int)rows_d;
			*cols=(int)cols_d;
			*n_taus=(int)n_taus_d;
			*n_steps=(int)n_steps_d;
			*n_measurments=(int)n_measurments_d;
			
			validator_count+=validate_int(rows_d) + validate_int(cols_d) + validate_int(n_taus_d) + 
							 validate_int(n_steps_d) + validate_int(n_measurments_d) +
							 validate_positive(*rows) + validate_positive(*cols) + 
							 validate_positive(*n_taus) + validate_positive(*n_steps) +
							 validate_positive(*n_measurments) + validate_positive(*t_min) + 
							 validate_positive(*t_max) + validate_bc_value(*bc) + (*n_measurments>*n_steps); 
							 
			log_section_error(errors_file,validator_count,"DATA");
			
			
		}
		//READ "ELEMENTS"
		if(n==15){
			errors_before=validator_count;
			for(i=0;i<(*rows);i++){
				for(j=0;j<(*cols);j++){	
					fscanf(input_file,"%d	%lf",&lattice[i][j],&magn_field[i][j]);
					validator_count += validate_spin_value(lattice[i][j]);
				}
			}
			log_section_error(errors_file,validator_count-errors_before,"ELEMENTS");
		}
		//READ "CORRELATION SPIN"
		if(n==33){
			errors_before=validator_count;
			fscanf(input_file,"%lf",&corr_spin_d);
			*corr_spin=(int)corr_spin_d;
			validator_count+=validate_int(corr_spin_d) + validate_spin_range(*rows,*cols,*corr_spin);
			log_section_error(errors_file,validator_count-errors_before,"CORRELATION SPIN");
		}
		//READ "RANDOM NUMBERS"
		if(n==48){
			errors_before=validator_count;
			for(i=0;i<(*n_steps);i++){
					fscanf(input_file,"%lf	%lf",&random_numbers[i][0],&random_numbers[i][1]);
					validator_count += validate_rnd_number(random_numbers[i][0]) + validate_rnd_number(random_numbers[i][1]) ;
				}
			log_section_error(errors_file,validator_count-errors_before,"RANDOM NUMBERS");	
		}
		
		
		}
		return validator_count;
	
		
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
//picks two random numbers from the rows and cols
//the probability to choose a site should be equal
//higher random number---->higer site index
// n=rows*cols  ,   i=0...n-1; 
// should be True:  i/n <= rnd <= (i+1)/n for the correct i
// or i <= rnd*n <= i+1
// and therefore i= lower integer part of rnd*n
void choose_spin(double random_num,int rows, int cols ,int* row ,int* col){
	int n=rows*cols;
	int site_number;	
	site_number = random_num ==1 ? n-1 : (int)(random_num*n);
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

//recalculate the average
//<A>n = <A>(n-1) * (n-1)/n  + An/n
//item_index is the index of the items. first one is 1 last is n
// [A1,A2,A3.......An]
double update_thermal_average(double thermal_average_n_minus_one,int item_index, double new_item){
	if(item_index==1){
		return new_item;
	}
	return thermal_average_n_minus_one*(item_index-1)/(double)(item_index) + (new_item/item_index);
}
double thermal_stdev_new(double average,double average_squared){
	double variance = average_squared - pow(average,2);
	double stdev=pow(variance ,0.5);
	if(stdev!=stdev){ //stdev is NaN,probably problem of double size precision
		return 0;
	}
	return stdev;
}


void update_correlation_averages(double corr_data[3][40],int lattice[40][40], int rows,int cols,int corr_spin,int item_index){
	int _c;
	int corr_spin_r,corr_spin_c;
	corr_spin=corr_spin-1;//my index is 0, there index start with 1
	cell_index_to_row_col(corr_spin,rows,cols,&corr_spin_r,&corr_spin_c);
	if(item_index==1){
			for(_c=0;_c<cols;_c++){
			corr_data[0][_c]= (double)lattice[corr_spin_r][corr_spin_c]*(double)lattice[corr_spin_r][_c]; //<SmSl>
			corr_data[1][_c]= (double)lattice[corr_spin_r][corr_spin_c]; //<Sl>
			corr_data[2][_c]= (double)lattice[corr_spin_r][_c]; //<Sm>
			}
	}
	else{
		for(_c=0;_c<cols;_c++){

			corr_data[0][_c]=corr_data[0][_c]*(item_index-1)/(double)(item_index) +
							(((double)lattice[corr_spin_r][corr_spin_c])*((double)(double)lattice[corr_spin_r][_c])/(double)item_index); //<SmSl>
			
			corr_data[1][_c]=corr_data[1][_c]*(item_index-1)/(double)(item_index) +
							((double)lattice[corr_spin_r][corr_spin_c]/(double)item_index); //<Sl>
							
			corr_data[2][_c]=corr_data[2][_c]*(item_index-1)/(double)(item_index) +
							((double)lattice[corr_spin_r][_c]/(double)item_index); //<Sm>
		}
	}
}

void print_G(FILE *output,double corr_data[3][40],int corr_spin,int rows,int cols){
	int _c;
	int corr_spin_r,corr_spin_c;
	corr_spin=corr_spin-1;//my index is 0, there index start with 1
	cell_index_to_row_col(corr_spin,rows,cols,&corr_spin_r,&corr_spin_c);
	double G;
	for(_c=0;_c<cols;_c++){
		G=corr_spin_c == _c ? 1 : corr_data[0][_c]-corr_data[1][_c]*corr_data[2][_c];
		fprintf(output,"%lf	",G);
	}
	
}

//runs the bussiness logic
int start(FILE *input,FILE *error_output,FILE *log_file,FILE *meas, FILE *direc, FILE *corr){
	//given variables
	int rows,cols,n_taus,n_steps,n_measurments,corr_spin;
	double J,t_min,t_max;
	char boundry_conditions;
	int original_lattice[40][40];
	double magn_field[40][40];
	double random_numbers[300000][2];
	logg("reading initial data",log_file);
	int errors=read_input_data(input,error_output,&rows,&cols,&J,&t_min,&t_max,&n_taus,&n_steps,&n_measurments,&boundry_conditions,original_lattice,magn_field,&corr_spin,random_numbers);
	if(errors>0){
		fprintf(error_output,"ERROR: we cannot proceed! you have %d errors on your input data\n",errors);
		return 1;
	}
	logg("finish reading initial data",log_file);
	
	//dynamic variables
	int step_counter; //in which step are we?
	int lattice[40][40];
	copy_arr(original_lattice,lattice); //duplicate the original lattice
	double tau;
	double current_energy,current_magnetization;
	double energy_average_n=0,magnetization_average_n=0;
	double energy_squared_average_n=0,magnetization_squared_average_n=0;
	double corr_data[3][40]; //first col is for <SlSm>, second is for <Sl>, third is for <Sm>
							   // each row corresponds to the interaction of corr_spin with someone of its row
	
	logg("prints files headers",log_file);
	//meas.magn direc.spin
	fprintf(meas,"%d	%lf	%d	%d	%d	%d\n",n_taus,J,1,rows,cols,n_measurments);
	fprintf(direc,"%d	%lf	%d	%d	%d	%d\n",n_taus,J,1,rows,cols,n_measurments);
	fprintf(meas,"H=\n");
	fprintf(direc,"H=\n");
	print_magn_field(meas,magn_field,rows,cols);
	print_magn_field(direc,magn_field,rows,cols);
	
	for(tau=t_min;tau <= t_max;tau = next_tau(tau,t_min,t_max,n_taus)){
		logg("running simulation for a specific temprature",log_file);
		reset_system(original_lattice,lattice);
		//meas.magn  direc.spin
		fprintf(meas, "%lf\n",tau);
		fprintf(direc, "%lf\n",tau);
		fprintf(corr, "%lf\n",tau);
		
		for(step_counter=0;step_counter<n_steps;step_counter++){
			
			current_energy=E(boundry_conditions,rows,cols,lattice,magn_field,J);
			current_magnetization=M(lattice,rows,cols);
			
			energy_average_n=update_thermal_average(energy_average_n,step_counter+1,current_energy);
			energy_squared_average_n=update_thermal_average(energy_squared_average_n,step_counter+1,pow(current_energy,2));

			magnetization_average_n=update_thermal_average(magnetization_average_n,step_counter+1,current_magnetization);
			magnetization_squared_average_n=update_thermal_average(magnetization_squared_average_n,step_counter+1,pow(current_magnetization,2));
			
			update_correlation_averages(corr_data,lattice,rows,cols,corr_spin,step_counter+1);
			
			if(need_to_print_measurment(step_counter,n_steps,n_measurments)){
				//meas.magn
				fprintf(meas,"%d	%lf	%lf	%lf	%lf\n",
									  step_index(step_counter,n_steps,n_measurments),
									  current_energy,
									  thermal_stdev_new(energy_average_n,energy_squared_average_n),
									  current_magnetization,
									  thermal_stdev_new(magnetization_average_n,magnetization_squared_average_n));
				//direc.spin
				fprintf(direc,"%d\n", step_index(step_counter,n_steps,n_measurments));
				print_lattice(direc,lattice,rows,cols);
				//corr.corr
				fprintf(corr,"%d	", step_index(step_counter,n_steps,n_measurments));
				print_G(corr,corr_data,corr_spin,rows,cols);
				fprintf(corr,"\n");
			}
			do_step(random_numbers[step_counter][0],random_numbers[step_counter][1],rows,cols,lattice,magn_field,J,tau,boundry_conditions);
		}
		//meas.magn
		fprintf(meas, "result:\n");
		fprintf(meas, "%lf	%lf	%lf	%lf \n",
									energy_average_n,
									thermal_stdev_new(energy_average_n,energy_squared_average_n),
									magnetization_average_n,
									thermal_stdev_new(magnetization_average_n,magnetization_squared_average_n));
		logg("done with that temprature",log_file);
	}
	logg("finish all simulations",log_file);
	return 0;
}
//returns 0 if program finished successfully, 1 otherwise
int main(int argc, char **argv)
{
	//prepare I/O files
	FILE *input_file = stdin; 
	FILE *errors_file = stderr;
	FILE *progress_file = stdout;
	FILE *meas_output_file = fopen("./meas.magn","w"); 
	FILE *direc_output_file = fopen("./direc.spin","w"); 
	FILE *corr_output_file = fopen("./corr.corr","w"); 
	
	int status=start(input_file, errors_file,progress_file, meas_output_file, direc_output_file, corr_output_file);
	
	if(status==0){
		logg("finished program successfully",progress_file);
	}
	else{
		logg("finished program with errors",progress_file);
	}
	
	return status;
	

}
