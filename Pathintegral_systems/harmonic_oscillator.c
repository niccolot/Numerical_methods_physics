#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define nlatt 100//lenght of the lattice

/*
Markov chain Montecarlo simulation for harmonic oscillator with lagrangian
L = 1/2m(dx/dt)^2 - 1/2m*omega^2(x^2)
using metropolis algorithm
*/

int npp[nlatt], nmm[nlatt];
double field[nlatt];

void geometry();//instantiate the 'step arrays' that encode the pbc
void initialize_lattice(int iflag);
void update_metropolis(double eta, double d_metro);//single metropolis step
void measures(double *obs0, double *obs1, double *obs2, double *obsk1, double *obsk2, int k);//returns observables to be measured


int main(){
	
	srand(time(NULL));
	int iflag = 1;//for the starting lattice, random or with all sites fixed to a point ('warm' or 'cold' initialization)
	int nstat = 1000000;//datapoints to be measured
	int nterm = 1000;//measures to be discarded in order to let the markov chain thermalize
	int idecorrel = nlatt;//lattice updates between 2 measures in order to decorrelate the data, one should at least update the whole lattice
	double eta = 0.2;
	double d_metro = 2*sqrt(eta);//optimized parameter for the metropolis step
	
	double obs0;
	double obs1;
	double obs2;
	double obsk1;
	double obsk2;
	double k = 2;//for the correlator <x_i*x_{i+k}>, given pbc ence discrete translation invariance, this is the same for all i
		
	geometry();
	initialize_lattice(iflag);
	
	FILE *fp_data;
	fp_data=fopen("data.txt","w");//one can save the measurements organized in columns

	//thermalization    
	for(int i_term = 1; i_term <= nterm; i_term++){
			update_metropolis(eta,d_metro);
	}
	
	//measures
	for(int iter = 1; iter <= nstat; iter++){
		//decorrelation
		for(int idec = 1; idec<= idecorrel; idec++){
			update_metropolis(eta,d_metro);
		}

		measures(&obs0,&obs1,&obs2,&obsk1,&obsk2,k);
		fprintf(fp_data,"%lf\t%lf\t%lf\t%lf\t%lf\n",obs0,obs1,obs2,obsk1,obsk2);
	}
	
	fclose(fp_data);
	
	return 0;
}


void geometry(){

    extern int npp[nlatt];
    extern int nmm[nlatt];

    for(int i=0;i<nlatt;i++){
        npp[i] = i+1;
        nmm[i] = i-1;
    }
	//pbc
    npp[nlatt-1] = 0;
    nmm[0] = nlatt-1;

    return;
}


void initialize_lattice(int iflag){

	double x;
	extern double field[nlatt];
	
	//case 1 is 'warm' start with random values given to the sites
	//case 2 is 'cold' start with sites initialli frozen in the arbitrary value of 10.0
	switch(iflag){
	
		case 1:{
			for(int i=0;i<nlatt;i++){
				x = rand()/(RAND_MAX+1.0);
				x = 1.0 - 2.*x;
				field[i]=x;
			}
		break;
		}
		
		case 0:{
			for(int i=0;i<nlatt;i++){
				field[i]=10.0;
			}
		break;
		}
		
		default:{
			for(int i=0;i<nlatt;i++){
				x = rand()/(RAND_MAX+1.0);
				x = 1.0 - 2.*x;
				field[i]=x;
			}
		break;
		}
	}
	return;
}

//metropolis step for the harmonic oscillator
void update_metropolis(double eta, double d_metro){

	extern double field[nlatt];
    extern int npp[nlatt];
    extern int nmm[nlatt];

	double c1 = 1.0/eta;
	double c2 = 1.0/eta + eta/2.0;
	double x;
	
	for(int i=0;i<nlatt;i++){
	
		int ip = npp[i];
		int im = nmm[i];
		
		double force = field[ip]+field[im];
		double phi = field[i];
		
		x = rand()/(RAND_MAX+1.0);
		
		double phi_try = phi + 2.0*d_metro*(0.5-x);
		
		double p_rat = c1*phi_try*force-c2*pow(phi_try,2);
		p_rat = p_rat - c1*phi*force+c2*pow(phi,2);
		
		double y = log(rand()/(RAND_MAX+1.0));
		
		if(y<p_rat){
			field[i] = phi_try;
		}
	}
	return;
}


void measures(double *obs0, double *obs1,double *obs2,double *obsk1, double *obsk2,int k){

/*
obs0 = <x>
obs1 = <x**2>
obs2 = <x_i*x_{i-1}> correlator, i is not important given the pbc ence discrete translation invariance
obsk1 = <x_i*x{i+k}>
obsk2 = <[(x_i)**2]*[(x_{i+k})**2]>, the correlator with squared inputs
*/
	
	extern double field[nlatt];
    extern int npp[nlatt];
    extern int nmm[nlatt];

	*obs0 = *obs1 = *obs2 = *obsk1  = *obsk2 = 0.0;
	
	for(int i=0;i<nlatt;i++){
	
		*obs0 = *obs0 + field[0];
		*obs1 = *obs1 + pow(field[i],2);
		*obs2 = *obs2 + pow(field[i]-field[npp[i]],2);
		*obsk1 = *obsk1 + field[i]*field[(i+k % nlatt)];
		*obsk2 = *obsk2 + pow(field[i],2)*pow(field[((i+k)%nlatt)],2);

	}
	
	*obs0 = *obs0/((double) nlatt);
	*obs1 = *obs1/((double) nlatt);
	*obs2 = *obs2/((double) nlatt);
	*obsk1 = *obsk1/((double) nlatt);
	*obsk2 = *obsk2/((double) nlatt);
	
	return;
}
