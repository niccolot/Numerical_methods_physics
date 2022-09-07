#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define nlatt 100

/*
SIMULAZIONE OSCILLATORE ARMONICO QUANTISTICO 
L = 1/2m(dx/dt)^2 - 1/2m*omega^2(x^2)
*/

int npp[nlatt], nmm[nlatt];
double field[nlatt];

void geometry();
void initialize_lattice(int iflag);
void update_metropolis(double eta, double d_metro);
void measures(double *obs0, double *obs1, double *obs2, double *obsk1, double *obsk2, int k);


int main(){
	
	srand(time(NULL));
	int iflag = 1;
	int nstat = 1000000;
	int nterm = 1000;//termalizzazione
	int idecorrel = nlatt;
	double eta = 0.2;
	double d_metro = 2*sqrt(eta);
	
	double obs0;
	double obs1;
	double obs2;
	double obsk1;
	double obsk2;
		
	geometry();
	initialize_lattice(iflag);
	
	FILE *fp;
	fp=fopen("hist.txt","w");

    
	for(int i_term = 1; i_term <= nterm; i_term++){
			update_metropolis(eta,d_metro);
	}
	
	    
	for(int iter = 1; iter <= nstat; iter++){
	
		for(int idec = 1; idec<= idecorrel; idec++){
			update_metropolis(eta,d_metro);
		}

		measures(&obs0,&obs1,&obs2,&obsk1,&obsk2,0);
		fprintf(fp,"%lf\n",obs0);
	}
	
	

	
	return 0;
}


void geometry(){

    extern int npp[nlatt];
    extern int nmm[nlatt];

    for(int i=0;i<nlatt;i++){
        npp[i] = i+1;
        nmm[i] = i-1;
    }

    npp[nlatt-1] = 0;
    nmm[0] = nlatt-1;

    return;
}


void initialize_lattice(int iflag){

	double x;
	extern double field[nlatt];
	
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
		
		double phi_prova = phi + 2.0*d_metro*(0.5-x);
		
		double p_rat = c1*phi_prova*force-c2*pow(phi_prova,2);
		p_rat = p_rat - c1*phi*force+c2*pow(phi,2);
		
		double y = log(rand()/(RAND_MAX+1.0));
		
		if(y<p_rat){
			field[i] = phi_prova;
		}
	}
	return;
}


void measures(double *obs0, double *obs1,double *obs2,double *obsk1, double *obsk2,int k){

/*obsk1 è il correlatore a 2 punti <X_i X_i+k>
obsk2 è il correlatore a due punti con le grandezze al quadrato <(X_i**2 X_i+k**2>
per i correlatori connessi va sottratto a obsk1 e obsk2 rispettivamente 
<x>**2, nullo se l' energia è pari e <x**2>**2, <x**2> si calcola già con obs1,
quindi resta solo da ELEVARALA AL QUADRATO DOPO AVER FATTO LA MEDIA*/
	
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









