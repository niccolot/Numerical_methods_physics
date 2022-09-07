#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define n_meas 131072

/*
dato array di dati a[n_meas], n_meas numero di dati, restituisce la varianza calcolata 
con algoritmo bootstrap s2k con binning di blocco 2**k
*/ 

/*se il file da analizzare contiene più colonne e si è interessati solo
ad alcune, il placeholder %*type ignora l' elemento di tipo 'type'*/

double bootbin(double a[], int k);

int main(){
	
	double a[n_meas];
	double stddev;
	int nlatt = 40;
	
	FILE *fp_err_N40;
	
	FILE *fp_dati_corr[nlatt];
    
	fp_dati_corr[0] = fopen("corrk0n40.txt","r");
	fp_dati_corr[1] = fopen("corrk1n40.txt","r");
	fp_dati_corr[2] = fopen("corrk2n40.txt","r");
	fp_dati_corr[3] = fopen("corrk3n40.txt","r");
	fp_dati_corr[4] = fopen("corrk4n40.txt","r");
	fp_dati_corr[5] = fopen("corrk5n40.txt","r");
	fp_dati_corr[6] = fopen("corrk6n40.txt","r");
	fp_dati_corr[7] = fopen("corrk7n40.txt","r");
	fp_dati_corr[8] = fopen("corrk8n40.txt","r");
	fp_dati_corr[9] = fopen("corrk9n40.txt","r");
	fp_dati_corr[10] = fopen("corrk10n40.txt","r");
	fp_dati_corr[11] = fopen("corrk11n40.txt","r");
	fp_dati_corr[12] = fopen("corrk12n40.txt","r");
	fp_dati_corr[13] = fopen("corrk13n40.txt","r");
	fp_dati_corr[14] = fopen("corrk14n40.txt","r");
	fp_dati_corr[15] = fopen("corrk15n40.txt","r");
	fp_dati_corr[16] = fopen("corrk16n40.txt","r");
	fp_dati_corr[17] = fopen("corrk17n40.txt","r");
	fp_dati_corr[18] = fopen("corrk18n40.txt","r");
	fp_dati_corr[19] = fopen("corrk19n40.txt","r");
	fp_dati_corr[20] = fopen("corrk20n40.txt","r");
	fp_dati_corr[21] = fopen("corrk21n40.txt","r");
	fp_dati_corr[22] = fopen("corrk22n40.txt","r");
	fp_dati_corr[23] = fopen("corrk23n40.txt","r");
	fp_dati_corr[24] = fopen("corrk24n40.txt","r");
	fp_dati_corr[25] = fopen("corrk25n40.txt","r");
	fp_dati_corr[26] = fopen("corrk26n40.txt","r");
	fp_dati_corr[27] = fopen("corrk27n40.txt","r");
	fp_dati_corr[28] = fopen("corrk28n40.txt","r");
	fp_dati_corr[29] = fopen("corrk29n40.txt","r");
	fp_dati_corr[30] = fopen("corrk30n40.txt","r");
	fp_dati_corr[31] = fopen("corrk31n40.txt","r");
	fp_dati_corr[32] = fopen("corrk32n40.txt","r");
	fp_dati_corr[33] = fopen("corrk33n40.txt","r");
	fp_dati_corr[34] = fopen("corrk34n40.txt","r");
	fp_dati_corr[35] = fopen("corrk35n40.txt","r");
	fp_dati_corr[36] = fopen("corrk36n40.txt","r");
	fp_dati_corr[37] = fopen("corrk37n40.txt","r");
	fp_dati_corr[38] = fopen("corrk38n40.txt","r");
	fp_dati_corr[39] = fopen("corrk39n40.txt","r");
	
	fp_err_N40 = fopen("err_N40.txt","w");
	
	for(int i=0;i<nlatt;i++){
	
		for(int j=0;j<n_meas;j++){
			fscanf(fp_dati_corr[i],"%*d%lf%*d%*f%*f",&a[j]);
		}
		stddev = sqrt(bootbin(a,10));
		fprintf(fp_err_N40,"%f\n",stddev);
	}


	return 0;

}

/*funzione per bootstrap con binning, restituisce la varianza s2k dei
dati con bin 2**k*/
//a[] array di dati da elaborare
//k potenza dimensione blocco, dim = 2**k
double bootbin(double a[], int k){
	
	int m = 100;//m numero ricampionamenti, in genere 50/100 bastano
	double b[n_meas];//conterrà il set di dati ricampionato
	double f[m];//conterrà i risultati di ogni ricampionamento
	
	for(int l=0;l<m;l++){
	
		for(int i=0;i<n_meas;i+=pow(2,k)){
		
			int x = (rand()%(n_meas));
			while(x >= (n_meas - (int)pow(2,k)+1)){
				x = (rand()%(n_meas));
			}
		
			b[i] = a[x];
			for(int j=1;j<=(int)pow(2,k);j++){
				b[i+j] = a[x+j];
			}
		}

		double mean = 0.0;
		for(int i=0;i<n_meas;i++){
			mean += b[i];
		}	
		
		f[l] = mean/(double)n_meas;
	}
	
	double dum1 = 0.0;
	double dum2 = 0.0;
	
	for(int i=0;i<m;i++){
		dum1 += pow(f[i],2);
		dum2 += f[i];
	}
	
	dum1 = dum1/(double)n_meas;
	dum2 = pow(dum2/(double)n_meas,2);
	
	double s2k = dum1-dum2;//varianza
	
	return s2k;
}







