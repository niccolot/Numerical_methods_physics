#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define n_meas 131072//to be a power of 2

/*
Given a[n_meas] array of n_meas size, it returns the VARIANCE with bootstrap algorithm

For correlated data it returns the 'binned' variation of the algo with 2**k size of the bin
*/

double bootbin(double a[], int k);

int main(){
	
	/*
	Here it is assumed one has a .txt file with float data arranged by columns called data.txt, if one has a file with multiple columns use the %*type 
	placeholder in order to ignore some columns
	
	At the end it saves the STANDARD DEVIATION in a separate file error.txt
	
	It is assumed the data.txt file has a power of 2 of datapoints in order to be divided in bins
	*/
	
	double stddev;
	FILE *fp_data;
	FILE *fp_errors;
	double k = 10; //2**k bin_size
	
	fp_data = fopen("data.txt", "r");
	fp_error = fopen("error.txt","w");
	
	
	for(int i=0;i<n_meas;i++){
		fscanf(fp_data,"%f",&a[i]);
	}
	
	stddev = sqrt(bootbin(a,k))
	fprintf(fp_error,"%f\n",stddev)

	return 0;

}

/*
It returns the variance with a 2**k binnig
*/
double bootbin(double a[], int k){
	
	int m = 100;//number of resampling, usually 50/100 are enough
	double b[n_meas];//it will contain the resampled data
	double f[m];//it will contain the observables one has to calculate the variance of
	
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
		
		/*
		Here it is assumed one as a series of n_meas x_i datapoints and one has to calculate the variance of the observable <x> = mean(x_i)
		
		For more complicated observables like the binder cumulant <x**4>/3<x**2>**2 one as to modify this point in order that f[l] contains the
		observable calculated for each resampling, 
		
		e.g. in the binder cumulant casa one needs two arrays b1[], b2[] during the resampling, one to store 
		the x**4 values and the other for the x**2, than you take the mean of the two and feed f[l] mean(b1)/3mean(b2) for each l = 1...m resampling 
		*/
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
	
	double s2k = dum1-dum2;//variance formula
	
	return s2k;
}
