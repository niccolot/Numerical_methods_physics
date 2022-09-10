#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#define nlatt 40//side lenth of the square lattice
#define nvol 1600//nlatt*nlatt

/*
Markov chain Montecarlo simulation for the 2d square lattice ising model
*/

double field[nlatt][nlatt];
int npp[nlatt], nmm[nlatt];

void initialize_lattice(int iflag);
void update_metropolis(double beta, double extfield);
void geometry();//instantiate the 'step arrays' that encode the pbc
void measures(double extfield, double *mag, double *mag2, double *mag4, double *ene, double *ene2);


int main(void){
	
	int numbeta = 56;
    srand(time(NULL));
	/*
	In order to plot the phase transition here we make a file for different thermodinamic beta values (beta ~ 1/temperature),
	later one can extract the various observables (energy, magnetization etc.) and plot some graphs to show the behaviour at
	different temperratures
	*/
    double beta[] = {0.360,0.370,0.380,0.385,0.390,0.395,0.400,
    				0.4025,0.405,0.4075,0.410,0.4125,0.4150,0.4175,0.420,
    				0.4225,0.4250,0.4275,0.430,0.4325,0.4350,0.4360,0.4375,
    				0.4380,0.4390,0.4395,0.4398,0.4399,0.440,0.4401,0.4402,0.4403,0.4406,0.4409,0.4412,
    				0.4425,0.4437,0.4450,0.4462,0.4475,0.4487,0.450,
    				0.4525,0.4550,0.4575,0.460,0.4625,0.4650,0.4675,0.470,
    				0.475,0.480,0.485,0.490,0.500,0.510};
    				
    double extfield = 0.0;

    double mag,mag2,mag4,ene,ene2;

    int iflag = 1;//for the starting lattice, random or with all sites fixed to a point ('warm' or 'cold' initialization)
    int nstat = 100000;//datapoints to be measured
	int nterm = 1000;//measures to be discarded in order to let the markov chain thermalize
    int i_decorrel = 100;//lattice updates between 2 measures in order to decorrelate the data, idelly one should at least update the whole lattice
						//here due to the very time comsuming execution of the program the decorrelation is only on 100 sites

    FILE *fp_data[numbeta];
    data[0] = fopen("ising_360L40.txt","w");
    data[1] = fopen("ising_370L40.txt","w");
    data[2] = fopen("ising_380L40.txt","w");
    data[3] = fopen("ising_385L40.txt","w");
    data[4] = fopen("ising_390L40.txt","w");
	data[5] = fopen("ising_395L40.txt","w");
	data[6] = fopen("ising_400L40.txt","w");
	data[7] = fopen("ising_4025L40.txt","w");
	data[8] = fopen("ising_405L40.txt","w");
	data[9] = fopen("ising_4075L40.txt","w");
	data[10] = fopen("ising_410L40.txt","w");
	data[11] = fopen("ising_4125L40.txt","w");
	data[12] = fopen("ising_4150L40.txt","w");
	data[13] = fopen("ising_4175L40.txt","w");
	data[14] = fopen("ising_420L40.txt","w");
	data[15] = fopen("ising_4225L40.txt","w");
	data[16] = fopen("ising_4250L40.txt","w");
	data[17] = fopen("ising_4275L40.txt","w");
	data[18] = fopen("ising_430L40.txt","w");
	data[19] = fopen("ising_4325L40.txt","w");
	data[20] = fopen("ising_4350L40.txt","w");
	data[21] = fopen("ising_4360L40.txt","w");
	data[22] = fopen("ising_4375L40.txt","w");
	data[23] = fopen("ising_4380L40.txt","w");
	data[24] = fopen("ising_4390L40.txt","w");
	data[25] = fopen("ising_4395L40.txt","w");
	data[26] = fopen("ising_4398L40.txt","w");
	data[27] = fopen("ising_4399L40.txt","w");
	data[28] = fopen("ising_440L40.txt","w");
	data[29] = fopen("ising_4401L40.txt","w");
	data[30] = fopen("ising_4402L40.txt","w");
	data[31] = fopen("ising_4403L40.txt","w");
	data[32] = fopen("ising_4406L40.txt","w");
	data[33] = fopen("ising_4409L40.txt","w");
	data[34] = fopen("ising_4412L40.txt","w");
	data[35] = fopen("ising_4425L40.txt","w");
	data[36] = fopen("ising_4437L40.txt","w");
	data[37] = fopen("ising_4450L40.txt","w");
	data[38] = fopen("ising_4462L40.txt","w");
	data[39] = fopen("ising_4475L40.txt","w");
	data[40] = fopen("ising_4487L40.txt","w");
	data[41] = fopen("ising_450L40.txt","w");
	data[42] = fopen("ising_4525L40.txt","w");
	data[43] = fopen("ising_4550L40.txt","w");
	data[44] = fopen("ising_4575L40.txt","w");
	data[45] = fopen("ising_460L40.txt","w");
	data[46] = fopen("ising_4625L40.txt","w");
	data[47] = fopen("ising_4650L40.txt","w");
	data[48] = fopen("ising_4675L40.txt","w");
	data[49] = fopen("ising_470L40.txt","w");
	data[50] = fopen("ising_475L40.txt","w");
	data[51] = fopen("ising_480L40.txt","w");
	data[52] = fopen("ising_485L40.txt","w");
	data[53] = fopen("ising_490L40.txt","w");
	data[54] = fopen("ising_500L40.txt","w");
	data[55] = fopen("ising_510L40.txt","w");
	
    geometry();
    initialize_lattice(iflag);
	
	//thermalization    
	for(int i_term = 1; i_term <= nterm; i_term++){
			update_metropolis(eta,d_metro);
	}

    for(int b=0;b<numbeta;b++){
    
    	for(int iter=1;iter<=nstat;iter++){
    	
        	for(int idec=1;idec<=i_decorrel;idec++){
        	update_metropolis(beta[b],extfield);
        	}

        	measures(extfield,&mag,&mag2,&mag4,&ene,&ene2);
        	fprintf(fp_data[b],"%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",iter,mag,mag2,mag4,ene,ene2);
    	}
    }
	
	for(int b=0;b<numbeta;b++){
		fclose(fp_dati[b]);
	}

    return 0;
}

void initialize_lattice(int iflag){

    double x;

    extern double lattice[nlatt][nlatt];

    switch(iflag){
    //case 1 is 'warm' start with random values given to the sites
	//case 2 is 'cold' start with sites initially frozen in the arbitrary value of 1.0
        case 1:{
            for(int i=0;i<nlatt;i++){
                for(int j=0;j<nlatt;j++){
                    x = rand()/(RAND_MAX+1.0);
                    field[i][j] = 1.0;
                    if(x<0.5) field[i][j] = -1.0;
                }
            }
        break;
        }

        case 0:{
            for(int i=0;i<nlatt;i++){
                for(int j=0;j<nlatt;j++){
                    field[i][j] = 1.0;
                }
            }
        break;
        }

        default:{
            for(int i=0;i<nlatt;i++){
                for(int j=0;j<nlatt;j++){
                    field[i][j] = 1.0;
                }
            }
        break;
        }
    }
    return;
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

void measures(double extfield, double *mag, double *mag2, double *mag4, double *ene, double *ene2){

    extern double lattice[nlatt][nlatt];
    extern int npp[nlatt];
    extern int nmm[nlatt];

    *mag = 0.0;
    *ene = 0.0;
    double force;
    int ip,im,jp,jm;

    for(int i=0;i<nlatt;i++){
        for(int j=0;j<nlatt;j++){

            *mag = *mag+field[i][j];

            int ip = npp[i];
            int im = nmm[i];
            int jp = npp[j];
            int jm = nmm[j];

            force = field[i][jp]+field[i][jm]+field[ip][j]+field[im][j]; 
            *ene = *ene - 0.5*force*field[i][j] - extfield*field[i][j];
        }
    }

    *mag = *mag/((double) nvol);
    //in a finite size lattice one uses |M|, the absolute value of the magnetization as order parameter
    if (*mag<0.){
    	*mag = -*mag;
    }
    *mag2 = pow(*mag,2);
    *mag4 = pow(*mag,4);
    *ene = *ene/((double) nvol);
    *ene2 = pow(*ene,2);

    return;
}

//single metropolis step
void update_metropolis(double beta, double extfield){

    extern double lattice[nlatt][nlatt];
    extern int npp[nlatt];
    extern int nmm[nlatt];

    for(int ivol=1;ivol<=nvol;ivol++){
        
        double x = rand()/(RAND_MAX+1.0);
        double y = rand()/(RAND_MAX+1.0);
        int i = (int) (x*nlatt);
        int j = (int) (y*nlatt);
 
        int ip = npp[i];
        int im = nmm[i];
        int jp = npp[j];
        int jm = nmm[j];
 
        double force = field[i][jp]+field[i][jm]+field[ip][j]+field[im][j];
        force = beta*(force+extfield);

        double phi = field[i][j];
        double p_rat = exp(-2.0*phi*force);

        double z = rand()/(RAND_MAX+1.0);

        if(z<p_rat){
            field[i][j] = -phi;
        } 
    }
    return;
}
