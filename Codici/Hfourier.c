#define RROT																				/*Refine the rotation of the snapshot --> -DRROT*/
#define A2PROF 																			/*Write a file 'A2_profile' for each snapshot --> -DA2PROF*/	
//#define MASS 																				/*Perform the mass analysis --> -DMASS*/
//#define BOXY 																				/*Performs an analisys of the bar structure --> -DBOXY*/
#define EQNPART																				/*Set an equal number of particle for each bin*/
#define SMOOTH 																				/*Perform a smoothing on the A2 and the phase differential profiles*/
#define DISP 																				/*Evaluate the velocity dispersion ratio within the bar*/
//#define PSNAP																				/*Print the snapshot in a tipsy-format with the galaxy translated and rotated to the centre of the galaxy*/
#define LIM
#define REFINE																				/*Refine the evaluation of the bar properties*/
//#define DEBUG																				/*Debugging key --> -DDEBUG*/

#ifdef _OPENMP
#include <omp.h>																			/*Including the openMP library*/
#endif
#include "thdf5.h"																			/*Including XDR-binary reader library*/
#include "tmisc.h"																			/*Including utility library*/
#include <unistd.h>
#ifdef DEBUG
#include <time.h>																			/*In debugging session this library marks time*/
#endif

#define kpc_to_cm 3.08567758130573e21														/*Centimeters in a kpc*/
#define msol_to_g 1.98892e33 																/*Mass of the sun in g*/
#define km_to_cm 1e5																		/*Centimeters in a km*/

#define NA 4294967295																		/*Not assigned*/

//----------------------------------------------------------------------------------------------------------------------------------------------

double dist_f2 (const double first[DIM], const double second[DIM])							/*Evaluating the 2-distance between a float and a double*/
{
	int i;
	double d = 0.;
	for (i=0; i<DIM-1; i++)
	{
		d += (first[i]-second[i])*(first[i]-second[i]);
	}
	
	return sqrt(d);
}

//----------------------------------------------------------------------------------------------------------------------------------------------

#ifdef EQNPART
int compare (const void *pa, const void *pb)												/*Function to compare the distance between particles*/
{
    const double *a = *(const double **)pa;
    const double *b = *(const double **)pb;
    return 1e4*(a[1]-b[1]);
}
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------

#ifdef REFINE
double avg_phase (double phase[], int first, int last)										/*Function to compute the average phase of the bar within an interval [first; last]*/
{
	double avg_phase = 0.;
	double tempval = phase[first];															/*Set the first temporary buffer to the first phase bin of the bar*/
	int j;
	const double ra = M_PI*0.5;

	for (j=first; j<=last; j++)																/*Cycle over the bar, till the peak*/
	{
		if (fabs(tempval-phase[j])>ra)														/*If the displacement between the current phase bin and the previous one is larger than 90°*/
		{	
			if (tempval>0)																	/*and the previous phase bin was positive*/
			{
				avg_phase += phase[j]+2*ra;													/*Sum up the current values shifting it*/
				tempval = phase[j]+2*ra;
			}
			else
			{
				avg_phase += phase[j]-2*ra;
				tempval = phase[j]-2*ra;
			}
		}
		else
		{
			avg_phase += phase[j];
			tempval = phase[j];
		}
	}
	avg_phase/=(last-first+1);																/*Evaluate the mean*/
	if (avg_phase<-ra) avg_phase+=(2*ra);													/*Shift it in order to bring it back to [-90°; 90°]*/
	else if (avg_phase>ra) avg_phase-=(2*ra);

	return avg_phase;
}
#endif

//--------------------------------------------------------------------------------------------------------------------------------------------

int main (int argc, char *argv[])															/*Version of 'fourier.c' for the Illustris TNG (the snapshots are in HDF5 format). Compile with '-I/usr/include/hdf5/serial -lhdf5_cpp -lhdf5 -lm'*/
{
#ifdef _OPENMP																				/*Conditional compilation for OpenMP*/
	printf("\nOpenMP support: %d\n", _OPENMP);
#else
	printf("\nSerial execution\n");
#endif

	char input[200];																		/*Definition of the strings for the input file name*/
	char file_name[100];

		const char centrefile[] = "offsets";
		const char out_path[] = ".";														/*Path of the output files*/

	const int idlength = 6;																	/*Number of digits of the galaxy id*/
	char *id = (char *)malloc(idlength*sizeof(char));										/*Galaxy ID*/

	long int id_buff;
	double displ_cm[DIM];
	double displ_vcm[DIM];

	int i, j, n;																			/*Definition of the counters*/
#ifndef DEBUG
	FILE *turin;																			/*Output file pointer for 'fourier.out'*/
	char outfou[200];
#ifdef A2PROF	
	char A2out[200];																		/*Definition of the strings for the 'A2_profile...' output file name (one for each snapshot)*/	
#endif		
#else
	char A2smooth[200];
#endif

	double in_to_kpc, in_to_msol, in_to_kmpersec;

	double ran;																				/*Radius of the area within which to start evaluating the centre of mass evaluation*/
#ifdef RROT		
	double zin;																				/*Bin half-legth for the translation and rotation routines*/
#endif
	double zmax;																			/*Cylindrical bin height*/
	double rad;																				/*Radius of the area within which to perform the analysis*/
	double radcm;																			/*Radius actually used to evaluate the centre of mass*/
	double tempdist;																		/*Temporary distance to evaluate the centre of mass*/
	double draft[DIM], cm[DIM], vcm[DIM];													/*Dynamic centre position, centre of mass position and velocity*/
	double mtot;																			/*Total mass to evaluate the fourier mass integral*/
	double L_x, L_y, L_z;																	/*Angular momentum components*/
	double sint, cost, sinf, cosf;															/*Trigonometric quantities (used to rotate the reference frame and to store the fourier integrals)*/
	double vtan, Krot, K;

		//const double Krot_min = 0.4;														/*Minimum Krot to trust the rotation procedure*/
		//cons sigma_ratio_max = 1;
		//const u_d_max = 1.3;

#ifdef EQNPART
	int selpart;																			/*Number of selected particles*/
	int npart;																				/*Particles per bin*/
	double **list = NULL;																	/*List of the selected particles with id and distance*/
	double firstP;																			/*Temporary distance of the first particle in the bin*/

		const int min_npart = 300;															/*Threshold number of particles for each bin*/
#endif

	float TotStarMass;
	float StarMass;
	double hmr;
	double r;																				/*Bin width*/
	
	double **bin = NULL;;																	/*Matrix containing 1. bin radius, 2. total bin mass, 3. total bin cos2th, 4. total bin sin2th, 5. cumulative A2, 6. cumulative th*/
	double maxA2;																			/*Value of the maximum of A2 for the diferential quantities*/
	unsigned int BAR;																		/*Flag to state where is the bar*/
	double floor;
	double Bmin;																			/*Minimum length for the bar*/
	unsigned int RANGE;																		/*This flag tells if the bar actually belongs to the range where the peak occurs*/
	int flag;
	const double ra = M_PI*0.5;
	int npeaks;
	double req;																				/*Maximum distance from the centre imposed to the bar inner edge in pkpc*/																				/*Number of peaks counter*/
	double maxdist_prel;																		/*Maximum distance from the centre allowed for a preliminar peak*/
	double bar_param[5];	 	
	int KERN_a2t, KERN_a2, KERN_ph;															/*Smoothing kernels for the A2 profiles and phase*/
	int round_temp_kern;
	double temp_kern;
	int nbin;																				/*Number of bins to perform the analysis*/

		const int nbin_max = 400;															/*Initial/Maximum number of bins*/
		const double ran0 = 10.;															/*Updating the value of the radius for the rotation and translation routines*/
		const double Dth = asin(0.15);														/*Setting the maximum phase dispacement*/
		const int nmaxp = 100;																/*Maximum number of expected peaks in each snapshot*/

	double deriv[nbin_max];
	double smooth_der[nbin_max];															/*smoothA2 derivative*/
	double smooth_ph[nbin_max];																/*smoothA2 punctual phase*/
	double smooth_A2[nbin_max];
	int peaks[nmaxp+1];																		/*List of peak indices. The '+1' is to count the first annulus after the last bar*/
	double R_th[nmaxp+1][2];																/*List of the R_th extents found (if any), for each A2 peak, along with the starting bin index of the bar*/
	double Phi[nmaxp+1];
	int bar_first_bin[nmaxp+1];
	int bar_last_bin[nmaxp+1];


#ifdef REFINE
	double famp, tempval;
		const int mi = 3;

	float vphys[DIM];
	double rxy;
	double H, A, rm, wm, hm, mr;
	double MeanV_z, MeanV_R, Sigma_z, Sigma_R;
	int N_bar_part;

	double barUP, barDOWN;
#endif

#if defined LIM
	double lim[4];																			/*Limits on the galactic radius where to look for the bar edge*/
#endif

	FILE *huor;																				/*Output file pointer for 'A2_profile..' or A2smooth...*/

//----------------------------------------------------------------

#ifdef MASS

#ifndef DEBUG
	FILE *hurin;																			/*Output file pointer for 'mass.out'*/
	char outmass[200];																		/*'mass.out' file name*/
#endif
	double rad1, rad2, rad3, rad4;															/*Definition of the radii within which to calculate the integral mass*/
	double mass1g, mass2g, mass3g, mass4g;													/*Definition of the variables to contain the mass inside 4 different radii GAS*/
	double mass1d, mass2d, mass3d, mass4d;													/*Definition of the variables to contain the mass inside 4 different radii DARK MATTER*/
	double mass1s, mass2s, mass3s, mass4s;													/*Definition of the variables to contain the mass inside 4 different radii STAR*/
	double massBg, massBd, massBs;															/*Definition of the variables to contain the mass inside the bar*/
#ifndef REFINE
	double rm, wm, hm;																		/*Limits for the region to analyse*/
#endif
#endif

//----------------------------------------------------------------

#ifdef BOXY
#ifndef DEBUG	
	FILE *galadriel;																		/*Boxyness-output file pointer per snapshot*/
	char boxyout[200];																		/*Definition of the strings for the 'boxy_profile...' output file name (one for each snapshot)*/
#endif

		const float xbin_mult = 1;															
		const float zbin_mult = 2;
		const float bmult = 0.9;															/*Find the z-oordinate of the 100*bmult% of mass enclosed*/

	int xbin, zbin;																			/*Number of bins to subdivide the bar over the x-axis*/
	double **density;																		/*Definition of the pointer to the struct where to save the bar particles*/
	double *median_d, *median_u;															/*Definition of the matrix to save the density map of the bar*/
	double sum, half, delta;																/*Partial integral of the mass and half of the total*/
	double rb, wb, hb;
#endif

//----------------------------------------------------------------

#if defined PSNAP || defined DEBUG
	char outdeb[200];																		/*Output snapshot for the debugging session*/
	char bar_ids[200];
#endif

#ifdef DEBUG
	struct timespec t_i, t_par, t_stop;														/*Stopwatch times*/
#endif

#if defined BOXY  || defined MASS || defined REFINE
	double C1, C2, C3;																		/*Constants to contrain the bar*/
#endif

//----------------------------------------------------------------

#ifdef DEBUG
	clock_gettime(CLOCK_MONOTONIC, &t_i);													/*Initial time*/
#endif

//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------

	if (argc == 1)																			/*If there is no input file*/
	{
		printf("\nNo input file\n\n");
		exit(1);
	}
	else if (argc == 2) strcpy(input, argv[1]);
	else
	{
		printf("\n\nThe syntax of the command is incorrect:\nUse ./Hfourier file.in");
		exit(1);
	}

	strcpy(file_name, input);
	strtok(input, "_");
	id = strtok(NULL, ".");																	/*Copy the address of a string with the name of the snapshot*/

//------------------------------------------------------------------------------------------Creating the file names---------------------------

#ifndef DEBUG
	sprintf(outfou, "%s/fourier.out", out_path);

#ifdef MASS
	sprintf(outmass, "%s/mass.%s", out_path, id); 
#endif
#ifdef A2PROF
	sprintf(A2out, "%s/A2_profile.%s", out_path, id);
#endif
#ifdef BOXY 
	sprintf(boxyout, "%s/boxy_profile.%s", out_path, id);
#endif

#ifdef PSNAP
	sprintf(outdeb, "%s/snap.%s", out_path, id);
#endif

#else
	sprintf(A2smooth, "%s/A2_smooth.%s", out_path, id);
	sprintf(outdeb, "%s/snap.%s", out_path, id);
#endif

//------------------------------------------------------------------------------------------Memory allocation outside the main cycle------------

	allmat(&bin, 8, nbin_max);																/*Allocating the memory for the matrix with radius, bin mass, cos2th, sin2th, A2, th, softening and number of particles*/

//----------------------------------------------------------------------------------------------------------------------------------------------

	if (hdf5read(file_name)==TRUE)	printf("\nHDF5 file '%s' successfully read\n", file_name);		/*Snapshot reading functions*/
	else exit (1);

	in_to_kpc = (header.UnitLength_in_cm/kpc_to_cm)*(header.time/header.HubbleParam);		/*Conversion factor from internal units to kpc*/
	in_to_msol = (header.UnitMass_in_g/msol_to_g)/header.HubbleParam;						/*Conversion factor from internal units to solar masses*/
	in_to_kmpersec = (header.UnitVelocity_in_cm_per_s/km_to_cm)*sqrt(header.time);			/*Conversion factor from internal unit to km/sec*/

	ran = ran0/in_to_kpc;																	/*Updating the value of the radius for the rotation and translation routines*/
		Bmin = 1.4*gsoft(header.time);
		floor = 1.4*gsoft(header.time);														/*Getting the softening as threshold*/

#ifdef DEBUG
	clock_gettime(CLOCK_MONOTONIC, &t_par);													/*Reading and initializing time elapsed*/
	printf("(Elapsed: %.2f seconds)\n", (float)((t_par.tv_sec-t_i.tv_sec)+(t_par.tv_nsec-t_i.tv_nsec)/1.0e9));
#endif

//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------Finding the centre of mass of the galaxy------------

	flag=0;
	if (access(centrefile, R_OK)==0)														/*If the centre file is present*/
	{
		printf("File '%s' found.\nLooking for the centre of mass from file...\n", centrefile);
		
		FILE *finwe = fopen (centrefile, "rb");												/*Open the file with the centres*/					
		if (finwe == NULL)
		{
			printf("File '%s' cannot be opened\n\n", centrefile);
			exit(1);
		}
		while(!feof(finwe))																	/*Reading the lines from the centre file file*/
		{
			n=1;
			if (n != fread(&id_buff, sizeof(long int), n, finwe))							/*Read the id*/
			{
				printf("The format of id in '%s' is not compatible\n\n", centrefile);
				break;
			}
			n=3;
			if (n != fread(displ_cm, sizeof(double), n, finwe))
			{
				printf("The format of positions in '%s' is not compatible\n\n", centrefile);
				break;
			}
			if (n != fread(displ_vcm, sizeof(double), n, finwe))
			{
				printf("The format of velocities in '%s' is not compatible\n\n", centrefile);
				break;
			}
			if (id_buff==atoi(id))
			{
				printf("Centre found in list: %lf %lf %lf\n", displ_cm[0], displ_cm[1], displ_cm[2]);
				flag=1;
				break;
			}
		}
	}

//------------------------------------------------------------------------------------------Finding the rough preliminary centre of mass--------

	if ((access(centrefile, R_OK)==0) && (flag==1))											/*If the displacements have been read in the file*/
	{
#pragma omp parallel for collapse(2)
		for (i=0; i<header.nsph; i++)														/*Shifting the gas particles to the centre of mass reference frame*/
		{
			for (j=0; j<DIM; j++)
			{
				gas_particles[i].pos[j] -= displ_cm[j]/in_to_kpc;
				gas_particles[i].vel[j] -= displ_vcm[j]/in_to_kmpersec;
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.ndark; i++)														/*Shifting the dark matter particles to the centre of mass reference frame*/																					
		{
			for (j=0; j<DIM; j++)
			{
				dark_particles[i].pos[j] -= displ_cm[j]/in_to_kpc;
				dark_particles[i].vel[j] -= displ_vcm[j]/in_to_kmpersec;
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.nstar; i++)														/*Shifting the star particles to the centre of mass reference frame*/
		{
			for (j=0; j<DIM; j++)
			{
				star_particles[i].pos[j] -= displ_cm[j]/in_to_kpc;
				star_particles[i].vel[j] -= displ_vcm[j]/in_to_kmpersec;
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.nsink; i++)														/*Shifting the sink particles to the centre of mass reference frame*/
		{
			for (j=0; j<DIM; j++)
			{
				sink_particles[i].pos[j] -= displ_cm[j]/in_to_kpc;
				sink_particles[i].vel[j] -= displ_vcm[j]/in_to_kmpersec;
			}
		}
		for (i=0; i<DIM; i++) cm[i] = 0.;													/*Setting the centre of mass in the origin*/
	}
	else 																					/*If there's no centre file or the centre wasn't found in it*/
	{
		if (access(centrefile, R_OK)==0) printf("\nNo centre of mass found in file '%s' for snapshot %s\nCalculating it...\n", centrefile, input);

		mtot = 0.;
		for (i=0; i<DIM; i++) draft[i] = 0.;
		for (i=0; i<DIM; i++) cm[i] = 0.;
		for (i=0; i<header.nsph; i++)
		{
			for (j=0; j<DIM; j++) draft[j] += gas_particles[i].pos[j]*gas_particles[i].mass;
			mtot += gas_particles[i].mass;
		}
		for (i=0; i<header.nstar; i++)
		{
			for (j=0; j<DIM; j++) draft[j] += star_particles[i].pos[j]*star_particles[i].mass;
			mtot += star_particles[i].mass;
		}
		for (i=0; i<header.nsink; i++)
		{
			for (j=0; j<DIM; j++) draft[j] += sink_particles[i].pos[j]*sink_particles[i].mass;
			mtot += sink_particles[i].mass;
		}	
	
		for (i=0; i<DIM; i++) draft[i] /= mtot;      	                                    /*Calculating the centre of mass position through the average*/

//------------------------------------------------------------------------------------------Using the minimum in the potential as draft to find the centre of mass
		
		tempdist=1.e10;																		/*Setting the temporary distance to an arbitrary high value*/
		i=0;
		do 																					/*Evaluate the mass centre around the dynamical centre*/
		{
			i++;
			radcm = ran/i;
			rcm(draft, radcm, cm, vcm);														/*Getting the centre of mass recursively around the dynamic centre*/
		} while (criterion(draft, cm, &tempdist, radcm, floor, 100));

#ifdef DEBUG 																				/*Do it only if the key DEBUG is defined*/
		printf("Position of the preliminar centre (x, y, z): (%f, %f, %f)\n", draft[0], draft[1], draft[2]);
		printf("Position of the centre of mass (x, y, z): (%f, %f, %f)\n", cm[0], cm[1], cm[2]);
#endif

//------------------------------------------------------------------------------------------Shifting to the centre of mass reference------------

#pragma omp parallel for collapse(2)
		for (i=0; i<header.nsph; i++)
		{
			for (j=0; j<DIM; j++)
			{
				gas_particles[i].pos[j] -= cm[j];
				gas_particles[i].vel[j] -= vcm[j];
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.ndark; i++)																									
		{
			for (j=0; j<DIM; j++)
			{
				dark_particles[i].pos[j] -= cm[j];
				dark_particles[i].vel[j] -= vcm[j];
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.nstar; i++)														/*Shifting the star particles to the centre of mass reference frame*/
		{
			for (j=0; j<DIM; j++)
			{
				star_particles[i].pos[j] -= cm[j];
				star_particles[i].vel[j] -= vcm[j];
			}
		}
#pragma omp parallel for collapse(2)	
		for (i=0; i<header.nsink; i++)														/*Shifting the sink particles to the centre of mass reference frame*/
		{
			for (j=0; j<DIM; j++)
			{
				sink_particles[i].pos[j] -= cm[j];
				sink_particles[i].vel[j] -= vcm[j];
			}
		}
		for (i=0; i<DIM; i++) cm[i] -= cm[i];												/*Setting the new centre of mass in the origin*/
	}

//------------------------------------------------------------------------------------------Finding the half-mass radius------------------------

	TotStarMass=0.;
	for (i=0; i<header.nstar; i++) 
		TotStarMass += star_particles[i].mass;
	
	hmr = 2.8*gsoft(header.time);
	StarMass=0.; 
	while (StarMass < 0.5*TotStarMass)
	{
		hmr += 0.7*gsoft(header.time);
		StarMass=0.;
		for (i=0; i<header.nstar; i++) 
			if (dist_f2(star_particles[i].pos, cm)<=hmr)
				StarMass += star_particles[i].mass;
	}
	if (StarMass==TotStarMass)
		hmr = 2.8*gsoft(header.time);

#ifdef DEBUG
	printf("Hmr = %e kpc\n", hmr*in_to_kpc);
#endif

//------------------------------------------------------------------------------------------Setting the remaining parameters--------------------

		rad = 4*hmr;																		/*Computing the radius of the galaxy area to analyse: four times the half mass radius*/
		zmax = max(2.8*gsoft(header.time), 0.18*hmr);										/*Computing the bin half height as 0.3/1.66*/
		ran = 3*hmr;																		/*Re-define ran*/
#ifdef RROT
		zin = max(2.8*gsoft(header.time), 0.2*hmr);											/*Updating the bin half height for the rotation and translation routines*/		
#endif

		req = max(2.8*gsoft(header.time), 0.5*hmr);											/*Maximum distance from the centre imposed to the bar inner edge in pkpc*/
		maxdist_prel = 1.5*hmr;																/*Maximum distance from the centre allowed for a preliminar peak*/
#ifdef LIM
		lim[0] = 1.4*gsoft(header.time);													/*Setting the minimum bar length*/
		lim[1] = rad;																		/*The limit is used to exclude a merger and to select the stellar disc only*/
		lim[2] = 0.1;																		/*Limits on the A2 profile, min*/
		lim[3] = 1.e5;																		/*...and max*/
#endif	

//------------------------------------------------------------------------------------------Evaluating the total angular momentum---------------

	L_x = L_y = L_z = 0.;
#pragma omp parallel for private(i) reduction(+:L_x), reduction(+:L_y), reduction(+:L_z)														
	for (i=0; i<header.nstar; i++)															/*Evaluating the total angular momentum of the star particles*/
	{
		if (dist(star_particles[i].pos, cm)<(ran*ran))
		{
			L_x += star_particles[i].mass*(star_particles[i].pos[1]*star_particles[i].vel[2]-star_particles[i].pos[2]*star_particles[i].vel[1]);
			L_y += star_particles[i].mass*(star_particles[i].pos[2]*star_particles[i].vel[0]-star_particles[i].pos[0]*star_particles[i].vel[2]);
			L_z += star_particles[i].mass*(star_particles[i].pos[0]*star_particles[i].vel[1]-star_particles[i].pos[1]*star_particles[i].vel[0]);
		}
	}

#ifdef DEBUG 																				/*Do it only if the key DEBUG is defined*/
	printf("\nOriginal angular momentum (Lx, Ly, Lz): (%e, %e, %e)\n", L_x, L_y, L_z);
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------

	sinf = L_y/sqrt((L_x*L_x)+(L_y*L_y));													/*Evaluating the projection of L on the reference frame axis*/
	cosf = L_x/sqrt((L_x*L_x)+(L_y*L_y));
	sint = sqrt((L_x*L_x)+(L_y*L_y))/sqrt((L_x*L_x)+(L_y*L_y)+(L_z*L_z));
	cost = L_z/sqrt((L_x*L_x)+(L_y*L_y)+(L_z*L_z));

//----------------------------------------------------------------------------------------------------------------------------------------------

#if defined MASS || defined DEBUG || defined PSNAP											/*Do it only if the key MASS or the key DEBUG is defined*/
#pragma omp parallel for private (i)
	for (i=0; i<header.nsph; i++)															/*Performing the inverse rotation on gas particles*/
	{	
		Drotate (gas_particles[i].pos, sinf, cosf, sint, cost);								/*Rotating the positions of the gas particles*/
		rotate (gas_particles[i].vel, sinf, cosf, sint, cost);								/*Rotating the velocities of the gas particles*/
	}
#pragma omp parallel for private (i)		
	for (i=0; i<header.ndark; i++)															/*Performing the inverse rotation on dark matter particles*/
	{	
		Drotate (dark_particles[i].pos, sinf, cosf, sint, cost);							/*Rotating the positions of the dark matter particles*/
		rotate (dark_particles[i].vel, sinf, cosf, sint, cost);								/*Rotating the velocities of the dark matter particles*/
	}
#pragma omp parallel for private (i)		
	for (i=0; i<header.nsink; i++)															/*Performing the inverse rotation on dark matter particles*/
	{	
		Drotate (sink_particles[i].pos, sinf, cosf, sint, cost);							/*Rotating the positions of the sink particles*/
		rotate (sink_particles[i].vel, sinf, cosf, sint, cost);								/*Rotating the velocities of the sink particles*/
	}	
#endif
#pragma omp parallel for private (i)		
	for (i=0; i<header.nstar; i++)															/*Performing the inverse rotation on star particles*/
	{	
		Drotate (star_particles[i].pos, sinf, cosf, sint, cost);							/*Rotating the positions of the star particles*/
		rotate (star_particles[i].vel, sinf, cosf, sint, cost);								/*Rotating the velocities of the star particles*/
	}

	Krot=K=0.;
#pragma omp parallel for private(i) reduction(+:Krot), reduction(+:K)
	for (i=0; i<header.nstar; i++)
	{
		if (dist(star_particles[i].pos, cm)<(ran*ran))
		{	
			vphys[0] = (sqrt(header.time)*star_particles[i].vel[0])+(0.1*header.time*star_particles[i].pos[0]);
			vphys[1] = (sqrt(header.time)*star_particles[i].vel[1])+(0.1*header.time*star_particles[i].pos[1]);
			vphys[2] = (sqrt(header.time)*star_particles[i].vel[2])+(0.1*header.time*star_particles[i].pos[2]);

			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);
			vtan = (vphys[1]*star_particles[i].pos[0]-vphys[0]*star_particles[i].pos[1])/rxy;
			Krot += star_particles[i].mass*(vtan*vtan);
			K += star_particles[i].mass*((vphys[0]*vphys[0])+(vphys[1]*vphys[1])+(vphys[2]*vphys[2]));
		}
	}
#ifdef DEBUG
	printf("Krot = %f\n", Krot/K);
#endif

#ifdef RROT																					/*Do it only if the key RROT is defined*/

	L_x = L_y = L_z = 0.;																	/*Resetting the angular momentum components*/
#pragma omp parallel for private(i) reduction(+:L_x), reduction(+:L_y), reduction(+:L_z)
	for (i=0; i<header.nstar; i++)															/*Evaluating again the total angular momentum of the star particles*/
	{
		if ((dist_f2(star_particles[i].pos, cm)<ran) && (star_particles[i].pos[2]*star_particles[i].pos[2]<zin*zin))
		{
			L_x += star_particles[i].mass*(star_particles[i].pos[1]*star_particles[i].vel[2]-star_particles[i].pos[2]*star_particles[i].vel[1]);
			L_y += star_particles[i].mass*(star_particles[i].pos[2]*star_particles[i].vel[0]-star_particles[i].pos[0]*star_particles[i].vel[2]);
			L_z += star_particles[i].mass*(star_particles[i].pos[0]*star_particles[i].vel[1]-star_particles[i].pos[1]*star_particles[i].vel[0]);
		}
	}

#ifdef DEBUG 																				/*Do it only if the key DEBUG is defined*/
	printf("RROT defined. Angular momentum (Lx, Ly, Lz) of the rotated system: (%e, %e, %e)\n", L_x, L_y, L_z);
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------

	sinf = L_y/sqrt((L_x*L_x)+(L_y*L_y));
	cosf = L_x/sqrt((L_x*L_x)+(L_y*L_y));
	sint = sqrt((L_x*L_x)+(L_y*L_y))/sqrt((L_x*L_x)+(L_y*L_y)+(L_z*L_z));
	cost = L_z/sqrt((L_x*L_x)+(L_y*L_y)+(L_z*L_z));

//----------------------------------------------------------------------------------------------------------------------------------------------

#if defined MASS || defined DEBUG || defined PSNAP											/*Do it only if the key MASS or the key DEBUG is defined*/
#pragma omp parallel for private (i)
	for (i=0; i<header.nsph; i++)
	{	
		Drotate (gas_particles[i].pos, sinf, cosf, sint, cost);
		rotate (gas_particles[i].vel, sinf, cosf, sint, cost);
	}	
#pragma omp parallel for private (i)	
	for (i=0; i<header.ndark; i++)
	{	
		Drotate (dark_particles[i].pos, sinf, cosf, sint, cost);
		rotate (dark_particles[i].vel, sinf, cosf, sint, cost);
	}
#pragma omp parallel for private (i)	
	for (i=0; i<header.nsink; i++)
	{	
		Drotate (sink_particles[i].pos, sinf, cosf, sint, cost);
		rotate (sink_particles[i].vel, sinf, cosf, sint, cost);
	}
#endif
#pragma omp parallel for private (i)
	for (i=0; i<header.nstar; i++)
	{	
		Drotate (star_particles[i].pos, sinf, cosf, sint, cost);
		rotate (star_particles[i].vel, sinf, cosf, sint, cost);
	}

#endif

	Krot=K=0.;
#pragma omp parallel for private(i) reduction(+:Krot), reduction(+:K)
	for (i=0; i<header.nstar; i++)
	{
		if (dist(star_particles[i].pos, cm)<(ran*ran))
		{	
			vphys[0] = (sqrt(header.time)*star_particles[i].vel[0])+(0.1*header.time*star_particles[i].pos[0]);
			vphys[1] = (sqrt(header.time)*star_particles[i].vel[1])+(0.1*header.time*star_particles[i].pos[1]);
			vphys[2] = (sqrt(header.time)*star_particles[i].vel[2])+(0.1*header.time*star_particles[i].pos[2]);

			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);
			vtan = (vphys[1]*star_particles[i].pos[0]-vphys[0]*star_particles[i].pos[1])/rxy;
			Krot += star_particles[i].mass*(vtan*vtan);
			K += star_particles[i].mass*((vphys[0]*vphys[0])+(vphys[1]*vphys[1])+(vphys[2]*vphys[2]));
		}
	}
#ifdef DEBUG
	printf("Krot = %f\n", Krot/K);
#endif

#if defined PSNAP || defined DEBUG
	if (thdf52tipsy(outdeb, gas_particles, dark_particles, star_particles, sink_particles)==TRUE)	/*Writing the rotated snapshot to debug*/
	{
		printf("\nFile '%s' successfully written\n", outdeb);
	}
#endif

//if (Krot<Krot_min)
//{
//	printf("Krot = % f is too low. Rotation likely failed. Code aborted.\n", Krot);
//	return 0;
//}
//else printf("System centered and rotated\n");

#ifdef DEBUG 																				/*Do it only if the key DEBUG is defined*/
	clock_gettime(CLOCK_MONOTONIC, &t_stop);												/*Translating and rotating time elapsed*/
	printf("(Elapsed: %.2f seconds)\n\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
	t_par = t_stop;
#endif

//------------------------------------------------------------------------------------------Performing the fourier analysis---------------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------

	for (i=0; i<nbin_max; i++) for (j=1; j<7; j++) bin[j][i] = 0.;

//------------------------------------------------------------------------------------------Fix the number of particles in each bin-------------

#ifdef EQNPART

	selpart = 0;
	for (i=0; i<header.nstar; i++)	if (((r = dist_f2(star_particles[i].pos, cm))<rad) && (r>0) && (star_particles[i].pos[2]*star_particles[i].pos[2]<zmax*zmax)) selpart ++; /*Counting the selected stellar particles*/
	allmat(&list, selpart, 2);																/*Allocatng the matrix to contain, for each selected particles, the id and the bidimentional distance from the centre of the galaxy*/
	j=0;
	for (i=0; i<header.nstar; i++)															/*Evaluating the differential quantities*/
	{
		if (((r = dist_f2(star_particles[i].pos, cm))<rad) && (r>0) && (star_particles[i].pos[2]*star_particles[i].pos[2]<zmax*zmax)) /*Isolating the central disc*/
		{
			list[j][0] = i;																	/*Getting the id of each particle in the disc*/
			list[j][1] = r;																	/*Getting also the corresponding 2-dimensional distance from the centre of the galaxy*/
			j++;
		}
	}
	qsort(list, selpart, sizeof(list[0]), compare);											/*Sorting the stellar particles from the nearest to the farest from the galactic centre*/
	npart = (int)(selpart/nbin_max);														/*Evaluating the number of particle for each bin*/
	if (selpart%npart!=0) npart++;															/*If there are more particles than npart*nbin_max, add one particle to each bin*/

	if (npart<min_npart) npart = min_npart;													/*Set a lower limit of particle in each bin*/

	for (i=0; i<selpart; i++)
	{
		if (i%npart==0) firstP = list[i][1];												/*Save the distance of the first particle in the bin*/
		if ((i%npart==npart-1) || (i==selpart-1)) bin[0][(int)(i/npart)] = 0.5*(list[i][1]+firstP); /*Save the distance of either the last particle in the bin or the last particle in the disc; then evaluate the mean distance of the bin*/

		bin[1][(int)(i/npart)] += star_particles[(int)(list[i][0])].mass;					/*Calculate the integral mass for each bin*/
		bin[2][(int)(i/npart)] += (2*((star_particles[(int)(list[i][0])].pos[0]*star_particles[(int)(list[i][0])].pos[0])/(list[i][1]*list[i][1]))-1)*star_particles[(int)(list[i][0])].mass; /*Evaluate the integral of m_j*cos(2*th) for each bin, where j is j-th particle*/
		bin[3][(int)(i/npart)] += (2*((star_particles[(int)(list[i][0])].pos[0]*star_particles[(int)(list[i][0])].pos[1])/(list[i][1]*list[i][1])))*star_particles[(int)(list[i][0])].mass; /*Evaluate the integral of m_j*sin(2*th) for each bin, where j is j-th particle*/
		bin[6][(int)(i/npart)] ++;															/*Counting the particle in each bin*/
	}

#ifdef DEBUG
	printf("\nSelected particles = %d\nParticle per bin = %d\n", selpart, npart);
	printf("Last populated bin with %d particles\n", (int)bin[6][nbin_max-(nbin_max-(int)(ceil(selpart/npart)))]);
	printf("Last %d bins are empty\n\n", (nbin_max-(int)(ceil(selpart/npart))));
#endif
	printf("%d bins used\n\n", (int)(ceil(selpart/npart)));

	nbin = (int)(ceil(selpart/npart));

#else

//------------------------------------------------------------------------------------------Fix the bin length----------------------------------		
	
	nbin = nbin_max;
	r = (rad/nbin);
	for (i=0; i<nbin; i++) bin[0][i] = (r*(0.5+i));											/*Evaluating the middle radius value for each bin*/

	for (i=0; i<header.nstar; i++)															/*Evaluating the differential quantities*/
	{
		if (((r = dist_f2(star_particles[i].pos, cm))<rad) && (r>0) && (star_particles[i].pos[2]*star_particles[i].pos[2]<zmax*zmax)) /*Isolating the central disc*/
		{
			bin[1][(int)((r/rad)*nbin)] += star_particles[i].mass;							/*Calculate the integral mass for each bin*/
			bin[2][(int)((r/rad)*nbin)] += (2*((star_particles[i].pos[0]*star_particles[i].pos[0])/(r*r))-1)*star_particles[i].mass; /*Evaluate the integral of m_j*cos(2*th) for each bin, where j is j-th particle*/
			bin[3][(int)((r/rad)*nbin)] += (2*((star_particles[i].pos[0]*star_particles[i].pos[1])/(r*r)))*star_particles[i].mass; /*Evaluate the integral of m_j*sin(2*th) for each bin, where j is j-th particle*/
			bin[6][(int)((r/rad)*nbin)] ++;													/*Counting the particle in each bin*/
		}
	}

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------	

#if defined A2PROF && !defined DEBUG														/*Don't write any output in debugging session and write the file 'A2_profile only if the key A2PROF is defined*/
	huor = fopen(A2out, "w");
	if (huor==NULL)																			/*Checking on the opening*/
	{
		printf("\nError in opening '%s' file\n\n", A2out);
		exit(1);
	}
	fprintf(huor, "Radius [kpc]\tA2\t\tA2tot\t\ttheta [deg]\ttheta_tot [deg]\n");			/*Print for each radius the profile of A2, its cumulative value and the angle*/
#endif

	mtot=sint=cost=maxA2=0.;																/*Resetting the variables*/

	for (i=0; i<nbin; i++)
	{
		mtot += bin[1][i];																	/*Compute the integral of the mass from 0 to x_i*/
		cost += bin[2][i];																	/*Compute the integral of m_j*cos(2*th) from 0 to x_i*/
		sint += bin[3][i];																	/*Compute the integral of m_j*sin(2*th) from 0 to x_i*/
		cosf = bin[2][i];
		sinf = bin[3][i];
		bin[2][i] = (sqrt(cosf*cosf+sinf*sinf))/bin[1][i];									/*Evaluate the module of the A2 coefficient for the differential quantities*/
		bin[3][i] = 0.5*atan2(sinf,cosf);													/*Angles in radians: since it returns an angle in [-PI/2, PI/2] (because atan2 returns in [-PI, PI] but it's divided for 2 because it comes from 2th) the case in which the bar is near the discontinuity point has to be considered. This angle is valid only in this reference frame. It cannot be computed any angle dispacement using these values*/
		bin[4][i] = (sqrt(cost*cost+sint*sint))/mtot;										/*Evaluate the module of the A2 coefficient for the integral quantities*/
		bin[5][i] = 0.5*atan2(sint,cost);													/*Angle of the integral overdensity*/

#if defined A2PROF && !defined DEBUG 														/*Don't write any output in debugging session and write the file 'A2_profile' only if the key A2PROF is defined*/
		fprintf(huor, "%lf\t%lf\t%lf\t%lf\t%lf\n", in_to_kpc*bin[0][i], bin[2][i], bin[4][i], bin[3][i]*(180./M_PI), bin[5][i]*(180./M_PI));
#endif
	}

#if defined A2PROF && !defined DEBUG
	if (fclose(huor)==0) printf("file '%s' successfully written\n", A2out);					/*Closing the file 'A2_profile'*/
	else
	{
		printf("\nCan't write on file '%s'\n\n", A2out);
		exit(1);
	}

#endif	

//------------------------------------------------------------------------------------------Finding the bar linked to the peaks-----------------

	temp_kern = pow(5, (nbin_max-nbin)/(nbin_max*0.75));									/*Function to evaluate the a2 smoothing kernel*/
	round_temp_kern = round(temp_kern);
	if ((round_temp_kern%2)!=0)
		KERN_a2t = round_temp_kern;
	else if ((temp_kern-round_temp_kern)>0)
		KERN_a2t = round_temp_kern+1;
	else
		KERN_a2t = round_temp_kern-1;

	KERN_a2 = KERN_a2t+2;
	KERN_ph = KERN_a2t+2;

#ifdef DEBUG
	printf("KERN_a2t = %d\nKERN_a2 = %d\nKERN_ph = %d\n\n", KERN_a2t, KERN_a2, KERN_ph);
#endif

	der_v(bin[4], bin[0], deriv, nbin);														/*Derivative of the A2 cumulative profile*/
	if (!smooth(deriv, smooth_der, nbin, KERN_a2t))											/*Smoothing out the derivative of the A2 cumulative profile*/
		return 0;
	if (!smooth_phase(bin[3], smooth_ph, nbin, KERN_ph))									/*Smoothing out the punctual phase in order to select a better phase-reference (The phase is now in degrees)*/
		return 0;
	if (!smooth(bin[2], smooth_A2, nbin, KERN_a2))
		return 0;		

#ifdef DEBUG
	huor = fopen(A2smooth, "w");
	if (huor==NULL)																			/*Checking on the opening*/
	{
		printf("\nError in opening '%s' file\n\n", A2smooth);
		exit(1);
	}
	fprintf(huor, "Radius [kpc]\tA2\t\tA2_smooth\tA2tot\t\tder\t\tder_smooth\tphase\t\tphase_smooth\n");
	for(i=0; i<nbin; i++) fprintf(huor, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", in_to_kpc*bin[0][i], bin[2][i], smooth_A2[i], bin[4][i], deriv[i], smooth_der[i], bin[3][i]*(180./M_PI), smooth_ph[i]*(180./M_PI));
	if (fclose(huor)==0) printf("file '%s' successfully written\n", A2smooth);
	else
	{
		printf("\nCan't write on file '%s'\n\n", A2smooth);
		exit(1);
	}
#endif

	i=j=0;
	flag=0;
	npeaks=0;																				/*Resetting the peak counter*/
	while (i<nbin)	 																		/*Do it for every radial bin*/
	{
		if (nmaxp==npeaks)
		{
			printf("\nToo many peaks in galaxy %s. Please, increment the buffer size 'nmaxp'\n\n", file_name);
			exit(0);
		}
		if (smooth_der[i]<0)																/*If the smoothA2 derivative is negative...*/
		{
			if (flag==1)																	/*...but, if it was positive, a peak is found:*/
			{
				peaks[npeaks] = i;
				npeaks++;
#ifdef DEBUG
				printf("peak n.%d found at %lf kpc with smoothA2 phase %lf\n", npeaks, in_to_kpc*bin[0][i], smooth_ph[i]*(180/M_PI));
#endif
			}
			flag=-1; 																		/*...set the flag to remember that, now, it's negative*/
		}
		else flag=1;
		i++;
	}

//----------------------------------------------------------------------------------------------------------------------------------------------	

	R_th[npeaks][1] = NA;																	/*Setting the fake first annulus of the last bar to an arbitrary large number*/

	for (i=0; i<npeaks; i++)
	{
		BAR=RANGE=0;																		/*Resetting the bar flag*/
		bar_first_bin[i]=nbin_max-1;														/*Setting the first bin of the bar to the end of the region*/
		Phi[i]=0;
		for (j=0; j<nbin; j++)																/*Getting the R_th*/
		{									
			if (fabs((fabs(smooth_ph[j])-ra))>=Dth)											/*If the annulus phase is far from the discontinuity -90°/+90°*/
			{				
				if (fabs(smooth_ph[peaks[i]]-smooth_ph[j])>=Dth)							/*If the phase displacement calculated is enough...*/
				{				
					if (((bin[0][j]-bin[0][bar_first_bin[i]])<Bmin) || (RANGE!=1)) 			/*...and if the bar minimum length has not been reached yet or the position of the A2 peak isn't part of the bar...*/
					{
						BAR=0; 																/*...reset the bar counter...*/
						Phi[i]=0.;															/*and the bar phase buffer*/
					}
					else break;																/*...otherwise break the cycle and save the extension*/
				}
				else 																		/*If the annulus phase is inside the tollerance...*/
				{
					if (BAR>0) 																/*and if the bar counter has already been updated one time*/
					{
						BAR++;																/*Increment the counter if a first annulus of the bar was already found*/
						Phi[i]+=bin[3][j];													/*Updating the bar phase*/
					}
					else if (bin[0][j]<=bin[0][peaks[i]])									/*If the bar counter has been updated not even once and the current position on the radius is BEFORE the peak position*/
					{
						BAR=1;																/*The first annulus of the bar must be inside the radius of the peak just found*/
						Phi[i]+=bin[3][j]; 													/*Update the bar phase also*/
						R_th[i][1]=bin[0][j];												/*Save the radius of the first bar annulus*/
						bar_first_bin[i] = j;												/*Save the first bin of the bar*/
					}
					if (peaks[i]==j) RANGE=1;
				}
			}
			else if (smooth_ph[j]>0)														/*This is to consider the discontinuity -90°/+90° when the annuls phase is positive*/
			{ 
				if ((fabs(smooth_ph[peaks[i]]-smooth_ph[j])>=Dth) && (smooth_ph[peaks[i]]-smooth_ph[j]+(2*ra)>=Dth))
				{
					if (((bin[0][j]-bin[0][bar_first_bin[i]])<Bmin) || (RANGE!=1))
					{
						BAR=0;
						Phi[i]=0.;
					}
					else break;
				}
				else
				{
					if (BAR>0)
					{
						BAR++;
						if (fabs(smooth_ph[j]-smooth_ph[peaks[i]])<fabs(smooth_ph[j]-smooth_ph[peaks[i]]-(2*ra))) Phi[i]+=bin[3][j];
						else Phi[i]+=bin[3][j]-(2*ra);
					}
					else if (bin[0][j]<=bin[0][peaks[i]])
					{
						BAR=1;
						if (fabs(smooth_ph[j]-smooth_ph[peaks[i]])<fabs(smooth_ph[j]-smooth_ph[peaks[i]]-(2*ra))) Phi[i]+=bin[3][j];
						else Phi[i]+=bin[3][j]-(2*ra);
						R_th[i][1]=bin[0][j];
						bar_first_bin[i] = j;			
					}
					if (peaks[i]==j) RANGE=1;
				}
			}
			else 																			/*This is to consider the discontinuity -90°/+90° when the annuls phase is negative*/
			{
				if ((fabs(smooth_ph[peaks[i]]-smooth_ph[j])>=Dth) && (smooth_ph[peaks[i]]-smooth_ph[j]-(2*ra)>=Dth))
				{ 
					if (((bin[0][j]-bin[0][bar_first_bin[i]])<Bmin) || (RANGE!=1))
					{
						BAR=0;
						Phi[i]=0.;
					}
					else break;
				}
				else
				{
					if (BAR>0)
					{
						BAR++;
						if (fabs(smooth_ph[j]-smooth_ph[peaks[i]])<fabs(smooth_ph[j]-smooth_ph[peaks[i]]+(2*ra))) Phi[i]+=bin[3][j];
						else Phi[i]+=bin[3][j]+(2*ra);
					}
					else if (bin[0][j]<=bin[0][peaks[i]])
					{
						BAR=1;
						if (fabs(smooth_ph[j]-smooth_ph[peaks[i]])<fabs(smooth_ph[j]-smooth_ph[peaks[i]]+(2*ra))) Phi[i]+=bin[3][j];
						else Phi[i]+=bin[3][j]+(2*ra);
						R_th[i][1]=bin[0][j];
						bar_first_bin[i] = j;												/*Save the first bin of the bar*/
					}
					if (peaks[i]==j) RANGE=1;
				}				
			}
		}
		if (j<nbin)
		{
			bar_last_bin[i] = j-1;
			R_th[i][0] = bin[0][j-1];
			Phi[i]/=BAR;
			if (Phi[i]<-ra) Phi[i]+=(2*ra);
			else if (Phi[i]>ra) Phi[i]-=(2*ra);
		}
		else R_th[i][0] = -1;
	}

//------------------------------------------------------------------------------------------Refine the properties of the bar--------------------

#ifdef REFINE

	for (i=0; i<npeaks; i++)
	{
		if (R_th[i][0] < 0) continue;														/*Avoid peaks that have been already rejected*/

		Phi[i] = avg_phase (smooth_ph, bar_first_bin[i], peaks[i]);							/*Evaluate the averaged phase between Rmin and Rpeak*/	

		famp = 0.;
		for (j=bar_first_bin[i]; j<=peaks[i]; j++)											/*Evaluate the stdv*/
		{
			tempval = fabs(Phi[i]-smooth_ph[j]);
			tempval = min(tempval, fabs(tempval-2*ra));
			famp += tempval*tempval;
		}
		famp = sqrt(famp/(peaks[i]-bar_first_bin[i]+1));									/*Get the amplitude of the fluctuations*/
#ifdef DEBUG		
		printf("Delta Phi = %lf\n", famp*(180./M_PI));
#endif

//----------------------------------------------------------------

		bar_last_bin[i] = peaks[i];
		for (j=peaks[i]+1; j<nbin; j++)														/*From the peak, toward the outskirt*/
		{
			tempval = fabs(Phi[i]-smooth_ph[j]);
			if (min(tempval, fabs(tempval-2*ra))<mi*famp)
				bar_last_bin[i] = j;														/*Update the bin of the outer bar limit*/
			else
				break;
		}
		R_th[i][0] = bin[0][bar_last_bin[i]];												/*Save the length of the bar*/

//----------------------------------------------------------------

		Phi[i] = avg_phase (smooth_ph, bar_first_bin[i], bar_last_bin[i]);					/*Evaluate the averaged phase between Rmin and Rmax*/

		famp = 0.;
		for (j=bar_first_bin[i]; j<=bar_last_bin[i]; j++)									/*Re-evaluate the stdv*/
		{
			tempval = fabs(Phi[i]-smooth_ph[j]);
			tempval = min(tempval, fabs(tempval-2*ra));
			famp += tempval*tempval;
		}
		famp = sqrt(famp/(bar_last_bin[i]-bar_first_bin[i]+1));								/*Get the amplitude of the fluctuations*/
#ifdef DEBUG		
		printf("Delta Phi = %lf\n", famp*(180./M_PI));
#endif

//----------------------------------------------------------------

		for (j=peaks[i]-1; j>=0; j--)														/*From the peak, toward the nucleus*/
		{
			tempval = fabs(Phi[i]-smooth_ph[j]);
			if (min(tempval, fabs(tempval-2*ra))<mi*famp)
				bar_first_bin[i] = j;														/*Update the bin of the inner bar limit*/
			else
				break;
		}
		R_th[i][1] = bin[0][bar_first_bin[i]];												/*Save the bar first bin position*/

//----------------------------------------------------------------

		Phi[i] = avg_phase (smooth_ph, bar_first_bin[i], bar_last_bin[i]);					/*Evaluate the averaged phase between Rmin and Rmax*/	
	}

#endif

#ifdef DEBUG
	clock_gettime(CLOCK_MONOTONIC, &t_stop);
	printf("(Elapsed: %.2f seconds)\n\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
	t_par = t_stop;
#endif

//------------------------------------------------------------------------------------------Apply the constraints to select the bar-------------

	bar_param[0] = NA;
	bar_param[1] = -1./in_to_kpc;															/*Set -1/in_to_kpc to print exactly -1*/
	bar_param[2] = -1.;
	bar_param[3] = -1.;
	bar_param[4] = -1./(180./M_PI);
	bar_param[5] = -1./in_to_kpc;


	j=0;
	for (i=0; i<npeaks; i++)
	{
	#ifdef LIM
		if ((R_th[i][0]>-1.) && 
		//(R_th[i+1][1]>R_th[i][0]) && 
			(R_th[i][1]<=req) && 
			(bin[0][peaks[i]]<maxdist_prel) && 
			(R_th[i][0]>=lim[0]) && 
			(R_th[i][0]<=lim[1]) && 
			(smooth_A2[peaks[i]]>=lim[2]) && 
			(smooth_A2[peaks[i]]<=lim[3])
			)
	#else
		if ((R_th[i][0]>-1.) 
			//&& //(R_th[i+1][1]>R_th[i][0]) && 
			//(R_th[i][1]<=req) && 
			//(bin[0][peaks[i]]<maxdist_prel)
			)
	#endif
		{
			maxA2 = 0.;																		/*Resetting the maximum of A2*/
			for (j=bar_first_bin[i]; j<=bar_last_bin[i]; j++)
				if (maxA2<smooth_A2[j]) maxA2 = smooth_A2[j]; 								/*Get the maximum of A2 in that range as the A2max of the current bar overdensity*/

			//if (bar_param[3]<maxA2)														/*This is to select the peak with the maximum of A2*/
			if (bar_param[0]>bin[0][peaks[i]])												/*This is to select the nearest peak to the centre*/
			{
				bar_param[0] = bin[0][peaks[i]];
				bar_param[1] = R_th[i][0];
				bar_param[2] = bin[4][peaks[i]];
				bar_param[3] = maxA2;
				bar_param[4] = Phi[i];
				bar_param[5] = R_th[i][1];													/*Bar minimum radius*/
				j=1;
			}
		}
	}																						/*End of the loop on the peaks*/
	if (bar_param[0]==NA) bar_param[0] = -1./in_to_kpc;

//------------------------------------------------------------------------------------------Evaluate the sigma ratio within the bar-------------	

#ifdef REFINE

	if (bar_param[0]>0)
	{
			A = bar_param[4];																/*Bar phase set to Rphi*/
			rm = bar_param[1];																/*Bar (half-)length*/
			wm = 2.8*gsoft(header.time);													/*Setting the bar full width*/
			hm = 2.8*gsoft(header.time);													/*Setting the bar full height*/
			mr = bar_param[5];																/*Minimum bar extent*/

		C1 = rm*rm;
		C2 = wm*wm*0.25;
		C3 = hm*hm*0.25;

		MeanV_z=MeanV_R=Sigma_z=Sigma_R = 0.;
		barUP=barDOWN=0.;
		N_bar_part = 0;

		#pragma omp parallel for private (i, vphys, rxy), reduction(+:MeanV_R), reduction(+:MeanV_z), reduction(+:N_bar_part), reduction(+:barUP), reduction(+:barDOWN)
		for (i=0; i<header.nstar; i++)
		{
			vphys[0] = (sqrt(header.time)*star_particles[i].vel[0])+(0.1*header.time*star_particles[i].pos[0]);
			vphys[1] = (sqrt(header.time)*star_particles[i].vel[1])+(0.1*header.time*star_particles[i].pos[1]);
			vphys[2] = (sqrt(header.time)*star_particles[i].vel[2])+(0.1*header.time*star_particles[i].pos[2]);

			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);	/*Getting for each particle the module of the vector projection on the xy plane*/
					
			if ((rxy>=mr) && 
				((star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))*(star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))<C1) && 
				((-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))*(-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))<C2) && 
				(star_particles[i].pos[2]*star_particles[i].pos[2]<C3))
			{
				MeanV_R += (vphys[0]*star_particles[i].pos[0]+vphys[1]*star_particles[i].pos[1])/rxy;
				MeanV_z += vphys[2];
				N_bar_part++;

				if (star_particles[i].pos[0]*cos(A)>-star_particles[i].pos[1]*sin(A)) barUP += star_particles[i].mass;
				else if (star_particles[i].pos[0]*cos(A)<-star_particles[i].pos[1]*sin(A)) barDOWN += star_particles[i].mass;				
			}
		}

		MeanV_R = MeanV_R/(double)N_bar_part;												/*Getting the mean radial velocity*/	
		MeanV_z = MeanV_z/(double)N_bar_part;												/*Getting the mean vertical velocity*/

#pragma omp parallel for private (i, vphys, rxy), reduction(+:Sigma_R), reduction(+:Sigma_z), reduction(+:N_bar_part)
		for (i=0; i<header.nstar; i++)
		{
			vphys[0] = (sqrt(header.time)*star_particles[i].vel[0])+(0.1*header.time*star_particles[i].pos[0]);
			vphys[1] = (sqrt(header.time)*star_particles[i].vel[1])+(0.1*header.time*star_particles[i].pos[1]);
			vphys[2] = (sqrt(header.time)*star_particles[i].vel[2])+(0.1*header.time*star_particles[i].pos[2]);

			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);

			if ((rxy>=mr) && ((star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))*(star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))<C1) && ((-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))*(-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))<C2) && (star_particles[i].pos[2]*star_particles[i].pos[2]<C3))
			{
				Sigma_R += (MeanV_R-(vphys[0]*star_particles[i].pos[0]+vphys[1]*star_particles[i].pos[1])/rxy)*(MeanV_R-(vphys[0]*star_particles[i].pos[0]+vphys[1]*star_particles[i].pos[1])/rxy);
				Sigma_z += (MeanV_z-vphys[2])*(MeanV_z-vphys[2]);	
			}
		}

		Sigma_R = sqrt(Sigma_R/(double)N_bar_part);
		Sigma_z = sqrt(Sigma_z/(double)N_bar_part);	

#ifdef DEBUG
		printf("Sigma_z/Sigma_R = %f\n", id, Sigma_z/Sigma_R);
		printf("barUP/barDOWN = %f\n", id, barUP/barDOWN);
#endif

	}
	else
	{
		Sigma_R = 1;
		Sigma_z = -1;

		barUP = -1.;
		barDOWN = 1;
	}

#endif

#ifdef DEBUG 																				/*Do it only if the key DEBUG is defined*/
	clock_gettime(CLOCK_MONOTONIC, &t_stop);												/*Fourier analysis time elapsed*/
	printf("Fourier analysis completed\n(Elapsed: %.2f seconds)\n\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
	t_par = t_stop;
#endif

//------------------------------------------------------------------------------------------Evaluating the Boxy-peanut shape--------------------	
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------

#ifdef BOXY

	if (bar_param[0]>0)
	{
			wb = 1./in_to_kpc;																/*Setting the bar full-width*/
			//rb = wid*bar_param[0];														/*Getting the bar full-length from the estimate of the position of A2tot*/
			rb = hmr/in_to_kpc;																/*Setting the maximum distance where to compute the median density*/
			hb = 1./in_to_kpc;																/*Setting the bar full-height*/

		xbin = xbin_mult*nbin;
		zbin = zbin_mult*nbin;
		allmat(&density, zbin, xbin);														/*Allocating the memory for the matrix of the density map: zbin x xbin*/
		median_d = (double *)malloc(xbin*sizeof(double));									/*Allocating the memory for the median lower profile*/
		median_u = (double *)malloc(xbin*sizeof(double));									/*Allocating the memory for the median upper profile*/

			C1 = (rb*cos(bar_param[4])+rb*tan(bar_param[4])*sin(bar_param[4]))*(rb*cos(bar_param[4])+rb*tan(bar_param[4])*sin(bar_param[4]));
			C2 = (wb*0.5*(sin(bar_param[4])*tan(bar_param[4])+cos(bar_param[4])))*(wb*0.5*(sin(bar_param[4])*tan(bar_param[4])+cos(bar_param[4])));
			C3 = (hb*hb*0.25);

#pragma omp parallel for collapse(2)
			for (j=0; j<xbin; j++)
			 	for (i=0; i<zbin; i++)  
			 		density[i][j] = 0.;														/*Resetting the density matrix*/

#pragma omp parallel for private(i) reduction(+:density[:zbin][0:xbin])
			for (i=0; i<header.nstar; i++)
				if (((star_particles[i].pos[0]+star_particles[i].pos[1]*tan(bar_param[4]))*(star_particles[i].pos[0]+star_particles[i].pos[1]*tan(bar_param[4]))<C1) && ((star_particles[i].pos[1]-star_particles[i].pos[0]*tan(bar_param[4]))*(star_particles[i].pos[1]-star_particles[i].pos[0]*tan(bar_param[4]))<C2) && ((star_particles[i].pos[2]*star_particles[i].pos[2])<C3) && (sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1])>0)) /*Isolating the bar*/
					density[(int)((((star_particles[i].pos[2])+hb)/hb)*0.5*zbin)][(int)((((star_particles[i].pos[0]*cos(bar_param[4])+star_particles[i].pos[1]*sin(bar_param[4]))+rb)/rb)*0.5*xbin)] += star_particles[i].mass;	/*Calculate the integral mass for each bin*/

			delta = 0.;																		/*Resetting the 'delta' value*/
#pragma omp parallel for private (j, half, i, sum), reduction(+:delta)
			for (j=0; j<xbin; j++)
			{
				half=0.;																	/*Resetting the total value*/
				for (i=zbin*0.5; i<zbin; i++) half += density[i][j];						/*Computing the integral of the mass in each column for the upper space*/
				half *= bmult;																/*Select the value of the mass distribution to show*/

				i=zbin*0.5;
				sum=0.;																		/*Resetting the partial integral*/
				while ((sum<half) && (i<zbin))												/*Evaluating the median in the upper space starting from the centre*/
				{
					sum += density[i][j];
					i++;
				}
				median_u[j] = (i-zbin*0.5)*(hb/zbin);										/*Getting the position of the selected bin*/

				half=0.;																	/*Resetting the total value*/
				for (i=zbin*0.5-(1-zbin%2); i>=0; i--) half += density[i][j];				/*Computing the integral of the mass in each column for the lower space*/
				half *= bmult;																/*Select the value of the mass distribution to show*/
						
				i=zbin*0.5-(1-zbin%2);
				sum=0.;																		/*Resetting the partial integral*/
				while ((sum<half) && (i>=0))												/*Evaluating the median in the lower space starting from the centre*/
				{
					sum += density[i][j];
					i--;
				}
				median_d[j] = (i-(zbin*0.5-(1-zbin%2)))*(hb/zbin);							/*Getting the position of the selected bin*/
				delta += 2*(median_u[j]+median_d[j])/(median_u[j]-median_d[j]);				/*Evaluating the mean value of the median z over the investigated range*/ 
			}

			delta /= xbin;

#ifndef DEBUG																				/*Don't write any output in debugging session*/
		galadriel = fopen(boxyout, "w");
		if (galadriel==NULL)																/*Checking on the opening*/
		{
			printf("\nError in opening '%s' file\n\n", boxyout);
			exit(1);
		}
		fprintf(galadriel, "x\t\tzmed_d\t\tzmed_u\t\tdelta\n");								/*Print for each snapshot the boxynes profile: the x coordinate and the median of the z-density (both for the lower and for the upper zones)*/
		for (i=0; i<xbin; i++) fprintf(galadriel, "%e\t%e\t%e\t%lf\n", (i-xbin*0.5)*(r/xbin)*in_to_kpc, median_u[i]*in_to_kpc, median_d[i]*in_to_kpc, 2*(median_u[i]+median_d[i])/(median_u[i]-median_d[i]));
		if (fclose(galadriel)==0) printf("file '%s' successfully written\n", boxyout);
		else
		{
			printf("\nCan't write on file '%s'\n\n", boxyout);								/*Closing the file 'boxy_profile..'*/
			exit(1);
		}
#else
		clock_gettime(CLOCK_MONOTONIC, &t_stop);											/*Boxyness calculation time elapsed*/
		printf("Boxyness analysis completed\n(Elapsed: %.2f seconds)\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
		t_par = t_stop;
#endif
	}
	else delta = -1;
#endif

//------------------------------------------------------------------------------------------Compute the mass inside given radii-----------------
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------

#ifdef MASS 																				/*Do it only if the ket MASS is defined*/

		rad1 = 20.0/in_to_kpc;																/*Updating the radii within which to calculate the integral mass [internal units]*/
		rad2 = 10.0/in_to_kpc;
		rad3 = 2.0/in_to_kpc;
		rad4 = 0.3/in_to_kpc;

		A = bar_param[4];																	/*Bar phase*/
		rm = bar_param[1];
		//wm = 0.2/in_to_kpc;																/*Setting the bar full width*/
		//hm = 0.2/in_to_kpc;																/*Setting the bar full height*/
		wm = 2.8*gsoft(header.time);
		hm = 2.8*gsoft(header.time);
		mr = bar_param[5];																	/*Minimum bar extent*/	
	
	rad1 *= rad1;
	rad2 *= rad2;
	rad3 *= rad3;
	rad4 *= rad4;

	C1 = rm*rm;
	C2 = wm*wm*0.25;
	C3 = hm*hm*0.25;
	mass1g=mass2g=mass3g=mass4g=mass1d=mass2d=mass3d=mass4d=mass1s=mass2s=mass3s=mass4s=massBg=massBd=massBs=0.;	/*Initializing the mass array*/

#pragma omp parallel for private (i), reduction(+:mass1g), reduction(+:mass2g), reduction(+:mass3g), reduction(+:mass4g), reduction(+:massBg)
	for (i=0; i<header.nsph; i++)															/*Evaluating the gaseous mass*/
	{
		if ((dist(gas_particles[i].pos, cm))<rad1)
		{
			mass1g += gas_particles[i].mass;
			if ((dist(gas_particles[i].pos, cm))<rad2)
			{
				mass2g += gas_particles[i].mass;
				if ((dist(gas_particles[i].pos, cm))<rad3)
				{
					mass3g += gas_particles[i].mass;
					if ((dist(gas_particles[i].pos, cm))<rad4) mass4g += gas_particles[i].mass;
				}
			}
		}
		rxy = sqrt(gas_particles[i].pos[0]*gas_particles[i].pos[0]+gas_particles[i].pos[1]*gas_particles[i].pos[1]);
		if ((rxy>=mr) &&
			((gas_particles[i].pos[0]*cos(A)+gas_particles[i].pos[1]*sin(A))*(gas_particles[i].pos[0]*cos(A)+gas_particles[i].pos[1]*sin(A))<C1) && 
			((-gas_particles[i].pos[0]*sin(A)+gas_particles[i].pos[1]*cos(A))*(-gas_particles[i].pos[0]*sin(A)+gas_particles[i].pos[1]*cos(A))<C2) && 
			(gas_particles[i].pos[2]*gas_particles[i].pos[2]<C3))		
				massBg += gas_particles[i].mass;
	}
#pragma omp parallel for private (i), reduction(+:mass1d), reduction(+:mass2d), reduction(+:mass3d), reduction(+:mass4d), reduction(+:massBd)
	for (i=0; i<header.ndark; i++)															/*Evaluating the dark matter mass*/
	{
		if ((dist(dark_particles[i].pos, cm))<rad1)
		{
			mass1d += dark_particles[i].mass;
			if ((dist(dark_particles[i].pos, cm))<rad2)
			{
				mass2d += dark_particles[i].mass;
				if ((dist(dark_particles[i].pos, cm))<rad3)
				{
					mass3d += dark_particles[i].mass;
					if ((dist(dark_particles[i].pos, cm))<rad4) mass4d += dark_particles[i].mass;
				}
			}
		}
		rxy = sqrt(dark_particles[i].pos[0]*dark_particles[i].pos[0]+dark_particles[i].pos[1]*dark_particles[i].pos[1]);
		if ((rxy>=mr) && 
			((dark_particles[i].pos[0]*cos(A)+dark_particles[i].pos[1]*sin(A))*(dark_particles[i].pos[0]*cos(A)+dark_particles[i].pos[1]*sin(A))<C1) && 
			((-dark_particles[i].pos[0]*sin(A)+dark_particles[i].pos[1]*cos(A))*(-dark_particles[i].pos[0]*sin(A)+dark_particles[i].pos[1]*cos(A))<C2) && 
			(dark_particles[i].pos[2]*dark_particles[i].pos[2]<C3))
				massBd += dark_particles[i].mass;
	}
#pragma omp parallel for private (i), reduction(+:mass1s), reduction(+:mass2s), reduction(+:mass3s), reduction(+:mass4s), reduction(+:massBs)
		for (i=0; i<header.nstar; i++)														/*Evaluating the stellar mass*/
		{
			if ((dist(star_particles[i].pos, cm)<rad1) && (star_particles[i].tform>=0))
			{
				mass1s += star_particles[i].mass;
				if ((dist(star_particles[i].pos, cm))<rad2)
				{
					mass2s += star_particles[i].mass;
					if ((dist(star_particles[i].pos, cm))<rad3)
					{
						mass3s += star_particles[i].mass;
						if ((dist(star_particles[i].pos, cm))<rad4) mass4s += star_particles[i].mass;
					}
				}
			}
			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);
			if ((rxy>=mr) && 
			((star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))*(star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))<C1) && 
			((-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))*(-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))<C2) && 
			(star_particles[i].pos[2]*star_particles[i].pos[2]<C3)) /*Both stars and Black holes*/		
				massBs += star_particles[i].mass;			
		}

#ifdef DEBUG
	clock_gettime(CLOCK_MONOTONIC, &t_stop);												/*Mass calculation time elapsed*/
	printf("Mass analysis completed\n(Elapsed: %.2f seconds)\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
	t_par = t_stop;
#endif

#endif

//------------------------------------------------------------------------------------------Printing bar particles if DEBUG key is active------

#ifdef PSNAP
	if (bar_param[0]>0)
	{
		sprintf(bar_ids, "%s/bar_ids.%s", out_path, id);

		A = bar_param[4];																	/*Bar phase set to Rphi*/
		rm = bar_param[1];																	/*Bar (half-)length*/
		wm = 1./in_to_kpc;																	/*Setting the bar full width*/
		hm = 1./in_to_kpc;																	/*Setting the bar full height*/
		mr = bar_param[5];																	/*Minimum bar extent*/
	
		C1 = rm*rm;
		C2 = wm*wm*0.25;
		C3 = hm*hm*0.25;
								
		turin = fopen(bar_ids, "wb");														/*Opening the output file*/
		if (turin==NULL)
		{
		printf("File '%s' cannot be opened\n\n", bar_ids);
		return FALSE;
		}	

		for (i=0; i<header.nstar; i++)
		{

			rxy = sqrt(star_particles[i].pos[0]*star_particles[i].pos[0]+star_particles[i].pos[1]*star_particles[i].pos[1]);	/*Getting for each particle the module of the vector projection on the xy plane*/

			if ((rxy>=mr) && ((star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))*(star_particles[i].pos[0]*cos(A)+star_particles[i].pos[1]*sin(A))<C1) && ((-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))*(-star_particles[i].pos[0]*sin(A)+star_particles[i].pos[1]*cos(A))<C2) && (star_particles[i].pos[2]*star_particles[i].pos[2]<C3))
			{
				fwrite(&star_particles[i].index, sizeof(unsigned long long int), 1, turin);		
			}
		}
		if (fclose(turin)==0) printf("file '%s' written\n", bar_ids);
		else
		{
			printf("\nCan't write on file '%s'\n\n", bar_ids);								/*Closing the file 'bar_ids.out'*/
			exit(1);
		}
	}
#endif

//-------------------------------------------------------------------------------------------Printing if the key DEBUG is not defined-----------															
//--------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------

#ifndef DEBUG 																				/*Don't write any output in debugging session*/

	if(access(outfou, F_OK) == -1) 															/*Check if the file doesn't exist*/
	{
		turin = fopen(outfou, "w");															/*Opening the final parameter file*/
		if (turin==NULL)																	/*Checking on the opening*/
		{
			printf("\nError in opening '%s' file\n\n", outfou);
			exit(1);
		}
#if !defined BOXY && !defined REFINE
		fprintf(turin, "Galaxy  z\t\tR(A2_max)\tR_Phi\t\tA2tot_max\tA2_max\t\ttheta_mean\tRmin\t\tkrot\n");
#endif
#if defined BOXY && !defined REFINE
		fprintf(turin, "Galaxy  z\t\tR(A2_max)\tR_Phi\t\tA2tot_max\tA2_max\t\ttheta_mean\tRmin\t\tkrot\t\t\tdelta\n");
#endif
#if defined REFINE && !defined BOXY
		fprintf(turin, "Galaxy  z\t\tR(A2_max)\tR_Phi\t\tA2tot_max\tA2_max\t\ttheta_mean\tRmin\t\tkrot\t\tSigma_z/Sigma_R\tUP/DOWN\n");	/*Printing the parameter file header*/
#endif
#if defined REFINE && defined BOXY
		fprintf(turin, "Galaxy  z\t\tR(A2_max)\tR_Phi\t\tA2tot_max\tA2_max\t\ttheta_mean\tRmin\t\tkrot\t\tSigma_z/Sigma_R\tUP/DOWN\t\tdelta\n");
#endif
	}
	else
	{
		turin = fopen(outfou, "a");															/*Opening the final parameter file*/
		if (turin==NULL)																	/*Checking on the opening*/
		{
			printf("\nError in opening '%s' file\n\n", outfou);
			exit(1);
		}
	}

#if !defined BOXY && !defined REFINE
	fprintf(turin, "%s\t%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, (1./header.time)-1, in_to_kpc*bar_param[0], in_to_kpc*bar_param[1], bar_param[2], bar_param[3], bar_param[4]*(180./M_PI), bar_param[5]*in_to_kpc, Krot/K);
#endif
#if defined BOXY && !defined REFINE
	fprintf(turin, "%s\t%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, (1./header.time)-1, in_to_kpc*bar_param[0], in_to_kpc*bar_param[1], bar_param[2], bar_param[3], bar_param[4]*(180./M_PI), bar_param[5]*in_to_kpc, Krot/K, delta);
#endif
#if defined REFINE && !defined BOXY
	fprintf(turin, "%s\t%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%f\t%f\n", id, (1./header.time)-1, in_to_kpc*bar_param[0], in_to_kpc*bar_param[1], bar_param[2], bar_param[3], bar_param[4]*(180./M_PI), bar_param[5]*in_to_kpc, Krot/K, Sigma_z/Sigma_R, barUP/barDOWN);
#endif
#if defined REFINE && defined BOXY
	fprintf(turin, "%s\t%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%f\t%f\t%lf\n", id, (1./header.time)-1, in_to_kpc*bar_param[0], in_to_kpc*bar_param[1], bar_param[2], bar_param[3], bar_param[4]*(180./M_PI), bar_param[5]*in_to_kpc, Krot/K, Sigma_z/Sigma_R, barUP/barDOWN, delta);
#endif
	if (fclose(turin)==0) printf("file '%s' updated\n", outfou);
	else
	{
		printf("\nCan't write on file '%s'\n\n", outfou);									/*Closing the file 'fourier.out'*/
		exit(1);
	}

#ifdef MASS 																				/*Do it only if the key MASS is defined*/

	if(access(outmass, F_OK) == -1) 														/*Check if the file doesn't exist*/
	{
		hurin = fopen(outmass, "w");														/*Opening the output file 'mass.out'*/
		if (hurin==NULL)																	/*Checking on the opening*/
		{
			printf("\nError in opening '%s' file\n\n", outmass);
			exit(1);
		}
		fprintf(hurin, "Galaxy  z\t\tMgas(r<20)   Mdark(r<20)  Mstar(r<20)  Mtot(r<20)   Mgas(r<10)   Mdark(r<10)  Mstar(r<10)  Mtot(r<10)   Mgas(r<2)    Mdark(r<2)   Mstar(r<2)   Mtot(r<2)    Mgas(r<0.3)  Mdark(r<0.3) Mstar(r<0.3) Mtot(r<0.3) Mgas(Bar) Mdark(Bar) Mstar(Bar)\n");		/*Printing the header*/
	}
	else
	{
		hurin = fopen(outmass, "a");														/*Opening the output file 'mass.out'*/
		if (hurin==NULL)																	/*Checking on the opening*/
		{
			printf("\nError in opening '%s' file\n\n", outmass);
			exit(1);
		}
	}

	fprintf(hurin, "%s\t%f\t%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n", id, (1./header.time)-1, mass1g*in_to_msol, mass1d*in_to_msol, mass1s*in_to_msol, (mass1g+mass1d+mass1s)*in_to_msol, mass2g*in_to_msol, mass2d*in_to_msol, mass2s*in_to_msol, (mass2g+mass2d+mass2s)*in_to_msol, mass3g*in_to_msol, mass3d*in_to_msol, mass3s*in_to_msol, (mass3g+mass3d+mass3s)*in_to_msol, mass4g*in_to_msol, mass4d*in_to_msol, mass4s*in_to_msol, (mass4g+mass4d+mass4s)*in_to_msol, massBg*in_to_msol, massBd*in_to_msol, massBs*in_to_msol);
	if (fclose(hurin)==0) printf("file '%s' updated\n", outmass);
	else
	{
		printf("\nCan't write on file '%s'\n\n", outmass);									 /*Closing the file 'mass.out'*/
		exit(1);
	}

#endif

#endif

//-------------------------------------------------------------------------------------------Printing if the key DEBUG is defined---------------															

#ifdef DEBUG 																				 /*Do it only if the key DEBUG is defined*/

	printf("\nGalaxy  z\t\tR(A2_max)\tR_Phi\t\tA2tot_max\tA2_max\t\ttheta_mean\tRmin\t\tk_rot\n");	/*Angles in degrees*/
	printf("%s\t%f\t%f\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\n", id, (1./header.time)-1, in_to_kpc*bar_param[0], in_to_kpc*bar_param[1], bar_param[2], bar_param[3], bar_param[4]*(180./M_PI), bar_param[5]*in_to_kpc, Krot/K);

#ifdef MASS
 
	printf("\nMgas(r<%.1fkpc) Mdark(r<%.1fkpc) Mstar(r<%.1fkpc) Mtot(r<%.1fkpc)\n", rad1*in_to_kpc, rad1*in_to_kpc, rad1*in_to_kpc, rad1*in_to_kpc);
	printf("%e\t%e\t%e\t%e\n", mass1g*in_to_msol, mass1d*in_to_msol, mass1s*in_to_msol, (mass1g+mass1d+mass1s)*in_to_msol);
	printf("Mgas(r<%.1fkpc) Mdark(r<%.1fkpc) Mstar(r<%.1fkpc) Mtot(r<%.1fkpc)\n", rad2*in_to_kpc, rad2*in_to_kpc, rad2*in_to_kpc, rad2*in_to_kpc);
	printf("%e\t%e\t%e\t%e\n", mass2g*in_to_msol, mass2d*in_to_msol, mass2s*in_to_msol, (mass2g+mass2d+mass2s)*in_to_msol);
	printf("Mgas(r<%.1fkpc) Mdark(r<%.1fkpc) Mstar(r<%.1fkpc) Mtot(r<%.1fkpc)\n", rad3*in_to_kpc, rad3*in_to_kpc, rad3*in_to_kpc, rad3*in_to_kpc);
	printf("%e\t%e\t%e\t%e\n", mass3g*in_to_msol, mass3d*in_to_msol, mass3s*in_to_msol, (mass3g+mass3d+mass3s)*in_to_msol);
	printf("Mgas(r<%.1fkpc) Mdark(r<%.1fkpc) Mstar(r<%.1fkpc) Mtot(r<%.1fkpc)\n", rad4*in_to_kpc, rad4*in_to_kpc, rad4*in_to_kpc, rad4*in_to_kpc);
	printf("%e\t%e\t%e\t%e\n", mass4g*in_to_msol, mass4d*in_to_msol, mass4s*in_to_msol, (mass4g+mass4d+mass4s)*in_to_msol);

#endif

#endif

//----------------------------------------------------------------------------------------------------------------------------------------------

#ifdef DEBUG
	clock_gettime(CLOCK_MONOTONIC, &t_stop);												/*Final time elapsed*/
	printf("(Elapsed: %.2f seconds)\n\n", (float)((t_stop.tv_sec-t_par.tv_sec)+(t_stop.tv_nsec-t_par.tv_nsec)/1.0e9));
	printf("(Total elapsed: %.2f seconds)\n\n", (float)((t_stop.tv_sec-t_i.tv_sec)+(t_stop.tv_nsec-t_i.tv_nsec)/1.0e9));	/*Total time elapsed*/
#endif

//----------------------------------------------------------------------------------------------------------------------------------------------

	if (bin!=NULL) freemat(bin, 6);
	if (list!=NULL) freemat(list, 6);
#ifdef BOXY																			
	if (bar_param[0]>0)
	{
		if (density!=NULL) freemat(density, zbin);
		if (median_d!=NULL) free(median_d);
		if (median_u!=NULL) free(median_u);
	}
#endif
	
	return 0;
}
