/*
	Name: Radar Simulation Program
	Date: Oct. 31, 2016
	Author: Devin Huyghebaert

	Modified: Jan. 17, 2017
		Devin Huyghebaert - Improved some IQ irregularity scatter stuff

	Goal: Create a program to modify IQ samples using a propagation scheme and a user generated ionosphere

	Things to do:
		- implement rf scattering from ground
		- take into account horizontal electron gradients (3-D scattering)
			- need to determine gradient of index in 3-D generated ionosphere
		- implement o- and x- modes due to magnetic field (Implement full Appleton-Hartree)
		- use IRI model to calculate the ionosphere properties

	Current State:
		- Ionosphere is created in function generate_ionosphere()
		- geometric propagation and 2-D index of refraction model implemented
		- current index of refraction calculation uses plasma frequency, does not include collisions or magnetic field
		- generates 2-D plot of ground range vs altitude of rf propagating ray
		- prints the inputs and outputs for the IQ data
		- uses geomag program from igrf to determine magnetic field at ray location (only works to 600 km)
		- scattering from ionosphere is done through a random sample of a gaussian given the spectral width, plasma velocity and power
			of reflection of the ionosphere

*/

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <string.h>

//generate struct for ionosphere
//electron_density, drift_velocity, power_refl, spectral_width x altitude x distance (1km resolution)
struct ionosphere_profile {
	long double e_density;
	double drift_vel;
	double power_refl;
	double spectral_width;
};

//create variable to store ionosphere properties (1000km x 2000km)
struct ionosphere_profile test_profile[1000][2000];

//generate structs for the radar path
//fwd_IQ,rev_IQ,angle_of_prop,altitude,ground_range, direction of propagation, latitude, longitude, aspect angle x sample distance (1km)
struct signal_path {
	long double complex fwd_IQ;
	long double complex rev_IQ;
	double angle_of_propagation;
	double altitude;	
	double ground_range;
	double propagation_direction;  //0 degrees is north, 90 degrees is east
	double latitude;
	double longitude;
	double aspect_angle;
};

//create a rf path of length 2000km
struct signal_path current_signal[2000];
struct signal_path future_signal[2000];

//verticle angle radar ray is propagating
float input_angle_of_propagation = 40.0;

//resolution between samples
float range_resolution = 1.0;

//input number of time steps
int input_time_steps = 2000;

//location of radar for magnetic data
float geo_lat_input = 60.0;
float geo_long_input = -110.0;
float year_input = 2016.2;
char igrf_input[] = "IGRF12.COF";

//horizontal angle of propagation (0 degrees is north)
float propagation_dir_input = 45.0;

//frequency of signal transmitted
float input_frequency = 10000000.0;

//input IQ values
long double complex in_IQ[2000];

//output IQ values
long double complex out_IQ[2000];

//variables to hold the reflection effects on IQ
float mag_reflection = 0.0;
float phase_reflection = 0.0;

//constants
#define ELECTRON_CHARGE 1.602E-19
#define ELECTRON_MASS 9.109E-31
#define PERM_FREE_SPACE 8.854E-12
#define EARTH_RADIUS 6371.5

//global variables to store the inclination of the magnetic field and the offset of the magnetic field from geomagnetic coordinates
double declination_global=0,inclination_global=0;

//--------------------------------------------------

//declare the multiple functions used
int generate_ionosphere(void);
double calculate_index_of_refraction(double input_altitude,double input_ground_range);
int get_mag_field(char file_in[],float year_in, float altitude_in, float latitude_in, float longitude_in);
int plot_ray_trace(void);
double complex IQ_modification(double input_altitude, double input_ground_range, double input_time, double time_factor);
double random_gaussian_generator(double mean, double standard_deviation);


//--------------------------------------------------

int main(void)
{

//int error_code = 0;

//variables to store calculated indices of refraction
double index_1;
double index_2;

//variable to store the dot product of the magnetic field and the rf ray
double mag_dot_ray;

//variables to store horizontal and verticle rf components
double vert_rf;
double horz_rf;
double vert_mag;
double horz_mag;

//variables to store altitude and ground range of rf signal
double temp_altitude;
double temp_ground_range;

//variable to store the forward path loss factor of signal
double fwd_path_loss_factor=0.0,rev_path_loss_factor=0.0;

//scatter factor
double scatter_factor = 0.0;

//irregularity scatter factor
double complex irregularity_scatter_factor = 0.0 + 0.0 * I;

//time between samples
double time_factor = 1.0*1000.0/300000000.0;

//generate an ionosphere
generate_ionosphere();

//set a test pulse
for (int i=0;i<10;i++)
	{
	in_IQ[i] = 1.0 + 1.0 * I;
	}

for (int i=40;i<50;i++)
	{
	in_IQ[i] = 1.0 - 1.0 * I;
	}

//convert the initial angle of propagation to radians
current_signal[0].angle_of_propagation = (input_angle_of_propagation*M_PI)/180.0;
current_signal[0].longitude = geo_long_input;
current_signal[0].latitude = geo_lat_input;
current_signal[0].propagation_direction = propagation_dir_input;

//create radar ray path
for (int i=1;i<2000;i++)
	{

	//Calculate the altitude of the signal
	temp_altitude = sqrt((range_resolution*range_resolution) + ((EARTH_RADIUS+current_signal[i-1].altitude)*(EARTH_RADIUS+current_signal[i-1].altitude)) - (2.0*range_resolution*(EARTH_RADIUS+current_signal[i-1].altitude)*cos(M_PI/2+current_signal[i-1].angle_of_propagation))) - EARTH_RADIUS;

	//calculate the ground range of the signal
	temp_ground_range = EARTH_RADIUS * acos((((EARTH_RADIUS+current_signal[i-1].altitude)*(EARTH_RADIUS+current_signal[i-1].altitude)) + ((EARTH_RADIUS+temp_altitude)*(EARTH_RADIUS+temp_altitude)) - (range_resolution*range_resolution))/(2.0*(EARTH_RADIUS+current_signal[i-1].altitude)*(EARTH_RADIUS+temp_altitude)));

	//determine the latitude of the propagating signal from geometric propagation
	current_signal[i].latitude = (180.0/M_PI)*asin( (sin(current_signal[i-1].latitude*M_PI/180.0)*cos(temp_ground_range/EARTH_RADIUS)) + (cos(current_signal[i-1].latitude*M_PI/180.0)*sin(temp_ground_range/EARTH_RADIUS)*cos(current_signal[i-1].propagation_direction*M_PI/180.0)));
	
	//determine the longitude of the propagating signal from geometric propagation
	current_signal[i].longitude = current_signal[i-1].longitude + ((180.0/M_PI)*atan2(cos(current_signal[i-1].latitude*M_PI/180.0)*sin(temp_ground_range/EARTH_RADIUS)*sin(current_signal[i-1].propagation_direction*M_PI/180.0),cos(temp_ground_range/EARTH_RADIUS) - sin(current_signal[i-1].latitude*M_PI/180.0)*sin(current_signal[i].latitude*M_PI/180.0)));

	//determine the new bearing of the signal from geometric propagation
	current_signal[i].propagation_direction = fmod(((180/M_PI)*(atan2(sin((current_signal[i-1].longitude*M_PI/180.0)-(current_signal[i].longitude*M_PI/180.0))*cos(current_signal[i-1].latitude*M_PI/180.0),cos(current_signal[i].latitude*M_PI/180.0)*sin(current_signal[i-1].latitude*M_PI/180.0) - (sin(current_signal[i].latitude*M_PI/180.0)*cos(current_signal[i-1].latitude*M_PI/180.0)*cos((current_signal[i-1].longitude*M_PI/180.0)-(current_signal[i].longitude*M_PI/180.0)))))+180.0),360.0);

	//record the new ground range and altitude of the propagating ray
	current_signal[i].altitude = temp_altitude;
	current_signal[i].ground_range = current_signal[i-1].ground_range + temp_ground_range;
	
	//calculate the previous index of refraction
	index_1 = calculate_index_of_refraction(current_signal[i-1].altitude,current_signal[i-1].ground_range);
	//calculate the index of refraction at the current point
	index_2 = calculate_index_of_refraction(current_signal[i].altitude,current_signal[i].ground_range);

	//determine the angle of propagation
	if (current_signal[i].altitude < 0.0) //ground scatter
		{
		current_signal[i].angle_of_propagation = -1.0 * current_signal[i-1].angle_of_propagation;
		current_signal[i].altitude = -1.0*(current_signal[i].altitude);
		}	
	else if ( ((index_1/index_2)*sin(M_PI/2.0-current_signal[i-1].angle_of_propagation)) > 1.0) //reflection from plasma
		{
		current_signal[i].angle_of_propagation = -1.0 * current_signal[i-1].angle_of_propagation;
		}
	else if (current_signal[i-1].angle_of_propagation < 0.0) //scattering up
		{
		current_signal[i].angle_of_propagation = -(M_PI/2.0) + asin((index_1/index_2)*sin(M_PI/2.0-current_signal[i-1].angle_of_propagation));
		}
	else //scattering down
		{
			current_signal[i].angle_of_propagation = (M_PI/2.0) - asin((index_1/index_2)*sin(M_PI/2.0-current_signal[i-1].angle_of_propagation));
		}

	//correct for Round Earth effects on angle of propagation
	current_signal[i].angle_of_propagation = current_signal[i].angle_of_propagation + acos((((EARTH_RADIUS+current_signal[i-1].altitude)*(EARTH_RADIUS+current_signal[i-1].altitude)) + ((EARTH_RADIUS+temp_altitude)*(EARTH_RADIUS+temp_altitude)) - (range_resolution*range_resolution))/(2.0*(EARTH_RADIUS+current_signal[i-1].altitude)*(EARTH_RADIUS+temp_altitude)));

	//calculate the magnetic field parameters
	if (current_signal[i].altitude < 599.0)
		{
		get_mag_field(igrf_input,year_input,(current_signal[i].altitude+EARTH_RADIUS),current_signal[i].latitude,current_signal[i].longitude);
		//determine the aspect angle of the signal
		vert_rf = sin(current_signal[i].angle_of_propagation);
		horz_rf = cos(current_signal[i].angle_of_propagation);
		vert_mag = sin(inclination_global*M_PI/180.0);
		horz_mag = cos(inclination_global*M_PI/180.0);
		//dot product of the magnetic field and the rf ray
		mag_dot_ray = (vert_rf*vert_mag)+(cos(current_signal[i].propagation_direction)*horz_rf*cos(declination_global*M_PI/180.0)*horz_mag)+(sin(current_signal[i].propagation_direction)*horz_rf*sin(declination_global*M_PI/180.0)*horz_mag);
		//calculate angle between the magnetic field and the rf ray
		current_signal[i].aspect_angle = acos(mag_dot_ray)*180.0/M_PI;	
		}
	//if the altitude is out of range of the magnetic field model	
	else
		{
		inclination_global = -1000.0;
		current_signal[i].aspect_angle = -1000.0;
		declination_global = -1000.0;
		}


//	printf("latitude of signal: %f\n", current_signal[i].latitude);
//	printf("longitude of signal: %f\n", current_signal[i].longitude);
//	printf("direction of propagation: %f\n", current_signal[i].propagation_direction);
//	printf("ratio of indices of refraction: %f\n",(index_1/index_2));
//	printf("sin of current signal: %f\n",sin(M_PI/2.0-current_signal[i-1].angle_of_propagation));
//	printf("asin of current signal: %f\n",asin((index_1/index_2)*sin(M_PI/2.0-current_signal[i-1].angle_of_propagation)));	
//	printf("current angle of prop: %f\n",current_signal[i].angle_of_propagation*180.0/M_PI);
//	printf("inclination is: %f\n", inclination_global);
	printf("aspect angle is: %f\n",current_signal[i].aspect_angle);
//	printf("declination is: %f\n", declination_global); 
	printf("-----------------------------------\n"); 


	//copy rf ray properties to future signal struct as well
	future_signal[i].angle_of_propagation = current_signal[i].angle_of_propagation;
	future_signal[i].altitude = current_signal[i].altitude;
	future_signal[i].ground_range = current_signal[i].ground_range;
	future_signal[i].propagation_direction = current_signal[i].propagation_direction;
	future_signal[i].latitude = current_signal[i].latitude;
	future_signal[i].longitude = current_signal[i].longitude;
	future_signal[i].aspect_angle = current_signal[i].aspect_angle;
	}

//determine how the rf signal is scattered in the ionosphere
for (int j=0;j<input_time_steps;j++)
	{

	//set the signal near the radar to be the input signal
	current_signal[0].fwd_IQ = in_IQ[j];
	//set the output signal to be the signal propagating to the radar
	out_IQ[j] = current_signal[0].rev_IQ;

	printf("input is: %.18lf + %.18lfi\n", creal(in_IQ[j]),cimag(in_IQ[j]));
	printf("output is: %.18lf + %.18lfi\n", creal(out_IQ[j]),cimag(out_IQ[j]));

	for (int i=0;i<1999;i++)
		{
		//determine the path loss for the signal
		fwd_path_loss_factor = (1.0/((i*1000.0+500.0)/((i-1.0)*1000.0+500.0)))*(1.0/((i*1000.0+500.0)/((i-1.0)*1000.0+500.0)));
		
		if (i > 0)
			{
			//new_fwd = prev_fwd*loss from propagation/spreading - loss from reflections
				//add in faraday rotations?
			future_signal[i].fwd_IQ = current_signal[i-1].fwd_IQ*(fwd_path_loss_factor);
			}

		//calculate the total path loss accumulated and double it for reverse (account for all path loss right at scatter location)
		rev_path_loss_factor = (1.0/((i*1000.0+500.0)/(500.0)))*(1.0/((i*1000.0+500.0)/(500.0)));

		//calculate the scatter factor based on the aspect angle (-10dB/degree)
		scatter_factor = 1.0/(powl(10,abs(current_signal[i].aspect_angle-90.0)));

		//determine the change in the IQ voltages due to irregularities in the ionosphere
		irregularity_scatter_factor = IQ_modification(current_signal[i].altitude,current_signal[i].ground_range,j*time_factor,time_factor);

		//calculate the power reflected by the ionosphere/ground
		future_signal[i].rev_IQ = current_signal[i+1].rev_IQ + current_signal[i].fwd_IQ * scatter_factor * rev_path_loss_factor * irregularity_scatter_factor;

		}

	//transfer future signal IQ to current signal IQ
	for (int i=0;i<2000;i++)
		{
		current_signal[i].fwd_IQ = future_signal[i].fwd_IQ;
		current_signal[i].rev_IQ = future_signal[i].rev_IQ;
		}

	}

//plot a trace of the rf path
plot_ray_trace();

return 0;

}

//--------------------------------------------------

int generate_ionosphere(void)
{

/*
	double drift_vel;
	double power_refl;
	double spectral_width;
*/
	
//update ionosphere model
for (int i=0;i<1000;i++)
	{

	for (int j=0;j<2000;j++)
		{
		//create a parabolic electron density profile for testing
		if ((i > 50) && (i<150))
			{
			test_profile[i][j].e_density = (double)(i-50)*(double)100000000.0;
			test_profile[i][j].drift_vel = 300.0;
			test_profile[i][j].power_refl = 1.0;
			test_profile[i][j].spectral_width = 10.0;
			}		
		if (i > 150)
			{
			test_profile[i][j].e_density = (-1.0)*(double)(i-100.0)*(double)(i-1000.0)*(double)100000.0;
			test_profile[i][j].drift_vel = 1000.0;
			test_profile[i][j].power_refl = 1.0;
			test_profile[i][j].spectral_width = 100.0;
			}		
		
		}
	}

return 0;

}

//---------------------------------------------

double calculate_index_of_refraction(double input_altitude,double input_ground_range) 
{

//variable for plasma freq
float plasma_freq = 0.0;

//create variable to store calculated index of refraction
double index_of_refraction;

	printf("ground range: %f altitude: %f \n",input_ground_range,input_altitude);

	if (input_ground_range < 0.0)
		{
		input_ground_range=0.0;
		}

	//determine the plasma frequency of the present and future path
	plasma_freq = sqrt(test_profile[(int)(input_altitude)][(int)(input_ground_range)].e_density)*ELECTRON_CHARGE/sqrt(PERM_FREE_SPACE*ELECTRON_MASS);

	//calculate the change in the index of refraction
	index_of_refraction = 1.0 - ((plasma_freq)*(plasma_freq)/((input_frequency)*(input_frequency)));

	printf("index of refraction: %f\n", index_of_refraction);

return index_of_refraction;

}

//------------------------------------------------------------------------------

double complex IQ_modification(double input_altitude, double input_ground_range, double input_time, double time_factor)
{

	double gaussian_sample = 0.0;

	double complex mod_factor;

	//generate gaussian with given doppler shift based on plasma parameters
	//test_profile[(int)(input_altitude)][(int)(input_ground_range)].e_density
	//test_profile[(int)(input_altitude)][(int)(input_ground_range)].drift_vel
	//test_profile[(int)(input_altitude)][(int)(input_ground_range)].power_refl
	//test_profile[(int)(input_altitude)][(int)(input_ground_range)].spectral_width

	gaussian_sample = random_gaussian_generator(test_profile[(int)(input_altitude)][(int)(input_ground_range)].drift_vel,test_profile[(int)(input_altitude)][(int)(input_ground_range)].spectral_width);

	mod_factor = test_profile[(int)(input_altitude)][(int)(input_ground_range)].power_refl*cexp(I*((input_time*test_profile[(int)(input_altitude)][(int)(input_ground_range)].drift_vel)+(gaussian_sample*time_factor)));

	return mod_factor;

	//do a weighted random sample of this gaussian to determine IQ transformation

}

//------------------------------------------------------------------------------

double random_gaussian_generator(double mean, double standard_deviation)
{
  float v1,v2,s;

  do {
    v1 = 2.0 * ((float) rand()/RAND_MAX) - 1;
    v2 = 2.0 * ((float) rand()/RAND_MAX) - 1;

    s = v1*v1 + v2*v2;
  } while ( s >= 1.0 );

  if (s == 0.0)
    return 0.0+mean;
  else
    return (standard_deviation*v1*sqrt(-2.0 * log(s) / s))+mean;
}

//------------------------------------------------------------------------------


//function to plot the ray trace (Unix only), requires gnuplot
int plot_ray_trace(void)
{
int NUM_COMMANDS = 7;

double xvals[2000];
double yvals[2000];

    char * commandsForGnuplot[] = {"set title \"Ray Trace\"", "set xlabel 'Ground Range (km)'",
	"set ylabel 'Altitude (km)'","set cbrange[6.0:13.0]",
	"set xrange [0.0:2000.0]","set yrange[0.0:1000.0]", "plot 'data2.temp' matrix with image, 'data.temp' lc 5"};

//, 'data.temp'"

	for (int i=0;i<2000;i++)
		{
		xvals[i] = current_signal[i].ground_range;
		yvals[i] = current_signal[i].altitude;
		}

    FILE * temp2 = fopen("data2.temp", "w");

	for (int i=0;i<1000;i++)
		{
		for (int j=0;j<2000;j++)
			{
    			fprintf(temp2, "%Lf ",log10l(test_profile[i][j].e_density)); //Write the data to a temporary file		
			}
		fprintf(temp2, "\n");
		}

    FILE * temp = fopen("data.temp", "w");
    /*Opens an interface that one can use to send commands as if they were typing into the
     *     gnuplot command line.  "The -persistent" keeps the plot open even after your
     *     C program terminates.
     */
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    int i;
    for (i=0; i < 2000; i++)
    {
    fprintf(temp, "%lf %lf \n", xvals[i], yvals[i]); //Write the data to a temporary file
    }

    for (i=0; i < NUM_COMMANDS; i++)
    {
    fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); //Send commands to gnuplot one by one.
    }
    return 0;
}
