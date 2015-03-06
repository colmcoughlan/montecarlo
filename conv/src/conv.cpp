#include <iostream>
#include <string>
#include <fftw3.h>
#include <omp.h>
#include <cmath>
extern "C"{
#include "quickfits.h"
}

	using namespace std;


int gen_gauss(double* matrix, int imsize, double cellsize, double bmaj, double bmin, double bpa, double* peak);
void arrange_ft(double* arr, int imsize);
int convolve(double* data, fftw_complex* response, int imsize, int pad_factor, double* output , fftw_plan& forward_transform, fftw_plan& backward_transform, double* double_buff, fftw_complex* complex_buff);
int ft_beam(double* beam, fftw_complex* ft_beam, int imsize, int pad_factor, fftw_plan& plan, double* double_buff, fftw_complex* complex_buff);

int main(int argc, char** argv)
{
	if(argc!=7)
	{
		cout<<"Convolution tool to convolved model map with specified beam (as,as,deg)"<<endl;
		cout<<"Residual map also added (set to 0 if none)"<<endl;
		cout<<"Useage: conv model_map residual_map bmaj bmin bpa output_map"<<endl;
		return(1);
	}

	string model_map_name;
	string residual_map_name;
	string output_map_name;

	fftw_complex *beam_ft;
	fftw_complex *complex_buff;

	double *model_map;
	double *residual_map;
	double *double_buff;
	bool do_residual = false;

	double bmaj_in , bmin_in , bpa_in;
	double bmaj_out , bmin_out , bpa_out;
	double null_double;

	int err;
	int pad_factor = 2;
	int stokes;
	double temp;
	int peak[2];

	double ra;
	double dec;
	int imsize , imsize2;
	double cell = 0;	// stored in degrees

	double centre_shift[2];	// where the peak of the source is on the map (x and y coords)
	double rotations[2];	// any rotation applied to the map

	double freq;	// stored in Hz
	double freq_delta;
	int ncc;

	char object[FLEN_VALUE];
	char observer[FLEN_VALUE];	// information about the source
	char telescope[FLEN_VALUE];
	double equinox;
	char date_obs[FLEN_VALUE];

	char history[] = "Convolution performed by conv (external C++ program).";

	// start up some of the fourier tranform settings and variables from FFTW


	fftw_init_threads();	// initialise fftw parallisation
	fftw_plan_with_nthreads(omp_get_max_threads());


	fftw_plan forward_transform, backward_transform;


	// get some information from the user

	model_map_name.assign(argv[1]);
	residual_map_name.assign(argv[2]);
	bmaj_out = atof(argv[3]);
	bmin_out = atof(argv[4]);
	bpa_out = atof(argv[5]);
	output_map_name.assign(argv[6]);
/*
	cout<<"Convolution Tool."<<endl;
	cout<<"Please enter the name of the model map"<<endl;
	cin>>model_map_name;
	cout<<"Please enter the name of the residual map"<<endl;
	cin>>residual_map_name;
	cout<<"Please enter the desired bmaj (as)"<<endl;
	cin>>bmaj_out;
	cout<<"Please enter the desired bmin (as)"<<endl;
	cin>>bmin_out;
	cout<<"Please enter the desired bpa (deg)"<<endl;
	cin>>bpa_out;
	cout<<"What output name would you like?"<<endl;
	cin>>output_map_name;
*/

	// convert bmaj, bmin to degrees

	bmaj_out = bmaj_out / (3600.0);
	bmin_out = bmin_out / (3600.0);

	if(residual_map_name.compare("0") != 0)
	{
		do_residual=true;
		cout<<"Residual map will be added."<<endl;
	}
	else
	{
		cout<<"No residuals will be added."<<endl;
	}

	// read in the information from the header of the RESIDUAL MAP if possible

	if(do_residual)
	{
		err = quickfits_read_map_header( residual_map_name.c_str() , &imsize , &cell , &ra , &dec , centre_shift , rotations , &freq , &freq_delta , &stokes , object , observer , telescope , &equinox , date_obs , &bmaj_in , &bmin_in , &bpa_in , &ncc , 0);
	}
	else
	{
		err = quickfits_read_map_header( model_map_name.c_str() , &imsize , &cell , &ra , &dec , centre_shift , rotations , &freq , &freq_delta , &stokes , object , observer , telescope , &equinox , date_obs , &bmaj_in , &bmin_in , &bpa_in , &ncc , 0);
	}
	
	imsize2 = imsize * imsize;

//	cout<<"Detected cellsize (in as) = "<<cell*3600<<endl;


	// assign memory

	ncc = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// allowing for padding in buffers and reducing memory needed for r2c and c2r transform by taking Herm. conjugacy into account

	model_map = new double[imsize2];
	beam_ft = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * ncc );	// saving memory with herm. conjugacy
	complex_buff = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * ncc );	// allowing room for padding in buffers
	double_buff = ( double* ) fftw_malloc( sizeof( double ) * pad_factor * pad_factor * imsize2 );
	
//	check if residual map has been provided

	if(do_residual)
	{
		residual_map = new double[imsize2];
	}

	// prepare FT plans



	forward_transform = fftw_plan_dft_r2c_2d(imsize * pad_factor , imsize * pad_factor , double_buff , complex_buff , FFTW_ESTIMATE );	// make FFT routines (just doing one, no optimisation)

	backward_transform = fftw_plan_dft_c2r_2d(imsize * pad_factor , imsize * pad_factor , complex_buff , double_buff , FFTW_ESTIMATE );	// r2c is always a forward transform etc

	// make ft of new beam
	
	temp = 0.0;
	for(int i=0;i<imsize;i++)
	{
		for(int j=0;j<imsize;j++)
		{
			if(model_map[i*imsize+j]>temp)
			{
				temp = model_map[i*imsize+j];
				peak[0] = i;
				peak[1] = j;
			}
		}
	}

	gen_gauss( model_map , imsize , cell , bmaj_out , bmin_out , bpa_out, peak);	// make restoring beam

	if( pad_factor == 1)
	{
		arrange_ft( model_map , imsize );	// arrange beam so that the convolution will work well
	}

	ft_beam( model_map , beam_ft , imsize , pad_factor  , forward_transform , double_buff , complex_buff);	// get ft of restoring beam

	// load in model map and convolve it with the new beam

	ncc = 0;
	quickfits_read_map( model_map_name.c_str() , model_map , imsize2 , &null_double , &null_double , &null_double , ncc, 0 );	// load model map
	
	if(do_residual)
	{
		quickfits_read_map( residual_map_name.c_str() , residual_map , imsize2 , &null_double , &null_double , &null_double , ncc, 0 );	// load residual map
	}

	convolve( model_map , beam_ft , imsize , pad_factor  , model_map , forward_transform , backward_transform , double_buff , complex_buff);	// convolve and overwrite model_map

	if(do_residual)
	{
		temp = bmaj_out * bmin_out * M_PI / ( 4.0 * log(2) * cell * cell );	// rescale residuals to new beam units
		temp = temp / ( bmaj_in * bmin_in * M_PI / ( 4.0 * log(2) * cell * cell ) );
	
		#pragma omp parallel for
		for(int i=0;i<imsize2;i++)
		{
			// rescale residuals and add back on
			model_map[i] = model_map[i] + residual_map[i] * temp;
		}
	}

	// write out answer

	err = quickfits_write_map( output_map_name.c_str() , model_map , imsize , cell , ra , dec , centre_shift , rotations , freq , freq_delta , stokes , object , observer , telescope , equinox , date_obs , history , bmaj_out , bmin_out , bpa_out , 0 , true);

	// clean up and exit

	fftw_destroy_plan(forward_transform);
	fftw_destroy_plan(backward_transform);
	fftw_cleanup_threads();

	fftw_free(beam_ft);
	fftw_free(double_buff);
	fftw_free(complex_buff);

	delete[] model_map;
	if(do_residual)
	{
		delete[] residual_map;
	}

	return(0);
}

int gen_gauss(double* matrix, int imsize, double cellsize, double bmaj, double bmin, double bpa, double* peak)
{
	int imsize2=imsize*imsize;
	int i,j,k;
	double x,y,a,b,c;
	int xorigin = peak[0];	// centre of distribution (pixels). Assumed centre is of the form (255,256) for a 512 map.
	int yorigin = peak[1];

	bpa = 90 - bpa; // convert from CCW measurement used in astronomy

	bpa*=(M_PI/180.0);	// convert to radiens
	bmaj*=((M_PI/180.0)/(2.354820045));	// convert to radiens from degrees and from FWHM to sigma (2*sqrt(2*log(2)))) = 2.354820045
	bmin*=((M_PI/180.0)/(2.354820045));

	cellsize*=M_PI/180.0;	// convert from degrees to radiens

	a=0.5*(pow(cos(bpa)/bmaj,2)+pow(sin(bpa)/bmin,2));
	b=0.25*sin(2.0*bpa)*(-1.0/pow(bmaj,2)+1.0/pow(bmin,2));
	c=0.5*(pow(sin(bpa)/bmaj,2)+pow(cos(bpa)/bmin,2));


	i=0;
	j=0;

	for(k=0;k<imsize2;k++)
	{
		x=double(i-xorigin)*cellsize;
		y=double(j-yorigin)*cellsize;
		matrix[k]=exp(-(a*pow(x,2)+2.0*b*x*y+c*pow(y,2)));	// 2d gaussian general form
		j++;
		if(j==imsize)
		{
			j=0;
			i++;
		}
	}

	return(0);
}

void arrange_ft(double* arr, int imsize)
{
	int i,j,k,l;
	int half=(imsize/2);	// 0.5*dimension

	double temp;


	#pragma omp parallel for collapse(2) private(k,l,temp)
	for(i=0;i<half;i++)	// swap values from q1 to q3, and q2 to q4
	{
		for(j=0;j<half;j++)
		{
			k=i*imsize+j;	// location in q1
			l=(i+half)*imsize+half+j;	// location in q3

			temp = arr[k];	// re
			arr[k] = arr[l];
			arr[l] = temp;


			k=i*imsize+j+half;	// location in q2
			l=(i+half)*imsize+j;	// location in q4

			temp = arr[k];	// re
			arr[k] = arr[l];
			arr[l] = temp;


		}
	}
}


/*
	"convolve"

	Convolves data with the repsonse function response and copies output to output.

	Inputs:

	data = data to be convolved
	response = response function (same size as data)
	imsize = size of data, response and output
	pad_factor = amount of zero padding (1 = no zero padding, 2 = imsize zero padding, etc.)
	forward_transform = fftw plan for forward transform
	backward_transform = fftw plan for backward transform
	double_buff = work space for fftw plans
	complex_buff = work space for fftw plans

	Outputs:

	output = the convolved data
	
*/

int convolve(double* data, fftw_complex* response, int imsize, int pad_factor, double* output , fftw_plan& forward_transform, fftw_plan& backward_transform, double* double_buff, fftw_complex* complex_buff)
{
	int i, j, k;

	int imsize_pad = pad_factor * imsize;
	double temp1 , temp2;



	#pragma omp parallel for collapse(2) private(k)	// initialise (padded) working array for forward transform
	for(i=0;i<imsize_pad;i++)
	{
		for(j=0;j<imsize_pad;j++)
		{
			k =  (i * imsize_pad) + j ;
			if( ( i < imsize ) && ( j <imsize ) )
			{
				double_buff[k] = data[ (i * imsize) + j ];
			}
			else
			{
					double_buff[k] = 0.0;
			}
		}
	}



	fftw_execute( forward_transform );	// forward transform



	j = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// this is the point where Hermitian conjugacy means the FFTW routine has stopped spitting out more data
									// no need to process the additional data, as the c2r transform assumes conjugacy too

	k = imsize_pad * imsize_pad;

	#pragma omp parallel for private(temp1,temp2)
	for( i=0 ; i < j ; i++ )	// convolution theorem
	{
		temp1 = ( complex_buff[i][0] * response[i][0] - complex_buff[i][1] * response[i][1] ) / double( k );
		temp2 = ( complex_buff[i][1] * response[i][0] + complex_buff[i][0] * response[i][1] ) / double( k );	// the k scales data_pad and response_pad correctly

		complex_buff[i][0] = temp1;
		complex_buff[i][1] = temp2;
	}


	fftw_execute( backward_transform );	// inverse transform


	if(pad_factor==1)	// k is the displacement necessary to read from the centre of the padded array
	{
		k = 0;
	}
	else
	{
		k = (pad_factor - 1) * imsize / 2;
	}


	#pragma omp parallel for collapse(2)	// copy data to result, reading from the "centre" of the padded output
	for(i=0;i<imsize;i++)
	{
		for(j=0;j<imsize;j++)
		{
			output[ (i * imsize) + j ] = double_buff[ ( (i + k) * imsize_pad) + j + k];
		}
	}

//	temp1 = average_region(output, 100, 100, 120, 120, imsize);
//	clip_edges(output , temp1 , imsize);



	return(0);
}

/*
	ft_beam

	ft_beam FFTs a beam (dirty or restoring) and saves the result in a complex array
	As the FFT is real to complex, the array does not need the same dimensions as the input due to the Hermitian conjugacy of the resulting FT frequencies
	This FFT is done outside the convolution function because the beam only needs to be FFTed once. This saves processing time.

	inputs:

	beam = the beam image to be FFTed (dirty or restoring)
	imsize = the length/width of the image
	pad_factor = amount of zero padding (1 = no zero padding, 2 = imsize zero padding, etc.)
	plan = fftw plan for forward transform
	double_buff = work space for fftw plans
	complex_buff = work space for fftw plans

	outputs:

	ft_beam = the FFT of the beam

*/

int ft_beam(double* beam, fftw_complex* ft_beam, int imsize, int pad_factor, fftw_plan& plan, double* double_buff, fftw_complex* complex_buff)
{
	int imsize_pad = pad_factor * imsize;
	int i, j, k;

	#pragma omp parallel for collapse(2) private(k)	// initialise (padded) working array for transform
	for(i=0;i<imsize_pad;i++)
	{
		for(j=0;j<imsize_pad;j++)
		{
			k =  (i * imsize_pad) + j ;
			if( ( i < imsize ) && ( j < imsize ) )
			{
				double_buff[k] = beam[ (i * imsize) + j ];
			}
			else
			{
				double_buff[k] = 0.0;
			}
		}
	}

//	pad_image( beam , double_buff , imsize, pad_factor);	// initialise (padded) working array for forward transform, using padding correctly

	fftw_execute(plan);	// forward transform

	j = pad_factor * imsize * ( pad_factor * imsize / 2 + 1 );	// only read in as much as necessary, taking the Herm. cong. into account

	for(i=0;i<j;i++)	// copy result to output array
	{
		ft_beam[i][0] = complex_buff[i][0];
		ft_beam[i][1] = complex_buff[i][1];
	}

	return(0);
}

