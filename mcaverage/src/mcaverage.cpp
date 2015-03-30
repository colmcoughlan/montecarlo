/*
 Copyright (c) 2014, Colm Coughlan
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// quick monte carlo averaging program. This is the important one!
// colm coughlan

#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cmath>
#include <stdlib.h>
#include <sstream>
#include <fitsio.h>
#include <stdio.h>
extern "C"{
#include "quickfits.h"
}
	using namespace std;



string int_to_str(int i);
double average_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize);
double rms_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize);

int main(int argc, char** argv)
{
	cout<<"Starting program."<<endl;

	if(argc!=23)
	{
		cout<<"Please use the correct number of arguements. Recorded "<<argc - 1<<", need 22."<<endl;
		cout<<"Specify shifting such that: model coords + shift = map coords."<<endl;
		cout<<"All other coords are interpreted in the model frame."<<endl;
		cout<<"This program takes arguments in the form"<<endl;
		cout<<"Filename of model map"<<endl;
		cout<<"Stem of the MC maps"<<endl;
		cout<<"Ending of the MC maps"<<endl;
		cout<<"number of files"<<endl;
		cout<<"Noise box blcx"<<endl;
		cout<<"Noise box blcy"<<endl;
		cout<<"Noise box trcx"<<endl;
		cout<<"Noise box trcy"<<endl;
		cout<<"Souce box 1 blcx"<<endl;
		cout<<"Souce box 1 blcy"<<endl;
		cout<<"Souce box 2 blcx"<<endl;
		cout<<"Souce box 2 blcy"<<endl;
		cout<<"Souce box 3 blcx"<<endl;
		cout<<"Souce box 3 blcy"<<endl;
		cout<<"Total source box blcx"<<endl;
		cout<<"Total source box blcy"<<endl;
		cout<<"Total source box trcx"<<endl;
		cout<<"Total source box trcy"<<endl;
		cout<<"Box size"<<endl;
		cout<<"X Shift"<<endl;
		cout<<"Y Shift"<<endl;
		cout<<"Output filename"<<endl;
		return(1);
	}

	double ra;
	double dec;
	int imsize;
	double cell;	// stored in degrees

	double centre_shift[2];	// where the peak of the source is on the map (x and y coords)
	double rotations[2];	// any rotation applied to the map

	double freq;	// stored in Hz
	double freq_delta;

	char object[FLEN_VALUE];
	char observer[FLEN_VALUE];	// information about the source
	char telescope[FLEN_VALUE];
	double equinox;
	char date_obs[FLEN_VALUE];

	char history[] = "MEM deconvolution performed by PMEM. See PEM logs for details.";

	int stokes;	// stokes parameter

	double bmaj;	// stored in degrees
	double bmin;	// stored in degrees
	double bpa;	// stored in degrees
	int ncc;	// number of clean components

	string modelfilename;
	string mcstem;
	string ending;
	string tfilename;
	string outname;
	stringstream shift_info;

	fstream fout;

	int i,j,k;
	int imsize2;
	int err;
	int nfiles;
	int kx , ky; // used to figure out shifted location if necessary
	int jx , jy;
	double pixels_per_beam;

	double temp;
	int shift[2];
	int box_size;
	int hbox_size;


	int blcx_noise, blcy_noise;
	int trcx_noise, trcy_noise;
	int blcx_sourcebox, blcy_sourcebox;
	int blcx_sourcebox2, blcy_sourcebox2;
	int blcx_sourcebox3, blcy_sourcebox3;
	int blcx_sourcebox4, blcy_sourcebox4;
	int trcx_sourcebox4, trcy_sourcebox4;
	double* null_double;

	double* map_rms;
	double* map_sum;
	double* point_value;
	double* point_value2;
	double* point_value3;

	double* model;
	double* data;
	double* errors;
	double* errors2;
	double* average;

	cout<<"Reading user input."<<endl;

	modelfilename.assign(argv[1]);		// read in the names of the model file, the stem of the Monte Carlo files, and the ending
	mcstem.assign(argv[2]);
	ending.assign(argv[3]);
	nfiles = atoi(argv[4]);			// number of Monte Carlo files
	blcx_noise = atoi(argv[5]);		// coords of a noise box
	blcy_noise = atoi(argv[6]);
	trcx_noise = atoi(argv[7]);
	trcy_noise = atoi(argv[8]);
	blcx_sourcebox = atoi(argv[9]);		// coords of the source box
	blcy_sourcebox = atoi(argv[10]);
	blcx_sourcebox2 = atoi(argv[11]);		// coords of the source box 2 (2nd component)
	blcy_sourcebox2 = atoi(argv[12]);
	blcx_sourcebox3 = atoi(argv[13]);		// coords of the source box 3 (3nd component)
	blcy_sourcebox3 = atoi(argv[14]);
	blcx_sourcebox4 = atoi(argv[15]);		// coords of the source box 4 (Entire source (for total flux))
	blcy_sourcebox4 = atoi(argv[16]);
	trcx_sourcebox4 = atoi(argv[17]);
	trcy_sourcebox4 = atoi(argv[18]);
	box_size = atoi(argv[19]);
	shift[0] = atoi(argv[20]);		// shift to apply between model and MC maps, if necessary (check carefully)
	shift[1] = atoi(argv[21]);
	outname.assign(argv[22]);		// output names

	cout<<"Loaded. Working on model file "<<modelfilename<<endl;
	cout<<"MC stem = "<<mcstem<<endl;
	cout<<"MC end = "<<ending<<endl;

	hbox_size = int(box_size / 2);


	// model coords + shift = map coords

	blcx_noise = blcx_noise + shift[0];
	blcy_noise = blcy_noise + shift[1];

	trcx_noise = trcx_noise + shift[0];
	trcy_noise = trcy_noise + shift[1];

	blcx_sourcebox = blcx_sourcebox  + shift[0];
	blcy_sourcebox = blcy_sourcebox  + shift[1];

	blcx_sourcebox2 = blcx_sourcebox2  + shift[0];
	blcy_sourcebox2 = blcy_sourcebox2  + shift[1];

	blcx_sourcebox3 = blcx_sourcebox3  + shift[0];
	blcy_sourcebox3 = blcy_sourcebox3  + shift[1];

	blcx_sourcebox4 = blcx_sourcebox4  + shift[0];
	blcy_sourcebox4 = blcy_sourcebox4  + shift[1];

	trcx_sourcebox4 = trcx_sourcebox4  + shift[0];
	trcy_sourcebox4 = trcy_sourcebox4  + shift[1];


	tfilename.assign(mcstem);
	tfilename.append(int_to_str(1));
	tfilename.append(ending);
	err = quickfits_read_map_header( tfilename.c_str() , &imsize , &cell , &ra , &dec , centre_shift , rotations , &freq , &freq_delta , &stokes , object , observer , telescope , &equinox , date_obs , &bmaj , &bmin , &bpa , &ncc , -1);
	if(err!=0)
	{
		cout<<endl<<"Error detected opening "<<tfilename<<", err = "<<err<<endl<<endl;
		return(1);
	}
	ncc = 0;	// don't read in any clean components, not necessary

	cout<<"Using cellsize "<<cell*M_PI/180.0<<" radiens and imsize "<<imsize<<" pixels."<<endl;

	imsize2=imsize*imsize;
	pixels_per_beam = bmaj * bmin * M_PI / ( 4.0 * log(2) * cell * cell );
	cout<<"Number of pixels per beam = "<<pixels_per_beam<<endl;



	model=new double[imsize2];
	data=new double[imsize2];
	errors=new double[imsize2];
	errors2=new double[imsize2];
	average=new double[imsize2];

	map_rms=new double[nfiles];
	map_sum=new double[nfiles];
	point_value=new double[nfiles];
	point_value2=new double[nfiles];
	point_value3=new double[nfiles];


	cout<<"Loading model fits file. Assuming cellsize etc. remain the same."<<endl;
	tfilename.assign(modelfilename);
	err = quickfits_read_map(tfilename.c_str() , model , imsize2 , null_double , null_double , null_double , ncc , -1);
	if(err!=0)
	{
		cout<<endl<<"Error detected opening "<<tfilename<<", err = "<<err<<endl<<endl;
		goto free_mem_exit;
	}


	cout<<"Summing up files."<<endl;

	for(k=0;k<imsize2;k++)
	{
		errors[k]=0.0;
		errors2[k]=0.0;
		average[k]=0.0;
	}

	for(i=0;i<nfiles;i++)
	{
		tfilename.assign(mcstem);
		tfilename.append(int_to_str(i+1));
		tfilename.append(ending);
		err=quickfits_read_map(tfilename.c_str() , data , imsize2 , null_double , null_double , null_double , ncc , -1 );
		if(err!=0)
		{
			cout<<endl<<"Error detected opening "<<tfilename<<", err = "<<err<<endl<<endl;
			goto free_mem_exit;
		}
		
		// record flux in region on map (not model)

		point_value[i] = average_region( data , blcx_sourcebox - hbox_size , blcy_sourcebox  - hbox_size ,  blcx_sourcebox + hbox_size , blcy_sourcebox  + hbox_size , imsize );	// find sum of flux
		point_value[i] = point_value[i] * double(pow( box_size , 2)) / pixels_per_beam;

		point_value2[i] = average_region( data , blcx_sourcebox2 - hbox_size , blcy_sourcebox2  - hbox_size ,  blcx_sourcebox2 + hbox_size , blcy_sourcebox2  + hbox_size , imsize );	// find sum of flux
		point_value2[i] = point_value2[i] * double(pow( box_size , 2)) / pixels_per_beam;

		point_value3[i] = average_region( data , blcx_sourcebox3 - hbox_size , blcy_sourcebox3  - hbox_size ,  blcx_sourcebox3 + hbox_size , blcy_sourcebox3  + hbox_size , imsize );	// find sum of flux
		point_value3[i] = point_value3[i] * double(pow( box_size , 2)) / pixels_per_beam;
		
		// record noise and total flux


		map_rms[i] = rms_region( data , blcx_noise , blcy_noise , trcx_noise , trcy_noise , imsize);
		map_sum[i] = average_region( data , blcx_sourcebox4 , blcy_sourcebox4 ,  trcx_sourcebox4 , trcy_sourcebox4 , imsize );	// find sum of flux
		map_sum[i] = map_sum[i] * (trcx_sourcebox4 - blcx_sourcebox4) * (trcy_sourcebox4 - blcy_sourcebox4);

		// generate total error map

		for(kx=0;kx<imsize;kx++)
		{
			jx = kx + shift[0];	// model coords + shift = data coords

			for(ky=0;ky<imsize;ky++)
			{
				jy = ky + shift[1];

				k = kx + imsize * ky;
				j = jx + imsize * jy;

				if( ( (jx>=0) && (jx<imsize) ) && ( (jy>=0) && (jy<imsize) ) )
				{
					average[k]+=data[j];	// update average map (note shift to model reference frame)
					
					temp=model[k]-data[j];
					errors[k]+=temp;
					errors2[k]+=(temp*temp);	// update error maps
				}
			}
		}
	}

	for(k=0;k<imsize2;k++)
	{
		average[k]/=nfiles;
		errors[k]/=nfiles;
		errors2[k]=sqrt(errors2[k]/nfiles);		// unbiased rms doesn't need n-1 factor because comparison map is not from same population 
		// note if you are making a thermal error map you should apply Bessel's correction (n-1)

	}


	if( (shift[0] != 0) || (shift[1] != 0) )	// Blank regions invalid due to shifting
	{
		for(kx=0;kx<imsize;kx++)
		{
			jx = kx + shift[0];

			for(ky=0;ky<imsize;ky++)
			{
				jy = ky + shift[1];

				if( ( (jx<0) && (jx>=imsize) ) && ( (jy<0) && (jy>=imsize) ) )
				{
					k = kx + imsize * ky;
					errors[k] = 666.0;
					errors2[k] = 666.0;
				}
			}
		}
	}

	shift_info<<outname<<"shift_"<<shift[0]<<"_"<<shift[1];

	tfilename.assign(shift_info.str());
	tfilename.append("_average_error.fits");
	err = quickfits_write_map( tfilename.c_str() , errors , imsize , cell , ra , dec , centre_shift , rotations , freq , freq_delta , stokes , object , observer , telescope , equinox , date_obs , history , bmaj , bmin , bpa , 0 , true);

	tfilename.assign(shift_info.str());
	tfilename.append("_rms_error.fits");
	err = quickfits_write_map( tfilename.c_str() , errors2 , imsize , cell , ra , dec , centre_shift , rotations , freq , freq_delta , stokes , object , observer , telescope , equinox , date_obs , history , bmaj , bmin , bpa , 0 , true);

	tfilename.assign(shift_info.str());
	tfilename.append("_average_map.fits");
	err = quickfits_write_map( tfilename.c_str() , average , imsize , cell , ra , dec , centre_shift , rotations , freq , freq_delta , stokes , object , observer , telescope , equinox , date_obs , history , bmaj , bmin , bpa , 0 , true);


	cout<<"Error maps generated. Also writing out distributions."<<endl;

	// rms


	tfilename.assign(shift_info.str());
	tfilename.append("_rmsdist.csv");
	fout.open(tfilename.c_str(), ios::out);

	temp=0.0;
	for(i=0;i<nfiles;i++)
	{
		temp+=map_rms[i];
		fout<<map_rms[i]<<endl;
	}
	temp/=double(nfiles);
	cout<<"Average rms of maps is "<<temp<<" Jy/Beam."<<endl;
	fout.close();
	shift_info.flush();



	// sum of flux in sourcebox



	tfilename.assign(shift_info.str());
	tfilename.append("_sumfluxdist.csv");
	fout.open(tfilename.c_str(), ios::out);

	temp=0.0;
	for(i=0;i<nfiles;i++)
	{
		temp+=map_sum[i];
		fout<<map_sum[i]/double(pixels_per_beam)<<endl;
	}
	temp/=double(nfiles*pixels_per_beam);
	cout<<"Average sum of flux in source region is "<<temp<<" Jy."<<endl;
	fout.close();
	shift_info.flush();



	// point 1



	tfilename.assign(shift_info.str());
	tfilename.append("_p1fluxdist.csv");
	fout.open(tfilename.c_str(), ios::out);

	temp=0.0;

	for(i=0;i<nfiles;i++)
	{
		temp+=point_value[i];
		fout<<point_value[i]<<endl;
	}
	temp/=double(nfiles);
	cout<<"Average value at point 1 = "<<temp<<" Jy."<<endl;
	fout.close();
	shift_info.flush();


	// point 2



	tfilename.assign(shift_info.str());
	tfilename.append("_p2fluxdist.csv");
	fout.open(tfilename.c_str(), ios::out);

	temp=0.0;

	for(i=0;i<nfiles;i++)
	{
		temp+=point_value2[i];
		fout<<point_value2[i]<<endl;
	}
	temp/=double(nfiles);
	cout<<"Average value at point 2 = "<<temp<<" Jy."<<endl;
	fout.close();
	shift_info.flush();

	// point 3



	tfilename.assign(shift_info.str());
	tfilename.append("_p3fluxdist.csv");
	fout.open(tfilename.c_str(), ios::out);

	temp=0.0;

	for(i=0;i<nfiles;i++)
	{
		temp+=point_value3[i];
		fout<<point_value3[i]<<endl;
	}
	temp/=double(nfiles);
	cout<<"Average value at point 3 = "<<temp<<" Jy."<<endl;
	fout.close();
	shift_info.flush();

	// Get the point values and total flux from the model map (undo shifting)

	blcx_sourcebox = blcx_sourcebox - shift[0];
	blcy_sourcebox = blcy_sourcebox - shift[1];

	blcx_sourcebox2 = blcx_sourcebox2 - shift[0];
	blcy_sourcebox2 = blcy_sourcebox2 - shift[1];

	blcx_sourcebox3 = blcx_sourcebox3 - shift[0];
	blcy_sourcebox3 = blcy_sourcebox3 - shift[1];

	blcx_sourcebox4 = blcx_sourcebox4 - shift[0];
	blcy_sourcebox4 = blcy_sourcebox4 - shift[1];
	trcx_sourcebox4 = trcx_sourcebox4 - shift[0];
	trcy_sourcebox4 = trcy_sourcebox4 - shift[1];

	point_value[0] = average_region( model , blcx_sourcebox - hbox_size , blcy_sourcebox  - hbox_size ,  blcx_sourcebox + hbox_size , blcy_sourcebox  + hbox_size , imsize );	// find sum of flux
	point_value[0] = point_value[0] * double(pow( box_size , 2)) / pixels_per_beam;
	cout<<"Model value at point 1 = "<<point_value[0]<<" Jy."<<endl;

	point_value2[0] = average_region( model , blcx_sourcebox2 - hbox_size , blcy_sourcebox2  - hbox_size ,  blcx_sourcebox2 + hbox_size , blcy_sourcebox2  + hbox_size , imsize );	// find sum of flux
	point_value2[0] = point_value2[0] * double(pow( box_size , 2)) / pixels_per_beam;
	cout<<"Model value at point 2 = "<<point_value2[0]<<" Jy."<<endl;

	point_value3[0] = average_region( model , blcx_sourcebox3 - hbox_size , blcy_sourcebox3  - hbox_size ,  blcx_sourcebox3 + hbox_size , blcy_sourcebox3  + hbox_size , imsize );	// find sum of flux
	point_value3[0] = point_value3[0] * double(pow( box_size , 2)) / pixels_per_beam;
	cout<<"Model value at point 3 = "<<point_value3[0]<<" Jy."<<endl;

	map_sum[0] = average_region( model , blcx_sourcebox4 , blcy_sourcebox4 ,  trcx_sourcebox4 , trcy_sourcebox4 , imsize );	// find sum of flux
	map_sum[0] = map_sum[0] * (trcx_sourcebox4 - blcx_sourcebox4) * (trcy_sourcebox4 - blcy_sourcebox4) / pixels_per_beam;
	cout<<"Model total flux = "<<map_sum[0]<<" Jy."<<endl;


	// clean up

	free_mem_exit:

	delete[] model;
	delete[] data;
	delete[] errors;
	delete[] errors2;
	delete[] average;

	delete[] map_rms;
	delete[] map_sum;
	delete[] point_value;
	delete[] point_value2;
	delete[] point_value3;


	cout<<"Program closing."<<endl;

	return(0);
}



string int_to_str(int i)
{
	stringstream ss;
	string str;

	ss<<i;

	str=ss.str();
	ss.flush();

	return(str);
}



/*
	average_region

	average_region takes the average of a given region in the map

	inputs:
		map = the map in question
		imsize = the size of the map (one side)
		blcx = the x coord of the bottom left corner
		blcy = the y coord of the bottom left corner
		trcx = the x coord of the top right corner
		trcy = the y coord of the top right corner

	outputs:
		on return = the average of the region
*/

double average_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize)
{
	int i,j;
	double sum = 0.0;

	trcx++;	// increment ends to make the for loop design a bit easier
	trcy++;

	#pragma omp parallel for collapse(2) reduction(+:sum)
	for(i=blcx;i<trcx;i++)
	{
		for(j=blcy;j<trcy;j++)
		{
			sum += map[j * imsize + i];	// add up over entire region, with openmp if possible
		}
	}

	i = trcx - blcx - 1;	// find number of pixels in region
	j = trcy - blcy - 1;

	i = i *j;

	if(i!=0)
	{
		return(sum/double(i));
	}
	else
	{
		return(map[blcy * imsize + blcx]);
	}
}

/*
	average_region

	rms_region takes the root mean square of a given region in the map

	inputs:
		map = the map in question
		imsize = the size of the map (one side)
		blcx = the x coord of the bottom left corner
		blcy = the y coord of the bottom left corner
		trcx = the x coord of the top right corner
		trcy = the y coord of the top right corner

	outputs:
		on return = the rms of the region
*/

double rms_region(double* map, int blcx, int blcy, int trcx, int trcy, int imsize)
{
	int i,j;
	double sum = 0.0;
	double mean;
	
	mean = average_region( map,  blcx,  blcy,  trcx,  trcy,  imsize );

	trcx++;	// increment ends to make the for loop design a bit easier
	trcy++;

	#pragma omp parallel for collapse(2) reduction(+:sum)
	for(i=blcx;i<trcx;i++)
	{
		for(j=blcy;j<trcy;j++)
		{
			sum += pow(map[j * imsize + i] - mean , 2.0);
		}
	}

	i = trcx - blcx - 1;	// find number of pixels in region
	j = trcy - blcy - 1;

	i = i *j;

	if(i!=0)
	{
		return(sqrt(sum/double(i)));
	}
	else
	{
		return(map[blcy * imsize + blcx] - mean);
	}
}


