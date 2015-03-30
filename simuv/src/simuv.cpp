/*
 Copyright (c) 2014, Colm Coughlan
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

// TEMPNOTE: Work on doing the transform centred on the peak. Make a grid of coords with (0,0) on peak

#include <cmath>
#include <iostream>
#include <fstream>
#include "string.h"
#include <sstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
extern "C"{
#include "quickfits.h"
}
	using namespace std;
	
	
#include <unistd.h>


string int_to_str(int i);
int dft_model(double* u, double* v, int nvis, double freq, double* if_array, int nif, int nchan, int central_chan, double chan_width, double* imodel, double* qmodel, double* umodel, double* vmodel, int imsize, double cellsize, double* uvblock, int blocksize, int peak[]);


// g++ -I/Users/admin/git/quickfits -o uvfill uvfill2.cpp -L/Users/admin/git/quickfits -lquickfits -lcfitsio -O3

int main()
{
	double ra;
	double dec;
	double freq;
	int nchan;
	int nif;
	int central_chan;
	double chan_width;
	char object[FLEN_VALUE];
	int blocksize;

	string uvfits;
	string modelname;
	string tempfilename;


	int err;
	int i,j,k;
	int nvis;
	int question;
	int nmaps;


	double temp;
	double stddev;
	double stddev_dterm;

	int imsize;
	int imsize2;
	double cellsize;

	double *imap;
	double *qmap;
	double *umap;
	double *vmap;
	
	double* u_array;
	double* v_array;
	double* uvblock_model;
	double* if_array;

	double* rdterm;
	double* ldterm;
	
	double shift[2];
	shift[0]=0.0;
	shift[1]=0.0;
	
	int peak[2];

	ifstream fin;
	ofstream fout;
	double* null_double;

	char outdata[]="simuv_model";


	cout<<"Welcome to simuv v1.0"<<endl;
	cout<<"This tool will create a simulated observation of model FITS files with the UV data from a real observation"<<endl;
	cout<<"Values in the FITS models are interpreted as fluxes in Jy."<<endl;
	cout<<"Please enter the imsize of the model map (pixels)"<<endl;
	cin>>imsize;
	imsize2 = imsize * imsize;
	cout<<"Please enter the cellsize of the model map (as)"<<endl;
	cin>>cellsize;

	imap=new double[imsize*imsize];
	qmap=new double[imsize*imsize];
	umap=new double[imsize*imsize];
	vmap=new double[imsize*imsize];
	cellsize*=(M_PI)/(180.0*3600);



	cout<<"Please enter the name of the I model map"<<endl;
	cin>>modelname;
	err = quickfits_read_map( modelname.c_str() , imap , imsize2 , null_double , null_double , null_double , 0 , 0 );
	if(err==0)
	{
		cout<<"Imap read successful"<<endl;
		temp = 0.0;
		for(i=0;i<imsize;i++)	// i is dec
		{
			for(j=0;j<imsize;j++)
			{
				if(imap[i*imsize+j]>temp)
				{
					temp = imap[i*imsize+j];
					peak[1] = i;
					peak[0] = j;
				}
			}
		}
		cout<<"Peak of "<<temp<<" detected at ("<<peak[1]<<","<<peak[0]<<")."<<endl;
	}
	else
	{
		cout<<"Problem reading model map."<<endl;
		return(1);
	}

	cout<<"Please enter the name of the Q model map"<<endl;
	cin>>modelname;
	err = quickfits_read_map( modelname.c_str() , qmap , imsize2 , null_double , null_double , null_double , 0 , 0 );
	if(err==0)
	{
		cout<<"Qmap read successful"<<endl;
	}
	else
	{
		cout<<"Problem reading model map."<<endl;
		return(1);
	}

	cout<<"Please enter the name of the U model map"<<endl;
	cin>>modelname;
	err = quickfits_read_map( modelname.c_str() , umap , imsize2 , null_double , null_double , null_double , 0 , 0 );
	if(err==0)
	{
		cout<<"Umap read successful"<<endl;
	}
	else
	{
		cout<<"Problem reading model map."<<endl;
		return(1);
	}
	
	cout<<"Please enter the name of the V model map"<<endl;
	cin>>modelname;
	err = quickfits_read_map( modelname.c_str() , vmap , imsize2 , null_double , null_double , null_double , 0 , 0 );
	if(err==0)
	{
		cout<<"Vmap read successful"<<endl;
	}
	else
	{
		cout<<"Problem reading model map."<<endl;
		return(1);
	}


	cout<<"Please enter the name of the FITS file with the UV information"<<endl;
	cin>>uvfits;

	

	cout<<"Attempting to open "<<uvfits<<endl;
	
	quickfits_read_uv_header(uvfits.c_str(),&ra,&dec,object,&freq,&nvis, &nchan, &central_chan, &chan_width, &nif);
	if(err!=0)
	{
		cout<<"Error attempting to open fits file."<<endl;
		cout<<"Program closing."<<endl;
		return(1);
	}
	
	blocksize = nvis*12*nif*nchan;
	u_array = new double[nvis];
	v_array = new double[nvis];
	uvblock_model = new double[blocksize];
	if_array = new double[nif];

	
	err=quickfits_read_uv_data(uvfits.c_str(),nvis,nchan,nif,u_array,v_array,uvblock_model, if_array);

	
	cout<<"\tRaw UV data and visibility read complete."<<endl;
	cout<<"\t"<<nvis<<" visibilites read."<<endl;
	cout<<"\t"<<nif<<" IF(s) with "<<nchan<<" channel(s) each detected."<<endl;


	err = quickfits_write_map( "model_imap.fits" , imap , imsize , cellsize*(180.0/M_PI) , ra , dec , shift , shift , freq , 0 , 1 , outdata , outdata , outdata , 0 , outdata , outdata , 0 , 0 , 0 , 0 , false);
	err = quickfits_write_map( "model_qmap.fits" , qmap , imsize , cellsize*(180.0/M_PI) , ra , dec , shift , shift , freq , 0 , 2 , outdata , outdata , outdata , 0 , outdata , outdata , 0 , 0 , 0 , 0 , false);
	err = quickfits_write_map( "model_umap.fits" , umap , imsize , cellsize*(180.0/M_PI) , ra , dec , shift , shift , freq , 0 , 3 , outdata , outdata , outdata , 0 , outdata , outdata , 0 , 0 , 0 , 0 , false);
	err = quickfits_write_map( "model_vmap.fits" , vmap , imsize , cellsize*(180.0/M_PI) , ra , dec , shift , shift , freq , 0 , 3 , outdata , outdata , outdata , 0 , outdata , outdata , 0 , 0 , 0 , 0 , false);


	cout<<"Performing DFT - this may take a moment."<<endl;
//  DFT model into uvblock_model
	err = dft_model(u_array, v_array, nvis, freq, if_array, nif, nchan, central_chan, chan_width, imap, qmap, umap, vmap, imsize, cellsize, uvblock_model,blocksize,peak);



	cout<<"DFT complete - writing results to fits file"<<endl;

	fin.open(uvfits.c_str(), fstream::binary);
	fout.open("output_no_noise.fits", fstream::trunc|fstream::binary);
	fout<<fin.rdbuf();
	fout.close();


	tempfilename.assign("output_no_noise.fits");
	err=quickfits_overwrite_uv_data(tempfilename.c_str(), nvis, nchan, nif, u_array, v_array, uvblock_model);

	cout<<"Would you like to add normally distributed noise to the DFT? (1 = yes, 0 = no)"<<endl;
	cin>>question;

	if(question==1)
	{
		const gsl_rng_type * T;
		gsl_rng * r;
		gsl_rng_env_setup();
		unsigned long int seed;
		double* uvblock_temp;
		uvblock_temp = new double[blocksize];
	     
		T = gsl_rng_default;
		r = gsl_rng_alloc (T);

		seed = time (NULL) * getpid();
		gsl_rng_set ( r, seed);


		rdterm=new double[2*10];
		ldterm=new double[2*10];
		
		for(i=2;i<blocksize;i+=3)	// copy model into temporary uvblock (only need weights)
		{
			uvblock_temp[i] = uvblock_model[i];
		}


		for(i=0;i<20;i++)	// format is (real, imag)*10
		{
			rdterm[i]=0.0;
			ldterm[i]=0.0;
		}

		cout<<"How many files would you like to make?"<<endl;
		cin>>nmaps;
		cout<<"What standard deviation would you like the noise distributed about zero to have?"<<endl;
		cin>>stddev;
		stddev/=sqrt(2.0);
		cout<<"Sample Noise values."<<endl;
		for(i=0;i<10;i++)
		{
			cout<<"\t "<<gsl_ran_gaussian(r,stddev)<<endl;
		}

		cout<<"What D term noise would you like to add to the antenna table?"<<endl;
		cout<<"Warning: Only tested for VLBA observations. Enter zero to ignore"<<endl;
		cin>>stddev_dterm;
		stddev_dterm/=sqrt(2.0);

		for(i=0;i<nmaps;i++)
		{
			fin.seekg(0, ios::beg);
			fin.clear();

			modelname.assign("output_");
			modelname.append(int_to_str(i+1));
			modelname.append(".fits");
			fout.open(modelname.c_str(), fstream::trunc|fstream::binary);
			fout<<fin.rdbuf();
			fout.close();

			// add noise to real and imaginary components at each freq, chan, rr, ll, lr, rl (same at each)

			#pragma omp parallel for
			for(j=0;j<blocksize;j+=3)
			{
				uvblock_temp[j] = uvblock_model[j] + gsl_ran_gaussian(r,stddev);
				uvblock_temp[j+1] = uvblock_model[j+1] + gsl_ran_gaussian(r,stddev);
			}

			err=quickfits_overwrite_uv_data(modelname.c_str(),nvis, nchan, nif, u_array, v_array, uvblock_temp);
			
			if(stddev_dterm != 0.0)
			{
				cout<<"Adding D-term noise"<<endl;
				for(j=0;j<20;j++)	// format is (real, imag)*10
				{
					rdterm[j]=gsl_ran_gaussian(r,stddev_dterm);
					ldterm[j]=gsl_ran_gaussian(r,stddev_dterm);
				}

				quickfits_replace_ant_info(modelname.c_str(),rdterm,ldterm);
			}
		}

		delete[] uvblock_temp;
		delete[] rdterm;
		delete[] ldterm;
		gsl_rng_free (r);
		
		cout<<"Process complete."<<endl<<"Closing Program"<<endl;
	}
	else
	{
		cout<<"Closing Program."<<endl;
	}

	delete[] imap;
	delete[] qmap;
	delete[] umap;
	delete[] vmap;
	delete[] u_array;
	delete[] v_array;
	delete[] uvblock_model;
	delete[] if_array;
	fin.close();

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


int dft_model(double* u, double* v, int nvis, double freq, double* if_array, int nif, int nchan, int central_chan, double chan_width, double* imodel, double* qmodel, double* umodel, double* vmodel, int imsize, double cellsize, double* uvblock, int blocksize, int peak[])
{
	int i,j,k, xctr, yctr,n;
	double u_curr, v_curr, temp1, temp2, curr_freq;
	int if_size, chan_size, skip;
	int llr,lli,rrr,rri,rlr,rli,lrr,lri;
	double twopi = 2.0*M_PI;
	double* grid_ra;
	double* grid_dec;
	
	grid_ra = new double[imsize];
	grid_dec = new double[imsize];
	
	for(i=0;i<imsize;i++)	// note ra should be reversed
	{
		grid_ra[i] = -(i-peak[0])*cellsize;
		grid_dec[i] = (i-peak[1])*cellsize;
	}


	chan_size = 12;	// each channel has 4 Stokes parameters with Real, Imag + weight
	if_size = nchan * chan_size;	// each IF has nchan channels

	#pragma omp parallel for
	for(i=0;i<blocksize;i+=3)	// zero everything but the weights
	{
		uvblock[i] = 0.0;
		uvblock[i+1] = 0.0;
	}

	for(i=0;i<nif;i++)
	{
		for(j=0;j<nchan;j++)
		{
			curr_freq = freq + if_array[i] + (j + 1 - central_chan)*chan_width;

			#pragma omp parallel for private(xctr, yctr, n, u_curr, v_curr, temp1, temp2, llr,lli,rrr,rri,rlr,rli,lrr,lri)
			for(k=0;k<nvis;k++)
			{
				u_curr = curr_freq * u[k];
				v_curr = curr_freq * v[k];
				llr = (k*if_size*nif) + i * if_size + j * chan_size;	// each visibility has nif_size * nif data points
				lli = llr +1;
				rrr = llr +3;
				rri = llr +4;	//weights are not changed
				rlr = llr +6;
				rli = llr +7;
				lrr = llr +9;
				lri = llr +10;
				n=0;

				for(xctr = 0; xctr<imsize; xctr++)	// note each thread only adds to one coord (safe)
				{
					for(yctr = 0; yctr<imsize; yctr++)
					{
						temp1 = twopi*((u_curr*grid_ra[yctr]) + (v_curr*grid_dec[xctr]) );
						temp2 = sin(temp1);
						temp1 = cos(temp1);

						uvblock[llr] += (imodel[n] - vmodel[n])*temp1; //ll real
						uvblock[lli] += (imodel[n] - vmodel[n])*temp2; //ll imag

						uvblock[rrr] += (imodel[n] + vmodel[n])*temp1; //rr real
						uvblock[rri] += (imodel[n] + vmodel[n])*temp2; //rr imag

						uvblock[rlr] += qmodel[n]*temp1 - umodel[n]*temp2; //rl real
						uvblock[rli] += umodel[n]*temp1 + qmodel[n]*temp2; //rl imag

						uvblock[lrr] += qmodel[n]*temp1 + umodel[n]*temp2; //lr real
						uvblock[lri] += qmodel[n]*temp2 - umodel[n]*temp1; //lr imag

						n++;
					}
				}
			}
		}
	}


	delete[] grid_ra;
	delete[] grid_dec;

	return(0);
}


