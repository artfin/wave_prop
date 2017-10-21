#include <iostream>
#include <vector>

#include <cmath>
#include <fstream>
#include <sstream>

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

#include <fftw3.h>

#define REAL_PART 0
#define IMAG_PART 1

using namespace std;

const double LIGHTSPEED = 3.0 * 1E10;
const long double PLANCKCONST = 6.62 * 1E-34;
const long double BOLTZCONST = 1.23 * 1E-23;
const double Temperature = 250.0;

void read_dipole( string filename, vector<double> &ts, vector<double> &ds )
{
	ifstream input( filename );
	
	string line;
	double temp;
	vector<double> v;

	if ( input.is_open() )
	{
		while ( getline( input, line ) )
		{
			stringstream iss( line );

			while ( iss >> temp )
			{
				v.push_back( temp );
			}

			ts.push_back( v[0] );
			ds.push_back( v[1] );

			v.clear();
		}
	
		input.close();
	} 
	else 
	{
		cout << "Unable to open the file." << endl;
	}
}

void save_spectrum( string filename, vector<double> &freqs, vector<double> &ints )
{
	ofstream output( filename );

	for ( int i = 0; i < freqs.size(); i++ )
	{
		output << freqs[i] << " " << ints[i] << endl;
	}

	output.close();
}

void plot_signal( Gnuplot &gp, vector<double> &ts, vector<double> &ds )
{
	vector< pair<double, double>> signal;
	for ( int i = 0; i < ts.size(); i++ )
	{
		signal.push_back( make_pair( ts[i], ds[i] ));	
	}

	gp << "set xrange [0:400]\n;";
		
	gp << "plot '-' with lines title 'signal'" << endl;
	gp.send1d( signal );

	gp.flush();
}

void show_vector( string name, vector<double> v )
{
	cout << "##################################" << endl;
	
	for ( int i = 0; i < v.size(); i++ )
	{
		cout << name << "[" << i << "] = " << v[i] << endl;
	}
	
	cout << "##################################" << endl;
}

vector<double> linspace( double min, double max, double npoints )
{
	double step = ( max - min ) / ( npoints - 1 ); 
	vector<double> res;
	for ( double temp = min; temp <= max; temp += step )
	{
		res.push_back( temp );
	}

	return res;
}

vector<double> sample_signal( vector<double> v )
{
	vector<double> res;
	for ( int i = 0; i < v.size(); i++ )
	{
		res.push_back( sin( 5 * v[i]) );
	}

	return res;
}

void multiply_vector( vector<double> &v, double factor )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		v[i] *= factor;
	}
}

void copy_to( vector<double> &v, fftw_complex* arr )
{
	for ( int i = 0; i < v.size(); i++ )
	{
		arr[i][0] = v[i];
		arr[i][1] = 0;
	}
}

vector<double> fft( vector<double> signal )
{
	int npoints = signal.size();

	fftw_complex *signal_time = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * npoints ); 

	fftw_complex *signal_freq = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * npoints );

	copy_to( signal, signal_time );

	//// +1 -- BACKWARD
	fftw_plan plan = fftw_plan_dft_1d( npoints, signal_time, signal_freq, +1, FFTW_ESTIMATE );

	fftw_execute( plan );
	
	fftw_destroy_plan( plan );

	double autocor;
	vector<double> ints;

	for ( int i = 0; i < npoints / 2; i++ )
	{
		autocor = signal_freq[i][REAL_PART] * signal_freq[i][REAL_PART] + signal_freq[i][IMAG_PART] * signal_freq[i][IMAG_PART];

		ints.push_back( autocor );
	}
	
	fftw_free( signal_time );
	fftw_free( signal_freq );

	return ints;
}

vector<double> binary_absorption( vector<double> &freq, vector<double> &ints )
{
	double temp;
	vector<double> res;

	for ( int i = 0; i < freq.size(); i++ )
	{
		temp = freq[i] * ints[i] * ( 1 - exp(- PLANCKCONST * LIGHTSPEED * freq[i] / (BOLTZCONST * Temperature) ) ); 
		res.push_back( temp );	
	}

	return res;
}

int main()
{
	// vectors of time and induced dipole moment
	vector<double> t, signal;
	read_dipole( "dipole_moment2.txt", t, signal );

	//vector<double> t = linspace( 0.0, npoints * sampling_time, npoints );
	//show_vector( "t", t );

	//vector<double> signal = sample_signal( t );
	//show_vector( "signal", signal );

	int npoints = t.size();
	double sampling_time = t[1] - t[0];

	vector<double> freqs = linspace( 0.0, 1.0 / ( 2.0 * sampling_time), npoints / 2.0 );
	multiply_vector( freqs, 2 * M_PI );
	
	const double ATU = 2.418884326505 * pow( 10, -17 );
	const double CMTOHZ = 3.335631 * pow( 10, -11 );
	multiply_vector( freqs, CMTOHZ / ATU );

	show_vector( "freqs", freqs );
	
	vector<double> ints = fft( signal );
	vector<double> bin_abs = binary_absorption( freqs, ints );

	show_vector( "ints", ints );

	Gnuplot gp;
	plot_signal( gp, freqs, ints );
	//plot_signal( gp, freqs, bin_abs );

	save_spectrum( "spectrum", freqs, ints );

	fftw_cleanup();

	return 0;
}
