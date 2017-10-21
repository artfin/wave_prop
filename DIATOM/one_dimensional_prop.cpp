#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>
#include <chrono>

#include <fstream>
#include <string>
#include <sstream>

#include <fftw3.h>

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

#include "ar_he_pes.h"
#include "ar_he_dip.h"

using namespace std;

typedef chrono::high_resolution_clock::time_point time_point;
typedef chrono::milliseconds milliseconds;

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

const complex<double> IMAG_I( 0, 1 );

double potential( double x )
{
	return ar_he_pot( x ); 
}

class Grid {
	
	double _min;
	double _max;
	int _npoints;

	double _d;

	double _step;
		
	vector<double> grid_coordinates;
	vector<double> grid_impulses;

	public:

		Grid( double min, double max, int npoints );
		~Grid(); 
		
		void create_coordinate_grid( void );
		void show_coordinate_grid( void );

		void fftfreq( void );
		void fftshift( void );
		void show_impulse_grid( void );

		vector<double> get_grid_impulses( void ) { return grid_impulses; };
		vector<double> get_grid_coordinates( void ) { return grid_coordinates; };
};	

void Grid::create_coordinate_grid( void )
{
	for ( int i = 0; i < this->_npoints; i++ )
	{
		grid_coordinates.push_back( this->_min + i * this->_step );
	}
}

void Grid::show_coordinate_grid( void )
{
	for ( int i = 0; i < grid_coordinates.size(); i++ )
	{
		cout << "grid_coordinates[" << i << "] = " << grid_coordinates[i] << endl;
	}
}

void Grid::show_impulse_grid( void )
{
	cout << "====================================" << endl;
	for ( int i = 0; i < grid_impulses.size(); i++ )
	{
		cout << "grid_impulses[" << i << "] = " << grid_impulses[i] << endl;
	}
	cout << "====================================" << endl;
}

Grid::Grid( double min, double max, int npoints )
{
	cout << "Grid is constructed; " << endl;
	_min = min;
	_max = max;
	_npoints = npoints;
	_step = ( _max - _min ) / ( _npoints - 1 );
	
	_d = ( _max - _min ) / _npoints;

	cout << "Parameters of grid are: " << endl;
	cout << ">> min = " << _min << endl;
	cout << ">> max = " << _max << endl;
	cout << ">> npoints = " << _npoints << endl;

	create_coordinate_grid();
	fftfreq();
	fftshift();
}

// return the Discrete Fourier Transform sample frequencies
void Grid::fftfreq()
{
	if ( this->_npoints % 2 == 0 )
	{
		int i;
		for ( i = 0; i < ( this->_npoints / 2 ); i++ )
		{
			grid_impulses.push_back( i / (this->_npoints * this->_d) );
		}
		for ( i = - this->_npoints / 2; i < 0; i++ )
		{
			grid_impulses.push_back( i / (this->_npoints * this->_d) );
		}	
	}
	else
	{
		int i;
		for ( i = 0; i < ( this->_npoints + 1) / 2; i++ )
		{
			grid_impulses.push_back( i / ( this->_d * this->_npoints ) );
		}
		for ( i = - (this->_npoints - 1) / 2; i < 0; i++ )
		{
			grid_impulses.push_back( i / ( this->_d * this->_npoints ) );
		}
	}
}

// shift the zero-frequency component to the center of the spectrum
void Grid::fftshift()
{
	int center = ( this->_npoints + 1 ) / 2 ;
	rotate( grid_impulses.begin(), grid_impulses.begin() + center, grid_impulses.end() );
}

Grid::~Grid( void )
{
	cout << "Grid is deleted" << endl;
}

class Wavepacket {

	double _width;

	double _q;
	double _p;
	double _j;

	complex<double> _xi;
	complex<double> _eta;

	double x_step;
	double p_step;

	double grid_length;

	double p_mean;

	vector<double> grid_coordinates;
	vector<double> grid_impulses;

	vector<double> grid_potential;
	vector<double> grid_kinetic_energy;

	vector<double> grid_dipole;

	public:
		vector< complex<double> > grid_wavefunction;
		
		Wavepacket( double q, double p, double j, double width, vector<double> grid_coordinates, vector<double> grid_impulses );
		~Wavepacket( void );

		double get_p_mean( void ) { return p_mean; }

		double calc_x_step( void );
		double calc_p_step( void );

		void calculate_potential( void );
		void calculate_kinetic_energy( void );
		void calculate_wavefunction( void );
		void calculate_dipole( void );

		void get_potential_grid( vector<double> &vec );

		void show_grid_potential( void );
		void show_grid_kinetic_energy( void );
		void show_grid_wavefunction( void );

		void normalize_wavefunction( void );
		void normalize_fft_result( fftw_complex* res, const int N );

		void copy_to( fftw_complex* arr );
		void copy_from (fftw_complex* arr );

		void fftshift();
		void ifftshift();

		double calculate_q_mean( void );
		double calculate_h_mean( void );
		double calculate_d_mean( void );

		void propagate( double dt );
		void propagate_kinetic_part( double dt );
		void propagate_potential_part( double dt );
		
		void plot( Gnuplot& gp, 
		   		   vector< pair<double, double>> &pot_pts,
				   int counter );
};

Wavepacket::Wavepacket( double q, double p, double j, double width, vector<double> _grid_coordinates, vector<double> _grid_impulses )
{
	cout << "Initializing wavepacket on x/p grid" << endl;
	_q = q;
	_p = p;
	_j = j;

	_width = width;

	grid_coordinates = _grid_coordinates;
	grid_impulses = _grid_impulses;
	
	_xi = 2 * _width * _q + IMAG_I * _p;
	_eta = - _width * pow( _q, 2 ) - IMAG_I * _q * _p;

	grid_length = _grid_coordinates.size();

	x_step = calc_x_step();
	p_step = calc_p_step();
}	

Wavepacket::~Wavepacket( void )
{
	cout << "Wavepacket is deleted" << endl;
}

double Wavepacket::calculate_q_mean( void )
{
	double q_mean = 0;

	q_mean += 0.5 * grid_coordinates[0] * pow( abs(grid_wavefunction[0]), 2 ) * this->x_step;
   	q_mean += 0.5 * grid_coordinates.end()[-1] * pow( abs(grid_wavefunction.end()[-1]), 2 ) * this->x_step;	

	for ( int i = 1; i < grid_coordinates.size() - 1; i++ )
	{
		q_mean += grid_coordinates[i] * pow( abs(grid_wavefunction[i]), 2 ) * this->x_step; 
	}

	return q_mean;
}

double Wavepacket::calculate_d_mean( void )
{
	double d_mean = 0;

	d_mean += 0.5 * grid_dipole[0] * pow( abs(grid_wavefunction[0]), 2 ) * this->x_step;
	d_mean += 0.5 * grid_dipole.end()[-1] * pow( abs(grid_wavefunction.end()[-1]), 2 ) * this->x_step;

	for ( int i = 1; i < grid_coordinates.size() - 1; i++ )
	{
		d_mean += grid_dipole[i] * pow( abs(grid_wavefunction[i]), 2 ) * this->x_step;
	}

	return d_mean;
}

double Wavepacket::calculate_h_mean( void )
{
}

double Wavepacket::calc_x_step( void )
{
	return grid_coordinates.end()[-1] - grid_coordinates.end()[-2];
}

double Wavepacket::calc_p_step( void )
{
	return grid_impulses.end()[-1] - grid_impulses.end()[-2];
}

void Wavepacket::calculate_dipole( void )
{
	double dip;
	for ( int i = 0; i < grid_coordinates.size(); i++ )
	{
		dip = ar_he_dip( grid_coordinates[i] );
	   	grid_dipole.push_back( dip );	
	}
}

void Wavepacket::calculate_potential( void )
{
	double pot;
	for ( int i = 0; i < grid_coordinates.size(); i++ )
	{
		pot = potential( grid_coordinates[i] ) + this->_j * ( this->_j + 1 ) / (2 * MU * pow( grid_coordinates[i], 2 ));
		grid_potential.push_back( pot );
	}
}

void Wavepacket::calculate_kinetic_energy( void )
{
	for ( int i = 0; i < grid_impulses.size(); i++ )
	{
		grid_kinetic_energy.push_back( grid_impulses[i] * grid_impulses[i] / 2.0 / MU );
	}
}

void Wavepacket::calculate_wavefunction( void )
{
	complex<double> i( 0, 1 );

	double x;
	for ( int i =0; i < grid_coordinates.size(); i++ )
	{
		x = grid_coordinates[i];
		grid_wavefunction.push_back( 
			exp(- this->_width * pow(x, 2) + this->_xi * x + this->_eta )
		);
	}	
}

void Wavepacket::show_grid_potential( void )
{
	cout << "====================================" << endl;
	for ( int i = 0; i < grid_potential.size(); i++ )
	{
		cout << "grid_potential[" << i << "] = " << grid_potential[i] << endl;
	}
	cout << "====================================" << endl;
}

void Wavepacket::show_grid_kinetic_energy( void )
{
	cout << "====================================" << endl;
	for ( int i = 0; i < grid_kinetic_energy.size(); i++ )
	{
		cout << "grid_kinetic_energy[" << i << "] = " << grid_kinetic_energy[i] << endl; 
	}
	cout << "====================================" << endl;
}

void Wavepacket::show_grid_wavefunction( void )
{
	cout << "====================================" << endl;
	for ( int i = 0; i < grid_wavefunction.size(); i++ )
	{
		cout << "grid_wavefunction[" << i << "] = " << real( grid_wavefunction[i] ) << " + i * " << imag( grid_wavefunction[i] ) << endl;
	}
	cout << "====================================" << endl;
}

// wtf is d ??
void Wavepacket::normalize_wavefunction( )
{
	double s = 0;
	for ( int i = 0; i < grid_wavefunction.size(); i++ )
	{
		s += pow( abs( grid_wavefunction[i] ), 2 );
	}
	s = sqrt(s * this->x_step);

	for ( int i = 0; i < grid_wavefunction.size(); i++ )
	{
		grid_wavefunction[i] /= s; 
	}
}

void Wavepacket::propagate_potential_part( double dt )
{
	// potential part	
	for ( int counter = 0; counter < this->grid_length; counter++ )
	{ 
		grid_wavefunction[counter] *= exp( - IMAG_I * grid_potential[counter] * dt );
	}
}

void Wavepacket::copy_to( fftw_complex* arr )
{
	for ( int counter = 0; counter < this->grid_length; counter++ )
	{
		arr[counter][0] = real( grid_wavefunction[counter] );
		arr[counter][1] = imag( grid_wavefunction[counter] );
	}
}

void Wavepacket::copy_from( fftw_complex* arr )
{
	for ( int counter = 0; counter < grid_length; counter++ )
	{
		grid_wavefunction[counter] = arr[counter][0] + IMAG_I * arr[counter][1]; 
	}
}

void Wavepacket::propagate_kinetic_part( double dt )
{
	fftw_complex *in = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * this->grid_length );
	fftw_complex *out = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * this->grid_length );
	
	copy_to( in );

	fftw_plan plan_backward = fftw_plan_dft_1d( this->grid_length, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );

	fftw_execute( plan_backward );

	// normalizing by factor N
	normalize_fft_result( out, this->grid_length );

	copy_from( out );

	fftw_destroy_plan( plan_backward );	
	fftw_free( in );
	fftw_free( out );

	ifftshift();

	this->p_mean = 0;
	this->p_mean += grid_impulses[0] * pow( abs(grid_wavefunction[0]), 2 ) * this->p_step;
	this->p_mean += grid_impulses.end()[-1] * pow( abs(grid_wavefunction.end()[-1]), 2 ) * this->p_step;

	// step in impulse space
	for ( int counter = 0; counter < this->grid_length; counter++ )
	{
		grid_wavefunction[counter] *= exp( - IMAG_I * grid_kinetic_energy[counter] * dt );

		if ( counter != 0 && counter != this->grid_length - 1 )
		{		
			this->p_mean += grid_impulses[counter] * pow( abs(grid_wavefunction[counter]), 2 ) * this->p_step;
		}	
	}

	fftshift();

	fftw_complex *in2 = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * this->grid_length );
	fftw_complex *out2 = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * this->grid_length );

	copy_to( in2 );

	fftw_plan plan_forward = fftw_plan_dft_1d( this->grid_length, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE );

	fftw_execute( plan_forward );

	copy_from( out2 );

	fftw_free( in2 );
	fftw_free( out2 );

	fftw_destroy_plan( plan_forward );
}


void Wavepacket::fftshift( void )
{
	const int N = grid_wavefunction.size();
	
	int center = ( N + 1 ) / 2;
	rotate( grid_wavefunction.begin(), grid_wavefunction.begin() + center, grid_wavefunction.end() );
}

// inverse of fftshift 
void Wavepacket::ifftshift( void )
{
	const int N = grid_wavefunction.size();

	int center = N - ( N + 1 ) / 2;
	rotate( grid_wavefunction.begin(), grid_wavefunction.begin() + center, grid_wavefunction.end() );
}

void Wavepacket::propagate( double dt )
{
	propagate_potential_part( dt );
	propagate_kinetic_part( dt );
}

void Wavepacket::normalize_fft_result( fftw_complex* res, const int N )
{
	for ( int counter = 0; counter < N; counter++ )
	{
		res[counter][0] /= N;
		res[counter][1] /= N;
	}
}

void Wavepacket::get_potential_grid( vector<double> &vec )
{
	cout << "inside get potential" << endl;
	for ( int i = 0; i < this->grid_length; i++ )
	{
		vec.push_back( grid_potential[i] );
	}
}

void Wavepacket::plot( Gnuplot& gp, 
			   		   vector< pair<double, double>> &pot_pts,
					   int counter )
{
	vector< pair<double, double> > wavepacket_pts;
	double q_mean = calculate_q_mean();

	for ( int pts_counter = 0; pts_counter < this->grid_length; pts_counter++ )
	{
		wavepacket_pts.push_back( make_pair( grid_coordinates[pts_counter], abs( grid_wavefunction[pts_counter] )) );

	}

	string plot_name = "plot" + static_cast<ostringstream*>( &(ostringstream() << counter) )->str() + ".png";

	gp << "set xrange [-1:20]\n;";
	gp << "set yrange [-0.05:1.5]\n";
			
	//gp << "set label 1 sprintf('q(mean) = " << q_mean << "') at 10,0.3\n";
	//gp << "set label 2 sprintf('p(mean) = " << p_mean << "') at 10,0.4\n"; 

	gp << "set terminal png size 1200,800 enhanced font 'Helvetica,20'\n";
	gp << "set output 'pics/" << plot_name << "'\n";

	gp << "plot '-' with lines title 'psi', '-' with points title 'potential'\n";
	gp.send1d( wavepacket_pts );
	gp.send1d( pot_pts );

	gp.flush();
}


int main ( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./... (int) nsteps" << endl;
		exit( 1 );
	}

	int nsteps = atoi( argv[1] ); 
	
	const double X_MIN = 1.0;
	const double X_MAX = 25.0;
	const int NPOINTS = 8192; 
	
	Grid grid( X_MIN, X_MAX, NPOINTS );

	vector<double> x = grid.get_grid_coordinates();
	vector<double> p = grid.get_grid_impulses();

	double q0 = 10.0;
	double p0 = -30.0;
	double j0 = 3.0;

	double width0 = 0.3;

	Wavepacket wp( q0, p0, j0, width0, x, p );   

	wp.calculate_potential();
	//wp.show_grid_potential();

	wp.calculate_dipole();

	wp.calculate_kinetic_energy();
	//wp.show_grid_kinetic_energy();

	wp.calculate_wavefunction();
	//wp.show_grid_wavefunction();

	wp.normalize_wavefunction( );
	//wp.show_grid_wavefunction();
	
	Gnuplot gp;

	double curr_time = 0.0;
	double dt = 1.0; // 0.05;
	
	const int block_size = 500;	
	vector< pair<double, double> > pot_pts;

	double d_mean;

	vector< time_point > blockTimes;
	blockTimes.push_back( chrono::high_resolution_clock::now() );

	milliseconds time_for_block;
	milliseconds time_for_blocks;

	int block_counter = 0;

	vector<double> potential_grid;
	wp.get_potential_grid ( potential_grid );

	ofstream file;
	file.open( "dipole_moment.txt" );

	for ( int pts_counter = 0; pts_counter < NPOINTS; pts_counter++ )
	{
		pot_pts.push_back( make_pair( x[pts_counter], potential_grid[pts_counter] ));
	}

	for ( int counter = 0; counter < nsteps; counter++ )
	{
		if ( counter % block_size == 0 && counter != 0 )
		{
			block_counter++;

			blockTimes.push_back( chrono::high_resolution_clock::now() );

			time_for_block = chrono::duration_cast<milliseconds>( blockTimes.end()[-1] - blockTimes.end()[-2] );
			time_for_blocks = chrono::duration_cast<milliseconds>( blockTimes.end()[-1] - blockTimes[0] );

			cout << endl;
			cout << "Block " << block_counter << " finished." << endl;
		   	cout << "Time for current block: " << time_for_block.count() / 1000.0 << " s; total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
	
			wp.plot( gp, pot_pts, block_counter );

			d_mean = wp.calculate_d_mean();
		
			// writing dipole to file
			file << curr_time << " " << d_mean << endl;
		}
		
		curr_time += dt;
		wp.propagate( dt );
	}

	blockTimes.push_back( chrono::high_resolution_clock::now() );
	cout << "Total time elapsed: " << chrono::duration_cast<milliseconds>( blockTimes.end()[-1] - blockTimes[0]).count() << " ms" << endl;

	file.close();

	// deallocating FFTW memory
	fftw_cleanup();
	
	return 0;	
}

