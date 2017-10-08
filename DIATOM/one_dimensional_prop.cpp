#include <iostream>
#include <algorithm>
#include <complex>
#include <vector>

#include <fftw3.h>

using namespace std;
	
const double MASS = 1.0;
const double OMEGA = 1.0;

double potential( double x )
{
	return 0.5 * OMEGA * pow(x, 2);
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
		
		double get_d( void ) { return _d; }

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

	complex<double> _xi;
	complex<double> _eta;

	vector< complex<double> > grid_wavefunction;

	vector<double> grid_coordinates;
	vector<double> grid_impulses;

	vector<double> grid_potential;
	vector<double> grid_kinetic_energy;

	public:
		Wavepacket( double q, double p, double width, vector<double> grid_coordinates, vector<double> grid_impulses );
		~Wavepacket( void );

		void calculate_potential( void );
		void calculate_kinetic_energy( void );
		void calculate_wavefunction( void );

		void show_grid_potential( void );
		void show_grid_kinetic_energy( void );
		void show_grid_wavefunction( void );

		void normalize_wavefunction( double d );
		void normalize_fft_result( fftw_complex* res, const int N );

		void copy_to( fftw_complex* arr );
		void copy_from (fftw_complex* arr );

		void fftshift();
		void ifftshift();

		void propagate( double dt );
		void propagate_kinetic_part( double dt );
		void propagate_potential_part( double dt );
};

Wavepacket::Wavepacket( double q, double p, double width, vector<double> _grid_coordinates, vector<double> _grid_impulses )
{
	complex<double> i (0, 1);
	
	cout << "Initializing wavepacket on x/p grid" << endl;
	_q = q;
	_p = p;
	_width = width;

	grid_coordinates = _grid_coordinates;
	grid_impulses = _grid_impulses;
	
	_xi = 2 * _width * _q + i * _p;
	_eta = - _width * pow( _q, 2 ) - i * _q * _p;
}	

Wavepacket::~Wavepacket( void )
{
	cout << "Wavepacket is deleted" << endl;
}

void Wavepacket::calculate_potential( void )
{
	for ( int i = 0; i < grid_coordinates.size(); i++ )
	{
		grid_potential.push_back( potential(grid_coordinates[i]) );
	}
}

void Wavepacket::calculate_kinetic_energy( void )
{
	for ( int i = 0; i < grid_impulses.size(); i++ )
	{
		grid_kinetic_energy.push_back( grid_impulses[i] * grid_impulses[i] / 2.0 / MASS );
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
void Wavepacket::normalize_wavefunction( double d )
{
	double s = 0;
	for ( int i = 0; i < grid_wavefunction.size(); i++ )
	{
		s += pow( abs( grid_wavefunction[i] ), 2 );
	}
	s = sqrt(s * d);

	for ( int i = 0; i < grid_wavefunction.size(); i++ )
	{
		grid_wavefunction[i] /= s; 
	}
}

void Wavepacket::propagate_potential_part( double dt )
{
	complex<double> i(0, 1);
	const int N = grid_wavefunction.size();

	// potential part	
	for ( int counter = 0; counter < N; counter++ )
	{
		grid_wavefunction[counter] *= exp( - i * grid_potential[counter] * dt );
	}
}

void Wavepacket::copy_to( fftw_complex* arr )
{
	for ( int counter = 0; counter < grid_wavefunction.size(); counter++ )
	{
		arr[counter][0] = real( grid_wavefunction[counter] );
		arr[counter][1] = imag( grid_wavefunction[counter] );
	}
}

void Wavepacket::copy_from( fftw_complex* arr )
{
	complex<double> i( 0, 1 );
	for ( int counter = 0; counter < grid_wavefunction.size(); counter++ )
	{
		grid_wavefunction[counter] = arr[counter][0] + i * arr[counter][1]; 
	}
}

void Wavepacket::propagate_kinetic_part( double dt )
{
	complex<double> i(0, 1);
	const int N = grid_wavefunction.size();

	fftw_complex *in = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * N );
	fftw_complex *out = ( fftw_complex* ) fftw_malloc( sizeof( fftw_complex ) * N );
	
	copy_to( in );
	for ( int i = 0; i < N; i++ )
	{
		cout << "in[" << i << "] = " << in[i][0] << " + i * " << in[i][1] << endl;
	}

	fftw_plan plan_backward = fftw_plan_dft_1d( N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE );

	fftw_execute( plan_backward );

	// normalizing by factor N
	normalize_fft_result( out, N );

	copy_from( out );

	ifftshift();

	for ( int counter = 0; counter < N; counter++ )
	{
			grid_wavefunction[counter] *= exp( - i * grid_kinetic_energy[counter] * dt );
	}	

	fftshift();

	cout << "after fftshift" << endl;
	show_grid_wavefunction();

	fftw_complex *in2 = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * N );
	fftw_complex *out2 = (fftw_complex*) fftw_malloc( sizeof( fftw_complex ) * N );

	copy_to( in2 );
	for ( int i = 0; i < N; i++ )
	{
			cout << "in2[" << i << "] = " << in2[i][0] << endl;
	}

	fftw_plan plan_forward = fftw_plan_dft_1d( N, in2, out2, FFTW_FORWARD, FFTW_ESTIMATE );

	fftw_execute( plan_forward );

	for ( int i = 0; i < N; i++ )
	{
			cout << "out2[" << i << "] = " << out2[i][0] << " + i * " << out2[i][1] << endl;
	}

	copy_from( out2 );

	show_grid_wavefunction();	
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
	for ( int i = 0; i < N; i++ )
	{
		res[i][0] /= N;
		res[i][1] /= N;
	}
}

int main ( int argc, char* argv[] )
{
	Grid grid( -10.0, 10.0, 10 );

	vector<double> x = grid.get_grid_coordinates();
	vector<double> p = grid.get_grid_impulses();

	double q0 = -2.0;
	double p0 =  1.0;
	double width0 = 0.5;

	double d = grid.get_d();

	Wavepacket wp( q0, p0, width0, x, p );   

	wp.calculate_potential();
	wp.show_grid_potential();

	wp.calculate_kinetic_energy();
	wp.show_grid_kinetic_energy();

	wp.calculate_wavefunction();
	wp.show_grid_wavefunction();

	wp.normalize_wavefunction( d );
	wp.show_grid_wavefunction();

	double dt = 0.001;
	wp.propagate( dt );

	return 0;	
}

