// #include "MultiYieldSurfaceMaterial.h"
// #include "RoundedMohrCoulomb_multi_surface.h"
#include "vonMises_multi_surface.h"
#include <fstream>
#include <string>
using namespace std;



int main(int argc, char const *argv[])
{
	double bulk_modulus = (double) atof(argv[1]);
	double scale_hardening = (double) atof(argv[2]);
	double max_strain_in = (double) atof(argv[3]);
	double strain_incr = (double) atof(argv[4]);
	int Nloop   = atoi(argv[5]);

	ofstream outfile;
	outfile.open("strain_stress.txt");

	int material_tag{1}  ;

	double K = bulk_modulus; //16750;
	double p_ratio = 0.15; 
	double G =  3*K*(1-2*p_ratio)/(2+2*p_ratio);  
	double E_in = 2 * G * (1 + p_ratio) ;
	double v_in = p_ratio ;
	double rho_in = 0.0 ;
	int NYS = 8 ;
	vector<double> radius{
		2.7, 2.74, 2.8
		, 2.82, 2.85, 
     	2.9, 3.0, 3.1

	};
	vector<double> HardingPara{
		5500, 4000, 2700
		, 2400, 1890, 1300, 
        915, 600
	};
	for(auto& item: HardingPara){
		item *= scale_hardening;
	}
	auto theMaterial= new vonMises_multi_surface(
		material_tag,
		E_in,
		v_in,
		rho_in,
		NYS,
		radius,
		HardingPara 
	);

	// double max_strain = 1E-3 ;
	// double incr_size = 1E-6;
	auto stress_ret = theMaterial->getStressTensor();
	auto strain_ret = theMaterial->getStrainTensor();
	auto Nactive = theMaterial->getNumActiveYS();

	outfile << strain_ret(0,1)  <<"\t" 
			<< stress_ret(0,1) <<"\t"
			<< Nactive << endl ;


	int Nsteps = max_strain_in/strain_incr ;
	// Loading 
	DTensor2 input_strain(3,3,0.) ;
	for (int i = 0; i < Nsteps; ++i)
	{
		cout<< "--------------------------------------------------" <<endl;
		cout<< "step " << i <<endl;
		input_strain *= 0. ;
		// input_strain(0,0) =  incr_size;
		input_strain(0,1) = strain_incr / 2.;
		input_strain(1,0) = strain_incr / 2.;

		theMaterial->setTrialStrainIncr(input_strain);

		theMaterial->commitState();
		stress_ret = theMaterial->getStressTensor();
		strain_ret = theMaterial->getStrainTensor();
		Nactive = theMaterial->getNumActiveYS();

		outfile << strain_ret(0,1)  <<"\t" 
				<< stress_ret(0,1) <<"\t"
				<< Nactive << endl ;
	}

	for (int loop = 0; loop < Nloop; ++loop)
	{
		// unloading 
		for (int i = 0; i < 2 * Nsteps; ++i)
		{
			cout<< "--------------------------------------------------" <<endl;
			cout<< "unloading step " << i <<endl;
			input_strain *= 0. ;
			// input_strain(0,0) = - strain_incr;
			input_strain(0,1) = - strain_incr / 2.;
			input_strain(1,0) = - strain_incr / 2.;

			theMaterial->setTrialStrainIncr(input_strain);

			theMaterial->commitState();
			stress_ret = theMaterial->getStressTensor();
			strain_ret = theMaterial->getStrainTensor();
			Nactive = theMaterial->getNumActiveYS();

			outfile << strain_ret(0,1)  <<"\t" 
					<< stress_ret(0,1) <<"\t"
					<< Nactive << endl ;
		}

		// reloading 
		for (int i = 0; i < 2 * Nsteps; ++i)
		{
			cout<< "--------------------------------------------------" <<endl;
			cout<< "reloading step " << i <<endl;
			input_strain *= 0. ;
			// input_strain(0,0) =  strain_incr;
			input_strain(0,1) = strain_incr / 2.;
			input_strain(1,0) = strain_incr / 2.;

			theMaterial->setTrialStrainIncr(input_strain);

			theMaterial->commitState();
			stress_ret = theMaterial->getStressTensor();
			strain_ret = theMaterial->getStrainTensor();
			Nactive = theMaterial->getNumActiveYS();

			outfile << strain_ret(0,1)  <<"\t" 
					<< stress_ret(0,1) <<"\t"
					<< Nactive << endl ;
		}
	}
	// Done the experiments and Clean
	delete theMaterial ;
	outfile.close();

	return 0;
}