// #include "MultiYieldSurfaceMaterial.h"
// #include "RoundedMohrCoulomb_multi_surface.h"
// #include "vonMises_multi_surface.h"
#include "DruckerPrager_multi_yield_surface.h"
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
	double initial_confinement = (double) atof(argv[6]);
	double reference_pressure= (double) atof(argv[7]);
	double modulus_n_in= (double) atof(argv[8]);
	double cohesion= (double) atof(argv[9]);
	// double RMC_shape_k= (double) atof(argv[10]);
	double dilation_angle_eta= (double) atof(argv[11]);
	double diletion_scale= (double) atof(argv[12]);

	ofstream outfile;
	outfile.open("strain_stress.txt");

	int material_tag{1}  ;

	double K = bulk_modulus; //16750;
	double p_ratio = 0.15; 
	double G =  3*K*(1-2*p_ratio)/(2+2*p_ratio);  
	double E_in = 2 * G * (1 + p_ratio) ;
	double v_in = p_ratio ;
	double rho_in = 0.0 ;
	int NYS = 4 ;
	// vector<double> radius{
	// 	2.7, 2.74, 2.8
	// 	, 2.82, 2.85, 
 //     	2.9, 3.0, 3.1

	// };
	vector<double> eta{
		0.1, 0.2, 0.3, 0.5
		// , 0.7, 0.9, 1.1, 1.3
		// , 1.5, 1.7, 1.9, 
		// 2.0, 2.2, 2.5, 
  //    	2.9, 3.3, 3.6
	};
	vector<double> HardingPara{
		5500, 4000, 2700, 2400
		// , 1890, 1300, 915, 600
	};
	for(auto& item: HardingPara){
		item *= scale_hardening;
	}

	double pp0_in = initial_confinement ;    // 1E5   ;
	double pa_in = reference_pressure ;    // 1E5    ;
	double modn_in = modulus_n_in ;    // 0.7    ;
	double pc_in = cohesion ;    // 0.0    ;
	// double kk_in = RMC_shape_k ;    // 1.0    ;
	double deta_in = dilation_angle_eta ;    // 1.5  ;
	double scal_in = diletion_scale ;    // 0.2  ;

	auto theMaterial= new DruckerPrager_multi_yield_surface(
		material_tag,
		E_in,
		v_in,
		rho_in,
		pp0_in,
		pa_in,
		modn_in,
		pc_in,
		// kk_in,
		deta_in,
		scal_in,
		NYS,
		eta,
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