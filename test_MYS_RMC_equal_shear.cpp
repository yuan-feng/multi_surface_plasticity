#include "MultiYieldSurfaceMaterial.h"
#include "RoundedMohrCoulomb_multi_surface.h"
#include "vonMises_multi_surface.h"
#include <fstream>
#include <string>
using namespace std;


// /////////////////////////////////////////////////////////////////////////////
// Pure Shear Increment
// dsigma(0,1) = E(0,1,0,1)*depsilon(0,1) 
//             + E(0,1,1,0)*depsilon(1,0)
// because other components in dsigma are zeroes. 
double stress2strain(double stress, DTensor4 const& E){
    double strain01{stress};
    strain01 /= (E(0,1,0,1) + E(0,1,1,0));
    return strain01;
}

int main(int argc, char const *argv[])
{
	double bulk_modulus = (double) atof(argv[1]);
	double scale_hardening = (double) atof(argv[2]);
	double max_stress_in = (double) atof(argv[3]);
	double strain_incr = (double) atof(argv[4]);
	int Nloop   = atoi(argv[5]);
	double initial_confinement = (double) atof(argv[6]);
	double reference_pressure= (double) atof(argv[7]);
	double modulus_n_in= (double) atof(argv[8]);
	double cohesion= (double) atof(argv[9]);
	double RMC_shape_k= (double) atof(argv[10]);
	double dilation_angle_eta= (double) atof(argv[11]);
	double diletion_scale= (double) atof(argv[12]);
	
	ofstream outfile;
	outfile.open("strain_stress.txt");

	int material_tag{1}  ;

	double K =  bulk_modulus ;
	double p_ratio = 0.15; 
	double G =  3*K*(1-2*p_ratio)/(2+2*p_ratio);  
	double E_in = 2 * G * (1 + p_ratio) ;
	double v_in = p_ratio ;
	double rho_in = 0.0 ;
	int NYS = 15 ;
	vector<double> eta{
		0.1, 0.2, 0.3, 0.5, 0.7,
		0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 
		2.0, 2.2, 2.5, 
     	2.9, 3.3, 3.6
	};
	vector<double> HardingPara{
		5500, 4000, 2700, 2400, 1890, 1300, 
        915, 600, 254, 167, 79,  65.4, 
        23, 2.2,  1.2
	};
	for(auto& item: HardingPara){
		item *= scale_hardening;
	}

	double pp0_in = initial_confinement ;    // 1E5   ;
	double pa_in = reference_pressure ;    // 1E5    ;
	double modn_in = modulus_n_in ;    // 0.7    ;
	double pc_in = cohesion ;    // 0.0    ;
	double kk_in = RMC_shape_k ;    // 1.0    ;
	double deta_in = dilation_angle_eta ;    // 1.5  ;
	double scal_in = diletion_scale ;    // 0.2  ;

	auto theMaterial= new RoundedMohrCoulomb_multi_surface(
		material_tag,
		E_in,
		v_in,
		rho_in,
		pp0_in,
		pa_in,
		modn_in,
		pc_in,
		kk_in,
		deta_in,
		scal_in,
		NYS,
		eta,
		HardingPara 
	);
	DTensor2 input_strain(3,3,0.) ;

	auto tStress = theMaterial->getStressTensor();
	auto tStrain = theMaterial->getStrainTensor();
	cout << "Stress and Strain State before Constitutive Testing\n";
	fprintf(stdout, "stress = \n");
	fprintf(stdout, "[%16.8f \t %16.8f \t %16.8f \n",   tStress(0, 0), tStress(0, 1), tStress(0, 2));
	fprintf(stdout, " %16.8f \t %16.8f \t %16.8f \n",   tStress(1, 0), tStress(1, 1), tStress(1, 2));
	fprintf(stdout, " %16.8f \t %16.8f \t %16.8f] \n",  tStress(2, 0), tStress(2, 1), tStress(2, 2));
	fprintf(stdout, "strain = \n");
	fprintf(stdout, "[%16.8f \t %16.8f \t %16.8f \n",   tStrain(0, 0), tStrain(0, 1), tStrain(0, 2));
	fprintf(stdout, " %16.8f \t %16.8f \t %16.8f \n",   tStrain(1, 0), tStrain(1, 1), tStrain(1, 2));
	fprintf(stdout, " %16.8f \t %16.8f \t %16.8f] \n",  tStrain(2, 0), tStrain(2, 1), tStrain(2, 2));

	// Loading 
	double max_stress = max_stress_in ;
	// double target_stress_incr = incr_stress_in ; 
	// int Nsteps = (int) (max_stress / target_stress_incr);
	auto tStiff  = theMaterial->getTangentTensor();
	auto stress_ret = theMaterial->getStressTensor();
	auto strain_ret = theMaterial->getStrainTensor();
	auto Nactive = theMaterial->getNumActiveYS();
	auto curStress = theMaterial->getStressTensor();
	double current_stress = curStress(0,1);

	double target_stress_range = 1E3 ;
	while(current_stress > (max_stress + target_stress_range) || 
			current_stress < (max_stress - target_stress_range))
	{
		input_strain *= 0.;
		input_strain(0,1) = strain_incr;
		input_strain(1,0) = strain_incr;
		curStress = theMaterial->getStressTensor();
		current_stress = curStress(0,1);
		theMaterial->revertToLastCommit();
	    theMaterial->setTrialStrainIncr(input_strain);
	    stress_ret = theMaterial->getStressTensor();
	    strain_ret = theMaterial->getStrainTensor();
	    Nactive = theMaterial->getNumActiveYS();
		theMaterial->commitState();

		double confinement = -1./3. * (stress_ret(0,0)+stress_ret(1,1)+stress_ret(2,2));
		outfile << strain_ret(0,1)  <<"\t" 
				<< stress_ret(0,1) <<"\t"
				<< Nactive << "\t"
				<< confinement << endl ;
	}

// // 	// unloading 
	for (int ll = 0; ll < Nloop; ++ll)
	{
		while( current_stress < -(max_stress + target_stress_range) || 
				current_stress > -(max_stress - target_stress_range))
		{
			input_strain *= 0.;
			input_strain(0,1) = - strain_incr;
			input_strain(1,0) = - strain_incr;
			curStress = theMaterial->getStressTensor();
			current_stress = curStress(0,1);
			theMaterial->revertToLastCommit();
		    theMaterial->setTrialStrainIncr(input_strain);
		    stress_ret = theMaterial->getStressTensor();
		    strain_ret = theMaterial->getStrainTensor();
		    Nactive = theMaterial->getNumActiveYS();
			theMaterial->commitState();

			double confinement = -1./3. * (stress_ret(0,0)+stress_ret(1,1)+stress_ret(2,2));
			outfile << strain_ret(0,1)  <<"\t" 
					<< stress_ret(0,1) <<"\t"
					<< Nactive << "\t"
					<< confinement << endl ;
		}

		// // reloading 
		while(current_stress > (max_stress + target_stress_range) || 
				current_stress < (max_stress - target_stress_range))
		{
			input_strain *= 0.;
			input_strain(0,1) = strain_incr;
			input_strain(1,0) = strain_incr;
			curStress = theMaterial->getStressTensor();
			current_stress = curStress(0,1);
			theMaterial->revertToLastCommit();
		    theMaterial->setTrialStrainIncr(input_strain);
		    stress_ret = theMaterial->getStressTensor();
		    strain_ret = theMaterial->getStrainTensor();
		    Nactive = theMaterial->getNumActiveYS();
			theMaterial->commitState();

			double confinement = -1./3. * (stress_ret(0,0)+stress_ret(1,1)+stress_ret(2,2));
			outfile << strain_ret(0,1)  <<"\t" 
					<< stress_ret(0,1) <<"\t"
					<< Nactive << "\t"
					<< confinement << endl ;
		}
	}


	// Done the experiments and Clean
	delete theMaterial ;
	outfile.close();

	return 0;
}