///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:
// CLASS:
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:
// PROGRAMMER:        Yuan Feng
// DATE:              Fri Mar 3 12:55:59 PST 2017
// UPDATE HISTORY:
//
//
/////////////////////////////////////////////////////////////////////////////


// #include <nDarray.h>
#include "MultiYieldSurfaceMaterial.h"
// #include <Information.h>
// #include <OPS_Globals.h>
//#include <DTensor2.h>
//#include <Vector.h>
//#include <stresst.h>
//#include <straint.h>
// #include <MaterialResponse.h>
#include <iostream>
// #include "../../ltensor/LTensor.h"
using namespace std;

#define BRENT_MAXITER 20
#define BRENT_TOLERANCE 1e-6

const DTensor2 MultiYieldSurfaceMaterial::ZeroStrain( 3, 3, 0. );
const DTensor2 MultiYieldSurfaceMaterial::ZeroStress( 3, 3, 0. );
const DTensor2 MultiYieldSurfaceMaterial::kronecker_delta(3, 3, "identity");
DTensor4 MultiYieldSurfaceMaterial::Ee( 3, 3, 3, 3, 0. );
DTensor4 MultiYieldSurfaceMaterial::Eep( 3, 3, 3, 3, 0. );
// DTensor2 vonMises_multi_surface::D(6, 6); 


DTensor2 MultiYieldSurfaceMaterial::errMatrix      (3, 3, 0.0);
DTensor1 MultiYieldSurfaceMaterial::errVector      (3,    0.0);
DTensor2 MultiYieldSurfaceMaterial::errTensor      (3, 3, 0.0);
DTensor4 MultiYieldSurfaceMaterial::errTensor4     (3, 3, 3, 3, 0.0);
DTensor2 MultiYieldSurfaceMaterial::errTensor2     (3, 3, 0.0);
DTensor2 MultiYieldSurfaceMaterial::errstresstensor(3, 3, 0.0);
DTensor2 MultiYieldSurfaceMaterial::errstraintensor(3, 3, 0.0);
// Vector  MultiYieldSurfaceMaterial::errVectorVector(0.0);
MultiYieldSurfaceMaterial_Constitutive_Integration_Method MultiYieldSurfaceMaterial::constitutive_integration_method(MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Not_Set);     //
double MultiYieldSurfaceMaterial::f_relative_tol(-1);
double MultiYieldSurfaceMaterial::stress_relative_tol(-1);
double MultiYieldSurfaceMaterial::allowed_subincrement(0.01);
int MultiYieldSurfaceMaterial::n_max_iterations(0);
// vector<double> MultiYieldSurfaceMaterial::errvec(1,0.);
MultiYieldSurfaceMaterial::MultiYieldSurfaceMaterial(
        int tag, 
        int classTag,
        double E_in,
        double v_in,
        double rho_in,
        int TNYS_in,
        vector<double> const& radius_in ,
        vector<double> const& HardingPara_in
        )
    : 
      // NDMaterialLT(tag, classTag),
      iterate_stress( 3, 3, 0. ),
      iterate_strain( 3, 3, 0. ),
      iterate_plastic_strain( 3, 3, 0. ),
      iterate_alpha_vec(TNYS_in+1),
      iterate_N_active(0),
      converge_commit_stress( 3, 3, 0. ),
      converge_commit_strain( 3, 3, 0. ),
      converge_commit_plastic_strain( 3, 3, 0. ),
      converge_commit_alpha_vec(TNYS_in+1),
      converge_commit_N_active(0),
      save_iter_stress( 3, 3, 0.),
      save_iter_strain( 3, 3, 0.),
      save_iter_plastic_strain( 3, 3, 0.),
      save_iter_alpha_vec(TNYS_in+1),
      save_iter_N_active(0),
      E( E_in ),
      v( v_in ),
      rho( rho_in ),
      TNYS(TNYS_in)
{
    // make yield_size and HardingPara 1-based instead of 0-based.
    yield_size = radius_in ;
    if(yield_size.size() == (uint) TNYS) { yield_size.insert(yield_size.begin(), 0) ;     }
    HardingPara = HardingPara_in;
    if(HardingPara.size() == (uint) TNYS) { HardingPara.insert(HardingPara.begin(), E/2./(1+v)) ;     }

    initial_E = E_in;

    for(auto& item: iterate_alpha_vec){item = ZeroStress;}
    for(auto& item: converge_commit_alpha_vec){item = ZeroStress;}
    for(auto& item: save_iter_alpha_vec){item  = ZeroStress ;}

    revertToStart();
    Eep = Ee;
}

MultiYieldSurfaceMaterial::MultiYieldSurfaceMaterial()
    : 
      // NDMaterialLT(0, 0),
      iterate_stress( 3, 3, 0. ),
      iterate_strain( 3, 3, 0. ),
      iterate_plastic_strain( 3, 3, 0. ),
      iterate_alpha_vec(0),
      iterate_N_active(0),
      converge_commit_stress( 3, 3, 0. ),
      converge_commit_strain( 3, 3, 0. ),
      converge_commit_plastic_strain( 3, 3, 0. ),
      converge_commit_alpha_vec(0),
      converge_commit_N_active(0),
      save_iter_stress( 3, 3, 0.),
      save_iter_strain( 3, 3, 0.),
      save_iter_plastic_strain( 3, 3, 0.),
      save_iter_alpha_vec(0),
      save_iter_N_active(0),
      E(0.),
      v(0.),
      rho(0.),
      TNYS(0)
{
    revertToStart();
    Eep = Ee;
}

MultiYieldSurfaceMaterial::~MultiYieldSurfaceMaterial()
{

}


//================================================================================
// Set the Total Strain, usually call by element
//================================================================================
int MultiYieldSurfaceMaterial::setTrialStrain( const DTensor2 &strain )
{
    DTensor2 strain_incre_(3, 3, 0.);
    strain_incre_ *= 0;
    iterate_strain(i, j) = strain(i,j);
    strain_incre_(i, j) = strain(i, j) - converge_commit_strain(i, j);

    return setTrialStrainIncr(strain_incre_);
}

// ================================================================================
// Set the Increment Strain, usually call by element
// ================================================================================
int MultiYieldSurfaceMaterial::setTrialStrainIncr( DTensor2 const& strain_increment ){
    iterate_strain(i, j) = converge_commit_strain(i, j) + strain_increment(i, j);
    return compute_stress(strain_increment);
}

// ================================================================================
// Compute the Stress from the Increment Strain
// ================================================================================
int MultiYieldSurfaceMaterial::compute_stress(DTensor2 const& strain_incr, int Nsubsteps  ){
    if (iterate_N_active > TNYS ){
        cerr<< "MultiYieldSurfaceMaterial::compute_stress " <<endl;
        cerr<< "Exceed the length of alpha_vec " <<endl;
        cerr<< "iterate_N_active " << iterate_N_active <<endl;
        cerr<< "alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "Total NYS" << TNYS <<endl;
    }
    
    for (int isub = 0; isub < Nsubsteps; ++isub)
    {
        DTensor2 strain_increment(3,3,0.);
        strain_increment(i,j) = strain_incr(i,j)/Nsubsteps;

        saveIterateState();
        // For unloading situation, always check the initial stiffness first.
        if(isub == 0){
            update_modulus(0, converge_commit_stress);
        }
        DTensor2 PredictStress(3,3,0.);
        PredictStress(i,j) = converge_commit_stress(i,j) + Ee(i,j,k,l) * strain_increment(k,l);
        double curr_radius = yield_size[iterate_N_active];
        DTensor2 curr_alpha = iterate_alpha_vec[iterate_N_active];
        if(iterate_N_active == 0){
            curr_radius = yield_size[1];
            curr_alpha = iterate_alpha_vec[1];
        }
        double curr_yf_val = yield_surface_val(PredictStress, curr_alpha, curr_radius) ; 
        iterate_stress(i,j) = PredictStress(i,j);

        // // // Debug Printing
        // cout << " strain_increment(0,1) : " << strain_increment(0,1) <<endl;
        // cout << " Ee(0,1,0,1)           : " << Ee(0,1,0,1) <<endl;
        // cout << " converge_commit_stress(0,1)    : " << converge_commit_stress(0,1) <<endl;
        // cout << " iterate_stress(0,1)   : " << iterate_stress(0,1) <<endl;
        // cout << " curr_alpha(0,1)       : " << curr_alpha(0,1) <<endl;

        double pred_yf_val{curr_yf_val};

        int iter = 0; 
        int Max_iter = 20;
        double relative_TOL = 1E-6;
        double absolute_TOL = 1E-4;
        // ============================================================================
        // This part is like the flowchart of the multi-yield-surface algorithm
        // ============================================================================
        double lambda{0.}, lambda1{0.}, lambda2{0.};
        // I. Elastic Situation:
        if (curr_yf_val < absolute_TOL)
        {
            // iterate_N_active = 0 ;
            // update_modulus(iterate_N_active);
            return 0;
        }
        // II. Plastic Situation:
        double next_yf_val{0.};

        // int outer_iter{0}, outer_Max_iter{10000};
        // while(outer_iter< outer_Max_iter)   // Actually, can be while (true) in logic, just to avoid potential halt.
        while(true)               
        {
            // Check "next_yf_val<TOL" is important. For overshooting next larger YS multiple times.
            if ( ( curr_yf_val < absolute_TOL || abs(curr_yf_val/pred_yf_val) < relative_TOL  )  
                && next_yf_val< absolute_TOL )
            {
                break;
            }
            else
            {
                if( iterate_N_active == 0 ) iterate_N_active = 1; 
                // Step 1. Return to current yield surface
                iter = 0 ;
                while( (abs(curr_yf_val/pred_yf_val) > relative_TOL && abs(curr_yf_val) > absolute_TOL ) 
                         && iter<Max_iter)
                {
                    if(iter == 0 && isub == (Nsubsteps-1)  )
                    // if( isub == (Nsubsteps-1)  )
                    {
                        double intersection_factor = zbrentstress(iterate_N_active, converge_commit_stress, PredictStress, 0.0, 1.0, BRENT_TOLERANCE );
                        // cout<< "intersection_factor = " << intersection_factor <<endl;
                        DTensor2 intersection_stress(3,3,0.);
                        intersection_stress(i, j) = converge_commit_stress(i, j) * (1 - intersection_factor) + PredictStress(i, j) * intersection_factor;

                        update_modulus(iterate_N_active, intersection_stress);
                        compute_elastoplastic_tangent(iterate_N_active, intersection_stress);
                        if(intersection_factor<0 || intersection_factor> 1)
                        {
                            compute_elastoplastic_tangent(iterate_N_active, converge_commit_stress);
                        }
                    }
                    // else
                    // {
                    update_stress(iterate_stress, lambda, iterate_N_active, iterate_stress);
                    // }

                    update_current_yield_surface(curr_yf_val, iterate_N_active, lambda, iterate_stress);
                    update_inner_yield_surfaces(iterate_N_active, iterate_stress);
                    iter++;
                }
                update_plastic_strain(iterate_plastic_strain, lambda, iterate_N_active, iterate_stress);
                
                // if( iter >= Max_iter)
                // {
                //     cerr<<"debug printing MultiYieldSurfaceMaterial::compute_stress() exceeds the max number of iteration."<<endl;
                //     cerr<<"debug printing iterate_N_active: " << iterate_N_active <<endl;
                //     cerr<<"debug printing curr_yf_val     : " << curr_yf_val <<endl;
                //     cerr<<"debug printing pred_yf_val     : " << pred_yf_val <<endl;
                // }

                // If not converge, Do subincrements. 
                if( iter >= Max_iter ){
                    backToLastIterateState();
                    if (compute_stress( strain_increment, Nsubsteps*10 ) == 0 && Nsubsteps < 10){
                        break;
                    }
                    else
                    // If subincrements still do not converge. Do elastic with the final elastic modulus.
                    {
                        iterate_N_active = TNYS;
                        update_modulus(iterate_N_active, converge_commit_stress);
                        iterate_stress(i,j) = converge_commit_stress(i,j) + Ee(i,j,k,l) * strain_increment(k,l);
                        update_failure_surface(iterate_stress);
                        compute_elastoplastic_tangent( iterate_N_active, iterate_stress, true);
                        saveIterateState();
                        return 0;
                    }
                }

                // Step 2. Check the next larger yield surface
                // If this is already the greatest yield surface. Keep going with the final elastic modulus.
                if( iterate_N_active  >= TNYS){
                    saveIterateState();
                    break;
                }
                else
                {
                    double next_radius = yield_size[iterate_N_active+1];
                    DTensor2 next_alpha = iterate_alpha_vec[iterate_N_active+1];
                    next_yf_val = yield_surface_val(iterate_stress, next_alpha, next_radius) ; 

                    // If NO overshooting the next larger yield surface
                    if( next_yf_val <= 0 )
                    {
                        break;
                    }
                    // If overshoot the next larger yield surface
                    else
                    {
                        // cout<< "next_yf_val before correct_update " << next_yf_val <<endl;
                        // cout<< "iterate_stress(0,1) before_correct: " << iterate_stress(0,1) <<endl;

                        correct_update_stress(iterate_stress, lambda1, lambda2, iterate_N_active );
                        correct_update_plastic_strain(iterate_plastic_strain, lambda1, lambda2, iterate_N_active, iterate_stress);

                        // cout<< "iterate_stress(0,1) after _correct: " << iterate_stress(0,1) <<endl;
                        // cout<< "MultiYieldSurfaceMaterial::compute_stress() \t";
                        // cout<< "iterate_N_active "<< iterate_N_active <<endl;

                        iterate_N_active++ ;
                        lambda = lambda2;
                        update_current_yield_surface(curr_yf_val, iterate_N_active, lambda, iterate_stress);
                        update_inner_yield_surfaces(iterate_N_active, iterate_stress);
                        update_modulus(iterate_N_active, iterate_stress);
                        // compute_elastoplastic_tangent(iterate_N_active, iterate_stress);

                        if( iterate_N_active  >= TNYS){
                            next_yf_val = 0. ;
                            continue;
                        }

                        // If overshoot the second larger yield surface, go back to while-loop.
                        next_radius = yield_size[iterate_N_active+1];
                        next_alpha = iterate_alpha_vec[iterate_N_active+1];
                        next_yf_val = yield_surface_val(iterate_stress, next_alpha, next_radius) ;

                        // cout<< "next_yf_val after correct_update " << next_yf_val <<endl;

                    }
                }
                update_modulus(iterate_N_active, iterate_stress);
            } // End if-else
            iter++; 
        } // End while-loop
        saveIterateState();
    }
    // ============================================================================
    // cout<< "N_active = " << N_active <<endl;
    return 0;

}

// ================================================================================
// After the yield, Update the Stress 
// ================================================================================
void MultiYieldSurfaceMaterial::update_stress(DTensor2& target_stress, double& lambda, int N_active_ys, DTensor2 const& normal_refer_stress){
    // df_dsigma 
    DTensor2 curr_nn = df_dsigma(N_active_ys, normal_refer_stress);

    // df_dalpha
    DTensor2 curr_xi = df_dalpha(N_active_ys, normal_refer_stress);

    // alpha_direction
    DTensor2 bar_alpha = alpha_bar(N_active_ys, normal_refer_stress);

    // yield_surface_val
    double curr_radius = yield_size[N_active_ys];
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    double yf_val = yield_surface_val(target_stress, curr_alpha, curr_radius);

    // get lambda
    double denominator = curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) - curr_xi(o,t) * bar_alpha(o,t);
    if (denominator == 0 ){
        cerr<< "MultiYieldSurfaceMaterial::update_stress()" <<endl;
        cerr<< "Error denominator == 0 " <<endl;
        cerr<< "LEFT  = curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) == " << curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) <<endl;
        cerr<< "RIGHT = curr_xi(o,t) * bar_alpha(o,t)             == " << curr_xi(o,t) * bar_alpha(o,t) <<endl;
    }
    lambda = yf_val / denominator ; 
    target_stress(i,j) = target_stress(i,j) - lambda * Ee(i,j,k,l) * curr_nn(k,l);
}

// ================================================================================
// After the yield, Update the current active yield surface (alpha)
// ================================================================================
void MultiYieldSurfaceMaterial::update_current_yield_surface(double& the_yf_val, int N_active_ys, double lambda, DTensor2 const& stress){
    // df_dsigma 
    DTensor2 curr_nn = df_dsigma(N_active_ys, stress);

    // rate of alpha
    DTensor2 bar_alpha = alpha_bar(N_active_ys, stress);
    // update the alpha
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys] ; 
    curr_alpha(i,j) = curr_alpha(i,j) + lambda * bar_alpha(i,j) ;
    iterate_alpha_vec[N_active_ys] = curr_alpha ; 

    // update the yield surface value 
    double curr_radius = yield_size[N_active_ys] ;
    the_yf_val = yield_surface_val(stress, curr_alpha, curr_radius);


    // // // Debug Printing
    // cout<< " update_current_yield_surface \n" ;
    // cout<< " Yield_surface Number"  << N_active_ys << "\t" ;
    // cout<< " Yield_surface Value "  << yield_surface_val(stress, curr_alpha, curr_radius);
    // cout<< " \n+++++++++++++++++++++++++++" <<endl;

}



// ================================================================================
// After the yield, Update the inner yield surfaces (alpha) 
// ================================================================================
void MultiYieldSurfaceMaterial::update_inner_yield_surfaces(int N_active_ys, DTensor2 const& stress){

    for (int inner = 0; inner < N_active_ys; ++inner)
    {
        double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
        DTensor2 DevStress(3,3,0.);
        DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
        DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys] ;
        double curr_radius = yield_size[N_active_ys] ;
        DTensor2 the_alpha = iterate_alpha_vec[inner] ;
        double the_radius = yield_size[inner] ;
        if(curr_radius == 0){
            cerr<< "MultiYieldSurfaceMaterial::update_inner_yield_surfaces() " <<endl;
            cerr<< "Error Denominator -- curr_radius == 0 " <<endl;
        }
        the_alpha(i,j) = DevStress(i,j) - the_radius / curr_radius * (DevStress(i,j) - curr_alpha(i,j));
        iterate_alpha_vec[inner] = the_alpha ; 
    }

    // // debug printing
    // for (int inner = 0; inner < N_active_ys; ++inner)
    // {
    //     cout<< " Yield_surface Number"  << inner + 1 << "\t" ;
    //     cout<< " Yield_surface Value "  << yield_surface_val(iterate_stress, iterate_alpha_vec[inner], yield_size[inner]);
    //     cout<< " \n";
    // }
    // cout<< " \n" <<endl;
}


void MultiYieldSurfaceMaterial::update_failure_surface(DTensor2 const& stress){
    // 
    DTensor2 last_alpha = iterate_alpha_vec[TNYS] ; 
    DTensor2 curr_nn = df_dsigma(TNYS, stress);
    double pp = - 1./3. * (stress(0,0) + stress(1,1) + stress(2,2));

    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
    // direction
    curr_nn(i,j) = curr_nn(i,j) / sqrt(curr_nn(k,l) * curr_nn(k,l)) ; 
    // magnitude
    double curr_radius = yield_size[TNYS] ;
    last_alpha(i,j) = DevStress(i,j) - sqrt(2./3.) * curr_radius * curr_nn(i,j) ; 

    iterate_alpha_vec[TNYS] = last_alpha ; 

    update_inner_yield_surfaces(TNYS, stress); 

    // double yf_val = yield_surface_val(stress, last_alpha, curr_radius);
    // cout<< "failure surface YF value " << yf_val <<endl;
}


// ================================================================================
// After the yield, Update the plastic strain
// ================================================================================
void MultiYieldSurfaceMaterial::update_plastic_strain(DTensor2& pstrain, double lambda, int N_active_ys, DTensor2 const& stress ){
    // df_dsigma and associate flow 
    DTensor2 curr_nn = df_dsigma(N_active_ys, stress);

    // lambda * plastic_flow
    pstrain(i,j) = pstrain(i,j) + lambda * curr_nn(i,j) ; 
}

// ================================================================================
// After the overshooting, Correct the overshooting stress
// ================================================================================
void MultiYieldSurfaceMaterial::correct_update_stress(DTensor2& stress, double& lambda1, double& lambda2, int N_active_ys){
    // df_dalpha
    DTensor2 curr_xi = df_dalpha(N_active_ys, stress);

    // alpha_direction
    DTensor2 bar_alpha = alpha_bar(N_active_ys, stress);

    // curr_H_prime
    double next_radius = yield_size[N_active_ys + 1];
    DTensor2 next_alpha = iterate_alpha_vec[N_active_ys + 1];
    double next_yf_val = yield_surface_val(stress, next_alpha, next_radius) ; 
    double curr_H_prime = - curr_xi(i,j) * bar_alpha(i,j) ; 

    // next_bar_alpha
    DTensor2 next_bar_alpha(3,3,0.);
    next_bar_alpha = alpha_bar(N_active_ys + 1,  stress);

    // next_H_prime
    DTensor2 next_xi = df_dalpha(N_active_ys + 1,  stress);
    double next_H_prime = - next_xi(i,j) * next_bar_alpha(i,j) ; 

    // curr_nn
    DTensor2 curr_nn = df_dsigma(N_active_ys,  stress);

    // next_H0
    DTensor2 next_nn = df_dsigma(N_active_ys + 1,  stress);
    double next_H0 = next_nn(i,j) * Ee(i,j,k,l) * next_nn(k,l) ; 

    // lambda1 
    if(curr_H_prime == 0){
        cerr<< "MultiYieldSurfaceMaterial::correct_update_stress() " <<endl;
        cerr<< "Error Denominator -- curr_H_prime == 0 " <<endl;
    }
    lambda1 = next_yf_val / curr_H_prime ; 

    // numerator for lambda2
    double numerator = next_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) ;
    // lambda2
    if( (next_H0 + next_H_prime) == 0){
        cerr<< "MultiYieldSurfaceMaterial::correct_update_stress() " <<endl;
        cerr<< "Error Denominator -- (next_H0 + next_H_prime) == 0 " <<endl;
    }
    lambda2 = next_yf_val * (1 + numerator/curr_H_prime) / (next_H0 + next_H_prime)  ;

    stress(i,j) = stress(i,j) 
                        + lambda1 * Ee(i,j,k,l) * curr_nn(k,l) 
                        - lambda2 * Ee(i,j,k,l) * next_nn(k,l) ;
}

// ================================================================================
// After the overshooting, Correct the plastic strain
// ================================================================================
void MultiYieldSurfaceMaterial::correct_update_plastic_strain(DTensor2& pstrain, double lambda1, double lambda2, int N_active_ys, DTensor2 const& stress ){
    // curr_nn
    DTensor2 curr_nn = df_dsigma(N_active_ys,  stress);

    // next_nn
    DTensor2 next_nn = df_dsigma(N_active_ys + 1,  stress);

    // lambda * plastic_flow
    pstrain(i,j) = pstrain(i,j) 
                                - lambda1 * curr_nn(i,j) 
                                + lambda2 * next_nn(i,j); 

}







//================================================================================
// Commit/Store the Trial state at the end of converge state.
// This function should be called by the element class and convergence-check class.
// This material class itself should never call it.
//================================================================================
int MultiYieldSurfaceMaterial::commitState( void ){
    converge_commit_stress          = iterate_stress;
    converge_commit_strain          = iterate_strain;
    converge_commit_plastic_strain  = iterate_plastic_strain ;
    converge_commit_N_active        = iterate_N_active;
    converge_commit_alpha_vec       = iterate_alpha_vec;
    return 0;
}
    
//================================================================================
// Go back to last Commit/stored state.
//================================================================================
int MultiYieldSurfaceMaterial::revertToLastCommit( void ){
    iterate_stress         = converge_commit_stress ;
    iterate_strain         = converge_commit_strain ;
    iterate_plastic_strain = converge_commit_plastic_strain ;
    iterate_N_active       = converge_commit_N_active ;
    iterate_alpha_vec      = converge_commit_alpha_vec ;
    return 0;
}

//================================================================================
// Commit/store the Trial state.
//================================================================================
int MultiYieldSurfaceMaterial::saveIterateState( void ){
    save_iter_stress          = iterate_stress;
    save_iter_strain          = iterate_strain;
    save_iter_plastic_strain  = iterate_plastic_strain ;
    save_iter_N_active        = iterate_N_active;
    save_iter_alpha_vec       = iterate_alpha_vec;
    return 0;
}
    
//================================================================================
// Go back to last Commit/stored state.
//================================================================================
int MultiYieldSurfaceMaterial::backToLastIterateState( void ){
    iterate_stress         = save_iter_stress ;
    iterate_strain         = save_iter_strain ;
    iterate_plastic_strain = save_iter_plastic_strain ;
    iterate_N_active       = save_iter_N_active ;
    iterate_alpha_vec      = save_iter_alpha_vec ;
    return 0;
}

//================================================================================
// Go back to all zeroes.
//================================================================================
int MultiYieldSurfaceMaterial::revertToStart( void ){
    iterate_stress                  *= 0.  ;
    iterate_strain                  *= 0.  ;
    iterate_plastic_strain          *= 0.  ;
    iterate_N_active                 = 0;
    converge_commit_stress          *= 0.  ;
    converge_commit_strain          *= 0.  ;
    converge_commit_plastic_strain  *= 0.  ;
    converge_commit_N_active         = 0;
    save_iter_stress                *= 0.  ;
    save_iter_strain                *= 0.  ;
    save_iter_plastic_strain        *= 0.  ;
    save_iter_N_active               = 0;
    for(auto& item: iterate_alpha_vec){          item *= 0. ;}
    for(auto& item: converge_commit_alpha_vec){  item *= 0. ;}
    for(auto& item: save_iter_alpha_vec){        item *= 0. ;}

    return 0;
}

// //================================================================================
// // Message passing for parallel
// //================================================================================
// int MultiYieldSurfaceMaterial::sendSelf( int commitTag, Channel &theChannel )
// {
//     cerr<<"MultiYieldSurfaceMaterial::sendSelf() is not implemented yet! " <<endl;
//     ID idData(1);
//     Vector vectorData(6);
//     Matrix a(3, 3);

//     idData(0) = this->getTag();

//     if (theChannel.sendID(0, commitTag, idData) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send idData\n";
//         return -1;
//     }

//     vectorData(0) = converge_commit_N_active;
//     vectorData(1) = E;
//     vectorData(2) = v;
//     vectorData(3) = rho;
//     vectorData(4) = TNYS;
//     vectorData(5) = initial_E;

//     if (theChannel.sendVector(0, commitTag, vectorData) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send vectorData\n";
//         return -1;
//     }

//     Vector multi_ys_data( 2 * TNYS + 2 ); 
//     for (int it = 0; it < TNYS+1; ++it){
//         multi_ys_data(it) = yield_size[it];
//     }
//     for (int it = TNYS+1; it < 2*TNYS+2; ++it){
//         multi_ys_data(it) = HardingPara[it];
//     }
//     if (theChannel.sendVector(0, commitTag, multi_ys_data) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send multi_ys_data\n";
//         return -1;
//     }


//     a.setData(converge_commit_stress.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send converge_commit_stress\n";
//         return -1;
//     }

//     a.setData(converge_commit_strain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send converge_commit_strain\n";
//         return -1;
//     }

//     a.setData(converge_commit_plastic_strain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send converge_commit_plastic_strain\n";
//         return -1;
//     }
    
//     for (int it = 0; it < TNYS+1 ; ++i)
//     {
//         a.setData(converge_commit_alpha_vec[it].data, 3, 3);
//         if (theChannel.sendMatrix(0, 0, a) < 0)
//         {
//             cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send converge_commit_alpha_vec\n";
//             return -1;
//         }
//     }

//     return 0;
// }

// //================================================================================
// // Message passing for parallel
// //================================================================================
// int MultiYieldSurfaceMaterial::receiveSelf( int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker )
// {
//     cerr<<"MultiYieldSurfaceMaterial::receiveSelf() is not implemented yet! " <<endl;
//     ID idData(1);
//     Vector vectorData(6);
//     Matrix a(3, 3);

//     if (theChannel.receiveID(0, commitTag, idData) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive idData\n";
//         return -1;
//     }
//     this->setTag(idData(0));

//     if (theChannel.receiveVector(0, commitTag, vectorData) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive vectorData\n";
//         return -1;
//     }
//     converge_commit_N_active = vectorData(0) ;
//     E                        = vectorData(1) ;
//     v                        = vectorData(2) ;
//     rho                      = vectorData(3) ;
//     TNYS                     = vectorData(4) ;
//     initial_E                = vectorData(5) ;


//     Vector multi_ys_data( 2 * TNYS + 2 ); 
//     if (theChannel.receiveVector(0, commitTag, multi_ys_data) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send multi_ys_data\n";
//         return -1;
//     }
//     for (int it = 0; it < TNYS+1; ++it){
//         yield_size[it]  = multi_ys_data(it)  ;
//     }
//     for (int it = TNYS+1; it < 2*TNYS+2; ++it){
//         HardingPara[it] = multi_ys_data(it);
//     }



//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             converge_commit_stress(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             converge_commit_strain(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "MultiYieldSurfaceMaterial::receiveSelf -- could not receive Elastic Constant Tensor\n";
//         return -1;
//     }

    

//     for (int it = 0; it < TNYS+1 ; ++i)
//     {
//         if (theChannel.sendMatrix(0, 0, a) < 0)
//         {
//             cerr << "MultiYieldSurfaceMaterial::sendSelf -- could not send converge_commit_alpha_vec\n";
//             return -1;
//         }
//         for (int ii = 0; ii < 3; ii++){
//             for (int jj = 0; jj < 3; jj++){
//                 (converge_commit_alpha_vec[it])(ii, jj) = a(ii, jj);
//             }
//         }
//     }


//     update_modulus(converge_commit_N_active);
//     return 0;
// }



void MultiYieldSurfaceMaterial::zeroStrain()
{
    iterate_strain *= 0. ; 
    iterate_plastic_strain *= 0. ; 
    converge_commit_strain *= 0. ; 
    converge_commit_plastic_strain *= 0. ; 
    save_iter_strain *= 0. ;
    save_iter_plastic_strain *= 0. ;
}


bool MultiYieldSurfaceMaterial::set_constitutive_integration_method(int method, double f_relative_tol, double stress_relative_tol, int n_max_iterations, double allowed_subincrement)
{
    if ( method == (int) MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Not_Set
            || method == (int) MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Forward_Euler
            || method == (int) MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Backward_Euler
            || method == (int) MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Forward_Euler_Subincrement
            || method == (int) MultiYieldSurfaceMaterial_Constitutive_Integration_Method::Backward_Euler_Subincrement
    )
    {
        MultiYieldSurfaceMaterial::constitutive_integration_method = (MultiYieldSurfaceMaterial_Constitutive_Integration_Method) method ;
        MultiYieldSurfaceMaterial::f_relative_tol = f_relative_tol ;
        MultiYieldSurfaceMaterial::stress_relative_tol = stress_relative_tol ;
        MultiYieldSurfaceMaterial::n_max_iterations = n_max_iterations ;
        MultiYieldSurfaceMaterial::allowed_subincrement = allowed_subincrement ;

        cout << "Setting set_constitutive_integration_method = " << method << endl;

        return true;
    }
    else
    {
        cerr << "MultiYieldSurfaceMaterial::set_constitutive_integration_method - Unknown constitutive_integration_method\n";
        return false;
    }
}


//================================================================================
double MultiYieldSurfaceMaterial::zbrentstress(
                    int num_active_ys,
                    const DTensor2& start_stress,
                    const DTensor2& end_stress,
                    double x1, double x2, double tol) 
{
    // using namespace ClassicElastoplasticityGlobals;
    double EPS = numeric_limits<double>::epsilon();

    int iter;
    double a = x1;
    double b = x2;
    double c = 0.0;
    double d = 0.0;
    double e = 0.0;
    double min1 = 0.0;
    double min2 = 0.0;
    double fc = 0.0;
    double p = 0.0;
    double q = 0.0;
    double r = 0.0;
    double s = 0.0;
    double tol1 = 0.0;
    double xm = 0.0;

    // double fa = func(start_stress, end_stress, *ptr_material_parameter, a);
    // double fb = func(start_stress, end_stress, *ptr_material_parameter, b);

    static DTensor2 sigma_a(3, 3, 0.0);
    static DTensor2 sigma_b(3, 3, 0.0);

    sigma_a(i, j) = start_stress(i, j) * (1 - a)  + end_stress(i, j) * a;
    sigma_b(i, j) = start_stress(i, j) * (1 - b)  + end_stress(i, j) * b;

    DTensor2 curr_alpha = iterate_alpha_vec[num_active_ys];
    double curr_sz = yield_size[num_active_ys];
    double fa = yield_surface_val(sigma_a, curr_alpha, curr_sz);
    double fb = yield_surface_val(sigma_b, curr_alpha, curr_sz);

    // cout << "   brent fa = " << fa << " fb = " << fb << endl;


    if ( (fb * fa) > 0.0)
    {
        // cout<< " fa = " << fa << endl;
        // cout<< " fb = " << fb << endl;
        // std::cout << "\a\n Root must be bracketed in ZBRENTstress " << std::endl;
        if( fabs(fa) < 1E-3){
            fa = - fabs(fa);
        }
        // return 1. ;
        // exit(1);
    }

    fc = fb;

    for ( iter = 1; iter <= BRENT_MAXITER; iter++ )
    {
        if ( (fb * fc) > 0.0)
        {
            c = a;
            fc = fa;
            e = d = b - a;
        }

        if ( fabs(fc) < fabs(fb) )
        {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol1 = 2.0 * EPS * fabs(b) + 0.5 * tol;
        xm = 0.5 * (c - b);

        if ( fabs(xm) <= tol1 || fb == 0.0 )
        {
            return b;
        }

        if ( fabs(e) >= tol1 && fabs(fa) > fabs(fb) )
        {
            s = fb / fa;

            if (a == c)
            {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            }
            else
            {
                q = fa / fc;
                r = fb / fc;
                p = s * ( 2.0 * xm * q * (q - r) - (b - a) * (r - 1.0) );
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (p > 0.0)
            {
                q = -q;
            }

            p = fabs(p);
            min1 = 3.0 * xm * q - fabs(tol1 * q);
            min2 = fabs(e * q);

            if (2.0 * p < (min1 < min2 ? min1 : min2))
            {
                e = d;
                d = p / q;
            }
            else
            {
                d = xm;
                e = d;
            }
        }
        else
        {
            d = xm;
            e = d;
        }

        a = b;
        fa = fb;

        if (fabs(d) > tol1)
        {
            b += d;
        }
        else
        {
            b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
        }

        // fb = func(start_stress, end_stress, *ptr_material_parameter, b);
        sigma_b(i, j) = start_stress(i, j) * (1 - b)  + end_stress(i, j) * b;
        fb = yield_surface_val(sigma_b, curr_alpha, curr_sz);
    }

    return 0.0;
}

//================================================================================
// The Getter
//================================================================================
double MultiYieldSurfaceMaterial::getE() const {return E;}
double MultiYieldSurfaceMaterial::getv() const {return v;}
double MultiYieldSurfaceMaterial::getRho() const {return rho;}
double MultiYieldSurfaceMaterial::getNumActiveYS() const {return iterate_N_active;}
double MultiYieldSurfaceMaterial::getTNYS() const {return TNYS;}
vector<double> const& MultiYieldSurfaceMaterial::getYieldSize() const {return yield_size;}
vector<double> const& MultiYieldSurfaceMaterial::getHardPara() const {return HardingPara;}
const DTensor2  &MultiYieldSurfaceMaterial::getStressTensor() const {return iterate_stress;}
const DTensor2  &MultiYieldSurfaceMaterial::getStrainTensor() const {return iterate_strain;}
const DTensor2  &MultiYieldSurfaceMaterial::getPlasticStrainTensor() const {return iterate_plastic_strain;}










/************************************************************************************
* Yuan @ Thu Mar 2 21:35:05 PST 2017                                                *
* Following Methods are Subclass Responsibility.                                    *
************************************************************************************/

void MultiYieldSurfaceMaterial::Print( ostream &s, int flag )
{
    s << "MultiYieldSurfaceMaterial::" << endl;
    // s << "\t Tag      : " << this->getTag() << endl;
    s << "\t E        : " << this->getE() << endl;
    s << "\t v        : " << this->getv() << endl;
    s << "\t Rho      : " << this->getRho() << endl;
    s << "\t TNYS     : " << this->getTNYS() << endl;
    // s << "\t Radius   : " << this->getYieldSize() << endl;
    // s << "\t HardPara : " << this->getHardPara() << endl;
}


void MultiYieldSurfaceMaterial::update_modulus(int num_active_ys, DTensor2 const& stress){
    cerr << "MultiYieldSurfaceMaterial::update_modulus -> Subclass responsability" << endl;
}

double MultiYieldSurfaceMaterial::yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double radius){
    cerr << "MultiYieldSurfaceMaterial::yield_surface_val -> Subclass responsability" << endl;
    return 0. ;
}

DTensor2 MultiYieldSurfaceMaterial::df_dsigma(int num_active_ys, DTensor2 const& stress){
    cerr << "MultiYieldSurfaceMaterial::df_dsigma -> Subclass responsability" << endl;
    return ZeroStress;
}

DTensor2 MultiYieldSurfaceMaterial::df_dalpha(int num_active_ys, DTensor2 const& stress){
    cerr << "MultiYieldSurfaceMaterial::df_dalpha -> Subclass responsability" << endl;
    return ZeroStress;
}

DTensor2 MultiYieldSurfaceMaterial::alpha_bar(int num_active_ys, DTensor2 const& stress){
    cerr << "MultiYieldSurfaceMaterial::alpha_bar -> Subclass responsability" << endl;
    return ZeroStress;
}

DTensor2 MultiYieldSurfaceMaterial::plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys){
    cerr << "MultiYieldSurfaceMaterial::plastic_flow_direct -> Subclass responsability" << endl;
    return ZeroStress;
}


MultiYieldSurfaceMaterial *MultiYieldSurfaceMaterial::getCopy(void)
{
    cerr << "MultiYieldSurfaceMaterial::getCopy -> Subclass responsability" << endl;
    return 0;
}

    
DTensor4 const& MultiYieldSurfaceMaterial::getTangentTensor(){
    cerr << "MultiYieldSurfaceMaterial::getTangentTensor -> Subclass responsability" << endl;
    return Ee;
}
void MultiYieldSurfaceMaterial::compute_elastoplastic_tangent(int N_active_ys, DTensor2 const& intersection_stress , bool elastic){
    cerr << "MultiYieldSurfaceMaterial::compute_elastoplastic_tangent -> Subclass responsability" << endl;
}


const char *MultiYieldSurfaceMaterial::getType(void) const
{
    cerr << "MultiYieldSurfaceMaterial::getType -> Subclass responsability" << endl;
    return 0;
}

