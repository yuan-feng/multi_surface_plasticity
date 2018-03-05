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

#include "RoundedMohrCoulomb_multi_surface.h"
// #include <Channel.h>
// #include <HDF5_Channel.h>
// #include <Matrix.h>
// #include <ID.h>
// #include <Vector.h>
#include <vector>

//================================================================================
// Constructor
//================================================================================
RoundedMohrCoulomb_multi_surface::RoundedMohrCoulomb_multi_surface( 
        int tag,
        double E_in,
        double v_in,
        double rho_in,
        double pp0_in,
        double pa_in,
        double modulus_n_in,
        double pc_in,
        double kk_in,
        double deta_in,
        double scal_in,
        int TNYS_in,
        vector<double> const& yield_surface_size_in ,
        vector<double> const& HardingPara_in
        )
    : 
    MultiYieldSurfaceMaterial( 
            tag, 
            0, //ND_TAG_RoundedMohrCoulomb_multi_surface, 
            E_in, 
            v_in, 
            rho_in, 
            TNYS_in, 
            yield_surface_size_in, 
            HardingPara_in),
      pp0( pp0_in ),
      pa ( pa_in ),
      modn ( modulus_n_in ),
      pc ( pc_in ),
      kk ( kk_in ),
      deta ( deta_in ),
      scal ( scal_in )
{
    for (int it = 0; it < 3; ++it){
        iterate_stress(it,it) = - pp0;
        converge_commit_stress(it,it) = - pp0;
        save_iter_stress(it,it) = - pp0;
    }
    update_modulus(0, iterate_stress);
    Eep = Ee;
}

//================================================================================
// Empty Constructor for parallel
//================================================================================
RoundedMohrCoulomb_multi_surface::RoundedMohrCoulomb_multi_surface( )
    : 
    MultiYieldSurfaceMaterial(),
    pp0(0. ),
    pa (0. ),
    modn (0. ),
    pc (0. ),
    kk ( 0. ),
    deta (0. ),
    scal (0. )
{
    // revertToStart();
    for (int it = 0; it < 3; ++it){
        iterate_stress(it,it) = - pp0;
        converge_commit_stress(it,it) = - pp0;
        save_iter_stress(it,it) = - pp0;
    }
    update_modulus(0, iterate_stress);
    Eep = Ee;

}


//================================================================================
// Destructor
//================================================================================
RoundedMohrCoulomb_multi_surface::~RoundedMohrCoulomb_multi_surface()
{

}


// Next Lower Level Equations
// ================================================================================
// Return the yield surface Value
// ================================================================================
double RoundedMohrCoulomb_multi_surface::yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double yield_sz){


    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 s(3,3,0.);
    s(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;

    s(i,j) -= pp * alpha(i,j) ;

    return sqrt( s(i,j) * s(i,j) ) - Rtheta(stress, alpha) * sqrt(2./3.) * yield_sz * (pp - pc);
}

// ================================================================================
// Return the normal to the yield surface w.r.t stress
// ================================================================================
DTensor2 RoundedMohrCoulomb_multi_surface::df_dsigma(int N_active_ys, DTensor2 const& stress){
    if (N_active_ys > TNYS )
    {
        cerr<< "RoundedMohrCoulomb_multi_surface::df_dsigma " <<endl;
        cerr<< "Exceed the length of iterate_alpha_vec " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
        cerr<< "Total NYS" << TNYS <<endl;
    }

    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    double curr_sz = yield_size[N_active_ys] ;
    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
    DevStress(i,j) -= pp * curr_alpha(i,j) ;

    static DTensor2 result(3,3,0.);
    result *= 0. ;

    double den = sqrt(DevStress(i,j) * DevStress(i,j));
    if (den == 0)
    {
        return result; //Elastic
    }
    else
    {
        result(i, j) =
            (
                (DevStress(i,j)) + curr_alpha(o, t) * kronecker_delta(i, j) * DevStress(o, t) / 3
            )
            / den;
    }
    result(i, j) += Rtheta(stress, curr_alpha) * sqrt(2./27.) * curr_sz * kronecker_delta(i, j);
    static DTensor2 dr_dsigma(3,3,0.);
    dR_dsigma(stress, curr_alpha, dr_dsigma);
    result(i, j) -= dr_dsigma(i,j) * sqrt(2./3.) * curr_sz * pp ;

    return result;
}

// ================================================================================
// Return the normal to the yield surface w.r.t alpha(backstress)
// ================================================================================
DTensor2 RoundedMohrCoulomb_multi_surface::df_dalpha(int N_active_ys, DTensor2 const& stress){
    if (N_active_ys > TNYS )
    {
        cerr<< "RoundedMohrCoulomb_multi_surface::df_dsigma " <<endl;
        cerr<< "Exceed the length of iterate_alpha_vec " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
        cerr<< "Total NYS" << TNYS <<endl;
    }

    double curr_sz = yield_size[N_active_ys] ;
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];

    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
    DevStress(i,j) -= pp * curr_alpha(i,j) ;

    static DTensor2 result(3,3,0.);
    result *= 0. ;

    double den = sqrt(DevStress(i,j) * DevStress(i,j));
    result(i,j) = - pp * (DevStress(i,j)) / den ;
    static DTensor2 dr_dalpha(3,3,0.);
    dR_dalpha(stress, curr_alpha, dr_dalpha);
    result(i,j) -= curr_sz*pp*dr_dalpha(i,j);

    return result ;
}


// ================================================================================
// Return the plastic flow direction
// ================================================================================
DTensor2 RoundedMohrCoulomb_multi_surface::plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys ){
    static DTensor2 mm(3,3,0.);
    // (1) The deviatoric plastic flow is associative. 
    mm = nn;

    // (2) The dilation component is different. 
    // double I1{0.}, J2{0.}, J3{0.};
    // calc_I1J2J3(stress, I1, J2, J3);
    // // ==============================================
    // double stress_ratio2 = J2 / pow(abs(I1/3.),2);
    // // ==============================================
    // // curr_alpha = alpha_vec[N_active];
    // // double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    // // DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
    // // DevStress(i,j) -= pp * curr_alpha(i,j) ;
    // // double J2bar = 1.5 * DevStress(i,j) * DevStress(i,j);
    // // stress_ratio2 = J2bar / pow(abs(I1/3.),2);
    // // ==============================================
    // double flow{0.};
    // flow = (stress_ratio2/pow(deta,2) - 1 ) / 
    //        (stress_ratio2/pow(deta,2) + 1 ) / 3.;
    // mm(0,0) = flow * scal ;
    // mm(1,1) = flow * scal ;
    // mm(2,2) = flow * scal ;

    // cout<< "stress_ratio2: " << stress_ratio2 <<"\t";
    // cout<< "pow(deta,2)  : " << pow(deta,2)   <<"\n";
    // cout<< "mm(0,1)      : " << mm(0,1)   <<"\n";
    // cout<< "plastic_flow : " << flow   <<"\n";
    return mm;

}

// ================================================================================
// Return the rate of alpha(backstress)
// ================================================================================
DTensor2 RoundedMohrCoulomb_multi_surface::alpha_bar(int N_active_ys, DTensor2 const& stress){

    DTensor2 curr_nn(3,3,0.);
    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    if (N_active_ys > (TNYS-1) )
    {
        // Direction
        curr_nn = df_dsigma(N_active_ys, stress);
        curr_nn(i,j) = curr_nn(i,j)/sqrt( curr_nn(k,l)*curr_nn(k,l) ) ;

        // Magnitude
        double hardening_rate_after_failure = HardingPara[TNYS]/pp  ;
        // cout<< "hardening_rate_after_failure = " << hardening_rate_after_failure <<endl;
        curr_nn(i,j) = hardening_rate_after_failure * curr_nn(i,j) ;
        return curr_nn ;
    }
    double curr_sz = yield_size[N_active_ys];
    double next_sz = yield_size[N_active_ys+1];
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    DTensor2 next_alpha = iterate_alpha_vec[N_active_ys+1];
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;

    DTensor2 direct(3,3,0.);
    if(curr_sz == 0){
        cerr<< "RoundedMohrCoulomb_multi_surface::alpha_bar " <<endl;
        cerr<< "curr_sz == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    direct(i,j) = next_sz/curr_sz * (DevStress(i,j) - pp * curr_alpha(i,j))
                  - (DevStress(i,j) - pp * next_alpha(i,j)) ;

    double denom = sqrt(  direct(i,j)*direct(i,j)  );
    if(denom == 0){
        cerr<< "RoundedMohrCoulomb_multi_surface::alpha_bar " <<endl;
        cerr<< "denom 1 == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    direct(i,j) = direct(i,j) / denom ;
    
    // Change the direct of alpha to the rate of alpha
    curr_nn = df_dsigma( N_active_ys , stress);
    double H_prime = HardingPara[N_active_ys] ; 
    denom = curr_nn(i,j) * direct(i,j); 
    if(denom == 0){
        cerr<< "RoundedMohrCoulomb_multi_surface::alpha_bar " <<endl;
        cerr<< "denom 2 == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    direct(i,j) = direct(i,j) * H_prime / denom / pp ;

    return direct ;
}


// void RoundedMohrCoulomb_multi_surface::dq_dsigma_ij(DTensor2 const& sigma, DTensor2 const& alpha, DTensor2 & result) // Stress derivative of deviatoric stress q
// {
//     static DTensor2 s(3, 3, 0.0);
//     s *= 0;
//     double p = -1./3. * (sigma(0,0)+sigma(1,1)+sigma(2,2)) ;
//     s(i,j) = sigma(i,j) + p * kronecker_delta(i,j) ;

//     s(i,j) -= p * alpha(i,j) ;
//     double s_minus_palpha2 = s(i,j)*s(i,j) ; 
//     double salpha_minus_palapha2  = s(i,j) * alpha(i,j) ; 


//     if (s_minus_palpha2 < machine_epsilon)
//     {
//         result *= 0;
//     }
//     else
//     {
//         result(i, j) = (3. * s(i,j) + salpha_minus_palapha2 * kronecker_delta(i,j)) 
//                         /
//                         sqrt(6. * s_minus_palpha2 ) ;
//     }

//     return;
// }

// void RoundedMohrCoulomb_multi_surface::dq_dalpha_ij(DTensor2 const& sigma, DTensor2 const& alpha, DTensor2 & result) // Stress derivative of deviatoric stress q
// {
//     static DTensor2 s(3, 3, 0.0);
//     s *= 0;
//     double p = -1./3. * (sigma(0,0)+sigma(1,1)+sigma(2,2)) ;
//     s(i,j) = sigma(i,j) + p * kronecker_delta(i,j) ;

//     s(i,j) -= p * alpha(i,j) ;
//     double s_minus_palpha2 = s(i,j)*s(i,j) ; 

//     if (s_minus_palpha2 < machine_epsilon)
//     {
//         result *= 0;
//     }
//     else
//     {
//         result(i, j) = (-1.) * sqrt(1.5) * s(i,j) * p
//                         /
//                         sqrt(s_minus_palpha2) ;
//     }

//     return;
// }


double RoundedMohrCoulomb_multi_surface::Rtheta(DTensor2 const& stress, DTensor2 const& alpha){

    double rtheta{0.};

    rtheta = 2 * kk / 
            (
                (1 + kk) - (1 - kk) * calc_sin3theta(stress, alpha)
            );
    return rtheta;
}


void RoundedMohrCoulomb_multi_surface::dR_dsigma(DTensor2 const& stress, DTensor2 const& alpha, DTensor2& ret){
    double dr_dtheta = dR_dtheta(stress, alpha);
    static DTensor2 dtheta_dsigma(3,3,0.);
    static DTensor2 sigmaBar(3,3,0.); 
    double pp = -1./3. * (stress(0,0) + stress(1,1) + stress(2,2));
    sigmaBar(i,j) = 3 * (stress(i,j) - pp * alpha(i,j)); //According to Prevost 1985. Eq(19).
    dtheta_dsigma_ij(sigmaBar, dtheta_dsigma);
    ret(i,j) = dr_dtheta * dtheta_dsigma(i,j);
}

double RoundedMohrCoulomb_multi_surface::dR_dtheta(DTensor2 const& stress, DTensor2 const& alpha){
    double sin3theta = calc_sin3theta(stress, alpha);
    double cos3theta = sqrt(1 - pow(sin3theta,2));
    double ret{0.};
    ret = 6 * kk * ( 1 - kk) * cos3theta 
          *
          pow(
                1 + kk - (1-kk)*sin3theta ,
                2
            ) ;
    return ret;
}

void RoundedMohrCoulomb_multi_surface::dR_dalpha(DTensor2 const& stress, DTensor2 const& alpha , DTensor2& ret){
    double dr_dtheta = dR_dtheta(stress, alpha);
    static DTensor2 dtheta_dalpha(3,3,0.);
    dtheta_dalpha_ij(stress, alpha, dtheta_dalpha);
    ret(i,j) = dr_dtheta * dtheta_dalpha(i,j);
}


void RoundedMohrCoulomb_multi_surface::dtheta_dsigma_ij(const DTensor2 & sigma, DTensor2 & ret) // Stress derivative of Lode angle
{
    double I1, J2D, J3D, q, theta;
    double trace;

    // ===============================================================
    //Following implementation in stresst.cpp
    // ===============================================================
    static DTensor2 s(3, 3, 0);
    static DTensor2 t(3, 3, 0);
    s *= 0;
    t *= 0;

    const DTensor2& I2 = kronecker_delta;
    calc_I1J2J3(sigma, I1, J2D, J3D) ;
    calc_pqtheta(sigma, trace, q, theta);
    theta = theta * M_PI / 180; /// getpqtheta returns in degrees

    double c3t = cos(3.0 * theta);
    double s3t = sin(3.0 * theta);

    double tempS = (3.0 / 2.0) * (c3t / (q * q * s3t));
    double tempT = (9.0 / 2.0) * (1.0 / (q * q * q * s3t));

    double p = -1./3. * (sigma(0,0)+sigma(1,1)+sigma(2,2)) ;
    s(i,j) = sigma(i,j) + p * kronecker_delta(i,j) ;

    t(i, j) = s(i, k) * s(k, j) - I2(i, j) * (J2D * (2.0 / 3.0));

    ret(i, j) = s(i, j) * tempS - t(i, j) * tempT;
}


void RoundedMohrCoulomb_multi_surface::dtheta_dalpha_ij (DTensor2 const& stress, DTensor2 const& alpha, DTensor2& ret){
    double J2bar{0.}, J3bar{0.}, I1{0.} ;
    static DTensor2 ss(3,3,0.); 
    double pp = -1./3. * (stress(0,0) + stress(1,1) + stress(2,2));
    ss(i,j) = stress(i,j) + pp * kronecker_delta(i,j);
    ss(i,j) = 3. * (ss(i,j) - pp * alpha(i,j) );
    calc_I1J2J3(ss, I1, J2bar, J3bar) ;
    J2bar *= 2. ; 
    J3bar *= 3. ;
    double for_square_root = 1 - 6 * pow(J3bar,2)/pow(J2bar,3);
    if(for_square_root < 0 ){
        cerr<< " RoundedMohrCoulomb_multi_surface::dtheta_dalpha_ij " <<endl;
        cerr<< " Numerical Error: (1 - 6 * pow(J3bar,2)/pow(J2bar,3)) < 0 " <<endl;
    }
    ret(i,j) = 1./3. * sqrt(for_square_root) * (-1) * sqrt(6.) *
                (
                    1./pow(J2bar,1.5) * (-27*pp*ss(i,j)) 
                    +
                    J3bar * (-1.5) * pow(J2bar, -2.5) * (-6 * pp * ss(i,j))
                );
}




// //================================================================================
// // Return the 6*6 Tangent matrix.
// //================================================================================
// DTensor2 const& RoundedMohrCoulomb_multi_surface::getTangent(void){
//     // =========================================================
//     // (1) Update the modulus based on active yield surface.
//     // =========================================================
//     double pp = -1./3. * (TrialStress(0,0) + TrialStress(1,1) + TrialStress(2,2));
//     double curr_h{0.};
//     if (N_active>0){
//         curr_h = HardingPara[N_active];
//     }else{
//     // Unloading. Get back the initial modulus.
//         curr_h = initial_E / 2. / (1 + v) ;
//     }
//     E = 2 * curr_h * (1 + v) ;
//     // =========================================================
//     // (2) Update the modulus based on pressure
//     // =========================================================
//     E = E * pow(abs(pp/pa),modn);


//     double mu2 = E / (1.0 + v);
//     double lam = v * mu2 / (1.0 - 2.0 * v);
//     double mu  = 0.50 * mu2;

//     mu2 += lam;

//     D(0, 0) = D(1, 1) = D(2, 2) = mu2;
//     D(0, 1) = D(1, 0) = lam;
//     D(0, 2) = D(2, 0) = lam;
//     D(1, 2) = D(2, 1) = lam;
//     D(3, 3) = mu;
//     D(4, 4) = mu;
//     D(5, 5) = mu;

//     return D;
// }

//================================================================================
// Return the 3*3*3*3 Tangent tensor.
//================================================================================
DTensor4 const& RoundedMohrCoulomb_multi_surface::getTangentTensor( void )
{
    // compute_tangent_tensor();
    return Ee;
}

//================================================================================
// Compute the Elastic Tangent Stiffness
//================================================================================
void RoundedMohrCoulomb_multi_surface::update_modulus(int N_active_ys, DTensor2 const& stress)
{
    Ee *= 0;
    // =========================================================
    double curr_h = HardingPara[N_active_ys];
    E = 2 * curr_h * (1 + v) ;
    double pp = -1./3. * (stress(0,0) + stress(1,1) + stress(2,2));
    E = E * pow(abs(pp/pa),modn);
    // =========================================================
    // cout<< "N_active_ys " << N_active_ys <<endl;
    // cout<< "curr_h        " << curr_h    <<endl;
    double lambda = ( v * E ) / ( ( 1 + v ) * ( 1 - 2 * v ) );
    double mu = E / ( 2 * ( 1 + v ) );
    Ee(i,j,k,l) = mu * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k))
                    + lambda * kronecker_delta(k,l)*kronecker_delta(i,j) ;
 
    // Reduce the shear stiffness AND Keep the bulk Modulus.
    // double K = initial_E / 3. / (1 - 2 * v) ; 
    // Ee(i,j,k,l) = K * kronecker_delta(i,j) * kronecker_delta(k,l) + 
    //                 mu * ( kronecker_delta(i,k) * kronecker_delta(j,l) 
    //                     + kronecker_delta(i,l) * kronecker_delta(j,k)
    //                     + (-2./3.) * kronecker_delta(i,j) * kronecker_delta(k,l) ) ; 

    // if(N_active_ys > 0){
    //     compute_elastoplastic_tangent();
    // }
}

void RoundedMohrCoulomb_multi_surface::compute_elastoplastic_tangent( int N_active_ys, DTensor2 const& intersection_stress , bool elastic){
    DTensor2 curr_xi = df_dalpha(N_active_ys, intersection_stress);
    DTensor2 bar_alpha = alpha_bar(N_active_ys, intersection_stress);
    DTensor2 curr_nn = df_dsigma(N_active_ys,  intersection_stress);
    DTensor2 curr_mm = plastic_flow_direct( curr_nn,  intersection_stress , N_active_ys);
    update_modulus(N_active_ys, intersection_stress);
    double denominator = curr_nn(i,j) * Ee(i,j,k,l) * curr_mm(k,l) - curr_xi(o,t) * bar_alpha(o,t);
    // if (denominator == 0 ){
    //     cerr<< "DruckerPrager_multi_yield_surface::compute_elastoplastic_tangent()" <<endl;
    //     cerr<< "Error denominator == 0 " <<endl;
    //     cerr<< "LEFT  = curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) == " << curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) <<endl;
    //     cerr<< "RIGHT = curr_xi(o,t) * bar_alpha(o,t)             == " << curr_xi(o,t) * bar_alpha(o,t) <<endl;
    // }
    Eep(i,j,k,l) = Ee(i,j,k,l) - Ee(i,j,o,t) * curr_mm(o,t) * curr_nn(x,y) * Ee(x,y,k,l) / denominator;

    // if(elastic == true ){
    //     Eep = Ee  ;
    // }

    // For the moment, use elastic matrix temporarily.
    Eep = Ee  ;

    // // =========NEW Test=======================
    // double curr_sz = yield_size[N_active_ys];
    // double GG = HardingPara[0];
    // denominator = (12. * GG + 6 * HardingPara[N_active_ys]  ) * pow(1.5 * curr_sz, 2.);
    // DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    // DTensor2 HH(3,3,0.);
    // double pp = - 1. /3. * (intersection_stress(0,0) + intersection_stress(1,1) + intersection_stress(2,2));
    // DTensor2 DevStress(3,3,0.);
    // DevStress(i,j) = intersection_stress(i,j) + pp * kronecker_delta(i,j);
    // HH(i,j) = 6. * GG * (DevStress(i,j) - curr_alpha(i,j)) ; 
    // Eep(i,j,k,l) = Ee(i,j,k,l) - HH(i,j) * HH(k,l) / denominator ; 
}



//================================================================================
RoundedMohrCoulomb_multi_surface *RoundedMohrCoulomb_multi_surface::getCopy( void ){
    RoundedMohrCoulomb_multi_surface *tmp = new RoundedMohrCoulomb_multi_surface( 
            0, //this->getTag(),
            this->getE(),
            this->getv(),
            this->getRho(),
            this->getpp0(),
            this->getpa(),
            this->getmodn(),
            this->getpc(),
            this->getkk(),
            this->getdeta(),
            this->getscal(),
            this->getTNYS(),
            this->getYieldSize(),
            this->getHardPara() 
       );


    return tmp;
}



// //================================================================================
// // Message passing for parallel
// //================================================================================
// int RoundedMohrCoulomb_multi_surface::sendSelf( int commitTag, Channel &theChannel )
// {
//     cerr<<"RoundedMohrCoulomb_multi_surface::sendSelf() is not implemented yet! " <<endl;
//     static ID idData(1);
//     static Vector vectorData(3);
//     static Matrix a(3, 3);

//     idData(0) = this->getTag();

//     if (theChannel.sendID(0, commitTag, idData) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::sendSelf -- could not send idData\n";
//         return -1;
//     }

//     vectorData(0) = E;
//     vectorData(1) = v;
//     vectorData(2) = rho;

//     if (theChannel.sendVector(0, commitTag, vectorData) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::sendSelf -- could not send vectorData\n";
//         return -1;
//     }

//     a.setData(CommitStress.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::sendSelf -- could not send CommitStress\n";
//         return -1;
//     }

//     a.setData(CommitStrain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::sendSelf -- could not send CommitStrain\n";
//         return -1;
//     }

//     a.setData(CommitPlasticStrain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::sendSelf -- could not send CommitPlasticStrain\n";
//         return -1;
//     }

//     return 0;
// }

// //================================================================================
// // Message passing for parallel
// //================================================================================
// int RoundedMohrCoulomb_multi_surface::receiveSelf( int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker )
// {
//     cerr<<"RoundedMohrCoulomb_multi_surface::receiveSelf() is not implemented yet! " <<endl;
//     static ID idData(1);
//     static Vector vectorData(3);
//     static Matrix a(3, 3);

//     if (theChannel.receiveID(0, commitTag, idData) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive idData\n";
//         return -1;
//     }
//     this->setTag(idData(0));

//     if (theChannel.receiveVector(0, commitTag, vectorData) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive vectorData\n";
//         return -1;
//     }
//     E    =   vectorData(0) ;
//     v    =   vectorData(1) ;
//     rho  =   vectorData(2) ;

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             CommitStress(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             CommitStrain(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "RoundedMohrCoulomb_multi_surface::receiveSelf -- could not receive Elastic Constant Tensor\n";
//         return -1;
//     }

//     compute_tangent_tensor();
//     return 0;
// }

//================================================================================
void RoundedMohrCoulomb_multi_surface::Print( ostream &s, int flag )
{
    s << "RoundedMohrCoulomb_multi_surface::" << endl;
    // s << "\tTag:    " << this->getTag() << endl;


}


int RoundedMohrCoulomb_multi_surface::calc_I1J2J3(DTensor2 const& mystress, double& I1, double& J2, double& J3) const
{
    // ------------------------------------------------------------
    // preliminary
    I1 = mystress(0, 0) + mystress(1, 1) + mystress(2, 2);
    const double sigma_m = I1 / 3.0;
    DTensor2 s = mystress;
    s(0, 0) -= sigma_m;
    s(1, 1) -= sigma_m;
    s(2, 2) -= sigma_m;
    // J2=0.5*s(i,j)*s(i,j)
    J2 = 0.5 * (
                  s(0, 0) * s(0, 0) +  s(0, 1) * s(0, 1)  + s(0, 2) * s(0, 2)
                  +   s(1, 0) * s(1, 0) +  s(1, 1) * s(1, 1)  + s(1, 2) * s(1, 2)
                  +   s(2, 0) * s(2, 0) +  s(2, 1) * s(2, 1)  + s(2, 2) * s(2, 2)   );
    // 3by3 Determinant: Refer to http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Tensors/Tensors.htm
    J3 = s(0, 0) * (s(1, 1) * s(2, 2) - s(1, 2) * s(2, 1))
      +  s(0, 1) * (s(1, 2) * s(2, 0) - s(1, 0) * s(2, 2))
      +  s(0, 2) * (s(1, 0) * s(2, 1) - s(2, 0) * s(1, 1));


    return 0;
}


// --------------------------------------------
// refer to https://en.wikipedia.org/wiki/Yield_surface
// p = 1/3 * I1
// q = sqrt(3* J2)
// cos(3*theta) = 3/2 * sqrt(3) * J3 / J2^(3/2)
// ------------------------------------------------------------
int RoundedMohrCoulomb_multi_surface::calc_pqtheta(DTensor2 const& sigma, double& p, double& q, double& theta) const
{
    double I1, J2, J3;
    calc_I1J2J3(sigma, I1, J2, J3);

    if (J2 < machine_epsilon)
    {
        J2 = 0;
        p = -I1 / 3;
        q = 0;
        theta = 0;
        return 0;
    }
    else
    {
        // (1) calculate p
        p = -I1 / 3;
        // (2) calculate q
        q = sqrt(3 * J2);
        // (3) calculate theta
        double cos3theta = 1.5 * sqrt(3) * J3 / pow(J2, 1.5);
        cos3theta = cos3theta < 1.0 ? cos3theta : 1.0;
        cos3theta = cos3theta > -1.0 ? cos3theta : -1.0;
        theta = acos(cos3theta) / 3.0 * 180 / M_PI; //theta is in degree.

        return 0;
    }
    return -1;
}

double RoundedMohrCoulomb_multi_surface::calc_sin3theta(DTensor2 const& stress, DTensor2 const& alpha){
    double sin3theta{0.}, J2bar{0.}, J3bar{0.}, I1{0.} ;
    static DTensor2 ss(3,3,0.); 
    double pp = -1./3. * (stress(0,0) + stress(1,1) + stress(2,2));
    ss(i,j) = stress(i,j) + pp * kronecker_delta(i,j);
    ss(i,j) = 3. * (ss(i,j) - pp * alpha(i,j) );
    calc_I1J2J3(ss, I1, J2bar, J3bar) ;
    J2bar *= 2. ; 
    J3bar *= 3. ;
    sin3theta = - sqrt(6.) * J3bar/pow(J2bar,1.5);
    if (sin3theta >1)
    {
        cerr<<" RoundedMohrCoulomb_multi_surface::calc_sin3theta " <<endl;
        cerr<<" Numerical Error: sin3theta >1 " <<endl;
    }
    return sin3theta;
}


//================================================================================
// The Getter
//================================================================================
double RoundedMohrCoulomb_multi_surface::getpp0() const {return pp0;}
double RoundedMohrCoulomb_multi_surface::getpa() const {return pa;}
double RoundedMohrCoulomb_multi_surface::getmodn() const {return modn;}
double RoundedMohrCoulomb_multi_surface::getpc() const {return pc;}
double RoundedMohrCoulomb_multi_surface::getkk() const {return kk;}
double RoundedMohrCoulomb_multi_surface::getdeta() const {return deta;}
double RoundedMohrCoulomb_multi_surface::getscal() const {return scal;}
