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

#include "vonMises_multi_surface.h"
// #include <Channel.h>
// #include <HDF5_Channel.h>
// #include <Matrix.h>
// #include <ID.h>
// #include <Vector.h>
#include <vector>

using namespace std;

#define BRENT_MAXITER 20
#define BRENT_TOLERANCE 1e-6
//================================================================================
// Constructor
//================================================================================
vonMises_multi_surface::vonMises_multi_surface( 
        int tag,
        double E_in,
        double v_in,
        double rho_in,
        int TNYS_in,
        vector<double> const& yield_surface_size_in ,
        vector<double> const& HardingPara_in
        )
    : 
    MultiYieldSurfaceMaterial( 
            tag, 
            0, //ND_TAG_vonMises_multi_surface, 
            E_in, 
            v_in, 
            rho_in, 
            TNYS_in, 
            yield_surface_size_in, 
            HardingPara_in)
{
}

//================================================================================
// Empty Constructor for parallel
//================================================================================
vonMises_multi_surface::vonMises_multi_surface( )
    : 
    MultiYieldSurfaceMaterial()
{
}


//================================================================================
// Destructor
//================================================================================
vonMises_multi_surface::~vonMises_multi_surface()
{

}


// ================================================================================
// Return the normal to the yield surface w.r.t stress
// ================================================================================
DTensor2 vonMises_multi_surface::df_dsigma(int N_active_ys, DTensor2 const& stress){
    if (N_active_ys > TNYS )
    {
        cerr<< "vonMises_multi_surface::df_dsigma " <<endl;
        cerr<< "Exceed the length of alpha_vec " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
        cerr<< "Total NYS " << TNYS <<endl;
    }
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;

    DTensor2 s_minus_alpha(3,3,0.);
    s_minus_alpha(i,j) = DevStress(i,j) - curr_alpha(i,j) ; 
    double denominato = sqrt(  s_minus_alpha(i,j)*s_minus_alpha(i,j)   )  ; 

    s_minus_alpha(i,j) = s_minus_alpha(i,j) / denominato; 

    return s_minus_alpha ;
}

// ================================================================================
// Return the normal to the yield surface w.r.t alpha(backstress)
// ================================================================================
DTensor2 vonMises_multi_surface::df_dalpha(int N_active_ys, DTensor2 const& stress){
    if (N_active_ys > TNYS )
    {
        cerr<< "vonMises_multi_surface::df_dalpha " <<endl;
        cerr<< "Exceed the length of iterate_alpha_vec " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
        cerr<< "Total NYS" << TNYS <<endl;
    }
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;

    DTensor2 s_minus_alpha(3,3,0.);
    s_minus_alpha(i,j) = DevStress(i,j) - curr_alpha(i,j) ; 
    double denominato = sqrt(  s_minus_alpha(i,j)*s_minus_alpha(i,j)   )  ; 
    // change s_minus_alpha to result df_dalpha. 
    s_minus_alpha(i,j) = -1.0 * s_minus_alpha(i,j) / denominato; 

    return s_minus_alpha ;
}

// ================================================================================
// Return the plastic flow direction
// ================================================================================
DTensor2 vonMises_multi_surface::plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys ){
    static DTensor2 mm(3,3,0.);
    mm *= 0. ;
    // (1) The plastic flow is associative. 
    mm = nn;

    return mm;

}

// ================================================================================
// Return the rate of alpha(backstress)
// ================================================================================
DTensor2 vonMises_multi_surface::alpha_bar(int N_active_ys, DTensor2 const& stress){
    DTensor2 curr_nn(3,3,0.);
    if (N_active_ys > (TNYS-1) )
    {
        curr_nn = df_dsigma(N_active_ys, stress);
        // after the failure surface, give a small stiffness. 
        double hardening_rate_after_failure = HardingPara[TNYS]/1E4 ;
        curr_nn(i,j) = hardening_rate_after_failure * curr_nn(i,j) ;
        return curr_nn ;
    }
    double curr_radius = yield_size[N_active_ys];
    double next_radius = yield_size[N_active_ys+1];
    DTensor2 curr_alpha = iterate_alpha_vec[N_active_ys];
    DTensor2 next_alpha = iterate_alpha_vec[N_active_ys+1];
    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 DevStress(3,3,0.);
    DevStress(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;
    DTensor2 direct(3,3,0.);
    if(curr_radius == 0){
        cerr<< "vonMises_multi_surface::alpha_bar " <<endl;
        cerr<< "curr_radius == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    // cout<< " N_active_ys      " << N_active_ys <<endl;
    // cout<< " next_radius " << next_radius <<endl;
    // cout<< " curr_radius " << curr_radius <<endl;

    direct(i,j) = next_radius/curr_radius * (DevStress(i,j) - curr_alpha(i,j))
                  - (DevStress(i,j) - next_alpha(i,j)) ;
    double denom = sqrt(  direct(i,j)*direct(i,j)  );
    if(denom == 0){
        cerr<< "vonMises_multi_surface::alpha_bar " <<endl;
        cerr<< "denom 1 == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    direct(i,j) = direct(i,j) / denom ;
    
    // Change the direct of alpha to the rate of alpha
    curr_nn = df_dsigma( N_active_ys, stress );
    double H_prime = HardingPara[N_active_ys] ; 
    denom = curr_nn(i,j) * direct(i,j); 
    if(denom == 0){
        cerr<< "vonMises_multi_surface::alpha_bar " <<endl;
        cerr<< "denom 2 == 0 " <<endl;
        cerr<< "N_active_ys " << N_active_ys <<endl;
        cerr<< "yield_size.size() " << yield_size.size() <<endl;
        cerr<< "iterate_alpha_vec.size() " << iterate_alpha_vec.size() <<endl;
    }
    direct(i,j) = direct(i,j) * H_prime / denom ;
    return direct ;
}

// ================================================================================
// Return the yield surface Value
// ================================================================================
double vonMises_multi_surface::yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double radius){

    double ret{0.} ; 

    double pp = -1./3. * (stress(0,0)+stress(1,1)+stress(2,2)) ;
    DTensor2 s(3,3,0.) ;
    s(i,j) = stress(i,j) + pp * kronecker_delta(i,j) ;

    DTensor2 s_minus_alpha(3,3,0.);
    s_minus_alpha(i,j) = s(i,j) - alpha(i,j) ; 
    double stress_state = sqrt(  s_minus_alpha(i,j)*s_minus_alpha(i,j)   )  ; 
    ret = stress_state - sqrt(2./3.) * radius ; 
    return ret;

}


// //================================================================================
// // Return the 6*6 Tangent matrix.
// //================================================================================
// DTensor2 const& vonMises_multi_surface::getTangent(void){
//     // double mu2 = E / (1.0 + v);
//     // double lam = v * mu2 / (1.0 - 2.0 * v);
//     // double mu  = 0.50 * mu2;

//     // mu2 += lam;

//     // D(0, 0) = D(1, 1) = D(2, 2) = mu2;
//     // D(0, 1) = D(1, 0) = lam;
//     // D(0, 2) = D(2, 0) = lam;
//     // D(1, 2) = D(2, 1) = lam;
//     // D(3, 3) = mu;
//     // D(4, 4) = mu;
//     // D(5, 5) = mu;
//     return D;
// }


//================================================================================
// Compute the Elastic Tangent Stiffness
//================================================================================
void vonMises_multi_surface::update_modulus( int N_active_ys, DTensor2 const& s )
{
    Ee *= 0;
    // =========================================================
    double curr_h = HardingPara[N_active_ys];
    E = 2 * curr_h * (1 + v) ;
    // =========================================================
    // cout<< "N_active_ys " << N_active_ys <<endl;
    // cout<< "curr_h        " << curr_h        <<endl;
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

DTensor4 const& vonMises_multi_surface::getTangentTensor()
{
    if (iterate_N_active != 0)
    {
        // compute_elastoplastic_tangent(iterate_N_active, iterate_stress);
    }
    return Eep;
}

void vonMises_multi_surface::compute_elastoplastic_tangent(int N_active_ys, DTensor2 const& intersection_stress, bool elastic){
    DTensor2 curr_xi = df_dalpha(N_active_ys, intersection_stress);
    DTensor2 bar_alpha = alpha_bar(N_active_ys, intersection_stress);
    DTensor2 curr_nn = df_dsigma(N_active_ys,  intersection_stress);
    update_modulus(N_active_ys, intersection_stress);
    double denominator = curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) - curr_xi(o,t) * bar_alpha(o,t);
    // if (denominator == 0 ){
    //     cerr<< "vonMises_multi_surface::compute_elastoplastic_tangent()" <<endl;
    //     cerr<< "Error denominator == 0 " <<endl;
    //     cerr<< "LEFT  = curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) == " << curr_nn(i,j) * Ee(i,j,k,l) * curr_nn(k,l) <<endl;
    //     cerr<< "RIGHT = curr_xi(o,t) * bar_alpha(o,t)             == " << curr_xi(o,t) * bar_alpha(o,t) <<endl;
    // }
    Eep(i,j,k,l) = Ee(i,j,k,l) - Ee(i,j,o,t) * curr_nn(o,t) * curr_nn(x,y) * Ee(x,y,k,l) / denominator;

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
vonMises_multi_surface *vonMises_multi_surface::getCopy( void ){
    vonMises_multi_surface *tmp = new vonMises_multi_surface( 
            1, //this->getTag(),
            this->getE(),
            this->getv(),
            this->getRho(),
            this->getTNYS(),
            this->getYieldSize(),
            this->getHardPara()
            );


    return tmp;
}

//================================================================================
const char *vonMises_multi_surface::getType( void ) const
{
    return "vonMises_multi_surface";
}

// //================================================================================
// // Message passing for parallel
// //================================================================================
// int vonMises_multi_surface::sendSelf( int commitTag, Channel &theChannel )
// {
//     cerr<<"vonMises_multi_surface::sendSelf() is not implemented yet! " <<endl;
//     ID idData(1);
//     Vector vectorData(6);
//     Matrix a(3, 3);

//     idData(0) = this->getTag();

//     if (theChannel.sendID(0, commitTag, idData) < 0)
//     {
//         cerr << "vonMises_multi_surface::sendSelf -- could not send idData\n";
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
//         cerr << "vonMises_multi_surface::sendSelf -- could not send vectorData\n";
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
//         cerr << "vonMises_multi_surface::sendSelf -- could not send multi_ys_data\n";
//         return -1;
//     }


//     a.setData(converge_commit_stress.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::sendSelf -- could not send converge_commit_stress\n";
//         return -1;
//     }

//     a.setData(converge_commit_strain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::sendSelf -- could not send converge_commit_strain\n";
//         return -1;
//     }

//     a.setData(converge_commit_plastic_strain.data, 3, 3);
//     if (theChannel.sendMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::sendSelf -- could not send converge_commit_plastic_strain\n";
//         return -1;
//     }
    
//     for (int it = 0; it < TNYS+1 ; ++i)
//     {
//         a.setData(converge_commit_alpha_vec[it].data, 3, 3);
//         if (theChannel.sendMatrix(0, 0, a) < 0)
//         {
//             cerr << "vonMises_multi_surface::sendSelf -- could not send converge_commit_alpha_vec\n";
//             return -1;
//         }
//     }

//     return 0;
// }

// //================================================================================
// // Message passing for parallel
// //================================================================================
// int vonMises_multi_surface::receiveSelf( int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker )
// {
//     cerr<<"vonMises_multi_surface::receiveSelf() is not implemented yet! " <<endl;
//     ID idData(1);
//     Vector vectorData(6);
//     Matrix a(3, 3);

//     if (theChannel.receiveID(0, commitTag, idData) < 0)
//     {
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive idData\n";
//         return -1;
//     }
//     this->setTag(idData(0));

//     if (theChannel.receiveVector(0, commitTag, vectorData) < 0)
//     {
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive vectorData\n";
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
//         cerr << "vonMises_multi_surface::sendSelf -- could not send multi_ys_data\n";
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
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             converge_commit_stress(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive Elastic Constant strain\n";
//         return -1;
//     }
//     for (int ii = 0; ii < 3; ii++)
//         for (int jj = 0; jj < 3; jj++)
//         {
//             converge_commit_strain(ii, jj) = a(ii, jj);
//         }

//     if (theChannel.receiveMatrix(0, 0, a) < 0)
//     {
//         cerr << "vonMises_multi_surface::receiveSelf -- could not receive Elastic Constant Tensor\n";
//         return -1;
//     }

    

//     for (int it = 0; it < TNYS+1 ; ++i)
//     {
//         if (theChannel.sendMatrix(0, 0, a) < 0)
//         {
//             cerr << "vonMises_multi_surface::sendSelf -- could not send converge_commit_alpha_vec\n";
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

//================================================================================
void vonMises_multi_surface::Print( ostream &s, int flag )
{
    s << "vonMises_multi_surface::" << endl;
    // s << "\t Tag      : " << this->getTag() << endl;
    s << "\t E        : " << this->getE() << endl;
    s << "\t v        : " << this->getv() << endl;
    s << "\t Rho      : " << this->getRho() << endl;
    s << "\t TNYS     : " << this->getTNYS() << endl;
    // s << "\t Radius   : " << this->getYieldSize() << endl;
    // s << "\t HardPara : " << this->getHardPara() << endl;
}
