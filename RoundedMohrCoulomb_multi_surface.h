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
/****************************************************************************************************************
* Reference:                                                                                                    *
* 1. In the article below, the Section (4). Hardening Rule is the key difference                                *
*    from the classical plasticity.                                                                             *
*                                                                                                               *
*     @article{Prevost1985,                                                                                     *
*     title = "A simple plasticity theory for frictional cohesionless soils",                                   *
*     doi = "http://dx.doi.org/10.1016/0261-7277(85)90030-0",                                                   *
*     author = "Jean H. Prevost",                                                                               *
*     napomena = { local eCM2155 ; 6Feb2016}}                                                                   *
*                                                                                                               *
*                                                                                                               *
* 2. In the book below, the Section 7.3.1 Algorithmic Set-up and Flow-chart are very helpful.                   *
*                                                                                                               *
*     @book{Prevost1989,                                                                                        *
*     title={DYNA1D: a computer program for nonlinear seismic site response analysis technical documentation},  *
*     author={Prevost, Jean H},                                                                                 *
*     year={1989},                                                                                              *
*     napomena = { local eCM2156 ; 6Feb2016}}                                                                   *
*****************************************************************************************************************/
#ifndef RoundedMohrCoulomb_multi_surface_H
#define RoundedMohrCoulomb_multi_surface_H

#include "MultiYieldSurfaceMaterial.h"
// #include <G3Globals.h>
// #include "../../ltensor/LTensor.h"
// #include <Channel.h>

#include <iostream>
#include "./ltensor/LTensor.h"
using namespace std;

class RoundedMohrCoulomb_multi_surface : public MultiYieldSurfaceMaterial
{
public:
    RoundedMohrCoulomb_multi_surface( 
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
        vector<double> const& eta ,
        vector<double> const& HardingPara 
    );                                                      // Constructor for DSL
    RoundedMohrCoulomb_multi_surface();                     // Empty Constructor for Parallel
    virtual ~RoundedMohrCoulomb_multi_surface( void );      // Destructor

    // Message Passing in Parallel 
    // virtual int sendSelf( int commitTag, Channel &theChannel );                                 
    // virtual int receiveSelf( int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker ); 


    RoundedMohrCoulomb_multi_surface *getCopy( void );
    virtual DTensor4 const& getTangentTensor( void );       
    virtual void compute_elastoplastic_tangent( int N_active_ys, DTensor2 const& intersection_stress , bool elastic=false );
    virtual void Print( ostream &s, int flag = 0 );
    virtual const char *getType( void ) const{return "RoundedMohrCoulomb_multi_surface";};
    virtual const char *getClassType( void ) const{return "RoundedMohrCoulomb_multi_surface";};
    virtual int getObjectSize(){return sizeof(*this);};

    double getpp0() const;                  // Return the material constant pp0
    double getpa() const;                   // Return the material constant pa
    double getmodn() const;                 // Return the material constant modn
    double getpc() const;                   // Return the material constant pc
    double getkk() const;                   // Return the material constant kk
    double getdeta() const;                 // Return the material constant deta
    double getscal() const;                 // Return the material constant scal


private:
    virtual void update_modulus( int num_active_ys, DTensor2 const& s );                                   // Compute the Tangent Stiffness
    virtual double yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double eta);           // Return the yield surface Value
    virtual DTensor2 df_dsigma(int num_active_ys, DTensor2 const& stress);                                 // Return the normal to the yield surface w.r.t stress
    virtual DTensor2 df_dalpha(int num_active_ys, DTensor2 const& stress);                                 // Return the normal to the yield surface w.r.t alpha(backstress)
    virtual DTensor2 plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys);     // Return the plastic flow
    virtual DTensor2 alpha_bar(int num_active_ys, DTensor2 const& stress);                                 // Return the rate of alpha(backstress)

    int calc_I1J2J3(DTensor2 const& mystress, double& I1, double& J2, double& J3) const;
    int calc_pqtheta(DTensor2 const& sigma, double& p, double& q, double& theta) const;
    double calc_sin3theta(DTensor2 const& stress, DTensor2 const& alpha);
    double Rtheta(DTensor2 const& sigma, DTensor2 const& alpha);
    void dR_dsigma(DTensor2 const& stress, DTensor2 const& alpha, DTensor2& ret);
    void dR_dalpha(DTensor2 const& stress, DTensor2 const& alpha , DTensor2& ret);
    double dR_dtheta(DTensor2 const& stress, DTensor2 const& alpha);
    void dtheta_dalpha_ij (DTensor2 const& stress, DTensor2 const& alpha, DTensor2& ret);
    void dtheta_dsigma_ij(const DTensor2 & sigma, DTensor2 & ret);
private:
  
    double pp0;                            // Initial Confinement
    double pa;                             // Reference pressure p_1
    double modn;                           // Exponential parameter n for modulus change w.r.t Reference pressure.
    double pc;                             // Represents the conhesion p_c.
    double kk;                             // Shape of Yield Surface. When k = 1, same to the Drucker-Prager.
    double deta;                           // Dilation Angle \bar{eta}
    double scal;                           // Dilation Scale 
   
};


#endif
