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

#ifndef vonMises_multi_surface_H
#define vonMises_multi_surface_H

#include "MultiYieldSurfaceMaterial.h"

// #include <G3Globals.h>
// #include <iostream>
// #include "../../ltensor/LTensor.h"
// #include <Channel.h>

#include "./ltensor/LTensor.h"
using namespace std;

class vonMises_multi_surface : public MultiYieldSurfaceMaterial
{
public:
    vonMises_multi_surface( 
        int tag,
        double E_in,
        double v_in,
        double rho_in,
        int TNYS_in,
        vector<double> const& radius ,
        vector<double> const& HardingPara 
    );                                                      // Constructor for DSL
    vonMises_multi_surface();                               // Empty Constructor for Parallel
    virtual ~vonMises_multi_surface( void );                // Destructor


    virtual vonMises_multi_surface *getCopy( void );
    virtual DTensor4 const& getTangentTensor( void );
    virtual void compute_elastoplastic_tangent(int N_active_ys, DTensor2 const& intersection_stress , bool elastic=false);
    virtual void Print( ostream &s, int flag = 0 );
    virtual const char *getType( void ) const;
    virtual const char *getClassType( void ) const{return "vonMises_multi_surface";};
    virtual int getObjectSize(){return sizeof(*this);};


private:
    virtual void update_modulus( int N_active_ys, DTensor2 const& s );               // Compute the Tangent Stiffness
    virtual double yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double radius);                    // Return the yield surface Value
    virtual DTensor2 df_dsigma(int N_active_ys, DTensor2 const& stress);     // Return the normal to the yield surface w.r.t stress
    virtual DTensor2 df_dalpha(int N_active_ys, DTensor2 const& stress);     // Return the normal to the yield surface w.r.t alpha(backstress)
    virtual DTensor2 plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys);     // Return the plastic flow
    virtual DTensor2 alpha_bar(int N_active_ys, DTensor2 const& stress);     // Return the rate of alpha(backstress)


};


#endif
