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

#ifndef MultiYieldSurfaceMaterial_h
#define MultiYieldSurfaceMaterial_h

#include "./ltensor/LTensor.h"

// #include <Material.h>
// #include "../../ltensor/LTensor.h"
// #include "../LTensor_Based/NDMaterialLT.h"
// class ID;
// class Information;
// // class Response;
// class Channel;
// class HDF5_Channel;


constexpr double machine_epsilon = std::numeric_limits<double>::epsilon();


using namespace std;
enum struct MultiYieldSurfaceMaterial_Constitutive_Integration_Method : int
{
    Not_Set,
    Forward_Euler,
    Backward_Euler,
    Backward_Euler_Subincrement,
    Forward_Euler_Subincrement
};

// class MultiYieldSurfaceMaterial : public NDMaterialLT
class MultiYieldSurfaceMaterial 
{
public:
    MultiYieldSurfaceMaterial(
        int tag, 
        int classTag, 
        double E_in,
        double v_in,
        double rho_in,
        int TNYS_in,
        vector<double> const& radius ,
        vector<double> const& HardingPara 
    );
    MultiYieldSurfaceMaterial();
    virtual ~MultiYieldSurfaceMaterial();

/************************************************************************************
* Yuan @ Thu Mar 2 21:35:05 PST 2017                                                *
* Following Methods are Implementated in this class.                                *
************************************************************************************/
    // Driver to run the constitutive integrater.
    virtual int setTrialStrain(DTensor2 const& other);      // Set the Total Strain, usually call by element
    virtual int setTrialStrainIncr(DTensor2 const& other);  // Set the Increment Strain, usually call by element
    virtual int compute_stress(DTensor2 const& strain_increment, int Nsubsteps = 1 );                               // Compute the Stress from the Increment Strain

    // For finite-element call
    // virtual DTensor2 const& getTangent (void);              // Return the 6*6 Tangent matrix.

    // State and Step Control
    // 1. Store and Roll-back Converge State for incremental strain approach.
    virtual int commitState( void );                        // Commit/store the Trial state.
    virtual int revertToLastCommit( void );                 // Go back to last Commit/stored state.
    virtual int revertToStart( void );                      // Go back to all zeroes.
    // 2. Store and Roll-back Local State for constutitive-level iteration.
    virtual int saveIterateState();
    virtual int backToLastIterateState();

    // Message Passing in Parallel 
    // virtual int sendSelf( int commitTag, Channel &theChannel );                                  // Message passing for parallel
    // virtual int receiveSelf( int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker );  // Message passing for parallel

    virtual void zeroStrain();
    
    // Getter
    virtual double getE()   const;                                      // Return the elastic_modulus
    virtual double getv()   const;                                      // Return the possion ratio
    virtual double getRho() const;                                      // Return the density
    virtual double getNumActiveYS(void) const;                          // Return the number of active Yield Surface
    virtual double getTNYS() const;                                     // Return the total number of Yield Surface
    virtual vector<double> const& getYieldSize() const;                    // Return the vector of Yield Surface Radius
    virtual vector<double> const& getHardPara() const;                  // Return the vector of Yield Surface Hardening Parameter
    virtual DTensor2 const& getStressTensor(void) const;                // Return the current commit stress
    virtual DTensor2 const& getStrainTensor(void) const;                // Return the current commit strain
    virtual DTensor2 const& getPlasticStrainTensor(void) const;         // Return the current plastic strain
    static bool set_constitutive_integration_method(int method, double f_relative_tol, double stress_relative_tol, int n_max_iterations, double allowed_subincrement);
    virtual const char *getClassType( void ) const{return "MultiYieldSurfaceMaterial";};
private:
    // Lower-Level helper Subroutine internal this class.
    // Action of the algorithms in the flow-control
    virtual void update_stress(DTensor2& target_stress, double& lambda, int N_active_ys, DTensor2 const& normal_refer_stress);                  // After the yield, Update the Stress 
    virtual void update_plastic_strain(DTensor2& pstrain, double lambda, int num_active_ys, DTensor2 const& stress );           // After the yield, Update the plastic strain
    virtual void update_current_yield_surface(double& curr_yf_val, int num_active_ys, double lambda, DTensor2 const& stress);   // After the yield, Update the current active yield surface (alpha)
    virtual void update_inner_yield_surfaces(int num_active_ys, DTensor2 const& stress);    // After the yield, Update the inner yield surface (alpha) 
    virtual void update_failure_surface(DTensor2 const& stress);
    virtual void correct_update_stress(DTensor2& stress, double& lambda1, double& lambda2, int num_active_ys);          // After the overshooting, Correct the overshooting stress
    virtual void correct_update_plastic_strain(DTensor2& pstrain, double lambda1, double lambda2, int num_active_ys, DTensor2 const& stress );  // After the overshooting, Correct the plastic strain
    virtual double zbrentstress(int num_active_ys,const DTensor2& start_stress,const DTensor2& end_stress,double x1, double x2, double tol) ;



/************************************************************************************
* Yuan @ Thu Mar 2 21:35:05 PST 2017                                                *
* Following Methods are Subclass Responsibility.                                    *
************************************************************************************/
public:
    virtual MultiYieldSurfaceMaterial *getCopy( void );
    virtual DTensor4 const& getTangentTensor( void );       // Return the 3*3*3*3 Tangent tensor.
    virtual void compute_elastoplastic_tangent(int N_active_ys, DTensor2 const& intersection_stress , bool elastic=false);             // Compute the Tangent Stiffness
    virtual void Print( ostream &s, int flag = 0 );
    virtual const char *getType( void ) const;
    virtual int getObjectSize(){return sizeof(*this);};

    
private:
    virtual void update_modulus( int num_active_ys, DTensor2 const& s );               // Compute the Tangent Stiffness
    virtual double yield_surface_val(DTensor2 const& stress, DTensor2 const& alpha, double radius);                    // Return the yield surface Value
    virtual DTensor2 df_dsigma(int num_active_ys, DTensor2 const& stress);     // Return the normal to the yield surface w.r.t stress
    virtual DTensor2 df_dalpha(int num_active_ys, DTensor2 const& stress);     // Return the normal to the yield surface w.r.t alpha(backstress)
    virtual DTensor2 plastic_flow_direct(DTensor2 const& nn, DTensor2 const& stress, int N_active_ys);     // Return the plastic flow
    virtual DTensor2 alpha_bar(int num_active_ys, DTensor2 const& stress);     // Return the rate of alpha(backstress)


    //virtual MultiYieldSurfaceMaterial *getCopy(const char *code);
    // virtual const char *getType(void) const;
    // virtual int getOrder(void) const
    // {
    //     return 1;
    // };
    // virtual int sendSelf(int commitTag, Channel &theChannel);
    // virtual int receiveSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    // // virtual Response *setResponse (const char **argv, int argc, Information &matInformation);
    // // virtual int getResponse (int responseID, Information &matInformation);
    // // virtual int CheckMesh(ofstream &);
    // // virtual void setVolume(double vol);
    // virtual void Print(ostream &s, int flag = 0);
    // virtual int getOutputSize()
    // {
    //     return 0; // Can be overloaded for material point output
    // }
    // virtual void zeroStrain();

/************************************************************************************
* Yuan @ Thu Mar 2 21:35:05 PST 2017                                                *
* Member Variables                                                                  *
************************************************************************************/
protected:
    static MultiYieldSurfaceMaterial_Constitutive_Integration_Method constitutive_integration_method;     //
    static double f_relative_tol;
    static double stress_relative_tol;
    static int n_max_iterations;
    static double allowed_subincrement;

    DTensor2         iterate_stress;                    // Iterative Stress State
    DTensor2         iterate_strain;                    // Iterative Strain State
    DTensor2         iterate_plastic_strain;            // Iterative Plastic Strain State
    vector<DTensor2> iterate_alpha_vec ;                // The vector of alpha(backstress)
    int              iterate_N_active;                  // Iterative Number of Active yield surface

    DTensor2         converge_commit_stress;            // Commit/Stored Stress State
    DTensor2         converge_commit_strain;            // Commit/Stored Strain State
    DTensor2         converge_commit_plastic_strain;    // Commit/Stored Plastic Strain State
    vector<DTensor2> converge_commit_alpha_vec ;        // The vector of alpha(backstress)
    int              converge_commit_N_active;          // Commit/Stored Number of Active yield surface

    DTensor2         save_iter_stress;                  // Commit/Stored Stress State
    DTensor2         save_iter_strain;                  // Commit/Stored Strain State
    DTensor2         save_iter_plastic_strain;          // Commit/Stored Plastic Strain State
    vector<DTensor2> save_iter_alpha_vec ;              // The vector of alpha(backstress)
    int              save_iter_N_active;                // Commit/Stored Number of Active yield surface

    double E;                              // Elastic Modulus
    double v;                              // Poisson Ratio
    double rho;                            // Density
    int TNYS;                              // Total Number of Yield Surface
    vector<double> yield_size;             // Radius of each Yield Surface
    vector<double> HardingPara;            // Hardening Parameter of each Yield Surface
    double initial_E;                      // Initial Young's Modulus

    static const  DTensor2 ZeroStrain;
    static const  DTensor2 ZeroStress;
    static DTensor4 Ee;                    // elastic constant: 3*3*3*3 tensor
    static DTensor4 Eep;                    // elastic constant: 3*3*3*3 tensor
    // static DTensor2 D;                     // elastic constant: 6*6 matrix
    static const DTensor2 kronecker_delta ;// Delta 


    Index < 'i' > i;                       // Dummy or Free Index for LTensor
    Index < 'j' > j;                       // Dummy or Free Index for LTensor
    Index < 'k' > k;                       // Dummy or Free Index for LTensor
    Index < 'l' > l;                       // Dummy or Free Index for LTensor
    Index < 'o' > o;                       // Dummy or Free Index for LTensor
    Index < 't' > t;                       // Dummy or Free Index for LTensor
    Index < 'x' > x;                       // Dummy or Free Index for LTensor
    Index < 'y' > y;                       // Dummy or Free Index for LTensor





// To be deleted.
    static DTensor2 errMatrix;
    static DTensor1 errVector;
    static DTensor2 errTensor;
    static DTensor4 errTensor4;
    static DTensor2 errTensor2;
    static DTensor2 errstresstensor;
    static DTensor2 errstraintensor;
    // static Vector errVectorVector;
    // static vector<double> errvec;
    double volume;



};


#endif
