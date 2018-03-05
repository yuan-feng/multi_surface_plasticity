
void
ClassicElastoplasticityGlobals::flatten_9tensor_to_6matrix(DTensor2 const& source, DTensor1 & target){
    using namespace ClassicElastoplasticityGlobals;
    // %   Reference: William S.Slaughter. The Linearized Theory of Elasticity. 2002.
    // %   Page 202. Section 5.3.1 Engineering Notation.
    target(0) = source(0,0);
    target(1) = source(1,1);
    target(2) = source(2,2);
    target(3) = source(1,2);
    target(4) = source(2,0);
    target(5) = source(0,1);
}

void
ClassicElastoplasticityGlobals::expand_6matrix_to_9tensor(DTensor1 const& source, DTensor2 & target){
    // %   Reference: William S.Slaughter. The Linearized Theory of Elasticity. 2002.
    // %   Page 202. Section 5.3.1 Engineering Notation.
    using namespace ClassicElastoplasticityGlobals;
    target(0,0) = source(0) ;
    target(1,1) = source(1) ;
    target(2,2) = source(2) ;
    target(1,2) = source(3) ;
    target(2,0) = source(4) ;
    target(0,1) = source(5) ;
    // symmetric
    target(2,1) = target(1,2);
    target(0,2) = target(2,0);
    target(1,0) = target(0,1);
}

void
ClassicElastoplasticityGlobals::flatten_81tensor_to_36matrix(DTensor4 const& source, DTensor2 & target){
    // %   Reference: William S.Slaughter. The Linearized Theory of Elasticity. 2002.
    // %   Page 202. Section 5.3.1 Engineering Notation.
    using namespace ClassicElastoplasticityGlobals;
    // % column 1
    target(0,0) = source(0,0,0,0);
    target(1,0) = source(1,1,0,0);
    target(2,0) = source(2,2,0,0);
    target(3,0) = source(1,2,0,0);
    target(4,0) = source(2,0,0,0);
    target(5,0) = source(0,1,0,0);
    // % column 2
    target(0,1) = source(0,0,1,1);
    target(1,1) = source(1,1,1,1);
    target(2,1) = source(2,2,1,1);
    target(3,1) = source(1,2,1,1);
    target(4,1) = source(2,0,1,1);
    target(5,1) = source(0,1,1,1);
    // % column 3
    target(0,2) = source(0,0,2,2);
    target(1,2) = source(1,1,2,2);
    target(2,2) = source(2,2,2,2);
    target(3,2) = source(1,2,2,2);
    target(4,2) = source(2,0,2,2);
    target(5,2) = source(0,1,2,2);
    // % column 4
    target(0,3) = source(0,0,1,2);
    target(1,3) = source(1,1,1,2);
    target(2,3) = source(2,2,1,2);
    target(3,3) = source(1,2,1,2);
    target(4,3) = source(2,0,1,2);
    target(5,3) = source(0,1,1,2);
    // % column 5
    target(0,4) = source(0,0,2,0);
    target(1,4) = source(1,1,2,0);
    target(2,4) = source(2,2,2,0);
    target(3,4) = source(1,2,2,0);
    target(4,4) = source(2,0,2,0);
    target(5,4) = source(0,1,2,0);
    // % column 6
    target(0,5) = source(0,0,0,1);
    target(1,5) = source(1,1,0,1);
    target(2,5) = source(2,2,0,1);
    target(3,5) = source(1,2,0,1);
    target(4,5) = source(2,0,0,1);
    target(5,5) = source(0,1,0,1);
}


void
ClassicElastoplasticityGlobals::expand_36matrix_to_81tensor(DTensor2 const& source, DTensor4 & target){
    // %   Reference: William S.Slaughter. The Linearized Theory of Elasticity. 2002.
    // %   Page 202. Section 5.3.1 Engineering Notation.
    using namespace ClassicElastoplasticityGlobals;
    target(0,0,0,0) = source(0,0);
    target(1,1,0,0) = source(1,0);
    target(2,2,0,0) = source(2,0);
    target(1,2,0,0) = source(3,0);
    target(2,0,0,0) = source(4,0);
    target(0,1,0,0) = source(5,0);
    // % column 2
    target(0,0,1,1) = source(0,1);
    target(1,1,1,1) = source(1,1);
    target(2,2,1,1) = source(2,1);
    target(1,2,1,1) = source(3,1);
    target(2,0,1,1) = source(4,1);
    target(0,1,1,1) = source(5,1);
    // % column 3
    target(0,0,2,2) = source(0,2);
    target(1,1,2,2) = source(1,2);
    target(2,2,2,2) = source(2,2);
    target(1,2,2,2) = source(3,2);
    target(2,0,2,2) = source(4,2);
    target(0,1,2,2) = source(5,2);
    // % column 4
    target(0,0,1,2) = source(0,3);
    target(1,1,1,2) = source(1,3);
    target(2,2,1,2) = source(2,3);
    target(1,2,1,2) = source(3,3);
    target(2,0,1,2) = source(4,3);
    target(0,1,1,2) = source(5,3);
    // % column 5
    target(0,0,2,0) = source(0,4);
    target(1,1,2,0) = source(1,4);
    target(2,2,2,0) = source(2,4);
    target(1,2,2,0) = source(3,4);
    target(2,0,2,0) = source(4,4);
    target(0,1,2,0) = source(5,4);
    // % column 6
    target(0,0,0,1) = source(0,5);
    target(1,1,0,1) = source(1,5);
    target(2,2,0,1) = source(2,5);
    target(1,2,0,1) = source(3,5);
    target(2,0,0,1) = source(4,5);
    target(0,1,0,1) = source(5,5);

    // % ===========================
    // % Minor Symmetric property
    // % [ block 1,  block 2,
    // %   block 3,  block 4 ]
    // % Each block is 3*3 matrix. 
    // % ===========================
    // % (1) symmetric in block 2.
    target(0,0,2,1) = target(0,0,1,2);
    target(1,1,2,1) = target(1,1,1,2);
    target(2,2,2,1) = target(2,2,1,2);
    target(0,0,0,2) = target(0,0,2,0);
    target(1,1,0,2) = target(1,1,2,0);
    target(2,2,0,2) = target(2,2,2,0);
    target(0,0,1,0) = target(0,0,0,1);
    target(1,1,1,0) = target(1,1,0,1);
    target(2,2,1,0) = target(2,2,0,1);
    // % (2) symmetric in block 3.
    target(2,1,0,0) = target(1,2,0,0);
    target(2,1,1,1) = target(1,2,1,1);
    target(2,1,2,2) = target(1,2,2,2);
    target(0,2,0,0) = target(2,0,0,0);
    target(0,2,1,1) = target(2,0,1,1);
    target(0,2,2,2) = target(2,0,2,2);
    target(1,0,0,0) = target(0,1,0,0);
    target(1,0,1,1) = target(0,1,1,1);
    target(1,0,2,2) = target(0,1,2,2);
    // % (3) symmetric in block 4: Change first 2 indices
    target(2,1,1,2) = target(1,2,1,2);
    target(0,2,1,2) = target(2,0,1,2);
    target(1,0,1,2) = target(0,1,1,2);
    target(2,1,2,0) = target(1,2,2,0);
    target(0,2,2,0) = target(2,0,2,0);
    target(1,0,2,0) = target(0,1,2,0);
    target(2,1,0,1) = target(1,2,0,1);
    target(0,2,0,1) = target(2,0,0,1);
    target(1,0,0,1) = target(0,1,0,1);
    // % (4) symmetric in block 4: Change last 2 indices
    target(1,2,2,1) = target(1,2,1,2);
    target(2,0,2,1) = target(2,0,1,2);
    target(0,1,2,1) = target(0,1,1,2);
    target(1,2,0,2) = target(1,2,2,0);
    target(2,0,0,2) = target(2,0,2,0);
    target(0,1,0,2) = target(0,1,2,0);
    target(1,2,1,0) = target(1,2,0,1);
    target(2,0,1,0) = target(2,0,0,1);
    target(0,1,1,0) = target(0,1,0,1);
    // % (5) symmetric in block 4: Change both
    target(2,1,2,1) = target(1,2,1,2);
    target(0,2,2,1) = target(2,0,1,2);
    target(1,0,2,1) = target(0,1,1,2);
    target(2,1,0,2) = target(1,2,2,0);
    target(0,2,0,2) = target(2,0,2,0);
    target(1,0,0,2) = target(0,1,2,0);
    target(2,1,1,0) = target(1,2,0,1);
    target(0,2,1,0) = target(2,0,0,1);
    target(1,0,1,0) = target(0,1,0,1);

}