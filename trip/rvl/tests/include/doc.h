/** \mainpage RVL Unit Tests 

This simple collection of unit tests exercises the RVL base classes, by combining them with LocalRVL concrete subclasses. Since the base classes are not instantiable, there is no other way to test them than simultaneously with concrete subclasses. 

<ul>
<li>ut1: RVL::RnArray, RVL::StdSpace, RVL::RnSpace, and RVL::Vector 
    constructors, also RVL::RVLRandomize</li>
<li>ut2: RVL::RnSpace comparison and inner product.
    constructors, also RVL::RVLAssignConst</li>
<li>ut3: direct test of RVL::StdProductDC behaviour (as opposed 
    to indirect tests via RVL::ProductSpace in UT4). Also tests
    value extraction (method getValue) and retention in 
    RVL::UnaryFunctionObjectScalarRedn class.</li>
<li>ut4: RVL::StdProductSpace, RVL::Components
<li>ut5: RVL::Functional, RVL::GradientTest, apply to square of norm squared</li>
<li>ut6: RVL::ASCIIReader and RVL::ASCIIWriter with RVL::RnArray</li>
<li>ut7: RVL::Operator, RVL::DerivTest, apply to x.*x.*x</li>
<li>ut8: RVL::RnArray<complex>, RVL::StdSpace, RVL::RnSpace<complex>,
and RVL::Vector constructors, also RVL::RVLRandomize for complex types </li>
<li>ut9: RVL::RnSpace<complex> comparison and inner product,
    constructors, also RVL::RVLAssignConst for complex scalars</li>
<li>ut10: direct test of RVL::ScalarFO1 etc., construction of RVL FOs from C functions </li>
<li>utll: use RVL::GenMat to exercise RVL::AdjointTest</li>
<li>ut12: detection of input/domain mismatch in RVL::LinearOp</li>
<li>ut13: test RVL::LinearOpFO by comparison with explicit construction (RVL::GenMat)</li>
<li>ut14: detection of input/domain mismatch in RVL::Operator</li>
</ul>

<hr>
<a href="../../../doc/html/index.html">RVL Home Page</a>

*/
