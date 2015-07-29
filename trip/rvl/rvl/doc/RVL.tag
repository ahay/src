<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>The Rice Vector Library</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>adjtest.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>adjtest_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>Users williamsymes Applications RSFSRC trip rvl rvl include adjtest hh bool</type>
      <name>AdjointTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ac76d7d4ea39f3eccaad9c8c571aab5d5</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;op, FunctionObject &amp;randomize, ostream &amp;str, int tol=100)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>blockop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>blockop_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
    <class kind="class">RVL::BlockOperator</class>
    <class kind="class">RVL::TensorOp</class>
    <class kind="class">RVL::BlockLinearOp</class>
    <class kind="class">RVL::TensorLinearOp</class>
    <class kind="class">RVL::InjectOp</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>data.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>data_8hh</filename>
    <includes id="utility_8hh" name="utility.hh" local="yes" imported="no">utility.hh</includes>
    <class kind="class">RVL::FunctionObject</class>
    <class kind="class">RVL::FunctionObjectConstEval</class>
    <class kind="class">RVL::ScalarRedn</class>
    <class kind="class">RVL::FunctionObjectScalarRedn</class>
    <class kind="class">RVL::DataContainer</class>
    <class kind="class">RVL::DataContainerFactory</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>derivtest.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>derivtest_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>bool</type>
      <name>DerivTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6859e81cb61ca706a5440de06bfb4dcf</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op, Vector&lt; Scalar &gt; const &amp;y, Vector&lt; Scalar &gt; const &amp;p, ostream &amp;str, int n=10, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmin=0.1, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmax=1.0, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType minrat=1.95)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>except.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>except_8hh</filename>
    <includes id="std__cpp__includes_8hh" name="std_cpp_includes.hh" local="yes" imported="no">std_cpp_includes.hh</includes>
    <class kind="class">RVL::RVLException</class>
    <namespace>RVL</namespace>
    <member kind="define">
      <type>#define</type>
      <name>BUFLEN</name>
      <anchorfile>except_8hh.html</anchorfile>
      <anchor>ad974fe981249f5e84fbf1683b012c9f8</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>functional.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>functional_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
    <class kind="class">RVL::Functional</class>
    <class kind="class">RVL::FunctionalProductDomain</class>
    <class kind="class">RVL::FunctionalEvaluation</class>
    <class kind="class">RVL::HessianEvaluation</class>
    <class kind="class">RVL::FunctionalProductDomainEvaluation</class>
    <class kind="class">RVL::HessianBlockEvaluation</class>
    <class kind="class">RVL::LinCombFunctional</class>
    <class kind="class">RVL::StdFOFunctional</class>
    <class kind="class">RVL::NullFunctional</class>
    <class kind="class">RVL::FcnlOpComp</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>gradtest.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>gradtest_8hh</filename>
    <includes id="functional_8hh" name="functional.hh" local="yes" imported="no">functional.hh</includes>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>Users williamsymes Applications RSFSRC trip rvl rvl include gradtest hh bool</type>
      <name>GradientTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8beb1ab0cae2fa0be8a91bec7bf90232</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, ostream &amp;str, int n=11, Scalar hmin=0.1, Scalar hmax=1.0, Scalar minrat=1.95)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>linalg.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>linalg_8hh</filename>
    <includes id="data_8hh" name="data.hh" local="yes" imported="no">data.hh</includes>
    <includes id="write_8hh" name="write.hh" local="yes" imported="no">write.hh</includes>
    <class kind="class">RVL::LinCombObject</class>
    <class kind="class">RVL::LinearAlgebraPackage</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>linop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>linop_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
  </compound>
  <compound kind="file">
    <name>linop_base.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>linop__base_8hh</filename>
    <includes id="space_8hh" name="space.hh" local="yes" imported="no">space.hh</includes>
    <includes id="write_8hh" name="write.hh" local="yes" imported="no">write.hh</includes>
    <class kind="class">RVL::LinearOp</class>
    <class kind="class">RVL::LinearOpFO</class>
    <class kind="class">RVL::Invertible</class>
    <class kind="class">RVL::LinearOpWithInverse</class>
    <class kind="class">RVL::AdjLinearOp</class>
    <class kind="class">RVL::NormalLinearOp</class>
    <class kind="class">RVL::ScaleOpFwd</class>
    <class kind="class">RVL::ScaleOpInv</class>
    <class kind="class">RVL::LinCombLinearOp</class>
    <class kind="class">RVL::CompLinearOp</class>
    <class kind="class">RVL::SymmetricBilinearOp</class>
    <class kind="class">RVL::LinearBilinearOp</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>ls.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>ls_8hh</filename>
    <includes id="op_8hh" name="op.hh" local="yes" imported="no">op.hh</includes>
    <includes id="functional_8hh" name="functional.hh" local="yes" imported="no">functional.hh</includes>
    <class kind="class">RVL::ShiftOperator</class>
    <class kind="class">RVL::ResidualOperator</class>
    <class kind="class">RVL::EuclideanForm</class>
    <class kind="class">RVL::QuadraticForm</class>
    <class kind="class">RVL::ShiftedQuadraticForm</class>
    <class kind="class">RVL::LeastSquaresFcnlGN</class>
    <class kind="class">RVL::StdLeastSquaresFcnlGN</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>op.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>op_8hh</filename>
    <includes id="space_8hh" name="space.hh" local="yes" imported="no">space.hh</includes>
    <includes id="linop__base_8hh" name="linop_base.hh" local="yes" imported="no">linop_base.hh</includes>
    <includes id="productspace_8hh" name="productspace.hh" local="yes" imported="no">productspace.hh</includes>
    <includes id="write_8hh" name="write.hh" local="yes" imported="no">write.hh</includes>
    <class kind="class">RVL::Operator</class>
    <class kind="class">RVL::OperatorProductDomain</class>
    <class kind="class">RVL::OperatorWithInvertibleDeriv</class>
    <class kind="class">RVL::OperatorEvaluation</class>
    <class kind="class">RVL::DerivEvaluation</class>
    <class kind="class">RVL::Deriv2Evaluation</class>
    <class kind="class">RVL::InvertibleDerivEvaluation</class>
    <class kind="class">RVL::OperatorProductDomainEvaluation</class>
    <class kind="class">RVL::PartialDerivEvaluation</class>
    <class kind="class">RVL::LNLOperator</class>
    <class kind="class">RVL::ANLOperator</class>
    <class kind="class">RVL::OpFO</class>
    <class kind="class">RVL::LinCombOperator</class>
    <class kind="class">RVL::LinearOpEvaluation</class>
    <class kind="class">RVL::LinearOpAdjEvaluation</class>
    <class kind="class">RVL::OpComp</class>
    <class kind="class">RVL::IdentityOp</class>
    <namespace>RVL</namespace>
    <member kind="define">
      <type>#define</type>
      <name>RVL_OPERATOR_NEW_ENABLED</name>
      <anchorfile>op_8hh.html</anchorfile>
      <anchor>ae02bb4645579c63b119f7502f0f6ac81</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>product.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>product_8hh</filename>
    <includes id="utility_8hh" name="utility.hh" local="yes" imported="no">utility.hh</includes>
    <includes id="except_8hh" name="except.hh" local="yes" imported="no">except.hh</includes>
    <class kind="class">RVL::Product</class>
    <class kind="class">RVL::ROProduct</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>productdata.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>productdata_8hh</filename>
    <includes id="data_8hh" name="data.hh" local="yes" imported="no">data.hh</includes>
    <includes id="product_8hh" name="product.hh" local="yes" imported="no">product.hh</includes>
    <class kind="class">RVL::BlockFunctionObject</class>
    <class kind="class">RVL::DiagonalFunctionObject</class>
    <class kind="class">RVL::ProductDataContainer</class>
    <class kind="class">RVL::StdProductDataContainer</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>productspace.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>productspace_8hh</filename>
    <includes id="space_8hh" name="space.hh" local="yes" imported="no">space.hh</includes>
    <includes id="productdata_8hh" name="productdata.hh" local="yes" imported="no">productdata.hh</includes>
    <class kind="class">RVL::ProductSpace</class>
    <class kind="class">RVL::StdProductSpace</class>
    <class kind="class">RVL::CartesianPowerSpace</class>
    <class kind="class">RVL::Components</class>
    <namespace>RVL</namespace>
  </compound>
  <compound kind="file">
    <name>scantest.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>scantest_8hh</filename>
    <includes id="functional_8hh" name="functional.hh" local="yes" imported="no">functional.hh</includes>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>void</type>
      <name>Scan</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a259968906b55eccbf5956b532f6b3a15</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, int n=11, Scalar hmin=-ScalarFieldTraits&lt; Scalar &gt;::One(), Scalar hmax=ScalarFieldTraits&lt; Scalar &gt;::One(), ostream &amp;str=cout)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>space.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>space_8hh</filename>
    <includes id="data_8hh" name="data.hh" local="yes" imported="no">data.hh</includes>
    <includes id="linalg_8hh" name="linalg.hh" local="yes" imported="no">linalg.hh</includes>
    <includes id="write_8hh" name="write.hh" local="yes" imported="no">write.hh</includes>
    <class kind="class">RVL::Space</class>
    <class kind="class">RVL::StdSpace</class>
    <class kind="class">RVL::SpaceDCF</class>
    <class kind="class">RVL::Vector</class>
    <class kind="class">RVL::WatchedVecRef</class>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>void</type>
      <name>SpaceTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a634631fe354f7675d74a33d340d80167</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;sp, Vector&lt; Scalar &gt; const &amp;v, std::string msg)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>std_cpp_includes.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>std__cpp__includes_8hh</filename>
    <member kind="define">
      <type>#define</type>
      <name>F77NAME</name>
      <anchorfile>std__cpp__includes_8hh.html</anchorfile>
      <anchor>a9b91c72ded7de2f0b393ea75de53929b</anchor>
      <arglist>(x)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>utility.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>utility_8hh</filename>
    <includes id="except_8hh" name="except.hh" local="yes" imported="no">except.hh</includes>
    <class kind="struct">RVL::ScalarFieldTraits</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; bool &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; int &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; long &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; unsigned int &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; unsigned long &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; float &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; double &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; std::complex&lt; T &gt; &gt;</class>
    <class kind="class">RVL::Writeable</class>
    <class kind="class">RVL::Oracle</class>
    <class kind="class">RVL::Factory</class>
    <namespace>RVL</namespace>
    <member kind="function">
      <type>int</type>
      <name>numeric_precision&lt; float &gt;</name>
      <anchorfile>utility_8hh.html</anchorfile>
      <anchor>aa8f442077effd206a3d9b158d46ff900</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>numeric_precision&lt; double &gt;</name>
      <anchorfile>utility_8hh.html</anchorfile>
      <anchor>a01062411a447523607f3999336243d40</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>testRealOnly</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>af258e947fc55349d6ab55bede76b4ae8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ProtectedDivision</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a3857f0e1097eff88c8c1a7e7feaa4c1e</anchor>
      <arglist>(real a, real b, real &amp;quot, real tol=ScalarFieldTraits&lt; real &gt;::AbsZero())</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>write.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/rvl/include/</path>
    <filename>write_8hh</filename>
    <includes id="utility_8hh" name="utility.hh" local="yes" imported="no">utility.hh</includes>
  </compound>
  <compound kind="namespace">
    <name>RVL</name>
    <filename>namespaceRVL.html</filename>
    <class kind="class">RVL::BlockOperator</class>
    <class kind="class">RVL::TensorOp</class>
    <class kind="class">RVL::BlockLinearOp</class>
    <class kind="class">RVL::TensorLinearOp</class>
    <class kind="class">RVL::InjectOp</class>
    <class kind="class">RVL::FunctionObject</class>
    <class kind="class">RVL::FunctionObjectConstEval</class>
    <class kind="class">RVL::ScalarRedn</class>
    <class kind="class">RVL::FunctionObjectScalarRedn</class>
    <class kind="class">RVL::DataContainer</class>
    <class kind="class">RVL::DataContainerFactory</class>
    <class kind="class">RVL::RVLException</class>
    <class kind="class">RVL::Functional</class>
    <class kind="class">RVL::FunctionalProductDomain</class>
    <class kind="class">RVL::FunctionalEvaluation</class>
    <class kind="class">RVL::HessianEvaluation</class>
    <class kind="class">RVL::FunctionalProductDomainEvaluation</class>
    <class kind="class">RVL::HessianBlockEvaluation</class>
    <class kind="class">RVL::LinCombFunctional</class>
    <class kind="class">RVL::StdFOFunctional</class>
    <class kind="class">RVL::NullFunctional</class>
    <class kind="class">RVL::FcnlOpComp</class>
    <class kind="class">RVL::LinCombObject</class>
    <class kind="class">RVL::LinearAlgebraPackage</class>
    <class kind="class">RVL::LinearOp</class>
    <class kind="class">RVL::LinearOpFO</class>
    <class kind="class">RVL::Invertible</class>
    <class kind="class">RVL::LinearOpWithInverse</class>
    <class kind="class">RVL::AdjLinearOp</class>
    <class kind="class">RVL::NormalLinearOp</class>
    <class kind="class">RVL::ScaleOpFwd</class>
    <class kind="class">RVL::ScaleOpInv</class>
    <class kind="class">RVL::LinCombLinearOp</class>
    <class kind="class">RVL::CompLinearOp</class>
    <class kind="class">RVL::SymmetricBilinearOp</class>
    <class kind="class">RVL::LinearBilinearOp</class>
    <class kind="class">RVL::ShiftOperator</class>
    <class kind="class">RVL::ResidualOperator</class>
    <class kind="class">RVL::EuclideanForm</class>
    <class kind="class">RVL::QuadraticForm</class>
    <class kind="class">RVL::ShiftedQuadraticForm</class>
    <class kind="class">RVL::LeastSquaresFcnlGN</class>
    <class kind="class">RVL::StdLeastSquaresFcnlGN</class>
    <class kind="class">RVL::Operator</class>
    <class kind="class">RVL::OperatorProductDomain</class>
    <class kind="class">RVL::OperatorWithInvertibleDeriv</class>
    <class kind="class">RVL::OperatorEvaluation</class>
    <class kind="class">RVL::DerivEvaluation</class>
    <class kind="class">RVL::Deriv2Evaluation</class>
    <class kind="class">RVL::InvertibleDerivEvaluation</class>
    <class kind="class">RVL::OperatorProductDomainEvaluation</class>
    <class kind="class">RVL::PartialDerivEvaluation</class>
    <class kind="class">RVL::LNLOperator</class>
    <class kind="class">RVL::ANLOperator</class>
    <class kind="class">RVL::OpFO</class>
    <class kind="class">RVL::LinCombOperator</class>
    <class kind="class">RVL::LinearOpEvaluation</class>
    <class kind="class">RVL::LinearOpAdjEvaluation</class>
    <class kind="class">RVL::OpComp</class>
    <class kind="class">RVL::IdentityOp</class>
    <class kind="class">RVL::Product</class>
    <class kind="class">RVL::ROProduct</class>
    <class kind="class">RVL::BlockFunctionObject</class>
    <class kind="class">RVL::DiagonalFunctionObject</class>
    <class kind="class">RVL::ProductDataContainer</class>
    <class kind="class">RVL::StdProductDataContainer</class>
    <class kind="class">RVL::ProductSpace</class>
    <class kind="class">RVL::StdProductSpace</class>
    <class kind="class">RVL::CartesianPowerSpace</class>
    <class kind="class">RVL::Components</class>
    <class kind="class">RVL::Space</class>
    <class kind="class">RVL::StdSpace</class>
    <class kind="class">RVL::SpaceDCF</class>
    <class kind="class">RVL::Vector</class>
    <class kind="class">RVL::WatchedVecRef</class>
    <class kind="struct">RVL::ScalarFieldTraits</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; bool &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; int &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; long &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; unsigned int &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; unsigned long &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; float &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; double &gt;</class>
    <class kind="struct">RVL::ScalarFieldTraits&lt; std::complex&lt; T &gt; &gt;</class>
    <class kind="class">RVL::Writeable</class>
    <class kind="class">RVL::Oracle</class>
    <class kind="class">RVL::Factory</class>
    <member kind="function">
      <type>Users williamsymes Applications RSFSRC trip rvl rvl include adjtest hh bool</type>
      <name>AdjointTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ac76d7d4ea39f3eccaad9c8c571aab5d5</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;op, FunctionObject &amp;randomize, ostream &amp;str, int tol=100)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>DerivTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a6859e81cb61ca706a5440de06bfb4dcf</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op, Vector&lt; Scalar &gt; const &amp;y, Vector&lt; Scalar &gt; const &amp;p, ostream &amp;str, int n=10, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmin=0.1, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType hmax=1.0, typename ScalarFieldTraits&lt; Scalar &gt;::AbsType minrat=1.95)</arglist>
    </member>
    <member kind="function">
      <type>Users williamsymes Applications RSFSRC trip rvl rvl include gradtest hh bool</type>
      <name>GradientTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8beb1ab0cae2fa0be8a91bec7bf90232</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, ostream &amp;str, int n=11, Scalar hmin=0.1, Scalar hmax=1.0, Scalar minrat=1.95)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Scan</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a259968906b55eccbf5956b532f6b3a15</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;p, int n=11, Scalar hmin=-ScalarFieldTraits&lt; Scalar &gt;::One(), Scalar hmax=ScalarFieldTraits&lt; Scalar &gt;::One(), ostream &amp;str=cout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SpaceTest</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a634631fe354f7675d74a33d340d80167</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;sp, Vector&lt; Scalar &gt; const &amp;v, std::string msg)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>testRealOnly</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>af258e947fc55349d6ab55bede76b4ae8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ProtectedDivision</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a3857f0e1097eff88c8c1a7e7feaa4c1e</anchor>
      <arglist>(real a, real b, real &amp;quot, real tol=ScalarFieldTraits&lt; real &gt;::AbsZero())</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BlockOperator</name>
    <filename>classRVL_1_1BlockOperator.html</filename>
    <templarg></templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>BlockOperator</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>afc4a84943d85f9d089ee81bc349f8429</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BlockOperator</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ae0dff0e047479cdbc45c73151822a715</anchor>
      <arglist>(const BlockOperator&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BlockOperator</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ab35cd69c507597c0ce5eb6965eef2684</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductRange</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ac5ed425ebf422f99eedbf831d4ece783</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a200ba4badfa32adbea844fa2c5a7edde</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponent</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a7603de185c99e4104da0567cb2887698</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;yi) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ad86ebaee1bc9e16f891b0d4962e3e98d</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponentDeriv</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ade65f4efb2a7d61371245950f1abae06</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dyi) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a56669f8c09c666dad7d3967e256b5326</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponentAdjDeriv</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ad7ccd9d115f8b07b0a559063b980c0c8</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dyi, Vector&lt; Scalar &gt; &amp;dx) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ae6f42cbf3c3569ed6f1383638679a4a4</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponentDeriv2</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a5427420e6f8f17259910e0b661a36806</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dx1, Vector&lt; Scalar &gt; &amp;dyi) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv2</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>adbc8c162c28ab288d2d30fef50ff3ce6</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dx1, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponentAdjDeriv2</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a90686e68d8af1e2210d3bb69ab9019fc</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dyi, Vector&lt; Scalar &gt; &amp;dx1) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ab7792896b761a8c1c66b051508492ff9</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual BlockOperator&lt; Scalar &gt; *</type>
      <name>cloneBlockOp</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ababae89ddfe0fa0ef422134a089a623f</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>ac61848def4778506dae45087a939dbc4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1BlockOperator.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TensorOp</name>
    <filename>classRVL_1_1TensorOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::BlockOperator</base>
    <member kind="function">
      <type></type>
      <name>TensorOp</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a056c0fe7834ae12858905645a10eb2f6</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;_op1, Operator&lt; Scalar &gt; const &amp;_op2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TensorOp</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>aab709534ce6540e6fd820841fbc17671</anchor>
      <arglist>(TensorOp&lt; Scalar &gt; const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~TensorOp</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a9500a9fabf83a22e268876e4f9b6e50e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a11da1e26115db0cf668b36cccb301a05</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ProductSpace&lt; Scalar &gt; const &amp;</type>
      <name>getProductRange</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a7ed8ceeea09952cf822358b04202b1c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>ace42b31ba750a89af2f25a3be5e58f29</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponent</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>af9e61817128abb63ac16a9a3947c003c</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;yi) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponentDeriv</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a1ab56e15c7d0b010b14efc07899d464b</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dyi) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponentAdjDeriv</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>aa0e5c278ce6c1ac60a591f5ba6ed8537</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dyi, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponentDeriv2</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>ab0e609a1ab6bd7ae4721f0a92ebc4ac0</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dx1, Vector&lt; Scalar &gt; &amp;dyi) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponentAdjDeriv2</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>ac56062bd394eb33404f96d929f7c09ce</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dyi, Vector&lt; Scalar &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>TensorOp&lt; Scalar &gt; *</type>
      <name>cloneTensorOp</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a570dc13892c72fe64a75237817e76d50</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>BlockOperator&lt; Scalar &gt; *</type>
      <name>cloneBlockOp</name>
      <anchorfile>classRVL_1_1TensorOp.html</anchorfile>
      <anchor>a472bc62e4d7c71bb9145cd1fefdb7f15</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BlockLinearOp</name>
    <filename>classRVL_1_1BlockLinearOp.html</filename>
    <templarg></templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>BlockLinearOp</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>af6457e59124aabc967bad8e56ab4eeee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BlockLinearOp</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a0b620baf158b311e2f845e0b62af7b2d</anchor>
      <arglist>(const BlockLinearOp&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BlockLinearOp</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a173a8dee125367c2d01c8ebdc3da0881</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductRange</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a673e924e837adab72552100545075973</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>af98b9a0b8035cb47bf8c9163e85445af</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponent</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>ab49e496bcef9bb57e52e19abf7d64993</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;yi) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>aca2e58b2518642bf0b44e2594f5aec6e</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyComponentAdj</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a06e481cccdd693f46a6c47fd64bd3b15</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;yi, Vector&lt; Scalar &gt; &amp;x) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a000e5e68b4b2a552c0b16a5f46e9033e</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual BlockLinearOp&lt; Scalar &gt; *</type>
      <name>cloneBlockLinearOp</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>acc4e42b8600e9061770f9fc337d9d15a</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a3dee4cc08fa8b5d74f0feb37da3a657d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1BlockLinearOp.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TensorLinearOp</name>
    <filename>classRVL_1_1TensorLinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::BlockLinearOp</base>
    <member kind="function">
      <type></type>
      <name>TensorLinearOp</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>afd7699404f75e6c93b42fc98aab65557</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;_op1, LinearOp&lt; Scalar &gt; const &amp;_op2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TensorLinearOp</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a5e7ee3ce5fae8412be37f11d822cc2c6</anchor>
      <arglist>(TensorLinearOp&lt; Scalar &gt; const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~TensorLinearOp</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>ab7afa22ee811f94e0f8de96fab966ae3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>ab50845eae7c6eccbf7db69d21ce09d9e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ProductSpace&lt; Scalar &gt; const &amp;</type>
      <name>getProductRange</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a6dffa1adb2fa9ac5750fce6325a8e73b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a64871ea4c726bbb46854e46e62cae1ea</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponent</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a1759d1f9a97fb343d3d63b1a9fb1057a</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;yi) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyComponentAdj</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a0050a05b42976678817d7e6988669557</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;yi, Vector&lt; Scalar &gt; &amp;x) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>TensorLinearOp&lt; Scalar &gt; *</type>
      <name>cloneTensorLinearOp</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>ae04e2bdcf0aa24505c933ab76f89760c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>BlockLinearOp&lt; Scalar &gt; *</type>
      <name>cloneBlockLinearOp</name>
      <anchorfile>classRVL_1_1TensorLinearOp.html</anchorfile>
      <anchor>a3aec1664019c8c80a941c48280ef8453</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::InjectOp</name>
    <filename>classRVL_1_1InjectOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>InjectOp</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a0d367af70e93734e1ea78a6f6b675ce9</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;_ref, int _icomp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InjectOp</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>ac264d8acddd5eea1d260422ae2bfd4b4</anchor>
      <arglist>(InjectOp&lt; Scalar &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InjectOp</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a4ef7a8330713f25218338673233a3535</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a434aa0b0bb8333d6fff7006cb6d45881</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a105ae19e48094f0f077e219f90d1c93a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a629fdc28a387bf1b9b4fcd4eb73771d4</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a2b10c09ea8fade33cb31228a64cc9015</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a23b943e14f65758101e51484719a2432</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; const &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>a49eb9b6b36b070790ef58c4b49789d5d</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; const &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1InjectOp.html</anchorfile>
      <anchor>ae2390b1a44cb862dd1d7181a987cca3f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionObject</name>
    <filename>classRVL_1_1FunctionObject.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>FunctionObject</name>
      <anchorfile>classRVL_1_1FunctionObject.html</anchorfile>
      <anchor>a97817e8039e53bf7b2539a4038d33dc2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionObject</name>
      <anchorfile>classRVL_1_1FunctionObject.html</anchorfile>
      <anchor>a1838034746c30e8a636fb447d9f6c400</anchor>
      <arglist>(const FunctionObject &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FunctionObject</name>
      <anchorfile>classRVL_1_1FunctionObject.html</anchorfile>
      <anchor>ab63b5924cf26a9c9f583fbc7afbac5d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1FunctionObject.html</anchorfile>
      <anchor>a1dce09a66af5a455b7a59b450fa4e9cc</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1FunctionObject.html</anchorfile>
      <anchor>a940b1573eb25dff31417027a4bbb096a</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionObjectConstEval</name>
    <filename>classRVL_1_1FunctionObjectConstEval.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FunctionObjectConstEval</name>
      <anchorfile>classRVL_1_1FunctionObjectConstEval.html</anchorfile>
      <anchor>a70df25837e1eddf198c46c1ee1220157</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1FunctionObjectConstEval.html</anchorfile>
      <anchor>ab68357de37f3879c3bd1ea527876ee42</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1FunctionObjectConstEval.html</anchorfile>
      <anchor>ae59eea7076f87bfd32b6b1a3d374a126</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarRedn</name>
    <filename>classRVL_1_1ScalarRedn.html</filename>
    <templarg>Scalar</templarg>
    <member kind="function">
      <type></type>
      <name>ScalarRedn</name>
      <anchorfile>classRVL_1_1ScalarRedn.html</anchorfile>
      <anchor>adb5d47839e7502f616134f1d70c14387</anchor>
      <arglist>(Scalar _val)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScalarRedn</name>
      <anchorfile>classRVL_1_1ScalarRedn.html</anchorfile>
      <anchor>aedb127e309a3e08c6809bd0c09be0fd2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1ScalarRedn.html</anchorfile>
      <anchor>a4f7c1b0ee18f7946a33bc47f056f31ba</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1ScalarRedn.html</anchorfile>
      <anchor>a5e5b4f87190298c5c46a2dc41bdb6e8b</anchor>
      <arglist>(Scalar _val)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Scalar</type>
      <name>getValue</name>
      <anchorfile>classRVL_1_1ScalarRedn.html</anchorfile>
      <anchor>ac9572ac3f03529538ddb49972ae80579</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionObjectScalarRedn</name>
    <filename>classRVL_1_1FunctionObjectScalarRedn.html</filename>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectConstEval</base>
    <base>ScalarRedn&lt; ValType &gt;</base>
    <member kind="function">
      <type></type>
      <name>FunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1FunctionObjectScalarRedn.html</anchorfile>
      <anchor>aa01cb662f42d282465228b3b63337f05</anchor>
      <arglist>(ValType val)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1FunctionObjectScalarRedn.html</anchorfile>
      <anchor>a476eafdd5ffe15a1839cfa387ca37673</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; ValType &gt; &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1FunctionObjectScalarRedn.html</anchorfile>
      <anchor>a54d5d3318fe03b0bf8ff0ccba0c8080f</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::DataContainer</name>
    <filename>classRVL_1_1DataContainer.html</filename>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>DataContainer</name>
      <anchorfile>classRVL_1_1DataContainer.html</anchorfile>
      <anchor>a95219c8a1c0d8db8ef4ada9f1aa96bbc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DataContainer</name>
      <anchorfile>classRVL_1_1DataContainer.html</anchorfile>
      <anchor>a315fa87a57cd485c5f1d35912973bc94</anchor>
      <arglist>(const DataContainer &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DataContainer</name>
      <anchorfile>classRVL_1_1DataContainer.html</anchorfile>
      <anchor>ada51681caa78e62b85374094739182e6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1DataContainer.html</anchorfile>
      <anchor>a9f5f9041dc642da1184ade0032e8975a</anchor>
      <arglist>(FunctionObject &amp;f, vector&lt; DataContainer const * &gt; &amp;x)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1DataContainer.html</anchorfile>
      <anchor>ac7ac7e190d3b13b8be9367014af806eb</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, vector&lt; DataContainer const * &gt; &amp;x) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::DataContainerFactory</name>
    <filename>classRVL_1_1DataContainerFactory.html</filename>
    <base>Factory&lt; DataContainer &gt;</base>
    <member kind="function">
      <type></type>
      <name>DataContainerFactory</name>
      <anchorfile>classRVL_1_1DataContainerFactory.html</anchorfile>
      <anchor>a7f8f85cceceab2d9399b13f3d87f61c4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DataContainerFactory</name>
      <anchorfile>classRVL_1_1DataContainerFactory.html</anchorfile>
      <anchor>ad20db453bb02f34fa731ac40bd7a4f1e</anchor>
      <arglist>(const DataContainerFactory &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DataContainerFactory</name>
      <anchorfile>classRVL_1_1DataContainerFactory.html</anchorfile>
      <anchor>a0beda6099308a690001c80a510da7d31</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1DataContainerFactory.html</anchorfile>
      <anchor>aa7b76163dfe794d071cebebcf66a2721</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1DataContainerFactory.html</anchorfile>
      <anchor>acf498bf68cc0cf6b4d180a9bb1316640</anchor>
      <arglist>(DataContainer const &amp;dc) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLException</name>
    <filename>classRVL_1_1RVLException.html</filename>
    <member kind="function">
      <type></type>
      <name>RVLException</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a05686e9caa6ffed25707bfdbf086aa19</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLException</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a7dbf620b8f0e39dc6ee4fd6cc51c5525</anchor>
      <arglist>(const RVLException &amp;s)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~RVLException</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a460c054c792f1697b4013074e9c6761e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const char *</type>
      <name>what</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a0f1c35d1513952979345d83dde186359</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a40f5dc6eb8c110b330bb54e41d773157</anchor>
      <arglist>(string str)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>aaa1cdcec78c0ad9962687f78da67e518</anchor>
      <arglist>(const char *str)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a271b34519ff955373632bce73cd0113a</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>af6a89879d7d26cffe27657eff56b6d31</anchor>
      <arglist>(unsigned int i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>add053b1c6e68678dddfe59f440482d05</anchor>
      <arglist>(long i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>ab811597847cd3e91ddc37670d587a185</anchor>
      <arglist>(unsigned long i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>ab0f5d5fb36b968723d490bbc5e142af8</anchor>
      <arglist>(short i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a88c57535be71d047d4674b213fdb5ef8</anchor>
      <arglist>(unsigned short i)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a8570c9fe9dccee272fca702a7d71b23c</anchor>
      <arglist>(double d)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>aa40ab83faa299e823e82470b13972fa0</anchor>
      <arglist>(float d)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>ada42c074373251b8ec337104860aa831</anchor>
      <arglist>(complex&lt; T &gt; d)</arglist>
    </member>
    <member kind="function">
      <type>RVLException &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a78d2eb21b648c5a09a86aa19a35987c9</anchor>
      <arglist>(char c)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RVLException.html</anchorfile>
      <anchor>a992ea89d56a73105f418078880dfaebf</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Functional</name>
    <filename>classRVL_1_1Functional.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Functional</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a9a1d5db140296abad62b4de120d14658</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Functional</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a8a9f31d0051d8f7eab73b3eccf58dee0</anchor>
      <arglist>(const Functional&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Functional</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a767bae3d132e34600c5b8cd69aa90abc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>af2a918a855eb1c4b758259307b1a3fbf</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a34053a686f35dfe2b3b28827711e34e0</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>ad020543801a0aad3cad5143d004d420b</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a80189ed26b9f11ba15ff66bbc0486b5c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>aaf198a0b43903e89bf26c344638e24ff</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_apply</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a6e458c93c28d76df3f9ef7be0626c6d4</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyGradient</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a3eab7bcadb29320ecacdbb054885fa11</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyHessian</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a3d2762ccb0fc4e5a4cb8de5a4d376788</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a3ac05b4ef2e69a6ae0a8b2697c9abe06</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_clone</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>ae3bd243aa8c9253d6921e6611a5f47f1</anchor>
      <arglist>(Functional&lt; Scalar &gt; const &amp;fref, Functional&lt; Scalar &gt; **f) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FunctionalEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a6c50312f72d5be03c622c8b1492762ed</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FcnlOpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>ad6486a678792a5664c5c427696939eed</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>LinCombFunctional&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Functional.html</anchorfile>
      <anchor>a150324d676f92713a2363c6b71052484</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionalProductDomain</name>
    <filename>classRVL_1_1FunctionalProductDomain.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>FunctionalProductDomain</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>ad595eefdd44847ecee831574512e4923</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionalProductDomain</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a4b1e9382a47a177f765e85449371e1e8</anchor>
      <arglist>(const FunctionalProductDomain&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FunctionalProductDomain</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a9135c7bcb11e16acb827cd05d1f3904b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a44c60106c7c79347dd5609775df69ba0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductDomain</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>ae257eaff1f05241fd86d2081c6e0516a</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyPartialGradient</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a09cea2f1cdf74a7dfea9750b137d9fd1</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a044a64954104d92c3cd896312a852ba3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyHessianBlock</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a37ff474bd624e17731f43e1b74e82fab</anchor>
      <arglist>(int i, int j, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dxi, Vector&lt; Scalar &gt; &amp;dxj) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a7aa4a6c2914596e55d6d6479020edbcf</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual FunctionalProductDomain&lt; Scalar &gt; *</type>
      <name>clonePD</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>aad83787814470e2c04939b59c2b8e21f</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>ab17616136a018374d4f5772ec4baa966</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FunctionalProductDomainEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1FunctionalProductDomain.html</anchorfile>
      <anchor>a5552522d2f0dcd96ee525e09d7e2dc87</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionalEvaluation</name>
    <filename>classRVL_1_1FunctionalEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>FunctionalEvaluation</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a950f55878463ce6f98549586edf58be1</anchor>
      <arglist>(const Functional&lt; Scalar &gt; &amp;_f, const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FunctionalEvaluation</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a25dd29f6000336fbcc973e186362c64d</anchor>
      <arglist>(const FunctionalEvaluation&lt; Scalar &gt; &amp;ev)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~FunctionalEvaluation</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a0d968e0c6a22a927d10ed61787cc5266</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a624100288faf5d1e51e655be5298ffda</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a32d0fdff8faff04c88840f73e4906abe</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>getPoint</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>aad73b2b6e6db1ce68f09f163138509eb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getPoint</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a51dcf69e290082d4638adb84eb41ff8f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getValue</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a52b5a26d54519c12453440dbecb4a066</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getGradient</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a8548812dc9914fbc3c892050cb379d7d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getGradientNorm</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a0a6d35819d346dd9b35a659bec6330e2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LinearOp&lt; Scalar &gt; const &amp;</type>
      <name>getHessian</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a0ea4ef7765c8fa229f0c9406dc3113e7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Functional&lt; Scalar &gt; const &amp;</type>
      <name>getFunctional</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a589dc2bb9a88cfb3783716107ffaba20</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a1df41fa497da83dc7630cef8a94c8604</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a46192c1fa63967188d77c14f283dec80</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductDomain</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>ad3669c1c248517db243a0b0262474ecc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const Vector&lt; Scalar &gt; &amp;</type>
      <name>getGradientBlock</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>aa745004b207244166c5b931747b3cee6</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessianBlock</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a5c7abe5ca459e5be4f9f0e5140fe44ab</anchor>
      <arglist>(int i, int j, const Vector&lt; Scalar &gt; &amp;dxi, Vector&lt; Scalar &gt; &amp;dxj) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>HessianEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a16e3d2dbbaac61d80b04fb3e03c09528</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>HessianBlockEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>a0734fbfe745826211ee788c473ae9131</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FcnlOpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1FunctionalEvaluation.html</anchorfile>
      <anchor>ad6486a678792a5664c5c427696939eed</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::HessianEvaluation</name>
    <filename>classRVL_1_1HessianEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>HessianEvaluation</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>acdbc986451a56a83d4dc8bbd127a43d4</anchor>
      <arglist>(FunctionalEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~HessianEvaluation</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>aca121bc051b335c591cb2674031d4ffd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>a092dc664cfefa1670fbb15c1065dfe41</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>ad37d8e2231cc9def37b1c8a928c82b2f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>a4c7ca441914d1710672a6a998d6af855</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>a48c1b9b3384f86505c09994e507cba3f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>a32f3aec98a668c6bb5fcd1fcc6106945</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1HessianEvaluation.html</anchorfile>
      <anchor>a17d670158d5880f713c2eb7ceccce001</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FunctionalProductDomainEvaluation</name>
    <filename>classRVL_1_1FunctionalProductDomainEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::FunctionalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>FunctionalProductDomainEvaluation</name>
      <anchorfile>classRVL_1_1FunctionalProductDomainEvaluation.html</anchorfile>
      <anchor>ae4c78ba11761ae046a9b7409d9835f3b</anchor>
      <arglist>(FunctionalProductDomain&lt; Scalar &gt; &amp;_f, const Vector&lt; Scalar &gt; &amp;_x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FunctionalProductDomainEvaluation</name>
      <anchorfile>classRVL_1_1FunctionalProductDomainEvaluation.html</anchorfile>
      <anchor>a94222b838ea6c628345559ac5608a472</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Vector&lt; Scalar &gt; &amp;</type>
      <name>getGradientBlock</name>
      <anchorfile>classRVL_1_1FunctionalProductDomainEvaluation.html</anchorfile>
      <anchor>a7ee8e9019e0c660be4ddf886e4f715e8</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>const LinearOp&lt; Scalar &gt; &amp;</type>
      <name>getHessianBlock</name>
      <anchorfile>classRVL_1_1FunctionalProductDomainEvaluation.html</anchorfile>
      <anchor>a22ee6d01109ee0bd91125d0ce856fc62</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1FunctionalProductDomainEvaluation.html</anchorfile>
      <anchor>a12d08c02185a7007bb6ba683693b5d6d</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::HessianBlockEvaluation</name>
    <filename>classRVL_1_1HessianBlockEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>~HessianBlockEvaluation</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a3cdaa91fb630b8127f6db1d0b1c5d28b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a152e88aaed606cf409cf22089bea175b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a01286e03976180cbb32a9fb6c99fd8b9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBlock</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a9843b140922b213f09c0348027d0d89e</anchor>
      <arglist>(int ii, int jj)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a1e2f5cc06ecd9ed972e88f276a635d10</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>HessianBlockEvaluation</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a121354dd26f5862f5fe3c375add613e6</anchor>
      <arglist>(FunctionalProductDomainEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>ac7749666d62c5e1df3bf82eb88b10b13</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>a730d7bcf82536d1f0b717ab6aef0da0e</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1HessianBlockEvaluation.html</anchorfile>
      <anchor>ad4c060ed6502ba0f32898b6eef7f6e85</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinCombFunctional</name>
    <filename>classRVL_1_1LinCombFunctional.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>LinCombFunctional</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>aacf57b50cc2738d40c0eeb731a6da4ac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombFunctional</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a6376e0cbc5e02981d05c1414c06c47d2</anchor>
      <arglist>(LinCombFunctional&lt; Scalar &gt; const &amp;fn)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombFunctional</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a887eee0ef379f91d3ead26fee587c790</anchor>
      <arglist>(Scalar a1, const Functional&lt; Scalar &gt; &amp;fn1, Scalar a2, const Functional&lt; Scalar &gt; &amp;fn2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LinCombFunctional</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a743db040b1ab31652642b21a0944654c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>ab67d72f871909fe2229d38a14193c847</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a5d240ff35214da94336b6bf561d603a5</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>aea953b642bf6a233f2d3dce1bb9b07f2</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a0dcc019ec869b397b934d80f507acace</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>a0dcdbfbb3f48f246e4260e97c7860cd3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinCombFunctional.html</anchorfile>
      <anchor>aef74e400f95362caf39ddb9be9a5df6e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::StdFOFunctional</name>
    <filename>classRVL_1_1StdFOFunctional.html</filename>
    <templarg>Scalar</templarg>
    <templarg>DataType</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>StdFOFunctional</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a4374eb5cc7861d80531a332caa745cdf</anchor>
      <arglist>(FunctionObjectScalarRedn&lt; Scalar &gt; &amp;f_, FunctionObject &amp;gradf_, FunctionObject &amp;hessf_, const Space&lt; Scalar &gt; &amp;dom_)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdFOFunctional</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a3e3ccd4eeb20a81f577437b1faf897d9</anchor>
      <arglist>(const StdFOFunctional&lt; Scalar, DataType &gt; &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StdFOFunctional</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a4254475123e3d0f8b742a5dc7a87a058</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>aca38993605d106047e4a3074511e7916</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a1c2bc64d5220243d7537a5b6d6349139</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a9fceeb31e2c94c86797efeab1a09f2d3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>ac98a41c13597d1e85b996b0a5a37324c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a9df1a28e39e21654f0c8cc2b10d53fe3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a7c9c2905e38a0195b8744fe02d6f6a50</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>StdFOFunctional</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>aafd060870b184c5c1cbe5f67f3f3b136</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>f</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a3d52f1b7d419c829898f6dc49fcb681c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObject &amp;</type>
      <name>gradf</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>ab7bb0f4117699f4b133463b384eab651</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionObject &amp;</type>
      <name>hessf</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>ad867e45363624cf5baf72b0f24e8d979</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>dom</name>
      <anchorfile>classRVL_1_1StdFOFunctional.html</anchorfile>
      <anchor>a8e130b4843dd33bb8d3d27e5efa5f9ba</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::NullFunctional</name>
    <filename>classRVL_1_1NullFunctional.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>NullFunctional</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>ac996c1a8539fec461bfd8edadeaaa038</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>NullFunctional</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a26c4e9d1aeb444f9f39a606a029fef30</anchor>
      <arglist>(const Functional&lt; Scalar &gt; &amp;_f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NullFunctional</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a3d73266238b0ad538767f95fd551a99a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a75530b56386acf2c81f3f24fdbe38087</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a3ac26518f44d32f2f31a836ac37c93c1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a64d4087e2e0dd29049fc8a26d59634cf</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a389f31b260609e415546879ae7293b69</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a8d54f256aa05abac045513721d2eae73</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1NullFunctional.html</anchorfile>
      <anchor>a4b71da1d1385ee9f01d713a749903e67</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::FcnlOpComp</name>
    <filename>classRVL_1_1FcnlOpComp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>FcnlOpComp</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a2250d8164a532e72cd26dcb177aa9e66</anchor>
      <arglist>(const Functional&lt; Scalar &gt; &amp;fref, const Operator&lt; Scalar &gt; &amp;opref)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>FcnlOpComp</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>ae60be64b370577a88d069c3f46e751f3</anchor>
      <arglist>(const FcnlOpComp&lt; Scalar &gt; &amp;c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~FcnlOpComp</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>adcf4437154a72c85e7f7aad29b109421</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a39443232288f4b453b6925c65d86aa4c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>aa4c816698e8d337d6968922402759148</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function">
      <type>OperatorEvaluation&lt; Scalar &gt; const &amp;</type>
      <name>getOpEval</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>aabe44f4fb641f980a7e9fbcd354acc41</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>FunctionalEvaluation&lt; Scalar &gt; const &amp;</type>
      <name>getFcnlEval</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a9424424156711c7ee52c6ac9d6d55b89</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a1fb1e6549e8eec70025407df85493aeb</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a4c7a67e1b4e8cfef52c269eac2641a97</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a0f2605ce23d70abacb36a70e657d9dca</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>a3e160076b34c24a2095cebae142c9f06</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1FcnlOpComp.html</anchorfile>
      <anchor>ac004b95cca4830185ae1480d6687fa7c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinCombObject</name>
    <filename>classRVL_1_1LinCombObject.html</filename>
    <templarg></templarg>
    <base>RVL::FunctionObject</base>
    <member kind="function">
      <type></type>
      <name>LinCombObject</name>
      <anchorfile>classRVL_1_1LinCombObject.html</anchorfile>
      <anchor>a27a21af25cbdfb02077b3c4554ef580f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombObject</name>
      <anchorfile>classRVL_1_1LinCombObject.html</anchorfile>
      <anchor>a147648c14c1ca42867bde5a1d7dbd587</anchor>
      <arglist>(const LinCombObject&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LinCombObject</name>
      <anchorfile>classRVL_1_1LinCombObject.html</anchorfile>
      <anchor>ac4d1a7647e3982587df4f2ce2963b9c9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>setScalar</name>
      <anchorfile>classRVL_1_1LinCombObject.html</anchorfile>
      <anchor>af2e51a9c9044fd235b979717bb8e0258</anchor>
      <arglist>(Scalar a, Scalar b)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearAlgebraPackage</name>
    <filename>classRVL_1_1LinearAlgebraPackage.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>LinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>a4ca4507a989433dd4b8c5cded6794cfb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>af0aee0bf1f48075b6e0d79ba51f64d1d</anchor>
      <arglist>(const LinearAlgebraPackage&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>ada6872bfad3620f27ba6fbff37da2397</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>a472736ef63b8af2df77f903e11a414bc</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual FunctionObject &amp;</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>a3ddd9b55c3cf404e20036a8cadc4aac2</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LinCombObject&lt; Scalar &gt; &amp;</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>a2823bbc6263c3d819dbca8843ead2cef</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1LinearAlgebraPackage.html</anchorfile>
      <anchor>a2aca4cdffa581aa37391ced7856fce75</anchor>
      <arglist>(LinearAlgebraPackage&lt; Scalar &gt; const &amp;lap) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearOp</name>
    <filename>classRVL_1_1LinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>LinearOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a2b8058841a73c28a666470ae6b218fcd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a9a4a09c7e0dc02570081f0f88bf61b22</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;Op)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LinearOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>ae21498ad3a8216e23366d36a778d8b47</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a396464766d28352f88d09bb48079709d</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyAdjOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a5a4410ad0d4c5fcd90ab97b80395d422</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>applyOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a60c5fbd3cebd4923ea05e48433c0f901</anchor>
      <arglist>(Scalar alpha, Vector&lt; Scalar &gt; const &amp;x, Scalar beta, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjOp</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>ae1a6ee7d752dc9ac5959ec15fc868877</anchor>
      <arglist>(Scalar alpha, Vector&lt; Scalar &gt; const &amp;x, Scalar beta, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void *</type>
      <name>operator new</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>ad414ce056d2bd75c89f4adf64a2316b1</anchor>
      <arglist>(size_t size)</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>af23a334ef408c110391414f90f23973e</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a16878b16fc8dbf5b310b926b5955feb0</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>a7526b619c1108c7168dee16a0d229a5e</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>aa06b0ff4cb2c1d06918d03171556df34</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1LinearOp.html</anchorfile>
      <anchor>ae03a4ad2efe9eb404bd5c7a3f52f09a0</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, Vector&lt; Scalar &gt; &amp;dx1) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearOpFO</name>
    <filename>classRVL_1_1LinearOpFO.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>LinearOpFO</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>af215b1ab5eb43bd876e00740486d8932</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_dom, const Space&lt; Scalar &gt; &amp;_rng, FunctionObject &amp;_fwdfo, FunctionObject &amp;_adjfo)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearOpFO</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>a9f136fafb69955bed0701c804e5e5220</anchor>
      <arglist>(const LinearOpFO&lt; Scalar &gt; &amp;l)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LinearOpFO</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>acd6b578a0f3b7904fb0b0d86a7fbd2f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>ab06384c5c75d8b90710cd492767f6558</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>a0095ad5d6b031586791cd29cc6629a2a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>af3ee0468e2ed3d4dac7fa6e89ffbf9d5</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>a5540ea4cd5dc591ff0899553a0f2865a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>a9c87f05cc3b2bba3216a56631e2dce8c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;Input, Vector&lt; Scalar &gt; &amp;Output) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1LinearOpFO.html</anchorfile>
      <anchor>aa0a2bc4bbbbbdfecef98c44b85031265</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;Input, Vector&lt; Scalar &gt; &amp;Output) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Invertible</name>
    <filename>classRVL_1_1Invertible.html</filename>
    <templarg>Scalar</templarg>
    <member kind="function">
      <type></type>
      <name>Invertible</name>
      <anchorfile>classRVL_1_1Invertible.html</anchorfile>
      <anchor>a03c92a5e8b751cf915904c91b73a6493</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Invertible</name>
      <anchorfile>classRVL_1_1Invertible.html</anchorfile>
      <anchor>ad37ca21e87f7d8a668277843149344c7</anchor>
      <arglist>(const Invertible&lt; Scalar &gt; &amp;Op)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Invertible</name>
      <anchorfile>classRVL_1_1Invertible.html</anchorfile>
      <anchor>a8cdce05920c9dd61598075e5e2ad5148</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyInv</name>
      <anchorfile>classRVL_1_1Invertible.html</anchorfile>
      <anchor>a8fb9fa04f48520e639c76c95ab1d7cba</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyInvAdj</name>
      <anchorfile>classRVL_1_1Invertible.html</anchorfile>
      <anchor>a1e83ff9859892dc409cfc78bbbafe775</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearOpWithInverse</name>
    <filename>classRVL_1_1LinearOpWithInverse.html</filename>
    <templarg></templarg>
    <base>RVL::LinearOp</base>
    <base>RVL::Invertible</base>
    <member kind="function">
      <type></type>
      <name>LinearOpWithInverse</name>
      <anchorfile>classRVL_1_1LinearOpWithInverse.html</anchorfile>
      <anchor>a9e65dffa6cdfea802ce31861239854b9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LinearOpWithInverse</name>
      <anchorfile>classRVL_1_1LinearOpWithInverse.html</anchorfile>
      <anchor>a21e6f41cab45c168cc06f267ac339b2c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyInvOp</name>
      <anchorfile>classRVL_1_1LinearOpWithInverse.html</anchorfile>
      <anchor>adb5e1dc744bb62645eed059b39660e44</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyInvAdjOp</name>
      <anchorfile>classRVL_1_1LinearOpWithInverse.html</anchorfile>
      <anchor>a7e968d60ae823907d1555e24f800ea02</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::AdjLinearOp</name>
    <filename>classRVL_1_1AdjLinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>AdjLinearOp</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>acdaee7cf6d57fd3c513ee55218ed3f2e</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;_op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AdjLinearOp</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>a69a8d90036a4e23e27edd5fd5f0ebcf2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>ae51707beea95d7cdb856dd1713c16d8b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>a0a15a6259b5e3814ebcc6c5e48e7c64d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>af69be7382b29b0d37ba8109bbdcbdd8b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>a475b9c156cf4f0105427be10d6ead3f4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>ae40d9359a1c5444769e4419bb7679334</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1AdjLinearOp.html</anchorfile>
      <anchor>ac7fca5b18663a626ec0f05e23642c2ec</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::NormalLinearOp</name>
    <filename>classRVL_1_1NormalLinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>NormalLinearOp</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>abeed282e662826af2b02c16f037254aa</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;_op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NormalLinearOp</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>add4118e506ea5374fd563ee3c3eabe5e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>a3c7d827efaa7b8d3b854eba561372020</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>a1744f705eef8f6d73c48e84a44aa29de</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>a44925812527a9af32cbff6c407071447</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>a0fd960fab6fd861c2f34a4a48f9cf877</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>ac142ec25e3be7698cfe573ded40ef6e3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1NormalLinearOp.html</anchorfile>
      <anchor>af55e22893c1a74f9848cd368a7299007</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScaleOpFwd</name>
    <filename>classRVL_1_1ScaleOpFwd.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>ScaleOpFwd</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>ad7d4fd7a1360fba545e2ab568a5e0d31</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_sp, Scalar _mu)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ScaleOpFwd</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a76c708f016db52d1db51ebf032ae9c81</anchor>
      <arglist>(const ScaleOpFwd&lt; Scalar &gt; &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ScaleOpFwd</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a21b72b638316f1a46002234834b06697</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>aa7391aef4b7b1703f2b1c79fe011e569</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a543ac3cc8843d86d81808a1430ae1a6b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getScale</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a0c80048924f64e9c81208d7ef20468e7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>ad749c93b99ad65fb3b7a4179507dcade</anchor>
      <arglist>(Scalar _mu)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a0d89b138da0f9205dd0426397baec851</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>aca95836a444b691b316f247386fd8748</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a60149a4e7cb1089187615db7c42b26fd</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1ScaleOpFwd.html</anchorfile>
      <anchor>a03f31b835d0d8a874fc193ab1435d512</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScaleOpInv</name>
    <filename>classRVL_1_1ScaleOpInv.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>ScaleOpInv</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>ae3d00f3f4f82a34edceba7f418dfd749</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;op, Scalar _mu)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ScaleOpInv</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a686537525909e08ec4496253aa40976e</anchor>
      <arglist>(const ScaleOpFwd&lt; Scalar &gt; &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ScaleOpInv</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a2e1f3692839e917e5bf9461fd945b637</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>ad086b7affca67f36124313f554216ada</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a3bb71d013d9f506a9ecfeb07f6f5cd80</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getScale</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a31f11d5097e20603d2bb21b28a38a837</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a4ba2eacedc7bca9b13aa3eb020979d1b</anchor>
      <arglist>(Scalar _mu)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a287e217e5aac6a7cf61b7ac23ef23eda</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>a1b68bf50f17bff55527b65c384e3dfab</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>af2e32867b68197df154e8cde6b499b4b</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1ScaleOpInv.html</anchorfile>
      <anchor>adc960326ce71598fef17d86ab5d33756</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinCombLinearOp</name>
    <filename>classRVL_1_1LinCombLinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>LinCombLinearOp</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a8c8a628aec9e1bd852f4a0ad57e2f683</anchor>
      <arglist>(Scalar _w1, LinearOp&lt; Scalar &gt; const &amp;_op1, Scalar _w2, LinearOp&lt; Scalar &gt; const &amp;_op2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombLinearOp</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a1889cc781fb70885675443ee540b2491</anchor>
      <arglist>(LinCombLinearOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LinCombLinearOp</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a806e2bb6237b5ed088324d6a463796ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>aee1fdd01a612a9cf5373cea0231ba5e5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a58edfe8c8ef02a93c97e13b3b8960858</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a5c86e780796b474533d319ae13884431</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a151f5fb4157f05f22241161d954a0ea1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>ac86668f1617ff93968a7c699badb675a</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1LinCombLinearOp.html</anchorfile>
      <anchor>a5980b5a1b7cadf3fdb55dd69925413fe</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::CompLinearOp</name>
    <filename>classRVL_1_1CompLinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>CompLinearOp</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>a65d18aafaadda4d3b39c89d4d91a4918</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;_op1, LinearOp&lt; Scalar &gt; const &amp;_op2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CompLinearOp</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>af43172244c7e3d5357338c9b6f7796eb</anchor>
      <arglist>(CompLinearOp const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~CompLinearOp</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>add47f805a8102c78fb1822d593a77a66</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>a54613e586c2ae1fbb4c7c6fe3310b99e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>ae93e9dbca4297b0d0702f313a84a1932</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>a54f7e122679fd05c39166362dd55ae69</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>aa25d6fce8036d38d901e24f7d3049dae</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>a8dc9a9bf3858629f071d11b21cfdbd57</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1CompLinearOp.html</anchorfile>
      <anchor>adb366c31b60043b607c3768d3e38919a</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::SymmetricBilinearOp</name>
    <filename>classRVL_1_1SymmetricBilinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type>void *</type>
      <name>operator new</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>a33f659c5a22bef075f44e501ded36fee</anchor>
      <arglist>(size_t size)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SymmetricBilinearOp</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>a074429e6a9e8a82eaf57f3c606615c7f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>ab686773c186c5f939bb4cb4e6ecdb5b3</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Space&lt; Scalar &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>a2fa838dd63860898079a4cf83d045dfe</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyOp</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>ab511093023dd4a25d82a33b31ef40bbd</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x0, Vector&lt; Scalar &gt; const &amp;x1, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>applyAdjOp</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>a73a831731331797e2aef22cc415740fe</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x0, Vector&lt; Scalar &gt; const &amp;y, Vector&lt; Scalar &gt; &amp;x1) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>a6efde3c1ef02652c31903a3a326a4a9c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x0, const Vector&lt; Scalar &gt; &amp;x1, Vector&lt; Scalar &gt; &amp;y) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>ab55612ba178c2b5e41774a670df3179a</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x0, const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;x1) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual SymmetricBilinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1SymmetricBilinearOp.html</anchorfile>
      <anchor>aa04bc41d59160d1f95baac7678e61938</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearBilinearOp</name>
    <filename>classRVL_1_1LinearBilinearOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>LinearBilinearOp</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>a72258c12204a797c8e6d0a82aacede8b</anchor>
      <arglist>(SymmetricBilinearOp&lt; Scalar &gt; const &amp;_blop, Vector&lt; Scalar &gt; const &amp;_x0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinearBilinearOp</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>a227a512aaea26ba4cd524e710f73215f</anchor>
      <arglist>(LinearBilinearOp&lt; Scalar &gt; const &amp;lbl)</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>a02a94777442a14f3ccb759c049fff133</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>ad1fe0e994e98a88eb78c8fd50e789a33</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>a7db49d209c7e7461cc33ec99e49d4056</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>a6c56925e84184a92107a6b4440228d4f</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x1, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>ac0a8d4e195c94a0fe423644b15c03b52</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;y, Vector&lt; Scalar &gt; &amp;x1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinearBilinearOp.html</anchorfile>
      <anchor>ae37eaa7a4c5463cddcf9678ef048e2bf</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ShiftOperator</name>
    <filename>classRVL_1_1ShiftOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>ShiftOperator</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>ab8324ff4591b118b0f523b5a8cd73815</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;dd)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ShiftOperator</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>ac6536dfb8c53b0ace7ff29140cd33a69</anchor>
      <arglist>(const ShiftOperator&lt; Scalar &gt; &amp;a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ShiftOperator</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>ab8672bbdfe1a6d0598d19122ab755882</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a45db13474381210e81f25ac06fdb01b6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a08166a61cf85c9f55233a02585e01d9e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a0873799601d48aeb103ad991c5400c0d</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>ac43085b15bc7404d4505f5293b0d1307</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a981955378028760f951dc18a3a86aace</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a6eb4d8ea314c00d8f88af978ad042168</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ShiftOperator.html</anchorfile>
      <anchor>a5dbd2827f4c6bf70f6ccbf6b9fb201ae</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ResidualOperator</name>
    <filename>classRVL_1_1ResidualOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>ResidualOperator</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a57f9aebc1c69145ad779929d5dc017b6</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;GG, Vector&lt; Scalar &gt; const &amp;dd)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ResidualOperator</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a97c1821e1e49d9177b819f4bc485ff3c</anchor>
      <arglist>(ResidualOperator&lt; Scalar &gt; const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ResidualOperator</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a341cb79642cbe6f093694dd4c2178ce4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>aea74f1cf36645eddb1b7b9ab4940ecfc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>addc31e13f5c68c30a7094ef2c214f683</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a913fa210608a638fc46c45cce871cd55</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a5c05c4cabcbbf39402af0d1d31967ab7</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a847b456dff7249bca65aec7d6efa9b54</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>ac047b9d317ebb6990b39b42412cf0f1a</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>a77a02c4107ba2cbc28038ae633113480</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ResidualOperator.html</anchorfile>
      <anchor>af7a72c38f913b28d66ad5c2201b004f9</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::EuclideanForm</name>
    <filename>classRVL_1_1EuclideanForm.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>EuclideanForm</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>a0eddafa554d05ff219764ae1c1f4305a</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>EuclideanForm</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>a910406ad8fb62a251dd57c0a73ce86f7</anchor>
      <arglist>(const EuclideanForm&lt; Scalar &gt; &amp;q)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~EuclideanForm</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>aa79d4e6bea4b6c1266db9cf0483e8f4f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>aa8f2b2fa864d547ba7c255b877651f4e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>a3d4bdb318e4857b50b2e5f5d1bb235fb</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>aaa81fb8e73c3ad857cf2fbd90466a733</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>ae0c0649317b69d7976bc20bfc6539a9d</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>a4bbbcfb846b2212fbfbf833deb0046b7</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;delx, Vector&lt; Scalar &gt; &amp;dely) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1EuclideanForm.html</anchorfile>
      <anchor>a4da89430a93f5152b72000459db9a255</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuadraticForm</name>
    <filename>classRVL_1_1QuadraticForm.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>QuadraticForm</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a96509257e114477e32deddf1bb848410</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;AA)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>QuadraticForm</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>ac427e65d0ec3a56051d1a34fb0bfdd2b</anchor>
      <arglist>(const QuadraticForm&lt; Scalar &gt; &amp;q)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~QuadraticForm</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a0273d5df98dc0a3bfb80ea90d1121201</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a690b038fefbbf05faa8865de45228982</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>aa4bfba02fe24ab77be099c582c1d8d2e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a9b3f8ddec6676f629244fce462216e85</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a48a2474834104c5fb7ecb4f6b367ae73</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>a8dcf94d8a4f29dd08eb264f2cbdd1269</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;delx, Vector&lt; Scalar &gt; &amp;dely) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1QuadraticForm.html</anchorfile>
      <anchor>aa8596d8a281c82c31153c730c0eb0a7a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ShiftedQuadraticForm</name>
    <filename>classRVL_1_1ShiftedQuadraticForm.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>ShiftedQuadraticForm</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a3262189d43192bd2e6e3271fdf0179f6</anchor>
      <arglist>(const LinearOp&lt; Scalar &gt; &amp;AA, const Vector&lt; Scalar &gt; &amp;bb)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ShiftedQuadraticForm</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>ac33d00d814f42c9ded13f19ce658bb1d</anchor>
      <arglist>(const ShiftedQuadraticForm&lt; Scalar &gt; &amp;q)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ShiftedQuadraticForm</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a53d98b07367d51de642227ba80d105ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a856dc0614f795eec625931ebfd0df4a7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a7ff05f22ee9677946adb8d619e898e04</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a4ef8c383381107a6c2ba9a2db10f6678</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a9f8e4277c8dd30879811e1d0c97d19fb</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>acde3866b57d877c38d896412333f03ca</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;delx, Vector&lt; Scalar &gt; &amp;dely) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ShiftedQuadraticForm.html</anchorfile>
      <anchor>a628b5adb0fa3a03bd7a45697bf42ff38</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LeastSquaresFcnlGN</name>
    <filename>classRVL_1_1LeastSquaresFcnlGN.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>LeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a84af2ef01d2d139b7ed1531c7a9bcd61</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a79e238171c7392232af15e7961585daf</anchor>
      <arglist>(const LeastSquaresFcnlGN&lt; Scalar &gt; &amp;J)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a91df0255f796484aaa277f48bc9d280f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>af96ac46b1d6ee1b22cda68b5a26e9b2b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a4ff8d178a391725ab6589eb719a755e8</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a98d5531f2ee82d2e12b010035b61f14b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a612105d3cb0ff2209b34aab42a6aed5c</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a1f412cd31a71700d71c57139e58526c8</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a1c386187482dfd96815c224d8a8c34f9</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; const &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LeastSquaresFcnlGN.html</anchorfile>
      <anchor>a9b590b5e319937bba8a59b4dd84ad376</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::StdLeastSquaresFcnlGN</name>
    <filename>classRVL_1_1StdLeastSquaresFcnlGN.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Functional</base>
    <member kind="function">
      <type></type>
      <name>StdLeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a9e2ef36f26f31b34ed5b69af33d14ae5</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;oper, Vector&lt; Scalar &gt; const &amp;d)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdLeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a6d21421bed909e40f4550b9d2edd0cc7</anchor>
      <arglist>(const StdLeastSquaresFcnlGN&lt; Scalar &gt; &amp;JJ)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~StdLeastSquaresFcnlGN</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a155724c79e18296112e801dc43fa5313</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a13a7f65b3e6d6ab41676a42dad08e2fd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>aa7e7a22ba1b10daceddafe6735eb8fa6</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a5b9528d01e164d29cf7b37120d3fb770</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a994f5e21df700d3cc493bee5a7568f4a</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Scalar &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyGradient</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a08ea6ac07e0a7143f6611a16b83db711</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; &amp;g) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyHessian</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a2f4bea6cb1bc44563a8118ed09e128f6</anchor>
      <arglist>(Vector&lt; Scalar &gt; const &amp;x, Vector&lt; Scalar &gt; const &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Functional&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1StdLeastSquaresFcnlGN.html</anchorfile>
      <anchor>a18f3a656cb8e7e4523cb3cacaf2a77f6</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Operator</name>
    <filename>classRVL_1_1Operator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Operator</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a0b3f633e87d980ee67dca067fe0b81b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Operator</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a02d9dcd1c4c3943ee8d5ca32c0c738c4</anchor>
      <arglist>(const Operator&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Operator</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>ae43230421922d5b9e2de7e8c5b09f720</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a890f5616c0e3e75622fcc8196da10ea1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a95473c08951af07ed6e22b3fac1594c8</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>getMaxStep</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a4e881086245bd8a90d5b4bd3bd11f7ee</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>ac34f0a59d15e187324c054791c77dba4</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a9f368828159c2c1aaa04cdb778227da6</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>aedc7d3bf09770bdbe19f8d002148fcf6</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv2</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a58b81e63417f9b8f9b43e2492108f060</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, Vector&lt; Scalar &gt; &amp;) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a754c7531e917d2b90e55d9279bf1e934</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, const Vector&lt; Scalar &gt; &amp;, Vector&lt; Scalar &gt; &amp;) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_apply</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a4892a83d6e23f42b9bc36bdaaf587946</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyDeriv</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>aac08a95a81bfb373df77de07989f1df4</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyAdjDeriv</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a76d703da5d5c1aca4f0c2841811b3ab1</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyDeriv2</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a2157463b20e85e97bf0f96b419ee137c</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dx1, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>export_applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a204ecbac383afef9c10b36c662a93d37</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>abfe2e784e45b79b152f83e8158563814</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>DerivEvaluation&lt; Scalar &gt; *</type>
      <name>createDerivEvaluation</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a282b929923d30222fa89fa373e46ee4c</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;opeval) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Deriv2Evaluation&lt; Scalar &gt; *</type>
      <name>createDeriv2Evaluation</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>aae7f0005ace42edb1ddd7ef78bafab28</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;opeval) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FcnlOpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>ad6486a678792a5664c5c427696939eed</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a7f134df83affca21fd956dd946edd476</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>LinCombOperator&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Operator.html</anchorfile>
      <anchor>a240bccef36a7d07e16bb73086fa6ca06</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OperatorProductDomain</name>
    <filename>classRVL_1_1OperatorProductDomain.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>OperatorProductDomain</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>af20eadf09b455a7c3c163e4ebb72b5d3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OperatorProductDomain</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a2fc1b1a873b7b50151430dccad1965c1</anchor>
      <arglist>(const OperatorProductDomain&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~OperatorProductDomain</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a8ceea65f63268ae4f0a0dec58e102d74</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductDomain</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a3b69102de2e8af8eeac9a9aa3a56162a</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>af0019409f15c5aaeb3f69c6b732c83eb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyPartialDeriv</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a9808ae276fff75a3cf9337f8879f130d</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dxi, Vector&lt; Scalar &gt; &amp;dy) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a90ad9561096a9e73819bedf92f9a4e74</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyAdjPartialDeriv</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a047e1be5f217984694cd6e8017b0c272</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dxi) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>aaf102e81a33871ca2609517efd103093</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual OperatorProductDomain&lt; Scalar &gt; *</type>
      <name>clonePD</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>ad922248173f9557581259a52cef28970</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>ae3eab1fd4a290c6d97b30101d9f72b98</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorProductDomain.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OperatorWithInvertibleDeriv</name>
    <filename>classRVL_1_1OperatorWithInvertibleDeriv.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>OperatorWithInvertibleDeriv</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>ae4a950d526464531e21da4d0ad9daaf6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~OperatorWithInvertibleDeriv</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>a5b7e5cd337b74d79e887594c08ccde82</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyInverseDeriv</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>a5e3b4ea0e7581cc6dad9e0acb950a7f8</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>applyAdjInverseDeriv</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>a3f5adba1a56467738d93065e971a7fa2</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>DerivEvaluation&lt; Scalar &gt; *</type>
      <name>createDerivEvaluation</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>abbf1ce1f2d76ba9137bc006b96e8f868</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;opeval) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>InvertibleDerivEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorWithInvertibleDeriv.html</anchorfile>
      <anchor>a0d2ae11e58f3849b5797e6c059085bbe</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OperatorEvaluation</name>
    <filename>classRVL_1_1OperatorEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>OperatorEvaluation</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>ad3fdb53ddd2a3def778b5389be63a534</anchor>
      <arglist>(const Operator&lt; Scalar &gt; &amp;_f, const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OperatorEvaluation</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a948f174e723edbe124bbd75fab18f297</anchor>
      <arglist>(const OperatorEvaluation&lt; Scalar &gt; &amp;ev)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~OperatorEvaluation</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>abcdd74c6db1d077b23a68c7f9d58f11f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a859cd3031c69045da8380d5d3837941d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>ae6a44501822f24b69a7679cb8634ab83</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>getPoint</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>aeb290adab7199f57d17fe2fda885c0c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getPoint</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a0473da05c74ae1f4837403297cb4a654</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Operator&lt; Scalar &gt; const &amp;</type>
      <name>getOp</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a946f38465874026f2b6da5078a4d879b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getValue</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>aa0db386d0c9729ed04360e7cf41cd2a3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LinearOp&lt; Scalar &gt; const &amp;</type>
      <name>getDeriv</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>abbf80703c2f39953bdacc287eeda4e6d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>SymmetricBilinearOp&lt; Scalar &gt; const &amp;</type>
      <name>getDeriv2</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>aba7600e4a57e275b1c7b0043154841b5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>af2f1f4591e08d41f101c82956514dfc3</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a973d23c2905d6a42b3d60898955053bc</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a57ce439300a791b9f18c3e76e50671b4</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>aa3ded79f3987634df5f2d4e160cf75e9</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x0, const Vector&lt; Scalar &gt; &amp;x1, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a9d2d4b987f09e8dc2cb41c1c8612b008</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x0, const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;x1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>const ProductSpace&lt; Scalar &gt; &amp;</type>
      <name>getProductDomain</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a98a9c91d30449c81fc16d00b7609f5d7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyPartialDeriv</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>af55bf8981f10c6a915b4cdd6c6ada4ab</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjPartialDeriv</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a70e35e1f6bff4742c3d86c2e6b3af251</anchor>
      <arglist>(int i, const Vector&lt; Scalar &gt; &amp;yin, Vector&lt; Scalar &gt; &amp;yout) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a7f134df83affca21fd956dd946edd476</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>FcnlOpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>ad6486a678792a5664c5c427696939eed</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>LinCombOperator&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a240bccef36a7d07e16bb73086fa6ca06</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>DerivEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a6d8f3d2d0b2d47ce1daa9a232730abf4</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>Deriv2Evaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a0ecf57ed0f4e5052a3e68df57a0f4077</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>PartialDerivEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorEvaluation.html</anchorfile>
      <anchor>a47ceb3dcc493b46f4ac2919be0516f87</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::DerivEvaluation</name>
    <filename>classRVL_1_1DerivEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>~DerivEvaluation</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a5b3930df9d9c117a317141d36405cbdd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>acc31cd7faa62f23e6a34f6b4250a6304</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a873db908443a135fcccf4c46b17f0a60</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>acd724e7e97d6add72e63d5faa223209e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>DerivEvaluation</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a0d02965c4612ddf28a8cf863562616ea</anchor>
      <arglist>(const DerivEvaluation&lt; Scalar &gt; &amp;d)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>DerivEvaluation</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a386fef1f8833599dc705a3e18e19493e</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>ae3aa176711a7b1e095ffb3d9a41de259</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; &amp;</type>
      <name>getOp</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a92f67922bd0796a6a85f03cb6c410e7b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>ab755d93c0076038d531543117bf766ad</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a9a21c83b50484a7efdef4004ac1dc594</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>OperatorEvaluation&lt; Scalar &gt; &amp;</type>
      <name>fx</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>af5d7055eccaf7369d292048d3d44446d</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>Operator&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1DerivEvaluation.html</anchorfile>
      <anchor>a6d0f7023ea1310002fd044e336e0dc5d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Deriv2Evaluation</name>
    <filename>classRVL_1_1Deriv2Evaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::SymmetricBilinearOp</base>
    <member kind="function">
      <type></type>
      <name>~Deriv2Evaluation</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>af46af1bd11750204acf1e1c6ee601de7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a9a6670acbe51715f521e68b59e168395</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a45fb07895ac99e5dfb29c0f143ff75f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>ad6cd024ff3b5a7c15e6c382d956a4fbb</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>Deriv2Evaluation</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a86d5ed50f60d27a35df6f1d0c48e61b7</anchor>
      <arglist>(const Deriv2Evaluation&lt; Scalar &gt; &amp;d)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>Deriv2Evaluation</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a5126760a1a1e47b0daa797ba454eaab2</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>SymmetricBilinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a676b90d3efceb315ee5d17944f5a3da2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; &amp;</type>
      <name>getOp</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>ae1217fee181a326dda3526e33124706d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a01b67b774f508e82242a0dac07dbf2e2</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y0, const Vector&lt; Scalar &gt; &amp;y1, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a895f364042b28e420f7098af077f63c3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y0, const Vector&lt; Scalar &gt; &amp;z, Vector&lt; Scalar &gt; &amp;y1) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>OperatorEvaluation&lt; Scalar &gt; &amp;</type>
      <name>fx</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a298e75b8dd0a6c5ab52d31fe13754a1c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a12541e06e6dac82bc05145a0fa99a64c</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>Operator&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Deriv2Evaluation.html</anchorfile>
      <anchor>a6d0f7023ea1310002fd044e336e0dc5d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::InvertibleDerivEvaluation</name>
    <filename>classRVL_1_1InvertibleDerivEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::DerivEvaluation</base>
    <base>RVL::Invertible</base>
    <member kind="function">
      <type></type>
      <name>InvertibleDerivEvaluation</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a564e991cc782a11977c95d43bd8d4e1d</anchor>
      <arglist>(OperatorEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>adafa0487a8ebe3c9e670af0ae8a67403</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>aa59a3800a491fa2a8cdfc271550ce275</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a46526936c0abaab3573ccb7d198bc929</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a11484e695ad49344282f64a079fc05f9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>InvertibleDerivEvaluation</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a24e9fbdefe03bc589984bef3a5afea4c</anchor>
      <arglist>(const InvertibleDerivEvaluation&lt; Scalar &gt; &amp;s)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyInv</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a0beffe75ee4bbf1a70b5f07625254548</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyInvAdj</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a980a69bdc0aa989dc04fb897f541589c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorWithInvertibleDeriv&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1InvertibleDerivEvaluation.html</anchorfile>
      <anchor>a96b18dad618dfcf5d84588f4ff6e90b9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OperatorProductDomainEvaluation</name>
    <filename>classRVL_1_1OperatorProductDomainEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::OperatorEvaluation</base>
    <member kind="function">
      <type></type>
      <name>OperatorProductDomainEvaluation</name>
      <anchorfile>classRVL_1_1OperatorProductDomainEvaluation.html</anchorfile>
      <anchor>aac13300604e2e7e7f9941ddc07b88a6d</anchor>
      <arglist>(OperatorProductDomain&lt; Scalar &gt; &amp;_f, const Vector&lt; Scalar &gt; &amp;_x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~OperatorProductDomainEvaluation</name>
      <anchorfile>classRVL_1_1OperatorProductDomainEvaluation.html</anchorfile>
      <anchor>a3704dc2b753c307e9fe4d5f320d479e0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const LinearOp&lt; Scalar &gt; &amp;</type>
      <name>getPartialDeriv</name>
      <anchorfile>classRVL_1_1OperatorProductDomainEvaluation.html</anchorfile>
      <anchor>a92cce3c20a0759adda0c2a54e7ee0532</anchor>
      <arglist>(int ic)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1OperatorProductDomainEvaluation.html</anchorfile>
      <anchor>a80727a9cd51dc4b4d2646b6e540867b5</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>PartialDerivEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1OperatorProductDomainEvaluation.html</anchorfile>
      <anchor>a47ceb3dcc493b46f4ac2919be0516f87</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PartialDerivEvaluation</name>
    <filename>classRVL_1_1PartialDerivEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::LinearOp</base>
    <member kind="function">
      <type></type>
      <name>~PartialDerivEvaluation</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>a2628dbc37a793e8a58dd8d003d806b7c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>a513558c7d32a05bcd2f93045768f9bb3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>a5b949886f11fe1821a59de99237f5abb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setBlock</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>a6dae9a503364d65d2973387c347fba1b</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>ac0392892bfe633068f41f6d24f669b1b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>PartialDerivEvaluation</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>a173080f1b8de4faeca6c1b04c961b5d4</anchor>
      <arglist>(OperatorProductDomainEvaluation&lt; Scalar &gt; &amp;_fx)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>aac67685e9e4891547de63882c5f73696</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>ac2ee444daeaa8c8d4646331c71c27f59</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>ae33435cab4726eb0c3ce7899c00b1aae</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OperatorProductDomainEvaluation&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1PartialDerivEvaluation.html</anchorfile>
      <anchor>acf433797cf913ad6a8ccff2c199e5f17</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LNLOperator</name>
    <filename>classRVL_1_1LNLOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>LNLOperator</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>ad0a51f6cb71c6e09d925c642b9cdb35c</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;LL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LNLOperator</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>ade3fb5fe7882fbe8ceaaaf9a1c05afc1</anchor>
      <arglist>(const LNLOperator&lt; Scalar &gt; &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LNLOperator</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>a38da963dff7c5fe03ca6a3ab38e6c8c2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>ada971f3e91e56f8e7f727982b994f6d9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>aeef77ede6eebb0df448cda97fa138e0f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>a998961c3361ac27762aaa10e438474a0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>a3404391fbc3408d654bf6c9d79c471bd</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>ad2bf4870de4ac5792c08999887615037</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>a4838c3707bc50eb8469789ceef6b30d3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1LNLOperator.html</anchorfile>
      <anchor>a46ceee87faddc50b7f02d59e929c324f</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ANLOperator</name>
    <filename>classRVL_1_1ANLOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>ANLOperator</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a0c5d712b30b3d51347b7c056caabdcf4</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;LL, Vector&lt; Scalar &gt; const &amp;dd)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ANLOperator</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a74a37b416265757099b83b7fcce87f8f</anchor>
      <arglist>(const ANLOperator&lt; Scalar &gt; &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ANLOperator</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a874e06bf14acc2a5bdf42f53715b73fa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a33a9c979c44464f531634105a1c9c8b3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a7327c13bfb25e456ad391a5fc41bb3da</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>acca5c765016a0dfd6fa6d4793fff067f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>aa42cde2c5c29c3d3a6d283629441d419</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>a3fb74a17a3662f3858d2a618f25d9a03</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>abc6472b40d317677f6a93d0a0f2198e9</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1ANLOperator.html</anchorfile>
      <anchor>adc4fa351d49af0a300ef56763772d371</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OpFO</name>
    <filename>classRVL_1_1OpFO.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>OpFO</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>acc4ba220c0bb572c31f9daed0975151f</anchor>
      <arglist>(const OpFO&lt; Scalar &gt; &amp;A)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpFO</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>ad28283cc1c296d35499972e951eaad29</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;_dom, Space&lt; Scalar &gt; const &amp;_rng, FunctionObject &amp;_f, FunctionObject &amp;_dff, FunctionObject &amp;_dfa)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpFO</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a3f4c63d82057a4a91781fa253deffa53</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;_dom, Space&lt; Scalar &gt; const &amp;_rng, FunctionObject &amp;_f, FunctionObject &amp;_dff, FunctionObject &amp;_dfa, std::vector&lt; RVL::Vector&lt; Scalar &gt; const * &gt; par)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~OpFO</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a69d8d5af31e5774181e238438bf58876</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a35171272d1ddf7a6d20e2661624fc26c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a25d95375c9695ddc672a4c7be91e7fef</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a720db3d0d08a5f41e04ec5ac7dd5b50e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a10dfcc0994fc03fa9ef698c1ae3ddbad</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>ad0c1f33df4fa60004595ecd81f7374d2</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>a2dd1aae75d2c0db890577bef6b9d6ee8</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1OpFO.html</anchorfile>
      <anchor>af2664a13c536daea29179e2a90511120</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;y, Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinCombOperator</name>
    <filename>classRVL_1_1LinCombOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a91c571ced80b5ab3c4d935002807a8a2</anchor>
      <arglist>(LinCombOperator&lt; Scalar &gt; const &amp;op)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a1661cd9c0368552b02c89f5ed25d4332</anchor>
      <arglist>(Scalar a1, Operator&lt; Scalar &gt; const &amp;op1, Scalar a2, Operator&lt; Scalar &gt; const &amp;op2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>ad66512050d5a4265b614011e561cca73</anchor>
      <arglist>(Scalar a1, Operator&lt; Scalar &gt; const &amp;op1, Scalar a2, Operator&lt; Scalar &gt; const &amp;op2, Scalar a3, Operator&lt; Scalar &gt; const &amp;op3)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>af1b380e29d4650cb0025c256ecf56096</anchor>
      <arglist>(Scalar a1, Operator&lt; Scalar &gt; const &amp;op1, Scalar a2, Operator&lt; Scalar &gt; const &amp;op2, Scalar a3, Operator&lt; Scalar &gt; const &amp;op3, Scalar a4, Operator&lt; Scalar &gt; const &amp;op4)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a7d8259a5623fbef2a53a217db1b0afb8</anchor>
      <arglist>(Scalar a1, Operator&lt; Scalar &gt; const &amp;op1, Scalar a2, Operator&lt; Scalar &gt; const &amp;op2, Scalar a3, Operator&lt; Scalar &gt; const &amp;op3, Scalar a4, Operator&lt; Scalar &gt; const &amp;op4, Scalar a5, Operator&lt; Scalar &gt; const &amp;op5)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LinCombOperator</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a4b5a4b49f0503dffe6273f1056396214</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>aaac847cd13634e7614fd1be9ef52b65d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>aa420e9c83e50c6d8ed5d6cc07107d6f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>abc6bc7f3d64d9e0c271f60cd18e5bc30</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a7b1e9525b0e7faafe8f67df2edf2de88</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>ad7925ce62240ec9d7a27de3968a28681</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>aed89889887e5616ca20d28eec39e4fb3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1LinCombOperator.html</anchorfile>
      <anchor>a99eb4fa802aea09cfab9ea0269030d30</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearOpEvaluation</name>
    <filename>classRVL_1_1LinearOpEvaluation.html</filename>
    <templarg>Scalar</templarg>
    <member kind="function">
      <type></type>
      <name>LinearOpEvaluation</name>
      <anchorfile>classRVL_1_1LinearOpEvaluation.html</anchorfile>
      <anchor>a8f5ccc076cdc0a719d3ebdec7740c138</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;_op, Vector&lt; Scalar &gt; const &amp;_ref)</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getValue</name>
      <anchorfile>classRVL_1_1LinearOpEvaluation.html</anchorfile>
      <anchor>ae61fe61d5f0a2e4434ea0836d137dd3f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1LinearOpEvaluation.html</anchorfile>
      <anchor>a7f134df83affca21fd956dd946edd476</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LinearOpAdjEvaluation</name>
    <filename>classRVL_1_1LinearOpAdjEvaluation.html</filename>
    <templarg></templarg>
    <member kind="function">
      <type></type>
      <name>LinearOpAdjEvaluation</name>
      <anchorfile>classRVL_1_1LinearOpAdjEvaluation.html</anchorfile>
      <anchor>a194c3be8a99b2fd454084f50e723431d</anchor>
      <arglist>(LinearOp&lt; Scalar &gt; const &amp;_op, Vector&lt; Scalar &gt; const &amp;_ref)</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>getValue</name>
      <anchorfile>classRVL_1_1LinearOpAdjEvaluation.html</anchorfile>
      <anchor>a9d3b44a81e6f40597d4b5c1b8b8c9fed</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>OpComp&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1LinearOpAdjEvaluation.html</anchorfile>
      <anchor>a7f134df83affca21fd956dd946edd476</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OpComp</name>
    <filename>classRVL_1_1OpComp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a6f4ac07a4b2be33c7be0d002e4f78271</anchor>
      <arglist>(OpComp&lt; Scalar &gt; const &amp;oc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>afffe49d8f761cabac1a4409c3e015e37</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op1ref, Operator&lt; Scalar &gt; const &amp;op2ref)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a6b1965f7a949a0121d6a8d78c00c8bdb</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op1ref, Operator&lt; Scalar &gt; const &amp;op2ref, Operator&lt; Scalar &gt; const &amp;op3ref)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a1d12fcc409f56dec7c82101b13e45afd</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op1ref, Operator&lt; Scalar &gt; const &amp;op2ref, Operator&lt; Scalar &gt; const &amp;op3ref, Operator&lt; Scalar &gt; const &amp;op4ref)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>ab0cd6adbcb339c3389649feb6d0969e2</anchor>
      <arglist>(Operator&lt; Scalar &gt; const &amp;op1ref, Operator&lt; Scalar &gt; const &amp;op2ref, Operator&lt; Scalar &gt; const &amp;op3ref, Operator&lt; Scalar &gt; const &amp;op4ref, Operator&lt; Scalar &gt; const &amp;op5ref)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~OpComp</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a6e11a76bf90eb14a7130b27a1aade48c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a102866dc19108b264e36f944d0aa07af</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>aa548d72e238d353f58e0294dab50f1bd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>ad3c30f056f615d4c63ab2fadaa89da0f</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a4524b4b66a5527617838130cc97479ef</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a355bea986379aabac27e5eb175bf0677</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a5942bbf6001ddfa62496a56fe938630e</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv2</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a8cbb243e47c7ba75a1d35bcad43906da</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dx1, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv2</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a8556e455a3020983ee81f6232a86a1c8</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx0, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx1) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1OpComp.html</anchorfile>
      <anchor>a1a1ae91d51b2d04f6140d9c2f6729337</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::IdentityOp</name>
    <filename>classRVL_1_1IdentityOp.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Operator</base>
    <member kind="function">
      <type></type>
      <name>IdentityOp</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>aa493ce7902b912c59476e5df7ee70b68</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IdentityOp</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>a7f5ad8d22f1ead6b40c53a18e4933bf0</anchor>
      <arglist>(const IdentityOp&lt; Scalar &gt; &amp;c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IdentityOp</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>ac99d3f8106d77cdbfa6ab47c1ce39ff7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>af33811cfb0234fe78f59172bd2d02ad8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>a3bcdc0b1cf799444c6e90500baad7d0f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>a39d9670ee8808b098500f7f4e931b9b7</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>a84024f450174caf00912a2dd930fb490</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;val) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>a9103cf6074b44bb50d45738d99c0a6b7</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>ac590cfe14759e4056eb88439676f492a</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1IdentityOp.html</anchorfile>
      <anchor>add3fd4f4767c70b1e803bd7962a1233c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Product</name>
    <filename>classRVL_1_1Product.html</filename>
    <templarg>T</templarg>
    <member kind="function">
      <type></type>
      <name>Product</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>ac4b13e433a79c2ba283c22d2ea673435</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Product</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>a44acd1f2cc14a99082e1a6f04afe4680</anchor>
      <arglist>(const Product&lt; T &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Product</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>ac01b5c03c2cb64311680503d9c5ca6dd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>a56f322601ca59b5672aef0f9610c8e7e</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual T &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>ab67d7988e2adc784cc47d8cd0033106d</anchor>
      <arglist>(size_t i)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual T const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1Product.html</anchorfile>
      <anchor>af75e8b3bfa8c4d29cdcd00b07dcc7add</anchor>
      <arglist>(size_t i) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ROProduct</name>
    <filename>classRVL_1_1ROProduct.html</filename>
    <templarg>T</templarg>
    <member kind="function">
      <type></type>
      <name>ROProduct</name>
      <anchorfile>classRVL_1_1ROProduct.html</anchorfile>
      <anchor>ac4a62c96ec4f73348711c5f12332265d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ROProduct</name>
      <anchorfile>classRVL_1_1ROProduct.html</anchorfile>
      <anchor>a053bbff854054fda6b11896c80fd6cee</anchor>
      <arglist>(const ROProduct&lt; T &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ROProduct</name>
      <anchorfile>classRVL_1_1ROProduct.html</anchorfile>
      <anchor>a558c756ed2154c433184aafe7e192769</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1ROProduct.html</anchorfile>
      <anchor>a4e4388b561867fdea258f3667496bfcd</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual T const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1ROProduct.html</anchorfile>
      <anchor>a03f031ec9f02324d673e26084e5555bf</anchor>
      <arglist>(size_t i) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BlockFunctionObject</name>
    <filename>classRVL_1_1BlockFunctionObject.html</filename>
    <base>RVL::FunctionObject</base>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1BlockFunctionObject.html</anchorfile>
      <anchor>aca3289db5f2ed34bb808c9bddde22327</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::DiagonalFunctionObject</name>
    <filename>classRVL_1_1DiagonalFunctionObject.html</filename>
    <base>RVL::BlockFunctionObject</base>
    <member kind="function">
      <type></type>
      <name>DiagonalFunctionObject</name>
      <anchorfile>classRVL_1_1DiagonalFunctionObject.html</anchorfile>
      <anchor>ab48d867b0440dfaec7fa7c3a697dd836</anchor>
      <arglist>(size_t n, FunctionObject &amp;f)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ProductDataContainer</name>
    <filename>classRVL_1_1ProductDataContainer.html</filename>
    <base>RVL::DataContainer</base>
    <base>Product&lt; DataContainer &gt;</base>
    <member kind="function">
      <type></type>
      <name>ProductDataContainer</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>ae00434473c08f7d96337b315204869ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ProductDataContainer</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>a17bd6f9f2f9f8dc31eb2b1f8c4d9b1e2</anchor>
      <arglist>(const ProductDataContainer &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ProductDataContainer</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>a3cfe3097fa136bd2100f9e98c5977fdd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>a6ca44dc5811dd7e7d9ecbbf76c8d3d71</anchor>
      <arglist>(FunctionObject &amp;f, std::vector&lt; DataContainer const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>af7a112ed5e2c054cd32b75c65f73b31c</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, std::vector&lt; DataContainer const * &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ProductDataContainer.html</anchorfile>
      <anchor>acc8aac67807962794d535b2b2a136292</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::StdProductDataContainer</name>
    <filename>classRVL_1_1StdProductDataContainer.html</filename>
    <base>RVL::ProductDataContainer</base>
    <member kind="function">
      <type></type>
      <name>StdProductDataContainer</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>a731567d0f187c9e829ada55eb78dff42</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StdProductDataContainer</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>a16dc7b820c4b5a18243afcfae71ca779</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>aae9711ef383fa5d0810b5c23a393e67c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataContainer &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>a52ff64e673c436cf01a0b57f99111e29</anchor>
      <arglist>(size_t i)</arglist>
    </member>
    <member kind="function">
      <type>DataContainer const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>a4805a88073936abbb4a6fe2c664d02fc</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>push</name>
      <anchorfile>classRVL_1_1StdProductDataContainer.html</anchorfile>
      <anchor>aed808b77845914b9ecfba98f242b4ade</anchor>
      <arglist>(DataContainerFactory const &amp;f)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ProductSpace</name>
    <filename>classRVL_1_1ProductSpace.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Space</base>
    <base>ROProduct&lt; Space&lt; Scalar &gt; &gt;</base>
    <member kind="function">
      <type></type>
      <name>ProductSpace</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>a4bd948e8bd819f95c7480b3e71786b3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ProductSpace</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>a0a172c0ae339dd42cba76daef59b0a4b</anchor>
      <arglist>(const ProductSpace&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ProductSpace</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>ab4fc3c7f6c42ce32337fe4aab578a2d2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>aa26d5288b21a49469a26ecfb69091082</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>a711563374b612866cc639b3337c0645b</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Scalar</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>ab024d0d01e7f4e47e74f49685019b383</anchor>
      <arglist>(DataContainer const &amp;x, DataContainer const &amp;y) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>af8a46583e839df53f79060241fca6c6b</anchor>
      <arglist>(DataContainer &amp;x) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>a17e233359b49a759e163fa8856859608</anchor>
      <arglist>(Scalar a, DataContainer const &amp;x, Scalar b, DataContainer &amp;y) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ProductSpace.html</anchorfile>
      <anchor>a029382f936645420b5fb67fb1e493740</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::StdProductSpace</name>
    <filename>classRVL_1_1StdProductSpace.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::ProductSpace</base>
    <member kind="function">
      <type></type>
      <name>StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>aa119de7e1a7d01d43883cec02785873c</anchor>
      <arglist>(const StdProductSpace&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>a21654137641b854819c353c5bf2c14f7</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>ab74f15ef4e1c162b147fcb289d8a5362</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;s1, Space&lt; Scalar &gt; const &amp;s2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>a1ea6b20f623e014215c769ed46b23594</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;s1, Space&lt; Scalar &gt; const &amp;s2, Space&lt; Scalar &gt; const &amp;s3)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>a0773749f4e71206eb37a91278dd7b245</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;s1, Space&lt; Scalar &gt; const &amp;s2, Space&lt; Scalar &gt; const &amp;s3, Space&lt; Scalar &gt; const &amp;s4)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~StdProductSpace</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>accf5fb8a658f5f54d3851a6c9bb7a8b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>ac71c6b9f81b0f4c8e31d28855e8016e0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>a48149ea122971dd8907cb1d0bc355ca7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1StdProductSpace.html</anchorfile>
      <anchor>a748a6eaa6d9191da4aa33d9dcc26e6f3</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::CartesianPowerSpace</name>
    <filename>classRVL_1_1CartesianPowerSpace.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::ProductSpace</base>
    <member kind="function">
      <type></type>
      <name>CartesianPowerSpace</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a334cdf96d577d12272d86a86c12a9d1f</anchor>
      <arglist>(size_t _size, const Space&lt; Scalar &gt; &amp;_subspc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CartesianPowerSpace</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>ae6ad3eceb175a630ce9554aa45aca714</anchor>
      <arglist>(const CartesianPowerSpace&lt; Scalar &gt; &amp;s)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~CartesianPowerSpace</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>aa6a764816d05c60c951af70f68ad1155</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>aada9accc69795cd408d655bfb7d237fc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>ab964f0a6dc45e4f7d0836c0ef9a7dfd1</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a242ff7c3f86837a5492c40c094978cb4</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>adc447ab90c5def6112a2bf4469c92871</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp) const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a9971fd293b22dea37f2f1ab7a85b8a57</anchor>
      <arglist>(DataContainer const &amp;x, DataContainer const &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>ad68a40a593ee3c6aa14dff2236e96976</anchor>
      <arglist>(DataContainer &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a32e3a87710666d12e766a375b4056ffc</anchor>
      <arglist>(Scalar a, DataContainer const &amp;x, Scalar b, DataContainer &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a78906659143e9a370c8d762023057eb7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a525f7bff9ebb1c9dc8b77c84973aaa19</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>CartesianPowerSpace</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>afcf51c71c05e4d3024520afc9580baa2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>size_t</type>
      <name>size</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a9cadeadfc654cbeb061eaa9f880b711a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>subspc</name>
      <anchorfile>classRVL_1_1CartesianPowerSpace.html</anchorfile>
      <anchor>a1196555cee171ff269e0b8239194e104</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Components</name>
    <filename>classRVL_1_1Components.html</filename>
    <templarg>Scalar</templarg>
    <base>Product&lt; Vector&lt; Scalar &gt; &gt;</base>
    <member kind="function">
      <type></type>
      <name>Components</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>a1d2b69609c4f961a05410e6ad3665384</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Components</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>a325e5a2142ef1d8979cbabc6a3811eb2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>ad2651c16113ffc197faad0628c1d6b02</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>a168e14b08cf26e3ff64deb1ccda2f420</anchor>
      <arglist>(size_t i)</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>a72948c58af1cd9d66da331287359cf27</anchor>
      <arglist>(size_t i) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Components.html</anchorfile>
      <anchor>ae51b50ed398be9452b3d8710ee1cbfbb</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Space</name>
    <filename>classRVL_1_1Space.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Space</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>ab075d8c78a35ee12496cf0eaf2fc1fbd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Space</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>ae08172aec780eb721d2ec845b2d0b1a2</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Space</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a6dc8de7f603bc1c6452f1603d88a0a21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>afe67eec7137aed9bebd5497b9f169690</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a161154b7546c0859b87d181c7e9ed94c</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp) const =0</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a6a173e432da517007dacfae5699bcb78</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp) const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a9ce71ddc4d9a1b18d1e9a72a736167dd</anchor>
      <arglist>(DataContainer const &amp;dc) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual Scalar</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a3ff95d3052e0f9574a1801cdc90d508d</anchor>
      <arglist>(DataContainer const &amp;x, DataContainer const &amp;y) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a9839e31d760283377c3732e76ba09320</anchor>
      <arglist>(DataContainer &amp;x) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a55c2104166e35ea65aa15b006fb5d0ac</anchor>
      <arglist>(Scalar a, DataContainer const &amp;x, Scalar b, DataContainer &amp;y) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>copy</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>acd5fd9d73ee620a095e119627bc04e1e</anchor>
      <arglist>(DataContainer &amp;tgt, DataContainer const &amp;src) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>negate</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a84dae7aea98c33795c00fcdf9488d310</anchor>
      <arglist>(DataContainer &amp;tgt) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>negate</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a84f99f1b60b047cad1cbc4951f23f2bd</anchor>
      <arglist>(DataContainer &amp;tgt, DataContainer const &amp;src) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>scale</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a3202683b11b8570e47b6ac095166727a</anchor>
      <arglist>(DataContainer &amp;tgt, Scalar c) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>scale</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>abf79d0897e63ddef4fd2c0b0f2b281da</anchor>
      <arglist>(DataContainer &amp;tgt, Scalar c, DataContainer const &amp;src) const </arglist>
    </member>
    <member kind="function">
      <type>NormRetType</type>
      <name>normsq</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>ac6a1538e2690f52d6284dbb70e73fda2</anchor>
      <arglist>(DataContainer const &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>NormRetType</type>
      <name>norm</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a202cdc68762df8c1e676685325598523</anchor>
      <arglist>(DataContainer const &amp;x) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void *</type>
      <name>operator new</name>
      <anchorfile>classRVL_1_1Space.html</anchorfile>
      <anchor>a164fffe58887d2d9e7a051a10defdd9e</anchor>
      <arglist>(size_t size)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::StdSpace</name>
    <filename>classRVL_1_1StdSpace.html</filename>
    <templarg>Scalar</templarg>
    <templarg>DataType</templarg>
    <base>RVL::Space</base>
    <member kind="function">
      <type></type>
      <name>StdSpace</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>ad3b6a11819db1216b03f4f261671f869</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StdSpace</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>a3850ae5c2b0a8047c25ee93938ab8e03</anchor>
      <arglist>(const StdSpace&lt; Scalar, DataType &gt; &amp;sp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~StdSpace</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>ab0f652bc715a996ea29322b5a8389e96</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual DataContainerFactory const &amp;</type>
      <name>getDCF</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>ae5456c2097f4d22bcb9c5a608fa7d75b</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LinearAlgebraPackage&lt; Scalar &gt; const &amp;</type>
      <name>getLAP</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>ad5dd64b00f9ad93d83fd5442a9c29862</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>acded1c8aed3b869bbf49e878f96f1933</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>a58ee0ebcba0ddeb4eb1333d52873d4ac</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>aed0386e0e4261d42192eb81509290d7f</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>a6707245b3fb56dacdb5deb2f5a4996e5</anchor>
      <arglist>(DataContainer const &amp;x, DataContainer const &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>a2f4ecb62724fa0bc1568999ebf65cfb0</anchor>
      <arglist>(DataContainer &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>abbcdf0a5e2150bf651344bab7a4a83f2</anchor>
      <arglist>(Scalar a, DataContainer const &amp;x, Scalar b, DataContainer &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1StdSpace.html</anchorfile>
      <anchor>a23a21b6f9db412e0ac367bc2f297f54b</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::SpaceDCF</name>
    <filename>classRVL_1_1SpaceDCF.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::DataContainerFactory</base>
    <member kind="function">
      <type></type>
      <name>SpaceDCF</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a5a214822acec8a34bb2444c72c0d3819</anchor>
      <arglist>(Space&lt; Scalar &gt; const &amp;_sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SpaceDCF</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a65a949c40f6137beba11a1d6ddb1f3f9</anchor>
      <arglist>(SpaceDCF&lt; Scalar &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SpaceDCF</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a5166719ab4cabec607cac83e1d601850</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>ae0c303dfb465531db66bc4370eb50570</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; Scalar &gt; const &amp;</type>
      <name>getSpace</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a8c3d4571df692cee48080b959395be82</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>ab01de3cafe8233a3f8de24a760e55a99</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a245848f15c17d254575329daf05c786b</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1SpaceDCF.html</anchorfile>
      <anchor>a7d34d792e2dedbb1ecfaac553fe9d696</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Vector</name>
    <filename>classRVL_1_1Vector.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::Writeable</base>
    <member kind="typedef">
      <type>ScalarFieldTraits&lt; Scalar &gt;::AbsType</type>
      <name>NormRetType</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>af412a3efeaccdf39a88058fc9b1bca0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Vector</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aba390cc21a61449dcfbc5ec0b8b17cab</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Vector</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a1a6fa97611c8bac08a4915b2781fecf1</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_sp, bool _initZero=false)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Vector</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aa2504422e2d6b149b6f115c64dd1871b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const Space&lt; Scalar &gt; &amp;</type>
      <name>getSpace</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a062ceae1a6acc51c283d373325c2fb4b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>inSpace</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ac58df0e43a970511131178a5ab63ca87</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp1) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>inSameSpace</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aac297d9a393c1e3a84ff886ed1061c60</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>getVersion</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a8794234ab495a8cde098501ff08bf731</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>incrementVersion</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ac2c1abc8e9349788db01610e14e103e8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ab1999798deea0f938074eb81a7ff600b</anchor>
      <arglist>(FunctionObject &amp;f, vector&lt; Vector&lt; Scalar &gt; const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aabed115ea0fdd55d2a2a1521accf984d</anchor>
      <arglist>(FunctionObject &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aa64e765d7dfb0ac9557c47976ad7df40</anchor>
      <arglist>(FunctionObject &amp;f, const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a9142e94d0ab020f6181f5fc22b070226</anchor>
      <arglist>(FunctionObject &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a6c53391072a9b2f1a87363d8676a50e7</anchor>
      <arglist>(FunctionObject &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;z)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ac14952300202a1d8a52f16f618644141</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, vector&lt; Vector&lt; Scalar &gt; const * &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>af7d7cee3a2e3f9ec6f2f9cf75dde57dd</anchor>
      <arglist>(FunctionObjectConstEval &amp;f) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a4804a5bbaf23a4d5ad8db31dc12e9e7a</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, const Vector&lt; Scalar &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aba42504a825edf6e8093f70139393fd2</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>afda19817c00ad4ce67e26f9069629a4b</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;y, const Vector&lt; Scalar &gt; &amp;z) const </arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a7878dc6cc0c9ba9c66d0dd3f79ebf512</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a5b8eea9d05c1a75c22cee9771a7cba11</anchor>
      <arglist>(Scalar a, const Vector&lt; Scalar &gt; &amp;x, Scalar b=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aa8c47c3554f45732f9ab60b5ff851e3f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copy</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ae2108e8b3a379298915886d167ed4a8c</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a08f8a25fa8bf30073df4b41da4d3c9a8</anchor>
      <arglist>(Scalar c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>scale</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a37360abe75b7da2e641e528634c79882</anchor>
      <arglist>(Scalar c, const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>negate</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ad16e9c4dfe0f0c7fa9e17be0d3f10bdc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>negate</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a495957292786b2621039f0725570a3aa</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>NormRetType</type>
      <name>normsq</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a78b03431942199c003c32530aa7236f6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>NormRetType</type>
      <name>norm</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>afcfd02f45d28fd3e6a7965bb17431093</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a139432993aec993a0523066cc3881798</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>DataContainer *</type>
      <name>getDataContainer</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>ae6f2d3610e9698e70209ab513190550d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>Vector</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aecb83c3ea69cbbaeed0237c0639ea202</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_sp, DataContainer *_d, unsigned int &amp;_verref, bool _own=false)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>Vector&lt; Scalar &gt; *</type>
      <name>build_from_kit</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a4ff05aeaeb0ae5977ee31736f2798b1a</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;_sp, DataContainer *_d, unsigned int &amp;_verref, bool _own=false)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type></type>
      <name>Vector</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>aa8c157d08043d2e1d47352b883e8520d</anchor>
      <arglist>(const Vector&lt; Scalar &gt; *v)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void *</type>
      <name>operator new</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a3758de47490fa3f96df10ac3a410accb</anchor>
      <arglist>(size_t size)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>unsigned int &amp;</type>
      <name>getVersionRef</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>a6542e71b1f9795f1f777ef6f674aca39</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>Components&lt; Scalar &gt;</name>
      <anchorfile>classRVL_1_1Vector.html</anchorfile>
      <anchor>af7ba062e8dae66249e66beb7b150e48e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::WatchedVecRef</name>
    <filename>classRVL_1_1WatchedVecRef.html</filename>
    <templarg>Scalar</templarg>
    <member kind="function">
      <type></type>
      <name>WatchedVecRef</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>a66b6a3cb0eaa35c6079ace4d2b425def</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;_x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>WatchedVecRef</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>a229b4715e84c7a68cff41f92ec59c906</anchor>
      <arglist>(const WatchedVecRef&lt; Scalar &gt; &amp;w)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~WatchedVecRef</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>ae4741991b2ab93786508a0aecad6082b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>a10e744b0c31aae48482cc24b2396111b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Vector&lt; Scalar &gt; const &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>a54bacb6950e7e6d50a3547fdbd0dddcb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>update</name>
      <anchorfile>classRVL_1_1WatchedVecRef.html</anchorfile>
      <anchor>ad6d3f1ad06a9216a57d7705089f6d7e6</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits</name>
    <filename>structRVL_1_1ScalarFieldTraits.html</filename>
    <templarg>Scalar</templarg>
    <member kind="typedef">
      <type>Scalar</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>a22d2441a7eeb8b50c3c62654cf2f7481</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Scalar</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>aa33696400cc692b57d5d902bad8ce6e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static Scalar</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>a3b320ab2ebe856fcd95e6f3afdcd3b12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static Scalar</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>a7005a43bca38fd826842bd54d6936d3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static Scalar</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>a46d2cad49c30e3744d9557b0b5101c71</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static Scalar</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits.html</anchorfile>
      <anchor>a7fcf5c63fb8a8c2d4bbc67b5cf0fa3f3</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; bool &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</filename>
    <member kind="typedef">
      <type>bool</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>aa4cdcf74a494cf08e4b57ef379f7af3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>bool</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>a694a9aabef473dec6effc71ca87fa237</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>a77a59175f99e31f8e3b432125bbff88a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>a547a11599b8e763cd70292f021ef42c5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>a21afab7f686d42a77dbafa416f80adbf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01bool_01_4.html</anchorfile>
      <anchor>a88a70f63796dd5a0ae176cd32c313f49</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; int &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</filename>
    <member kind="typedef">
      <type>int</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>a82d63fa5166a8f1e5b6acaee22dc9ea5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>abfccabf110c48094b8172e2a05c3c28a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>ae7273b9bf67132a6d0cfbae9cbdc7db3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>a85f0ea32d3d96b8cfa5109b803e37b61</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>a13e4aaf2c0bd9f1cbd23421ab07bacb3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static int</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01int_01_4.html</anchorfile>
      <anchor>ac460cd0eaaba3d9822675b28ca56587d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; long &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</filename>
    <member kind="typedef">
      <type>long</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>a3038b7cff807d56c9c69bb8823a2de7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>long</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>aa53fe2c5ba8289a4727e88761bc64248</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static long</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>a3fba4dcb244d00388c9c4c9bc58fc8c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static long</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>ad13ee11772323910f3b929ac9bbe37a4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static long</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>a0970be87579fa55e5e8331c2e0d0fbff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static long</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01long_01_4.html</anchorfile>
      <anchor>a1872b1877ee3df560c297ce4ce471cca</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; unsigned int &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</filename>
    <member kind="typedef">
      <type>unsigned int</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>afaa1609d3738b8f2d7e33958d26257ff</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>unsigned int</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>af2a4992bdcbeed1c36c0216015c4d295</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>aadc0ddb0a910358998102e67d5f31c12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>a38bb344b38c0bd7c74c889666c08800e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>a7545843755bbdf6187cd32a9bfe69837</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned int</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01int_01_4.html</anchorfile>
      <anchor>ada7f7114f7a8bc0e3d0fe84a80cab29e</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; unsigned long &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</filename>
    <member kind="typedef">
      <type>unsigned long</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>a9d8f236619e441ddf68be0751c07d4af</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>unsigned long</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>ae1f2a3e696ee5ed82e634ff5535c01c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned long</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>aefeb4d8716d40d98959192f75476ada8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned long</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>af122dd4910b64ee27034189bf5ca8ec1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned long</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>ae988cc584d48e4d70f7bcbef025d29bb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static unsigned long</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01unsigned_01long_01_4.html</anchorfile>
      <anchor>a9919eb029776bbaeebfcb2c9fb23eda7</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; float &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</filename>
    <member kind="typedef">
      <type>float</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>a74ecb186827dc255eb8a894e3b6398cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>float</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>a0006f4f7c079d4175d84a0e1d06fb4e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static float</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>afedbdb4c89b3033b55a836853c200c41</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static float</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>a4aff022e633c0db39982d68363bbf9cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static float</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>a84909d697feb3b68b99cbdd9be61058f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static float</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01float_01_4.html</anchorfile>
      <anchor>a89ebb44a9012a23246134b28e9c3001b</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; double &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</filename>
    <member kind="typedef">
      <type>double</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>a37467539bf4cdbd68c8e45353c5da610</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>double</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>aa018714ac14cc0325187eb5c58d9e7d1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>a734feb67e6fefbd606dc40cb9ddea977</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>ab8430b3a39965cd46962b84bff63acbd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>a15bf6b39d98dc648790743b7cb947c25</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static double</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01double_01_4.html</anchorfile>
      <anchor>a6bb52be7ee18be6932032f076ad77935</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>RVL::ScalarFieldTraits&lt; std::complex&lt; T &gt; &gt;</name>
    <filename>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</filename>
    <templarg></templarg>
    <member kind="typedef">
      <type>complex&lt; T &gt;</type>
      <name>ScalarType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>a8b02ba925842f33ee5a7d68a0bc9cbdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>T</type>
      <name>AbsType</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>a459efbffb732636c0a982ef42cb7f3b9</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" static="yes">
      <type>static complex&lt; T &gt;</type>
      <name>Zero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>a2c669b6405ba4609fb7407ae2147f437</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static complex&lt; T &gt;</type>
      <name>One</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>a9b0dc8c4c9d5a6ee0c13d56716b6c8af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static T</type>
      <name>AbsZero</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>ae4a2ca9e3285dad3707eb18564a27a34</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" static="yes">
      <type>static T</type>
      <name>AbsOne</name>
      <anchorfile>structRVL_1_1ScalarFieldTraits_3_01std_1_1complex_3_01T_01_4_01_4.html</anchorfile>
      <anchor>a6d0ee28349dd3ef4748aa973f8da6141</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Writeable</name>
    <filename>classRVL_1_1Writeable.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Writeable.html</anchorfile>
      <anchor>ac4947a1a807aaebc1218651dfb14bbdf</anchor>
      <arglist>(ostream &amp;str) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Writeable</name>
      <anchorfile>classRVL_1_1Writeable.html</anchorfile>
      <anchor>a76fad523a09eb3320b8e495df7b48c7f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1Writeable.html</anchorfile>
      <anchor>aaa091cfe9ecf5a0f8b9c519d3102669c</anchor>
      <arglist>(RVLException &amp;e) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Oracle</name>
    <filename>classRVL_1_1Oracle.html</filename>
    <templarg></templarg>
    <base>RVL::Writeable</base>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>isFeasible</name>
      <anchorfile>classRVL_1_1Oracle.html</anchorfile>
      <anchor>a5e3e0aec08eea0c6a1b3dcae35a15e35</anchor>
      <arglist>(T const &amp;) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::Factory</name>
    <filename>classRVL_1_1Factory.html</filename>
    <templarg>T</templarg>
    <base>RVL::Writeable</base>
    <member kind="function">
      <type></type>
      <name>Factory</name>
      <anchorfile>classRVL_1_1Factory.html</anchorfile>
      <anchor>a145d706482ab0eaf40643e6b078a09f2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Factory</name>
      <anchorfile>classRVL_1_1Factory.html</anchorfile>
      <anchor>aa8c00eee4a75fce8467abe278e634269</anchor>
      <arglist>(Factory&lt; T &gt; const &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Factory</name>
      <anchorfile>classRVL_1_1Factory.html</anchorfile>
      <anchor>a8d352720d6d151d40b9a2262a55c26d5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual T *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1Factory.html</anchorfile>
      <anchor>adf2eaee4852e76630cd54cc1f4528ee3</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
</tagfile>
