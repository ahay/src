<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="page">
    <name>index</name>
    <title>Local RVL - a simple realization of the RVL data management classes, based on data containers which explose a scalar array</title>
    <filename>index</filename>
  </compound>
  <compound kind="file">
    <name>contentpackage.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>contentpackage_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <includes id="dejavu_8hh" name="dejavu.hh" local="yes" imported="no">dejavu.hh</includes>
    <class kind="class">RVL::ContentPackage</class>
    <class kind="class">RVL::PackageContainer</class>
    <class kind="class">RVL::PackageContainerFactory</class>
    <class kind="class">RVL::SingleDataContainer</class>
    <class kind="class">RVL::SingleDataContainerFactory</class>
    <member kind="function">
      <type>size_t</type>
      <name>getDataSize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aae08ca24133e3e01b5d8d43d2c037698</anchor>
      <arglist>(MetaType const &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getMetaSize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a577425785ce4e2b8cbaf9f33fe6f3428</anchor>
      <arglist>(MetaType const &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a209d727550a397da5144fbe6c5f46057</anchor>
      <arglist>(MetaType const &amp;mt, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>MetaType *</type>
      <name>deserialize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a4b2931b8e8d08f1abab289666023ce60</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>char *</type>
      <name>serialize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ac67d81369d5b7a7528bed31f91539a4d</anchor>
      <arglist>(ContentPackage&lt; DataType, MetaType &gt; const &amp;cp, size_t &amp;len)</arglist>
    </member>
    <member kind="function">
      <type>ContentPackage&lt; DataType, MetaType &gt; *</type>
      <name>deserialize</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a3d5ab8b89b83f6fbaa64710bcf72ca22</anchor>
      <arglist>(char *cbuf, size_t len)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>writeMeta</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a1d4ee5a7219c3b79a5244f4cef7ef411</anchor>
      <arglist>(MetaType const &amp;md, ostream &amp;e)</arglist>
    </member>
    <member kind="function">
      <type>DataType *</type>
      <name>newData</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a9026e8d302057dfdd80789ffd49f2d8e</anchor>
      <arglist>(MetaType &amp;md)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>deleteData</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ad86f7f66bc97c5df72a68aae976b08b8</anchor>
      <arglist>(DataType **d, MetaType **md)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>copyCP</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a9a2b30c8ad9d818fdf8da64697c13fcd</anchor>
      <arglist>(Source const &amp;in, Target &amp;out)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>dejavu.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>dejavu_8hh</filename>
    <member kind="function">
      <type>void</type>
      <name>dejavu</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>a8e1a196dc2b525039e713a499a7045df</anchor>
      <arglist>(size_t *i, std::vector&lt; T * &gt; vec)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>doc.h</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>doc_8h</filename>
  </compound>
  <compound kind="file">
    <name>fcnfo.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>fcnfo_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <class kind="class">RVL::ScalarFO1</class>
    <class kind="class">RVL::ScalarFO2</class>
    <class kind="class">RVL::ScalarFO3</class>
    <class kind="class">RVL::ScalarFO4</class>
    <class kind="class">RVL::ScalarFO5</class>
    <class kind="class">RVL::ScalarFO6</class>
    <member kind="function">
      <type>float</type>
      <name>test1</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>aee4220b7ea65aee283e484f1d1a3b7cd</anchor>
      <arglist>(float x)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>test2</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ad96676fb439f0c10b01274b5a1a75b4f</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>ScalarFOSanity</name>
      <anchorfile>namespaceRVL.html</anchorfile>
      <anchor>ad4845e5d7fd15ee231c564905a54d4bc</anchor>
      <arglist>(int na, LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>functions.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>functions_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <class kind="class">RVL::RVLCopy</class>
    <class kind="class">RVL::RVLCopy&lt; float &gt;</class>
    <class kind="class">RVL::RVLCopy&lt; double &gt;</class>
    <class kind="class">RVL::RVLScale</class>
    <class kind="class">RVL::RVLMax</class>
    <class kind="class">RVL::RVLMin</class>
    <class kind="class">RVL::RVLL2innerProd</class>
    <class kind="class">RVL::RVLL2innerProd&lt; complex&lt; Scalar &gt; &gt;</class>
    <class kind="class">RVL::RVLAddAccumulate</class>
    <class kind="class">RVL::RVLAssignConst</class>
    <class kind="class">RVL::RVLRandomize</class>
    <class kind="class">RVL::RVLRandomize&lt; complex&lt; Scalar &gt; &gt;</class>
    <class kind="class">RVL::RVLRandomize&lt; int &gt;</class>
    <class kind="class">RVL::ASCIIReader</class>
    <class kind="class">RVL::ASCIIWriter</class>
    <class kind="class">RVL::BinaryReader</class>
    <class kind="class">RVL::BinaryWriter</class>
    <class kind="class">RVL::RVLBoxMaxStep</class>
    <class kind="class">RVL::ElementwiseMultiply</class>
    <class kind="class">RVL::ElementwiseDivision</class>
    <class kind="class">RVL::ElementwiseSqrtAbs</class>
    <class kind="class">RVL::RVLScalarLogistic</class>
    <class kind="class">RVL::RVLScalarLogisticInverse</class>
    <class kind="class">RVL::RVLScalarLogisticDeriv</class>
    <class kind="class">RVL::RVLVectorLogistic</class>
  </compound>
  <compound kind="file">
    <name>local.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>local_8hh</filename>
    <includes id="localdata_8hh" name="localdata.hh" local="yes" imported="no">localdata.hh</includes>
    <includes id="localevaluation_8hh" name="localevaluation.hh" local="yes" imported="no">localevaluation.hh</includes>
    <includes id="localreduction_8hh" name="localreduction.hh" local="yes" imported="no">localreduction.hh</includes>
  </compound>
  <compound kind="file">
    <name>localdata.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>localdata_8hh</filename>
    <class kind="class">RVL::LocalDataContainer</class>
    <class kind="class">RVL::LocalDataContainerFactory</class>
    <class kind="class">RVL::LocalDataContainerSection</class>
  </compound>
  <compound kind="file">
    <name>localevaluation.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>localevaluation_8hh</filename>
    <includes id="localdata_8hh" name="localdata.hh" local="yes" imported="no">localdata.hh</includes>
    <class kind="class">RVL::LocalEvaluation</class>
    <class kind="class">RVL::LocalFunctionObject</class>
    <class kind="class">RVL::UnaryLocalEvaluation</class>
    <class kind="class">RVL::UnaryLocalFunctionObject</class>
    <class kind="class">RVL::BinaryLocalEvaluation</class>
    <class kind="class">RVL::BinaryLocalFunctionObject</class>
    <class kind="class">RVL::TernaryLocalEvaluation</class>
    <class kind="class">RVL::TernaryLocalFunctionObject</class>
    <class kind="class">RVL::QuaternaryLocalEvaluation</class>
    <class kind="class">RVL::QuaternaryLocalFunctionObject</class>
  </compound>
  <compound kind="file">
    <name>locallinalg.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>locallinalg_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <includes id="functions_8hh" name="functions.hh" local="yes" imported="no">functions.hh</includes>
    <class kind="class">RVL::LocalLinearAlgebraPackage</class>
    <class kind="class">RVL::RVLLinCombObject</class>
    <class kind="class">RVL::RVLLinearAlgebraPackage</class>
  </compound>
  <compound kind="file">
    <name>localproduct.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>localproduct_8hh</filename>
    <includes id="localdata_8hh" name="localdata.hh" local="yes" imported="no">localdata.hh</includes>
    <class kind="class">RVL::ProductLocalDataContainer</class>
    <class kind="class">RVL::ProductDataContainerLDC</class>
    <class kind="class">RVL::PartitionedLocalDataContainer</class>
  </compound>
  <compound kind="file">
    <name>localreduction.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>localreduction_8hh</filename>
    <includes id="localdata_8hh" name="localdata.hh" local="yes" imported="no">localdata.hh</includes>
    <class kind="class">RVL::LocalConstEval</class>
    <class kind="class">RVL::UnaryLocalConstEval</class>
    <class kind="class">RVL::BinaryLocalConstEval</class>
    <class kind="class">RVL::TernaryLocalConstEval</class>
    <class kind="class">RVL::QuaternaryLocalConstEval</class>
    <class kind="class">RVL::UnaryLocalFunctionObjectConstEval</class>
    <class kind="class">RVL::BinaryLocalFunctionObjectConstEval</class>
    <class kind="class">RVL::TernaryLocalFunctionObjectConstEval</class>
    <class kind="class">RVL::QuaternaryLocalFunctionObjectConstEval</class>
    <class kind="class">RVL::UnaryLocalFunctionObjectScalarRedn</class>
    <class kind="class">RVL::BinaryLocalFunctionObjectScalarRedn</class>
    <class kind="class">RVL::TernaryLocalFunctionObjectScalarRedn</class>
    <class kind="class">RVL::QuaternaryLocalFunctionObjectScalarRedn</class>
  </compound>
  <compound kind="file">
    <name>localspace.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>localspace_8hh</filename>
    <includes id="localdata_8hh" name="localdata.hh" local="yes" imported="no">localdata.hh</includes>
    <class kind="class">RVL::LocalSpace</class>
    <class kind="class">RVL::LocalVector</class>
  </compound>
  <compound kind="file">
    <name>pcspace.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>pcspace_8hh</filename>
    <includes id="contentpackage_8hh" name="contentpackage.hh" local="yes" imported="no">contentpackage.hh</includes>
    <includes id="locallinalg_8hh" name="locallinalg.hh" local="yes" imported="no">locallinalg.hh</includes>
    <class kind="class">RVL::PackageContainerSpace</class>
  </compound>
  <compound kind="file">
    <name>polyop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>polyop_8hh</filename>
    <includes id="functions_8hh" name="functions.hh" local="yes" imported="no">functions.hh</includes>
    <class kind="class">RVL::PolynomialOperator</class>
  </compound>
  <compound kind="file">
    <name>rn.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>rn_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <class kind="class">RVL::RnArray</class>
    <class kind="class">RVL::RnDataContainerFactory</class>
  </compound>
  <compound kind="file">
    <name>rnmat.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>rnmat_8hh</filename>
    <includes id="local_8hh" name="local.hh" local="yes" imported="no">local.hh</includes>
    <includes id="rnspace_8hh" name="rnspace.hh" local="yes" imported="no">rnspace.hh</includes>
    <class kind="class">RVL::matvec</class>
    <class kind="class">RVL::fmatvec</class>
    <class kind="class">RVL::amatvec</class>
    <class kind="class">RVL::GenMat</class>
    <class kind="class">RVL::SymMat</class>
  </compound>
  <compound kind="file">
    <name>rnop.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>rnop_8hh</filename>
    <includes id="rnmat_8hh" name="rnmat.hh" local="yes" imported="no">rnmat.hh</includes>
    <class kind="class">RVL::CFunction</class>
    <class kind="class">RVL::CJacobian</class>
    <class kind="class">RVL::OpWithGenMatDeriv</class>
    <class kind="class">RVL::GenOp</class>
  </compound>
  <compound kind="file">
    <name>rnspace.hh</name>
    <path>/Users/williamsymes/Applications/RSFSRC/trip/rvl/local/include/</path>
    <filename>rnspace_8hh</filename>
    <includes id="functions_8hh" name="functions.hh" local="yes" imported="no">functions.hh</includes>
    <includes id="rn_8hh" name="rn.hh" local="yes" imported="no">rn.hh</includes>
    <includes id="localspace_8hh" name="localspace.hh" local="yes" imported="no">localspace.hh</includes>
    <includes id="locallinalg_8hh" name="locallinalg.hh" local="yes" imported="no">locallinalg.hh</includes>
    <class kind="class">RVL::RnSpace</class>
  </compound>
  <compound kind="class">
    <name>RVL::ContentPackage</name>
    <filename>classRVL_1_1ContentPackage.html</filename>
    <templarg>DataType</templarg>
    <templarg>MetaType</templarg>
    <base>RVL::LocalDataContainer</base>
    <member kind="function">
      <type></type>
      <name>ContentPackage</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>acf65cf7806cee935803f5e4c8c690aed</anchor>
      <arglist>(const ContentPackage&lt; DataType, MetaType &gt; &amp;bin)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ContentPackage</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a0b3e7ad0c09da94824b5d15952e5bda2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ContentPackage</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>abf6ea358bae9fffe575497cc07eb66eb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>initialize</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>af894d8d758dfc04a1eaddc21712a245e</anchor>
      <arglist>(MetaType const &amp;_md, LocalDataContainer&lt; DataType &gt; *p=NULL, size_t start=0)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isInitialized</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a80cee4b94e4dfa679218b36a847fb69e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ContentPackage&lt; DataType, MetaType &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>ae6ffa27eb49fa3bf735690de113ab98d</anchor>
      <arglist>(ContentPackage&lt; DataType, MetaType &gt; const &amp;src)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a6b1183c7af22c8e21963eac59d448795</anchor>
      <arglist>(ContentPackage&lt; DataType, MetaType &gt; const &amp;a) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a5411f1117b4dfe46869f0651855c6fd5</anchor>
      <arglist>(ContentPackage&lt; DataType, MetaType &gt; const &amp;a) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a9bf7dbc2fcb8f7061033ee1e4b604274</anchor>
      <arglist>(ContentPackage&lt; DataType, MetaType &gt; const &amp;a) const </arglist>
    </member>
    <member kind="function">
      <type>MetaType &amp;</type>
      <name>getMetadata</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a625e7fbbce32852e332fb15968a8aa90</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MetaType const &amp;</type>
      <name>getMetadata</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a0570f27578b0718a20c164bf97ee02c0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a238614ea3d4d13a286a4354f2e26b3f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataType *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>ae535cf04bfedb726ed5ec208dd0de964</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataType const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a50793f87a4417e7725b9626138164679</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1ContentPackage.html</anchorfile>
      <anchor>a9bda6686967a1f6067f85b89579299c8</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PackageContainer</name>
    <filename>classRVL_1_1PackageContainer.html</filename>
    <templarg>DataType</templarg>
    <templarg>MetaType</templarg>
    <base>RVL::DataContainer</base>
    <member kind="function">
      <type></type>
      <name>PackageContainer</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a0ed22c86a752c945808f81ab70fb9c21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PackageContainer</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a5325b486603ad98484da2a04fa6d546a</anchor>
      <arglist>(PackageContainer&lt; DataType, MetaType &gt; const &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PackageContainer</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a48c1fdcb0dd21477f4348d8eab3885d2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a7b5b768ddd6e32e315ff6d81d85fa66d</anchor>
      <arglist>(FunctionObject &amp;f, std::vector&lt; DataContainer const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>aedb3a333bd42bdb0fa5c265524c2b3bd</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, vector&lt; DataContainer const * &gt; &amp;x) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual ContentPackage&lt; DataType, MetaType &gt; &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a0dac23046fdcf39fa4f359c16ff00464</anchor>
      <arglist>(bool &amp;more)=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual ContentPackage&lt; DataType, MetaType &gt; const &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>ac6d8d70eacbf56ed310331e7a1cbc216</anchor>
      <arglist>(bool &amp;more) const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>put</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>a28f230972688def2ff6ad9c72142de2d</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual void</type>
      <name>reset</name>
      <anchorfile>classRVL_1_1PackageContainer.html</anchorfile>
      <anchor>ae72b2568e4e6ae76ad5127e0972de4e3</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PackageContainerFactory</name>
    <filename>classRVL_1_1PackageContainerFactory.html</filename>
    <templarg>DataType</templarg>
    <templarg>MetaType</templarg>
    <base>RVL::DataContainerFactory</base>
    <member kind="function" virtualness="pure">
      <type>virtual PackageContainerFactory&lt; DataType, MetaType &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>ac277f5a623e337482ad006c9abfc7d9e</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>a47d80dc88c01bff703f226535b799f4e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>afc5c92ac423809cda7290ef52a89fcf0</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>abf7046f4661cdb9a4c723a5ceb86c16a</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>a36649ad34c0661cd57d0c543816eb126</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual PackageContainer&lt; DataType, MetaType &gt; *</type>
      <name>buildPC</name>
      <anchorfile>classRVL_1_1PackageContainerFactory.html</anchorfile>
      <anchor>ab6e42929318ce89e3e5e88a685ab93dd</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::SingleDataContainer</name>
    <filename>classRVL_1_1SingleDataContainer.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>PackageContainer&lt; Datatype, Metatype &gt;</base>
    <member kind="function">
      <type></type>
      <name>SingleDataContainer</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a0367af29213723346a4a2e70002d821d</anchor>
      <arglist>(int _nrep=1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SingleDataContainer</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>ac6869d573f0d8bea4b6222fe7dcb98d5</anchor>
      <arglist>(SingleDataContainer&lt; Datatype, Metatype &gt; const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SingleDataContainer</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a685e114222275ca70de21f70ea0dd417</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a2acc7b5f3b46e9db34310ce7bc6c6037</anchor>
      <arglist>(Metatype const &amp;g)</arglist>
    </member>
    <member kind="function">
      <type>Metatype const &amp;</type>
      <name>getMetadata</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a87a69078fe7e1982f753d29a0ae894e7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a02287f5255396d86b6a4c004a55ecf21</anchor>
      <arglist>(ostream &amp;e) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>ContentPackage&lt; Datatype, Metatype &gt; &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a3513a191d3b6f14a294ed8b9a0e75617</anchor>
      <arglist>(bool &amp;more)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>ContentPackage&lt; Datatype, Metatype &gt; const &amp;</type>
      <name>get</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>ada3224864d5f7644087a815a2ae46e63</anchor>
      <arglist>(bool &amp;more) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>put</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a59c7eed77d44f0eff889e95f066cffb6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>reset</name>
      <anchorfile>classRVL_1_1SingleDataContainer.html</anchorfile>
      <anchor>a92076a772bf7308b593e128647f24732</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::SingleDataContainerFactory</name>
    <filename>classRVL_1_1SingleDataContainerFactory.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>PackageContainerFactory&lt; Datatype, Metatype &gt;</base>
    <member kind="function">
      <type></type>
      <name>SingleDataContainerFactory</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a3affea9a854ea8a77503c8211cd13734</anchor>
      <arglist>(int _nrep=1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SingleDataContainerFactory</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a36d6a8154918751f3c8bf6c428ba9a7d</anchor>
      <arglist>(SingleDataContainerFactory&lt; Datatype, Metatype &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SingleDataContainerFactory</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>ae9769ec2c0cac52e31b50ca2ac48539a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>PackageContainerFactory&lt; Datatype, Metatype &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a37449d7034e242c821bc4cfaf096ef53</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>initialize</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>ac0df74e79037af561b64d15dd71d85c8</anchor>
      <arglist>(Metatype const &amp;_g)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a62720416d48969e33d47f53bbcf581d4</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a4a022497871c84fbae0b2f7f4c368321</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>ae262567d5b1482be9e53be500b054973</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>PackageContainer&lt; Datatype, Metatype &gt; *</type>
      <name>buildPC</name>
      <anchorfile>classRVL_1_1SingleDataContainerFactory.html</anchorfile>
      <anchor>a06f904e4aee70a72df9d9da78905fc61</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO1</name>
    <filename>classRVL_1_1ScalarFO1.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO1.html</anchorfile>
      <anchor>a7bcdf1687bbc86652569078081dc02f6</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO1.html</anchorfile>
      <anchor>af16c3a91ffbdc9f2380a4f5b734d21f5</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO2</name>
    <filename>classRVL_1_1ScalarFO2.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO2.html</anchorfile>
      <anchor>a8f601c80fa63ab4e5b34fea7fb3dcfce</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO2.html</anchorfile>
      <anchor>a82f1e344e8eb4647f60190ab428e3b43</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO3</name>
    <filename>classRVL_1_1ScalarFO3.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO3.html</anchorfile>
      <anchor>a7d79785fbcd662141c93886f0c9d6788</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO3.html</anchorfile>
      <anchor>ae50a9b4a2cf41dc019f4ba88bf265246</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO4</name>
    <filename>classRVL_1_1ScalarFO4.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO4.html</anchorfile>
      <anchor>aa0c5d829f94cd9c0b3b9f9ace98ba860</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO4.html</anchorfile>
      <anchor>adfab25e38aa5cdb8f25317a887836203</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO5</name>
    <filename>classRVL_1_1ScalarFO5.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO5.html</anchorfile>
      <anchor>a7e68a8e346c681668e682f44eeac883c</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO5.html</anchorfile>
      <anchor>aeffcff0cb0d4f3dd4ad9333f84d9914e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ScalarFO6</name>
    <filename>classRVL_1_1ScalarFO6.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>LocalFunctionObject&lt; T &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ScalarFO6.html</anchorfile>
      <anchor>a84b125f00703dbcf00e055b0d0641c36</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;target, vector&lt; LocalDataContainer&lt; T &gt; const * &gt; &amp;sources)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ScalarFO6.html</anchorfile>
      <anchor>ad032981cdbeca8ca84d30108b52dbbc0</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLCopy</name>
    <filename>classRVL_1_1RVLCopy.html</filename>
    <templarg>Scalar</templarg>
    <base>BinaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy.html</anchorfile>
      <anchor>ad532b7738ac0eab3bce138eee178ac09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy.html</anchorfile>
      <anchor>a891ed691af2faf33b86101e975b67a9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLCopy.html</anchorfile>
      <anchor>a02231552013394056c0d73e99fc6f747</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLCopy.html</anchorfile>
      <anchor>a41cd273f90f9ab833752027d2cd63001</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLCopy&lt; float &gt;</name>
    <filename>classRVL_1_1RVLCopy_3_01float_01_4.html</filename>
    <base>BinaryLocalFunctionObject&lt; float &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01float_01_4.html</anchorfile>
      <anchor>ac34fdc24196d2d72efac66fc4a5dfdd0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01float_01_4.html</anchorfile>
      <anchor>a1fb3b9fd647bab58e93e173131d8aeaa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01float_01_4.html</anchorfile>
      <anchor>ac0b2e6a6e1aa6bef8666048dcc4cf119</anchor>
      <arglist>(LocalDataContainer&lt; float &gt; &amp;x, LocalDataContainer&lt; float &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01float_01_4.html</anchorfile>
      <anchor>a91fc4647e68e9b9f1feda4e85f29ffaa</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLCopy&lt; double &gt;</name>
    <filename>classRVL_1_1RVLCopy_3_01double_01_4.html</filename>
    <base>BinaryLocalFunctionObject&lt; double &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01double_01_4.html</anchorfile>
      <anchor>aeab8bf29e55499c85e5f25be8ce6c22a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLCopy</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01double_01_4.html</anchorfile>
      <anchor>a3775c206baa3036e92647d72fd2b97be</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01double_01_4.html</anchorfile>
      <anchor>ab4e50686449da4cb38731afec6563e3b</anchor>
      <arglist>(LocalDataContainer&lt; double &gt; &amp;x, LocalDataContainer&lt; double &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLCopy_3_01double_01_4.html</anchorfile>
      <anchor>a85c1d1bfa5dfab14fd218cf36527225f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLScale</name>
    <filename>classRVL_1_1RVLScale.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLScale</name>
      <anchorfile>classRVL_1_1RVLScale.html</anchorfile>
      <anchor>adf5e38af6b78186d51f8667e8f144d0f</anchor>
      <arglist>(Scalar _c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLScale</name>
      <anchorfile>classRVL_1_1RVLScale.html</anchorfile>
      <anchor>a1b82edc46b992b7f7df88f839d37d63b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLScale.html</anchorfile>
      <anchor>a7fc1baadb86913c2ff9ed57326014ee8</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLScale.html</anchorfile>
      <anchor>a9b22b2368ab27368b8a227b13954ce56</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLMax</name>
    <filename>classRVL_1_1RVLMax.html</filename>
    <templarg>Scalar</templarg>
    <base>UnaryLocalFunctionObjectScalarRedn&lt; Scalar, Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLMax</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>aa626aae5b0e7e403be9229782b8371e7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLMax</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>a207f86b2a5b2f9c27419adb5a3b3a988</anchor>
      <arglist>(Scalar _res)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLMax</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>a0fd30ffbbd9167f66cbd03594736f02d</anchor>
      <arglist>(const RVLMax&lt; Scalar &gt; &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLMax</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>ad6a4c6c29640c444c4f77bce49f223ef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>a69cbd02617c8efd84a0f7a27b2f2af91</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>af21a8b469534c3faf4a2ed1c436ff79b</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLMax.html</anchorfile>
      <anchor>a2b274f29508f154d4b51a5bf98805a59</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLMin</name>
    <filename>classRVL_1_1RVLMin.html</filename>
    <templarg>Scalar</templarg>
    <base>UnaryLocalFunctionObjectScalarRedn&lt; Scalar, Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLMin</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a9a633548e8e4853fd97c2018455de587</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLMin</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>af0833f2931e0043ddd8e6d9ae49dfb80</anchor>
      <arglist>(Scalar _res)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLMin</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a669a0ff9ffea833c3633fde032f17cc9</anchor>
      <arglist>(const RVLMin&lt; Scalar &gt; &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLMin</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a8fdfa3aea460f656cd939a62b8c6720c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a36d5436be6a620e1b8e5b3833085e0f1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a5aadab07a72651c58530471c9577e751</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLMin.html</anchorfile>
      <anchor>a124e3673ee80572acf9c356f47c155bb</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLL2innerProd</name>
    <filename>classRVL_1_1RVLL2innerProd.html</filename>
    <templarg>Scalar</templarg>
    <base>BinaryLocalFunctionObjectScalarRedn&lt; Scalar, Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>a7a3dde4c569285e19fc754c70e74e8d6</anchor>
      <arglist>(Scalar _scale=ScalarFieldTraits&lt; Scalar &gt;::One(), Scalar _init=ScalarFieldTraits&lt; Scalar &gt;::Zero())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>a0e3f152ddfe1e9f97add2b5ad752d0ff</anchor>
      <arglist>(const RVLL2innerProd&lt; Scalar &gt; &amp;ipfo)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>a2063f03e11cd4770c4ebd20e9f6cce76</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>ab409486895b62ff7ba93fdb884c622a3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>a423f5afe87426f2de8d3ed4af7134bcf</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;v, LocalDataContainer&lt; Scalar &gt; const &amp;w)</arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getScale</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>ac1e7ab7ad420c083d73f6e70e33d5506</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>a763a952838b8775205d96398b13f5b7d</anchor>
      <arglist>(Scalar newscale)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLL2innerProd.html</anchorfile>
      <anchor>ac9796710df53558f709e7acf64500e7d</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLL2innerProd&lt; complex&lt; Scalar &gt; &gt;</name>
    <filename>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObjectScalarRedn&lt; complex&lt; Scalar &gt;, complex&lt; Scalar &gt; &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>ae0af39e533f260ff42ec65396616f222</anchor>
      <arglist>(Scalar _scale=ScalarFieldTraits&lt; Scalar &gt;::One(), complex&lt; Scalar &gt; _init=complex&lt; Scalar &gt;(ScalarFieldTraits&lt; Scalar &gt;::Zero()))</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a22a6a36993e804788e1461a3e4385f8d</anchor>
      <arglist>(const RVLL2innerProd&lt; complex&lt; Scalar &gt; &gt; &amp;ipfo)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLL2innerProd</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a9d963e361427a7097afe4a797dfe0e62</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setValue</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>ae18eca6675eff696b26a5158bc1d6cae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>afd85b0c6f8e2b97b6e974290ac5acfcf</anchor>
      <arglist>(LocalDataContainer&lt; complex&lt; Scalar &gt; &gt; const &amp;v, LocalDataContainer&lt; complex&lt; Scalar &gt; &gt; const &amp;w)</arglist>
    </member>
    <member kind="function">
      <type>Scalar</type>
      <name>getScale</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>aae3efdd167d16dbfd8a7115c252932f1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a204548b5f64881068a36eb3c35f13f1e</anchor>
      <arglist>(Scalar newscale)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLL2innerProd_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a63b0aaf66fdd1e2fd46d90aeab7442fc</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLAddAccumulate</name>
    <filename>classRVL_1_1RVLAddAccumulate.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLAddAccumulate</name>
      <anchorfile>classRVL_1_1RVLAddAccumulate.html</anchorfile>
      <anchor>ac398ac60a091fb088ca94e2d5cc60fef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLAddAccumulate</name>
      <anchorfile>classRVL_1_1RVLAddAccumulate.html</anchorfile>
      <anchor>a2f135558df62af7157ae45b87577532f</anchor>
      <arglist>(const RVLAddAccumulate&lt; Scalar &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~RVLAddAccumulate</name>
      <anchorfile>classRVL_1_1RVLAddAccumulate.html</anchorfile>
      <anchor>a6da8e12cc9f1a92403ac3fa3ba677473</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLAddAccumulate.html</anchorfile>
      <anchor>ac3f2a635ab603fe627d3fef7ce951fb6</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;v, LocalDataContainer&lt; Scalar &gt; const &amp;w)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLAddAccumulate.html</anchorfile>
      <anchor>ab75a1d3d50ab355f1b3fa50ab63c920e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLAssignConst</name>
    <filename>classRVL_1_1RVLAssignConst.html</filename>
    <templarg>Scalar</templarg>
    <base>UnaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLAssignConst</name>
      <anchorfile>classRVL_1_1RVLAssignConst.html</anchorfile>
      <anchor>ac0f5b4de366851df967883f2b6890244</anchor>
      <arglist>(Scalar _c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLAssignConst</name>
      <anchorfile>classRVL_1_1RVLAssignConst.html</anchorfile>
      <anchor>a46592795f43a457b3cb3770cac6e02ca</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLAssignConst.html</anchorfile>
      <anchor>a87a9d046a103ae42f9da4d932fdffcf2</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLAssignConst.html</anchorfile>
      <anchor>a69ce02a0d28ce3baf1e6f0d99ebfc69a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLRandomize</name>
    <filename>classRVL_1_1RVLRandomize.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>a87cf5aee7f20de0beaca1bb6f03d7a33</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>aefd5350867097f03eea5426cac09b4a7</anchor>
      <arglist>(long seed, Scalar _a=0, Scalar _b=1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>ade81d6a40550523b3006b265eeea12c7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>readsData</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>a5cf645b08da2e3a7e04bba7abe05402f</anchor>
      <arglist>(size_t i=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>aa169284cfefaf6f5eabd0b70e57a748c</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLRandomize.html</anchorfile>
      <anchor>aa0475630d4c12c14825b9c0be00f6db7</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLRandomize&lt; complex&lt; Scalar &gt; &gt;</name>
    <filename>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObject&lt; complex&lt; Scalar &gt; &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>aee3a5c322bcf1967cdd40be4e6229dc9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>aeeb8a61639ce8f05ca78f6a26d642567</anchor>
      <arglist>(long seed, Scalar _a=ScalarFieldTraits&lt; Scalar &gt;::Zero(), Scalar _b=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a94716829fa817bcda0e0841319e41478</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>aa2e40cf22c5eab1db0841ecbe0660854</anchor>
      <arglist>(LocalDataContainer&lt; complex&lt; Scalar &gt; &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01complex_3_01Scalar_01_4_01_4.html</anchorfile>
      <anchor>a7937e21f5aeb11e74913ff077ee34173</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLRandomize&lt; int &gt;</name>
    <filename>classRVL_1_1RVLRandomize_3_01int_01_4.html</filename>
    <base>UnaryLocalFunctionObject&lt; int &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01int_01_4.html</anchorfile>
      <anchor>ac0a98a6a58771e6fef96f016848d7762</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLRandomize</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01int_01_4.html</anchorfile>
      <anchor>a6390e6193d5e1ccabea5dcdcab996ee8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>readsData</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01int_01_4.html</anchorfile>
      <anchor>a8a032748cdcaa5fbe5101e84ed0be2ab</anchor>
      <arglist>(size_t i=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01int_01_4.html</anchorfile>
      <anchor>a03781d5d94e1f8db99f38e3e82d23724</anchor>
      <arglist>(LocalDataContainer&lt; int &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLRandomize_3_01int_01_4.html</anchorfile>
      <anchor>a5d2f40ed60019342f4a18e6d70b78659</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASCIIReader</name>
    <filename>classRVL_1_1ASCIIReader.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>ASCIIReader</name>
      <anchorfile>classRVL_1_1ASCIIReader.html</anchorfile>
      <anchor>a54dd54a266fc0ed3fe07c9452464f0a0</anchor>
      <arglist>(string fname)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASCIIReader</name>
      <anchorfile>classRVL_1_1ASCIIReader.html</anchorfile>
      <anchor>a8170c4c736bccdad5908e6e265ee4d23</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASCIIReader.html</anchorfile>
      <anchor>a962d5b55ff7f59f4758cc0287d9df2c6</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASCIIReader.html</anchorfile>
      <anchor>a2b56f507d5ecec2e38b70c6cf0dce847</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ASCIIWriter</name>
    <filename>classRVL_1_1ASCIIWriter.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObjectConstEval&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>ASCIIWriter</name>
      <anchorfile>classRVL_1_1ASCIIWriter.html</anchorfile>
      <anchor>aa641c277c63c9e8a47d0096f64837f24</anchor>
      <arglist>(string fname)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ASCIIWriter</name>
      <anchorfile>classRVL_1_1ASCIIWriter.html</anchorfile>
      <anchor>a5014bd54a9d3331026fbecf3cf2f762b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ASCIIWriter.html</anchorfile>
      <anchor>a7570b5641786473283a665d89ee88a27</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ASCIIWriter.html</anchorfile>
      <anchor>a7c34e40bb7340f45b8802f2931ddf48e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryReader</name>
    <filename>classRVL_1_1BinaryReader.html</filename>
    <templarg></templarg>
    <base>UnaryLocalFunctionObjectConstEval&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>BinaryReader</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>ae7944247d6fa6b46c36ea57f4996d66f</anchor>
      <arglist>(char const *fname, long _first=0L)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BinaryReader</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>aaca44a51b2798ca6c21f76efb1f96346</anchor>
      <arglist>(const string fname, long _first=0L)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BinaryReader</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>a6afd21e1d632f592a838a8aad085a161</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>seek</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>a3ff930b7b2a2e668b24a7f0fb891a5a2</anchor>
      <arglist>(size_t firstword)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>a87b4a77c2cb003936b97f277510f59cf</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1BinaryReader.html</anchorfile>
      <anchor>a3dd1205a436c17b7458f6667e957792f</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryWriter</name>
    <filename>classRVL_1_1BinaryWriter.html</filename>
    <templarg></templarg>
    <base>UnaryLocalConstEval&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>BinaryWriter</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>af9989120ea0c6f7b2857a70549d77d49</anchor>
      <arglist>(char const *fname, long _first=0L)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BinaryWriter</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>a919fc7a7a2d5ebac7b085af0e7ce0c51</anchor>
      <arglist>(const string fname, long _first=0L)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BinaryWriter</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>ab527b603786ffa048711975442f47eac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>seek</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>a01d3b6f31b340fa29cc3a3cfe35f05cd</anchor>
      <arglist>(size_t firstword)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>a84c9adbc7d34ab68bf736886387ab9f7</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1BinaryWriter.html</anchorfile>
      <anchor>ad51d70709b9b78c0ae363a6b2d04e78a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLBoxMaxStep</name>
    <filename>classRVL_1_1RVLBoxMaxStep.html</filename>
    <templarg>Scalar</templarg>
    <base>QuaternaryLocalFunctionObjectScalarRedn&lt; Scalar, Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLBoxMaxStep</name>
      <anchorfile>classRVL_1_1RVLBoxMaxStep.html</anchorfile>
      <anchor>a8f29d83319bd31d50f1f91e095832b47</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLBoxMaxStep</name>
      <anchorfile>classRVL_1_1RVLBoxMaxStep.html</anchorfile>
      <anchor>a0d76597cc4a6565c3f3fc611e8b93c73</anchor>
      <arglist>(const RVLBoxMaxStep&lt; Scalar &gt; &amp;b)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLBoxMaxStep</name>
      <anchorfile>classRVL_1_1RVLBoxMaxStep.html</anchorfile>
      <anchor>af1c272ac78445632929e27500b12946f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLBoxMaxStep.html</anchorfile>
      <anchor>ad3dc79cd9edf60ebfad53b0de2418006</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; const &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;dx, LocalDataContainer&lt; Scalar &gt; const &amp;xmin, LocalDataContainer&lt; Scalar &gt; const &amp;xmax)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLBoxMaxStep.html</anchorfile>
      <anchor>a671638688f92578c9f3e6b3d7d5a3610</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ElementwiseMultiply</name>
    <filename>classRVL_1_1ElementwiseMultiply.html</filename>
    <templarg>Scalar</templarg>
    <base>TernaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>ElementwiseMultiply</name>
      <anchorfile>classRVL_1_1ElementwiseMultiply.html</anchorfile>
      <anchor>af89e2c568f2877bcbbe1d431c218c396</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ElementwiseMultiply</name>
      <anchorfile>classRVL_1_1ElementwiseMultiply.html</anchorfile>
      <anchor>a681236d09acb8807ac494c2b53b25763</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ElementwiseMultiply.html</anchorfile>
      <anchor>a03ad58725cd1a7584f8cefd67cc42f0f</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;u, LocalDataContainer&lt; Scalar &gt; const &amp;v, LocalDataContainer&lt; Scalar &gt; const &amp;w)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ElementwiseMultiply.html</anchorfile>
      <anchor>a4079db21525727f16e0b97e3d1e405f4</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ElementwiseDivision</name>
    <filename>classRVL_1_1ElementwiseDivision.html</filename>
    <templarg>Scalar</templarg>
    <base>TernaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>ElementwiseDivision</name>
      <anchorfile>classRVL_1_1ElementwiseDivision.html</anchorfile>
      <anchor>a7ae11616508981d5372a8c042ca76dd1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ElementwiseDivision</name>
      <anchorfile>classRVL_1_1ElementwiseDivision.html</anchorfile>
      <anchor>af83f38592f13d8848b2395be0fddce15</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ElementwiseDivision.html</anchorfile>
      <anchor>a6986ef1204e6a58c04c1aba5962afe27</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;u, LocalDataContainer&lt; Scalar &gt; const &amp;v, LocalDataContainer&lt; Scalar &gt; const &amp;w)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ElementwiseDivision.html</anchorfile>
      <anchor>ad4613982b2b2a61f7e1e4673cd07c65e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ElementwiseSqrtAbs</name>
    <filename>classRVL_1_1ElementwiseSqrtAbs.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>ElementwiseSqrtAbs</name>
      <anchorfile>classRVL_1_1ElementwiseSqrtAbs.html</anchorfile>
      <anchor>a7af958896d528852f27e63c43890b367</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ElementwiseSqrtAbs</name>
      <anchorfile>classRVL_1_1ElementwiseSqrtAbs.html</anchorfile>
      <anchor>a21d5e88c9bcc2734f88991caab26882f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1ElementwiseSqrtAbs.html</anchorfile>
      <anchor>a2ecedf5bfd68a7acd0a92f034176903f</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;u, LocalDataContainer&lt; Scalar &gt; const &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1ElementwiseSqrtAbs.html</anchorfile>
      <anchor>a3df4ca15a7bf85bb57083c5355a60ede</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLScalarLogistic</name>
    <filename>classRVL_1_1RVLScalarLogistic.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLScalarLogistic</name>
      <anchorfile>classRVL_1_1RVLScalarLogistic.html</anchorfile>
      <anchor>a95b3a83a124f172ec2428eb965a9c56a</anchor>
      <arglist>(Scalar fmin=ScalarFieldTraits&lt; Scalar &gt;::Zero(), Scalar fmax=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLScalarLogistic</name>
      <anchorfile>classRVL_1_1RVLScalarLogistic.html</anchorfile>
      <anchor>aeb880d4a91a49b35f4eb061f312631f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLScalarLogistic.html</anchorfile>
      <anchor>a86b7c9c753a5818b56c93a33cad65261</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLScalarLogistic.html</anchorfile>
      <anchor>aabb5ca9a5fe2138047eb99fbe3e00cc4</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLScalarLogisticInverse</name>
    <filename>classRVL_1_1RVLScalarLogisticInverse.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLScalarLogisticInverse</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticInverse.html</anchorfile>
      <anchor>a7918e615794b0ecb77e6e346ec3256dc</anchor>
      <arglist>(Scalar fmin=ScalarFieldTraits&lt; Scalar &gt;::Zero(), Scalar fmax=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLScalarLogisticInverse</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticInverse.html</anchorfile>
      <anchor>a677226a58df8d5d63c328a83f53957f7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticInverse.html</anchorfile>
      <anchor>a28fab136552326eb4ba21733ffdd051e</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;y)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticInverse.html</anchorfile>
      <anchor>ad5e46ce79709785455bf86c99a86bab6</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLScalarLogisticDeriv</name>
    <filename>classRVL_1_1RVLScalarLogisticDeriv.html</filename>
    <templarg></templarg>
    <base>TernaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLScalarLogisticDeriv</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticDeriv.html</anchorfile>
      <anchor>ab0195ef02c3e42d158ae94903e3960ee</anchor>
      <arglist>(Scalar fmin=ScalarFieldTraits&lt; Scalar &gt;::Zero(), Scalar fmax=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLScalarLogisticDeriv</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticDeriv.html</anchorfile>
      <anchor>a20d911bf1f5dd7ccf630a3b01277a893</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticDeriv.html</anchorfile>
      <anchor>a216ae9a1f78763688c5d1d9f6175eddb</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;y, LocalDataContainer&lt; Scalar &gt; const &amp;dy)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLScalarLogisticDeriv.html</anchorfile>
      <anchor>a56605a85ba365289546c1ae6247fbc6a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLVectorLogistic</name>
    <filename>classRVL_1_1RVLVectorLogistic.html</filename>
    <templarg></templarg>
    <base>QuaternaryLocalFunctionObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLVectorLogistic</name>
      <anchorfile>classRVL_1_1RVLVectorLogistic.html</anchorfile>
      <anchor>a754ea6e08bb58342bebf57854a699836</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLVectorLogistic</name>
      <anchorfile>classRVL_1_1RVLVectorLogistic.html</anchorfile>
      <anchor>aeb410d654a7b6ed9167e69ed26304d19</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLVectorLogistic.html</anchorfile>
      <anchor>aced03127ebf88a91a782ef30cac6a40a</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;x, LocalDataContainer&lt; Scalar &gt; const &amp;y, LocalDataContainer&lt; Scalar &gt; const &amp;lb, LocalDataContainer&lt; Scalar &gt; const &amp;ub)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLVectorLogistic.html</anchorfile>
      <anchor>a5e97aa4cecc373ce959b7006420a5dd0</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalDataContainer</name>
    <filename>classRVL_1_1LocalDataContainer.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::DataContainer</base>
    <member kind="function">
      <type></type>
      <name>LocalDataContainer</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>a50007ccc38be7a55791ef082d60ed485</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalDataContainer</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>a0e9230d60ee37df9f033a36722aafa45</anchor>
      <arglist>(const LocalDataContainer&lt; DataType &gt; &amp;D)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalDataContainer</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>add1dde143495df1a25bcf4d025dedb50</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>a016914da0c8ca944d662a317d1acaae6</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual DataType *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>abc9d2536f05ebb3aaf40ca2efa27c3da</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual DataType const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>aac9e044916c20784a54e884a2b9f41f5</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>a5234408b6a1fc3b77dc8fa2efa63d4be</anchor>
      <arglist>(FunctionObject &amp;f, vector&lt; DataContainer const * &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1LocalDataContainer.html</anchorfile>
      <anchor>a750a729af1fd4523a2a43d9a100fdd01</anchor>
      <arglist>(FunctionObjectConstEval &amp;f, vector&lt; DataContainer const * &gt; &amp;x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalDataContainerFactory</name>
    <filename>classRVL_1_1LocalDataContainerFactory.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::DataContainerFactory</base>
    <member kind="function">
      <type></type>
      <name>LocalDataContainerFactory</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>a408c510edf671fc6290e1351f84fb4fd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalDataContainerFactory</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>afd34200192a23658d87f3f0028a593b9</anchor>
      <arglist>(LocalDataContainerFactory&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalDataContainerFactory</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>aee2ad4176ece8fd3275f22b88f08340a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LocalDataContainer&lt; DataType &gt; *</type>
      <name>buildLocal</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>aca9e6e9b64fde36261f6679860b0dd3e</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>build</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>a31a86c432aeb62e35a40f16d8b7095cb</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1LocalDataContainerFactory.html</anchorfile>
      <anchor>a5dd8e9469e7dcd2f8b52a7d680e3e720</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalDataContainerSection</name>
    <filename>classRVL_1_1LocalDataContainerSection.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalDataContainer</base>
    <member kind="function">
      <type></type>
      <name>LocalDataContainerSection</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>aebbb4acef92665484ba57e55ec7b83af</anchor>
      <arglist>(const LocalDataContainerSection&lt; DataType &gt; &amp;D)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalDataContainerSection</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>a26f4fa0a3cff08cdf9085ee2e72a0fa9</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;_src, size_t _begin, size_t _length)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalDataContainerSection</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>aff820119d918cbb919a4f1e048a194a4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>a9719d2b553c18d9e500a2e5ee525e31b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataType *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>a70c7623ecd1c16784fdeb84280d85fd7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataType const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>a2bfd273143356a8878502c532297257f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LocalDataContainerSection.html</anchorfile>
      <anchor>a43fe9f35103027956c222402f0e36d7e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalEvaluation</name>
    <filename>classRVL_1_1LocalEvaluation.html</filename>
    <templarg>DataType</templarg>
    <member kind="function">
      <type></type>
      <name>LocalEvaluation</name>
      <anchorfile>classRVL_1_1LocalEvaluation.html</anchorfile>
      <anchor>ae31a08b1d6d0fbcbbd21063fee609961</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalEvaluation</name>
      <anchorfile>classRVL_1_1LocalEvaluation.html</anchorfile>
      <anchor>a9da81bcba25d63116edc689c0e0badf8</anchor>
      <arglist>(const LocalEvaluation&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalEvaluation</name>
      <anchorfile>classRVL_1_1LocalEvaluation.html</anchorfile>
      <anchor>a8992f0e64d06b8fbd7c75bba7c7b271f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1LocalEvaluation.html</anchorfile>
      <anchor>a840ca3c7daed94d66092230c1b1e28eb</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalFunctionObject</name>
    <filename>classRVL_1_1LocalFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObject</base>
    <base>RVL::LocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>LocalFunctionObject</name>
      <anchorfile>classRVL_1_1LocalFunctionObject.html</anchorfile>
      <anchor>aaa565e6f3cc2ddfdc16e3ee230271583</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalFunctionObject</name>
      <anchorfile>classRVL_1_1LocalFunctionObject.html</anchorfile>
      <anchor>a9188763bc1c3051428de144b1de169f7</anchor>
      <arglist>(const LocalFunctionObject&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalFunctionObject</name>
      <anchorfile>classRVL_1_1LocalFunctionObject.html</anchorfile>
      <anchor>aed6013bd64ade6a3d6acb99a110fe593</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::UnaryLocalEvaluation</name>
    <filename>classRVL_1_1UnaryLocalEvaluation.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>UnaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1UnaryLocalEvaluation.html</anchorfile>
      <anchor>a6b6ee22a221ad59d1d5e931c8b8fa00c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UnaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1UnaryLocalEvaluation.html</anchorfile>
      <anchor>a9900a8605b9fa07659e49136175ff65b</anchor>
      <arglist>(const UnaryLocalEvaluation&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~UnaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1UnaryLocalEvaluation.html</anchorfile>
      <anchor>ada36ca1036fd73a6b66d87cf7752ae8f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1UnaryLocalEvaluation.html</anchorfile>
      <anchor>a89f2f24b0691fc2b54da73bd6bd58186</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1UnaryLocalEvaluation.html</anchorfile>
      <anchor>a16a90992e5d0a78a8abc5991adbb3578</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::UnaryLocalFunctionObject</name>
    <filename>classRVL_1_1UnaryLocalFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObject</base>
    <base>RVL::UnaryLocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>UnaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObject.html</anchorfile>
      <anchor>a269a84a6a0c86d41379c93e1f481bfac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UnaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObject.html</anchorfile>
      <anchor>acefba77af909e098c33885e1299cd235</anchor>
      <arglist>(const UnaryLocalFunctionObject&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~UnaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObject.html</anchorfile>
      <anchor>a38d8e2814374b2817d5424ead24e2c63</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryLocalEvaluation</name>
    <filename>classRVL_1_1BinaryLocalEvaluation.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>BinaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1BinaryLocalEvaluation.html</anchorfile>
      <anchor>a809cb970b1b796d8ea5da801889f44f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BinaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1BinaryLocalEvaluation.html</anchorfile>
      <anchor>a5bc2cec4c291725acfa5e4e4fe402ae1</anchor>
      <arglist>(const BinaryLocalEvaluation&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BinaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1BinaryLocalEvaluation.html</anchorfile>
      <anchor>a2e514625302241e5854ff79179c0885f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryLocalEvaluation.html</anchorfile>
      <anchor>a33353217c50da621ca2790ce14d97d7c</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, LocalDataContainer&lt; DataType &gt; const &amp;source)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryLocalEvaluation.html</anchorfile>
      <anchor>aa9e560eb51eca672ed886fa00f1a1ef5</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryLocalFunctionObject</name>
    <filename>classRVL_1_1BinaryLocalFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObject</base>
    <base>RVL::BinaryLocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>BinaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObject.html</anchorfile>
      <anchor>a3de975cecaafa78392575e7b3ca057dc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BinaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObject.html</anchorfile>
      <anchor>afca47f0edf752dcd9fff23d26cc95d47</anchor>
      <arglist>(const BinaryLocalFunctionObject&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BinaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObject.html</anchorfile>
      <anchor>a6d379465dc09822f88ae2dcc41c74063</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TernaryLocalEvaluation</name>
    <filename>classRVL_1_1TernaryLocalEvaluation.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>TernaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1TernaryLocalEvaluation.html</anchorfile>
      <anchor>a05681de885f7dd52ce4f3f72046a27f6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TernaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1TernaryLocalEvaluation.html</anchorfile>
      <anchor>a94c5f1735628d42252765063987d1dc3</anchor>
      <arglist>(const TernaryLocalEvaluation&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TernaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1TernaryLocalEvaluation.html</anchorfile>
      <anchor>a4bc0b7aa30f92e5d643c4faa41bfafe6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1TernaryLocalEvaluation.html</anchorfile>
      <anchor>a08bf135485f16cfd81eea7b2affcd3c4</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, LocalDataContainer&lt; DataType &gt; const &amp;source1, LocalDataContainer&lt; DataType &gt; const &amp;source2)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1TernaryLocalEvaluation.html</anchorfile>
      <anchor>a8f371c7fbcfb7232ad5c9dc6a9968dac</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TernaryLocalFunctionObject</name>
    <filename>classRVL_1_1TernaryLocalFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObject</base>
    <base>RVL::TernaryLocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>TernaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObject.html</anchorfile>
      <anchor>a8e3dd22003fbff442bfbd066115c510f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TernaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObject.html</anchorfile>
      <anchor>a20701e37b04f649c5f008ff68864fe06</anchor>
      <arglist>(const TernaryLocalFunctionObject&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TernaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObject.html</anchorfile>
      <anchor>a69c4f6869a69ef68c2545d2fa327e410</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuaternaryLocalEvaluation</name>
    <filename>classRVL_1_1QuaternaryLocalEvaluation.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1QuaternaryLocalEvaluation.html</anchorfile>
      <anchor>afe794a86b0acd851027382a7f52b07ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1QuaternaryLocalEvaluation.html</anchorfile>
      <anchor>a6a98a4b0e0b8211c2317fb46b08eaf58</anchor>
      <arglist>(const QuaternaryLocalEvaluation&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~QuaternaryLocalEvaluation</name>
      <anchorfile>classRVL_1_1QuaternaryLocalEvaluation.html</anchorfile>
      <anchor>a30c421a78f6db96f9ead92abf60db880</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1QuaternaryLocalEvaluation.html</anchorfile>
      <anchor>afb277ddc40bc62e68c16d9748d8a0ff1</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, LocalDataContainer&lt; DataType &gt; const &amp;source1, LocalDataContainer&lt; DataType &gt; const &amp;source2, LocalDataContainer&lt; DataType &gt; const &amp;source3)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1QuaternaryLocalEvaluation.html</anchorfile>
      <anchor>a9513223b52e567b372289e30172e3cb7</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;target, vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuaternaryLocalFunctionObject</name>
    <filename>classRVL_1_1QuaternaryLocalFunctionObject.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObject</base>
    <base>RVL::QuaternaryLocalEvaluation</base>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObject.html</anchorfile>
      <anchor>ad6b2486742d88b6f8243100c7edaadc4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObject.html</anchorfile>
      <anchor>ae0217d46150a35721a081826bb70206a</anchor>
      <arglist>(const QuaternaryLocalFunctionObject&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~QuaternaryLocalFunctionObject</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObject.html</anchorfile>
      <anchor>a5e8b4038ea51f64e00d4fdf2848505fa</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalLinearAlgebraPackage</name>
    <filename>classRVL_1_1LocalLinearAlgebraPackage.html</filename>
    <templarg>DataType</templarg>
    <templarg>Scalar</templarg>
    <base>RVL::LinearAlgebraPackage</base>
    <member kind="function">
      <type></type>
      <name>LocalLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>ab95c6086516d63321f55147633f185c1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>a03da312fbf0a3606f939374045c5857a</anchor>
      <arglist>(const LocalLinearAlgebraPackage&lt; Scalar, DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>a2f7450ff17dfd4ba07a0a6577f19456b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual BinaryLocalFunctionObjectScalarRedn&lt; DataType, Scalar &gt; &amp;</type>
      <name>localinner</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>ab964812a426d93169a2e5fd8e6e93a89</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>FunctionObjectScalarRedn&lt; Scalar &gt; &amp;</type>
      <name>inner</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>a1259a4fc2dbdd706f35a4371e0feacb6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual UnaryLocalFunctionObject&lt; DataType &gt; &amp;</type>
      <name>localzero</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>ac16263b88b665bf40b072834b650b4e4</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>FunctionObject &amp;</type>
      <name>zero</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>a47789f2d3417bfae9225cc694ad5b04b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>a2cb7d741cf365f79badc536fafdc2f4e</anchor>
      <arglist>(LinearAlgebraPackage&lt; Scalar &gt; const &amp;lap) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>ae2a513af6bb3646404723ab181b47203</anchor>
      <arglist>(RVLException &amp;e) const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1LocalLinearAlgebraPackage.html</anchorfile>
      <anchor>ad7786130cf00eaae761c8715428f858e</anchor>
      <arglist>(ostream &amp;str) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLLinCombObject</name>
    <filename>classRVL_1_1RVLLinCombObject.html</filename>
    <templarg>Scalar</templarg>
    <base>BinaryLocalEvaluation&lt; Scalar &gt;</base>
    <base>LinCombObject&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLLinCombObject</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>abb33170c08b92d02dbc63c00fa25c7b2</anchor>
      <arglist>(Scalar ain=ScalarFieldTraits&lt; Scalar &gt;::One(), Scalar bin=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLLinCombObject</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>a8562b206b1274a617a36985468454d71</anchor>
      <arglist>(const RVLLinCombObject&lt; Scalar &gt; &amp;lc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLLinCombObject</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>a106c77d9fa8313dfc1a039e2f805ec82</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScalar</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>ad5eaf18fc39f708ba683ffeb96d9b109</anchor>
      <arglist>(Scalar ain, Scalar bin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>a1829e074e59eb06a39204dfc418aac79</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;u, LocalDataContainer&lt; Scalar &gt; const &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1RVLLinCombObject.html</anchorfile>
      <anchor>a6ff7e25ab9b48328f1feb02629a42441</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RVLLinearAlgebraPackage</name>
    <filename>classRVL_1_1RVLLinearAlgebraPackage.html</filename>
    <templarg>Scalar</templarg>
    <base>LocalLinearAlgebraPackage&lt; Scalar, Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RVLLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a15d0dd2e5ec8d47a06d31a57060664aa</anchor>
      <arglist>(Scalar ipscale=ScalarFieldTraits&lt; Scalar &gt;::One())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RVLLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a4b8db137a13f465adf38e1865d29dbdd</anchor>
      <arglist>(const RVLLinearAlgebraPackage&lt; Scalar &gt; &amp;p)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RVLLinearAlgebraPackage</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a06679762ddff37f516b276f3e85c44a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BinaryLocalFunctionObjectScalarRedn&lt; Scalar, Scalar &gt; &amp;</type>
      <name>localinner</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a5c6bc284fdd2b1c2b1ea53e4816f3d69</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>UnaryLocalFunctionObject&lt; Scalar &gt; &amp;</type>
      <name>localzero</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a2c0db4f4e71f58b4bf6759309c1109b0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LinCombObject&lt; Scalar &gt; &amp;</type>
      <name>linComb</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a87803c2530b54418b36c1afe131fe856</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>ab92ad379cb9dfc943c5342a98c5cefdd</anchor>
      <arglist>(LinearAlgebraPackage&lt; Scalar &gt; const &amp;lap) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setScale</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a2bf23f7cc963e0d8604dbeb3c7955bca</anchor>
      <arglist>(Scalar newscale)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a8f93bd4498fcae0ed7df31967af38ff4</anchor>
      <arglist>(RVLException &amp;e) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RVLLinearAlgebraPackage.html</anchorfile>
      <anchor>a465bf0c7ef03870700749ae52215a136</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ProductLocalDataContainer</name>
    <filename>classRVL_1_1ProductLocalDataContainer.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalDataContainer</base>
    <member kind="function">
      <type></type>
      <name>ProductLocalDataContainer</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>aceafe9b0665cdee87ed27e1255906f7d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ProductLocalDataContainer</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>a3c338a2da1a065ce14d7c12b2480d4db</anchor>
      <arglist>(ProductLocalDataContainer&lt; DataType &gt; &amp;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~ProductLocalDataContainer</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>aadc0fe75d18b71e2c0d6e4d49e07f7e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>getNumberOfComponents</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>a8ab6e2269cabb01615af1e78b81d0e33</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LocalDataContainer&lt; DataType &gt; &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>a13f64190c5aacbe4b6f9137da907f5db</anchor>
      <arglist>(int i)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual LocalDataContainer&lt; DataType &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1ProductLocalDataContainer.html</anchorfile>
      <anchor>a7888859e01f6cf2830989e8b1c0da948</anchor>
      <arglist>(int i) const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::ProductDataContainerLDC</name>
    <filename>classRVL_1_1ProductDataContainerLDC.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::ProductDataContainer</base>
    <member kind="function">
      <type></type>
      <name>ProductDataContainerLDC</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a35da6b70fabf212248af02b601b81fd7</anchor>
      <arglist>(ProductLocalDataContainer&lt; DataType &gt; &amp;_pldc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ProductDataContainerLDC</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a9ebe0fb62aa6084e23e4c752738d3ae7</anchor>
      <arglist>(const ProductDataContainerLDC&lt; DataType &gt; &amp;pdldc)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ProductDataContainerLDC</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a9185bec3fdee5b2d323a47d8d418b9d4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a2c1647f5721dd85c5969b230e0ed6566</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual DataContainer &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a7b3b86e0f819a5b66fa1271d624567d2</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual DataContainer const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1ProductDataContainerLDC.html</anchorfile>
      <anchor>a48e68c2831064b23f3b9e14be789a2fb</anchor>
      <arglist>(int i) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PartitionedLocalDataContainer</name>
    <filename>classRVL_1_1PartitionedLocalDataContainer.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::ProductLocalDataContainer</base>
    <member kind="function">
      <type></type>
      <name>PartitionedLocalDataContainer</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a932e80c16b9759d87fd31337589bb0ad</anchor>
      <arglist>(const PartitionedLocalDataContainer&lt; DataType &gt; &amp;p)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PartitionedLocalDataContainer</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a86c35cc306c72184ebc26909c24f9dd7</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; &amp;_ldc, vector&lt; int &gt; _k)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PartitionedLocalDataContainer</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a4d67c3092af652d1202c2d0d78d81027</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a24dc2a50cf291f9e7d1fcad809fdd5fc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataType *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a5bcb68c6175f5a94b19e1988b4a17bc8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataType const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>ad8152e907bb6f8de065cf44a3a0ac1c1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNumberOfComponents</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>ae8a06836156d4242adf4bb09129c5da2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>LocalDataContainer&lt; DataType &gt; &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>aad00cfb69d96a407afdd03e9a19636e1</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>LocalDataContainer&lt; DataType &gt; const &amp;</type>
      <name>operator[]</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>aa93947472984770cfa5cd6ddadf288c3</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1PartitionedLocalDataContainer.html</anchorfile>
      <anchor>a161657c6d552baa30d54e04cddf02ddd</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalConstEval</name>
    <filename>classRVL_1_1LocalConstEval.html</filename>
    <templarg>DataType</templarg>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalConstEval</name>
      <anchorfile>classRVL_1_1LocalConstEval.html</anchorfile>
      <anchor>abea79a029d25010706f6abea1960bd71</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1LocalConstEval.html</anchorfile>
      <anchor>a9aba99ad1e24e67d5d8c54bbcd5b7181</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::UnaryLocalConstEval</name>
    <filename>classRVL_1_1UnaryLocalConstEval.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~UnaryLocalConstEval</name>
      <anchorfile>classRVL_1_1UnaryLocalConstEval.html</anchorfile>
      <anchor>a2ca7945bd8ecbe21b8e4c4133e75d3b7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1UnaryLocalConstEval.html</anchorfile>
      <anchor>a17a3e03ff0a95edd2c798e3275050c70</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; const &amp;source)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1UnaryLocalConstEval.html</anchorfile>
      <anchor>a37fe5e89551ac73aeaa200aa554a66d1</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryLocalConstEval</name>
    <filename>classRVL_1_1BinaryLocalConstEval.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BinaryLocalConstEval</name>
      <anchorfile>classRVL_1_1BinaryLocalConstEval.html</anchorfile>
      <anchor>a733d56bef12b0a91444f1a9a5b94211f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryLocalConstEval.html</anchorfile>
      <anchor>a0f376c0cf1c6f2d9e3b4f1159038f8f9</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; const &amp;source1, LocalDataContainer&lt; DataType &gt; const &amp;source2)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1BinaryLocalConstEval.html</anchorfile>
      <anchor>ae9a9058a4fe8605a86d5245e22d3c9ff</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TernaryLocalConstEval</name>
    <filename>classRVL_1_1TernaryLocalConstEval.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::LocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TernaryLocalConstEval</name>
      <anchorfile>classRVL_1_1TernaryLocalConstEval.html</anchorfile>
      <anchor>a346802449d7b79a90b38e8e4caf707ce</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1TernaryLocalConstEval.html</anchorfile>
      <anchor>aa3fd6c7d2180a3f78760809ef3616957</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; const &amp;source1, LocalDataContainer&lt; DataType &gt; const &amp;source2, LocalDataContainer&lt; DataType &gt; const &amp;source3)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1TernaryLocalConstEval.html</anchorfile>
      <anchor>ac500e8cc47def380bdf7886b488adc05</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuaternaryLocalConstEval</name>
    <filename>classRVL_1_1QuaternaryLocalConstEval.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::LocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~QuaternaryLocalConstEval</name>
      <anchorfile>classRVL_1_1QuaternaryLocalConstEval.html</anchorfile>
      <anchor>a301c7573ec75bbf4784c813a83d5fea6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1QuaternaryLocalConstEval.html</anchorfile>
      <anchor>a74ae2d79bc6f7d9d458ee2d18e11acef</anchor>
      <arglist>(LocalDataContainer&lt; DataType &gt; const &amp;source1, LocalDataContainer&lt; DataType &gt; const &amp;source2, LocalDataContainer&lt; DataType &gt; const &amp;source3, LocalDataContainer&lt; DataType &gt; const &amp;source4)=0</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1QuaternaryLocalConstEval.html</anchorfile>
      <anchor>a2f376d33f6b474ba776efa0b1b2f62b0</anchor>
      <arglist>(vector&lt; LocalDataContainer&lt; DataType &gt; const * &gt; &amp;sources)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::UnaryLocalFunctionObjectConstEval</name>
    <filename>classRVL_1_1UnaryLocalFunctionObjectConstEval.html</filename>
    <templarg>DataType</templarg>
    <base>RVL::FunctionObjectConstEval</base>
    <base>RVL::UnaryLocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~UnaryLocalFunctionObjectConstEval</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObjectConstEval.html</anchorfile>
      <anchor>aec04ce81e488127339b8e355ab4f3c8b</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryLocalFunctionObjectConstEval</name>
    <filename>classRVL_1_1BinaryLocalFunctionObjectConstEval.html</filename>
    <templarg></templarg>
    <base>RVL::FunctionObjectConstEval</base>
    <base>RVL::BinaryLocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BinaryLocalFunctionObjectConstEval</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObjectConstEval.html</anchorfile>
      <anchor>a1564c23bf2333afb506aaddd413e9652</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TernaryLocalFunctionObjectConstEval</name>
    <filename>classRVL_1_1TernaryLocalFunctionObjectConstEval.html</filename>
    <templarg></templarg>
    <base>RVL::FunctionObjectConstEval</base>
    <base>TernaryLocalConstEval&lt; DataType &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TernaryLocalFunctionObjectConstEval</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObjectConstEval.html</anchorfile>
      <anchor>a2b36ecd9c7d718271625604f65d61190</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuaternaryLocalFunctionObjectConstEval</name>
    <filename>classRVL_1_1QuaternaryLocalFunctionObjectConstEval.html</filename>
    <templarg></templarg>
    <base>RVL::FunctionObjectConstEval</base>
    <base>RVL::QuaternaryLocalConstEval</base>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~QuaternaryLocalFunctionObjectConstEval</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObjectConstEval.html</anchorfile>
      <anchor>a10e94488778198378b3f4f095c1318a2</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::UnaryLocalFunctionObjectScalarRedn</name>
    <filename>classRVL_1_1UnaryLocalFunctionObjectScalarRedn.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectScalarRedn</base>
    <base>RVL::UnaryLocalConstEval</base>
    <member kind="function">
      <type></type>
      <name>UnaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a58e6497153a1e95c07b9d38ccffc97da</anchor>
      <arglist>(ValType val)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UnaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a6326dc0f1a5930dfac84d731b8f38cb8</anchor>
      <arglist>(UnaryLocalFunctionObjectScalarRedn&lt; DataType, ValType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~UnaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1UnaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a66b786098b201fdc7cf98b6e44a1f962</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::BinaryLocalFunctionObjectScalarRedn</name>
    <filename>classRVL_1_1BinaryLocalFunctionObjectScalarRedn.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectScalarRedn</base>
    <base>RVL::BinaryLocalConstEval</base>
    <member kind="function">
      <type></type>
      <name>BinaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>ae6dc4db514d2954440dd244fee076283</anchor>
      <arglist>(ValType val)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BinaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a9d8264c4a5497c2210428f99b542b178</anchor>
      <arglist>(BinaryLocalFunctionObjectScalarRedn&lt; DataType, ValType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BinaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1BinaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a981d1c5b8f25cac54e6d773333e40928</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::TernaryLocalFunctionObjectScalarRedn</name>
    <filename>classRVL_1_1TernaryLocalFunctionObjectScalarRedn.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectScalarRedn</base>
    <base>TernaryLocalConstEval&lt; DataType &gt;</base>
    <member kind="function">
      <type></type>
      <name>TernaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>aca493f4682ab04d99f1fac646c87ea71</anchor>
      <arglist>(ValType val)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TernaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a42b9842666b82f3c3ce0536334c4605e</anchor>
      <arglist>(TernaryLocalFunctionObjectScalarRedn&lt; DataType, ValType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TernaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1TernaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a04d2c6bb4da6ee851fc91d4dcd7669fb</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::QuaternaryLocalFunctionObjectScalarRedn</name>
    <filename>classRVL_1_1QuaternaryLocalFunctionObjectScalarRedn.html</filename>
    <templarg>DataType</templarg>
    <templarg>ValType</templarg>
    <base>RVL::FunctionObjectScalarRedn</base>
    <base>RVL::QuaternaryLocalConstEval</base>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a755df12d4aed565566e3f809690841fe</anchor>
      <arglist>(ValType val)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>QuaternaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>aa72b864c336250d3612e59a72b5dade8</anchor>
      <arglist>(QuaternaryLocalFunctionObjectScalarRedn&lt; DataType, ValType &gt; const &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~QuaternaryLocalFunctionObjectScalarRedn</name>
      <anchorfile>classRVL_1_1QuaternaryLocalFunctionObjectScalarRedn.html</anchorfile>
      <anchor>a2df2ae7f0101fc756afc0eb2f4fa96f4</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalSpace</name>
    <filename>classRVL_1_1LocalSpace.html</filename>
    <templarg>Scalar</templarg>
    <templarg>DataType</templarg>
    <base>RVL::StdSpace</base>
    <member kind="function">
      <type></type>
      <name>LocalSpace</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>a4e19b83445c3f4cf9a1fe03b0e320d92</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalSpace</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>ae7aa32a3e64b12d972beb547f99acff1</anchor>
      <arglist>(const LocalSpace&lt; Scalar, DataType &gt; &amp;sp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~LocalSpace</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>a0ac43084a4661ccb3e6ff569b5756d00</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>LocalDataContainer&lt; DataType &gt; *</type>
      <name>buildLocalDataContainer</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>ab98a162ed5d8cfc00591e1816469a628</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataContainer *</type>
      <name>buildDataContainer</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>adb25114090e63e8bc9b8ee9fcce70f7e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>a8ce81ea5e5f0505b911c6d4413704dd4</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual LocalDataContainerFactory&lt; DataType &gt; &amp;</type>
      <name>getLDCF</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>a5ba6939bfd9e3155955f00311dce9559</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>DataContainerFactory &amp;</type>
      <name>getDCF</name>
      <anchorfile>classRVL_1_1LocalSpace.html</anchorfile>
      <anchor>aebc1072935876b1867390d2e2dd6acb4</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::LocalVector</name>
    <filename>classRVL_1_1LocalVector.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>RVL::Vector</base>
    <member kind="function">
      <type></type>
      <name>LocalVector</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>a886d2ec20c6d29345f1b5d727e5826ab</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;v)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>LocalVector</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>a4daaadf326bc16f0be71c60c8de6e478</anchor>
      <arglist>(const Space&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LocalVector</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>a5821f8dfabf59164b2d56cc6a1224913</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>ad5b47e3cf65c8c12ec454835e0354468</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DataType *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>aa8413f6b184b903313227fae328f87b8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>DataType const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1LocalVector.html</anchorfile>
      <anchor>a53646004ccaaf909977866c18fcda67e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PackageContainerSpace</name>
    <filename>classRVL_1_1PackageContainerSpace.html</filename>
    <templarg></templarg>
    <templarg></templarg>
    <base>StdSpace&lt; DataType &gt;</base>
    <member kind="function">
      <type></type>
      <name>PackageContainerSpace</name>
      <anchorfile>classRVL_1_1PackageContainerSpace.html</anchorfile>
      <anchor>a56f7c96f964cc1f434a95e2af611e811</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PackageContainerSpace</name>
      <anchorfile>classRVL_1_1PackageContainerSpace.html</anchorfile>
      <anchor>ae4f31ff68b099e237c33bc1bb2a5b5f0</anchor>
      <arglist>(const PackageContainerSpace &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~PackageContainerSpace</name>
      <anchorfile>classRVL_1_1PackageContainerSpace.html</anchorfile>
      <anchor>a34ef772223e3d4d132705c5acf246a08</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual PackageContainerFactory&lt; DataType, MetaType &gt; const &amp;</type>
      <name>getPCF</name>
      <anchorfile>classRVL_1_1PackageContainerSpace.html</anchorfile>
      <anchor>a5743dba60604d4d4c0a62ced0cf2e008</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function">
      <type>DataContainerFactory const &amp;</type>
      <name>getDCF</name>
      <anchorfile>classRVL_1_1PackageContainerSpace.html</anchorfile>
      <anchor>a7255cda462d703af4216f7e41dd0084c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::PolynomialOperator</name>
    <filename>classRVL_1_1PolynomialOperator.html</filename>
    <templarg>Scalar</templarg>
    <base>RVL::OperatorWithInvertibleDeriv</base>
    <member kind="function">
      <type></type>
      <name>PolynomialOperator</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a06cd78e1d3ccf603a9c9c5a1d42b1e41</anchor>
      <arglist>(const std::valarray&lt; Scalar &gt; &amp;_coef, Space&lt; Scalar &gt; &amp;_spc)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>PolynomialOperator</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a41569e7d72b02ecf07b20db859d3d108</anchor>
      <arglist>(const PolynomialOperator&lt; Scalar &gt; &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PolynomialOperator</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a08a201eaf6ca979f7d6d60be44af6c0c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a59f973deea7b5cb9771f4069dcb8449b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Space&lt; Scalar &gt; &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a6fd6f7bc31f9761e46d120587510ae81</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>af581795b7ca6832c3bebf426e0de3fad</anchor>
      <arglist>(RVLException &amp;e) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a47e52861f45383fae72a51acffd8c23c</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a197d28489d5f22b36610c12909628614</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, Vector&lt; Scalar &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>ae4f9d5e3a17e7ced60d9ea72317933c3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a5f17cb85379e61f2b91728939e0fc023</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyInverseDeriv</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a2c1d3bbc8294f20ec30467e2e3a8859f</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dy, Vector&lt; Scalar &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdjInverseDeriv</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a35d4750f5249167cb3cb37a69886f1f3</anchor>
      <arglist>(const Vector&lt; Scalar &gt; &amp;x, const Vector&lt; Scalar &gt; &amp;dx, Vector&lt; Scalar &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Operator&lt; Scalar &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a8a972977902fe12ee07755a927de7a57</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Space&lt; Scalar &gt; &amp;</type>
      <name>spc</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>added5b7dc2d53d814b82dbb15ded4e74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::valarray&lt; Scalar &gt;</type>
      <name>coef</name>
      <anchorfile>classRVL_1_1PolynomialOperator.html</anchorfile>
      <anchor>a56ca812aed7168be0855cf15fdf343cd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RnArray</name>
    <filename>classRVL_1_1RnArray.html</filename>
    <templarg>Scalar</templarg>
    <base>LocalDataContainer&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RnArray</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a499be6d3f5b681d3abfe29d8c1b3be80</anchor>
      <arglist>(const RnArray&lt; Scalar &gt; &amp;x)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RnArray</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a03adf1c0bdc31fd2503c05883a802604</anchor>
      <arglist>(size_t _n)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RnArray</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>ae337ac800f3d9421d5809c14c9a553d2</anchor>
      <arglist>(LocalDataContainer&lt; Scalar &gt; &amp;rn, size_t len, size_t _start=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RnArray</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a910178223354679311168f9c8504785f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a6ae5d48f2b941f2cf7c65fc4b04fcc5c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Scalar *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>afbaffcc8405345ef0489f49cbe286f1e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Scalar const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a438c25ae88feacd055a4db46da4892e8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a7c2a4734c65bd6d39a1fb78f0d5927c9</anchor>
      <arglist>(RVLException &amp;str) const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnArray.html</anchorfile>
      <anchor>a043212e44bab5e5223662e7566243cc4</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RnDataContainerFactory</name>
    <filename>classRVL_1_1RnDataContainerFactory.html</filename>
    <templarg>Scalar</templarg>
    <base>LocalDataContainerFactory&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RnDataContainerFactory</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>a020adf21838ee2e03bdbb8c58b7ed2b7</anchor>
      <arglist>(size_t _n=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RnDataContainerFactory</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>aa19ec41d3db5329e2f312e2f0e6d62a1</anchor>
      <arglist>(const RnDataContainerFactory&lt; Scalar &gt; &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~RnDataContainerFactory</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>a16ef41d0869f05443902036d0a972862</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual LocalDataContainer&lt; Scalar &gt; *</type>
      <name>buildLocal</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>acde1219c2365c5608c7ca3d2853d27df</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>compare</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>aaa2f9d7d4c61ac65a4aa7a479e7aa56e</anchor>
      <arglist>(DataContainerFactory const &amp;dcf) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>a8ef563e64e627790ca361d11d10e9268</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function">
      <type>size_t</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>ab2cab7be33b827f5a475c3eefbcb04d8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>a53fe415b1d9e01cea0f24f686f3fa595</anchor>
      <arglist>(RVLException &amp;e) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>aeac1caca4f546bfecf413d112f2336be</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>redim</name>
      <anchorfile>classRVL_1_1RnDataContainerFactory.html</anchorfile>
      <anchor>a5866482ce7f6372f7c930250100b0de1</anchor>
      <arglist>(size_t _n)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::matvec</name>
    <filename>classRVL_1_1matvec.html</filename>
    <templarg></templarg>
    <base>BinaryLocalFunctionObject&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>matvec</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>ab1eb85ea8720aeb679f03a3bb1f5fe3b</anchor>
      <arglist>(int _rows, int _cols)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>matvec</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a5bbe437174db65466369d804d17c51e2</anchor>
      <arglist>(matvec&lt; T &gt; const *m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~matvec</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a4588bcdf2f762bdaf9f1f62182dc8c81</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a9774b31ee6be5e2a96fda6abb359a1a2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a5d2e33910ef74ae42d9acc2914c76fea</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a0690f354d079b16b806fb1ece7557347</anchor>
      <arglist>(UnaryLocalFunctionObject&lt; T &gt; &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>getElement</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a1ac6eefd57180a2ddb26c05813b6b518</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>T const &amp;</type>
      <name>getElement</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>ac2573151c24b6a744e7f0fcbebe92b5d</anchor>
      <arglist>(int i, int j) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setAdj</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>abe7c305d5a1a6ad8f1a3708e3214e975</anchor>
      <arglist>(bool flag)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>getAdj</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>aaaf1afd5aae7325074dbb1848c330ce7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>a7803ddefb1ac495138e8447e3771772f</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;y, LocalDataContainer&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1matvec.html</anchorfile>
      <anchor>af7d989d31990f9156dcf85bdc3b78d92</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::fmatvec</name>
    <filename>classRVL_1_1fmatvec.html</filename>
    <templarg>T</templarg>
    <base>BinaryLocalFunctionObject&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>fmatvec</name>
      <anchorfile>classRVL_1_1fmatvec.html</anchorfile>
      <anchor>a8d8a1e44f0b71c72cb2d272f9de89db0</anchor>
      <arglist>(matvec&lt; T &gt; &amp;_m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>fmatvec</name>
      <anchorfile>classRVL_1_1fmatvec.html</anchorfile>
      <anchor>ae7dd557d0d2d8566c67ef52da7eeb22e</anchor>
      <arglist>(fmatvec&lt; T &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~fmatvec</name>
      <anchorfile>classRVL_1_1fmatvec.html</anchorfile>
      <anchor>aa381b017798e03989347474ea9682e2d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1fmatvec.html</anchorfile>
      <anchor>a624b7dad291bf5c1cdac9ab0e10f7f55</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;y, LocalDataContainer&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1fmatvec.html</anchorfile>
      <anchor>a0309626d846388ff28b3947f2819bd7e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::amatvec</name>
    <filename>classRVL_1_1amatvec.html</filename>
    <templarg>T</templarg>
    <base>BinaryLocalFunctionObject&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>amatvec</name>
      <anchorfile>classRVL_1_1amatvec.html</anchorfile>
      <anchor>a7b638a1d7bc7f8d8b31e9fa657b696e0</anchor>
      <arglist>(matvec&lt; T &gt; &amp;_m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>amatvec</name>
      <anchorfile>classRVL_1_1amatvec.html</anchorfile>
      <anchor>ab956e90440a712108abe1589fd1fe625</anchor>
      <arglist>(amatvec&lt; T &gt; const &amp;f)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~amatvec</name>
      <anchorfile>classRVL_1_1amatvec.html</anchorfile>
      <anchor>a8be0e46eb1b5af469018bec24f01831a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1amatvec.html</anchorfile>
      <anchor>ab8b0dd518e889e31c9e8246667c820b9</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;y, LocalDataContainer&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1amatvec.html</anchorfile>
      <anchor>a5c94862953945dc490f11ac144b47fb8</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GenMat</name>
    <filename>classRVL_1_1GenMat.html</filename>
    <templarg>T</templarg>
    <base>LinearOp&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>GenMat</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a5c823ee9f00ba8a8552fa3120cacd54d</anchor>
      <arglist>(RnSpace&lt; T &gt; const &amp;_dom, RnSpace&lt; T &gt; const &amp;_rng)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GenMat</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a9ad3304598ffbef2233274288591e796</anchor>
      <arglist>(GenMat&lt; T &gt; const &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GenMat</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>abc8dd1218c70d1c92c399799223c6eff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a70f329038ee4009cf40f1827ffa91669</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>abbd0f2cf8838d005424a0320b519cd71</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNRows</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>afbb2ed7c101b4936e5ac319ccf5ff7aa</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getNCols</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a4c3e573f8663d23ae9ca209ec5b78b9e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a26a5a385465aaeffb1faff5558f54538</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T const *</type>
      <name>getData</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>aa8d0cba183aa67b1f6e43f4351f37943</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>T const &amp;</type>
      <name>getElement</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>abb828b750e47afa94ac5b073d5f7c715</anchor>
      <arglist>(int i, int j) const </arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>getElement</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>adac3722ba054ce4505a517019449699c</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>setElement</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>aa0fbfadee0f9f52bd59903897d6c8da6</anchor>
      <arglist>(int i, int j, T e)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a8d868603f4c3f712a3945c3481971730</anchor>
      <arglist>(UnaryLocalFunctionObject&lt; T &gt; &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>ab41173e8e5f6f18b6618ecc87240362f</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual LinearOp&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a7e52da125675d118276478c3b477faca</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a044faef7137c2a8fe1adec7d9b99d447</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>applyAdj</name>
      <anchorfile>classRVL_1_1GenMat.html</anchorfile>
      <anchor>a238b85431dcd54e78e8d0a3266941172</anchor>
      <arglist>(Vector&lt; T &gt; const &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::SymMat</name>
    <filename>classRVL_1_1SymMat.html</filename>
    <templarg></templarg>
    <base>RVL::GenMat</base>
    <member kind="function">
      <type></type>
      <name>SymMat</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>a06f54c90ea96557de65d390ee6628b38</anchor>
      <arglist>(RnSpace&lt; T &gt; const &amp;dom)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>SymMat</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>aa4720fbc79144deb522e5a2562951810</anchor>
      <arglist>(SymMat&lt; T &gt; const &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SymMat</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>ae4b25c45cb3585ea1c7086cf9f1e22f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>setElement</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>a6e3f43fa216fd6dd4233c0510fd7f221</anchor>
      <arglist>(int i, int j, T e)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>aa40569d937b1b4c73fafc10d2f3f540e</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearOp&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1SymMat.html</anchorfile>
      <anchor>a018c2f820f00bb67c3545b543d42005c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::CFunction</name>
    <filename>classRVL_1_1CFunction.html</filename>
    <templarg></templarg>
    <templarg>f</templarg>
    <base>BinaryLocalFunctionObject&lt; T &gt;</base>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1CFunction.html</anchorfile>
      <anchor>abdbd0721033f3d2e2e1452a50a47c4c2</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; &amp;y, LocalDataContainer&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1CFunction.html</anchorfile>
      <anchor>a78c57c5a238dccac65429af11834bb6d</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::CJacobian</name>
    <filename>classRVL_1_1CJacobian.html</filename>
    <templarg>T</templarg>
    <templarg>df</templarg>
    <base>UnaryLocalFunctionObjectConstEval&lt; T &gt;</base>
    <member kind="function">
      <type></type>
      <name>CJacobian</name>
      <anchorfile>classRVL_1_1CJacobian.html</anchorfile>
      <anchor>a15cc5b128feaac00abc363f4aaa373b6</anchor>
      <arglist>(GenMat&lt; T &gt; &amp;_a)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator()</name>
      <anchorfile>classRVL_1_1CJacobian.html</anchorfile>
      <anchor>add02a7a266e99d21e08b4cdc6de3c640</anchor>
      <arglist>(LocalDataContainer&lt; T &gt; const &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getName</name>
      <anchorfile>classRVL_1_1CJacobian.html</anchorfile>
      <anchor>a29971aef5758588763c0a9ddc965904e</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::OpWithGenMatDeriv</name>
    <filename>classRVL_1_1OpWithGenMatDeriv.html</filename>
    <templarg></templarg>
    <base>Operator&lt; T &gt;</base>
    <member kind="function" virtualness="pure">
      <type>virtual GenMat&lt; T &gt; const &amp;</type>
      <name>getGenMat</name>
      <anchorfile>classRVL_1_1OpWithGenMatDeriv.html</anchorfile>
      <anchor>a628fedd40c4e6f4134811aada07813bd</anchor>
      <arglist>() const =0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::GenOp</name>
    <filename>classRVL_1_1GenOp.html</filename>
    <templarg>T</templarg>
    <templarg>f</templarg>
    <templarg>df</templarg>
    <base>RVL::OpWithGenMatDeriv</base>
    <member kind="function">
      <type></type>
      <name>GenOp</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a34b4867c7259f5f6980607465fda9281</anchor>
      <arglist>(RnSpace&lt; T &gt; const &amp;_dom, RnSpace&lt; T &gt; const &amp;_rng)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GenOp</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>aaf39edf5ec8d5f8faa762c0df81a3424</anchor>
      <arglist>(GenOp&lt; T, f, df &gt; const &amp;a)</arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getDomain</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a34d5fc3dbe352cbc570788e31b13bf8e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Space&lt; T &gt; const &amp;</type>
      <name>getRange</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a4e6bf37de8c1c0c851f0f41ebd9e41e0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GenMat&lt; T &gt; const &amp;</type>
      <name>getGenMat</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>ae1f7014455dad1d405d7a2500fbc7161</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a50d266061e81911533ad77f80537a08a</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>apply</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>ae80e1b4f8d8503bc045ba4ea7fe1798b</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, Vector&lt; T &gt; &amp;y) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyDeriv</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>abec91bdc79192bf227605627215d93d5</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, const Vector&lt; T &gt; &amp;dx, Vector&lt; T &gt; &amp;dy) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual void</type>
      <name>applyAdjDeriv</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a72b4110daf82cdcbe7dfba79e3f65d6e</anchor>
      <arglist>(const Vector&lt; T &gt; &amp;x, const Vector&lt; T &gt; &amp;dy, Vector&lt; T &gt; &amp;dx) const </arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual Operator&lt; T &gt; *</type>
      <name>clone</name>
      <anchorfile>classRVL_1_1GenOp.html</anchorfile>
      <anchor>a96a7ed35239d3f750b93f88ad5f6173c</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>RVL::RnSpace</name>
    <filename>classRVL_1_1RnSpace.html</filename>
    <templarg>Scalar</templarg>
    <base>LocalSpace&lt; Scalar &gt;</base>
    <member kind="function">
      <type></type>
      <name>RnSpace</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>ad05c976c745e4b6e9dfa1f9d9a849441</anchor>
      <arglist>(int n=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>RnSpace</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a0ec5300cce881ca6cb231d1e0201eb9d</anchor>
      <arglist>(const RnSpace&lt; Scalar &gt; &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>isCompatible</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a085b63e06b6531aecc2f358a8e583cc1</anchor>
      <arglist>(DataContainer const &amp;dc) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~RnSpace</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a5907efec74016bf94ba8da3d9cadd855</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>getSize</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>acf72200492af3253b90bd9e84a770814</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a6eef721f3fd48bfeccb5a19441dbf84f</anchor>
      <arglist>(RVLException &amp;str) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual ostream &amp;</type>
      <name>write</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a6fbf2944b03969522adcd272af8914c9</anchor>
      <arglist>(ostream &amp;str) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>redim</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a1b518a51045f21a4886e8903e9315890</anchor>
      <arglist>(int newdim)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LocalDataContainerFactory&lt; Scalar &gt; &amp;</type>
      <name>getLDCF</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a5bc3259205e6cf0535e39a5aaf064989</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>LinearAlgebraPackage&lt; Scalar &gt; &amp;</type>
      <name>getLAP</name>
      <anchorfile>classRVL_1_1RnSpace.html</anchorfile>
      <anchor>a5ff240a270d7aa3cb13f7ea74b4a10fe</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
</tagfile>
