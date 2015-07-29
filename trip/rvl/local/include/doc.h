/** \mainpage Local RVL - a simple realization of the RVL data management classes, based on data containers which explose a scalar array

    Main components:

    <ul>
    <li>Local data package: 
    <ul>
    <li>RVL::LocalDataContainer - core abstraction, subtype of RVL::DataContainer, templated on Scalar type, with two essential pure virtual public methods to be implemented in subclasses:
    <ul>
    <li>RVL::LocalDataContainer::getSize returns the number of words in the data array (like the size function of the standard library containers)</li>
    <li> RVL::LocalDataContainer::getData returns a pointer to Scalar</li>
    </ul>
    In principle Scalar is very little restricted, and could be a complex type itself. In its intended use, Scalar is a number type with value semantics.
    <p>
    The eval methods inherited from DataContainer are implemented. This Element (Visitor design pattern) admits only two types of Visitor, namely RVL::LocalEvaluation and RVL::LocalConstEval (described in next bullet).
    </li>
    <li>RVL::LocalSpace - RVL::Space subtype whose buildDataContainer method returns a RVL::LocalDataContainer</li>
    <li>RVL::LocalVector - RVL::Vector subtype specializing vector in RVL::LocalSpace, offering direct access to instance array data through getSize and getData methods which delegate to the methods of the same name of the RVL::LocalDataContainer it owns.</li>
    </ul>
    </li>

    <li>Local function object package: function objects which may visit LocalDataContainers

    <ul>
    <li>RVL::LocalEvaluation - mixin class with operator() method taking std::vector of arguments </li>
    <li>RVL::LocalFunctionObject - child class of RVL::FunctionObject, RVL::LocalEvaluation</li>
    <li>LocalFunctionObject hierarchy - specializations of the above with 1 - 4 arguments:
    <ul>
    <li>RVL::UnaryLocalFunctionObject</li>
    <li>RVL::BinaryLocalFunctionObject</li>
    <li>RVL::TernaryLocalFunctionObject</li>
    <li>RVL::QuaternaryLocalFunctionObject</li>
    </ul>
    </li>
    <li>RVL::LocalConstEval - mixin class with operator() method taking std::vector of const arguments</li>
    <li>LocalFunctionObjectConstEval hierarchy - child classes of RVL::FunctionObjectConstEval, RVL::LocalConstEval, 1 - 4 arguments
    <ul>
    <li>RVL::UnaryLocalFunctionObjectConstEval</li>
    <li>RVL::BinaryLocalFunctionObjectConstEval</li>
    <li>RVL::TernaryLocalFunctionObjectConstEval</li>
    <li>RVL::QuaternaryLocalFunctionObjectConstEval</li>
    </ul>
    </li>
    <li>LocalFunctionObjectScalarRedn hierarchy - child classes of RVL::FunctionObjectScalarRedn, RVL::LocalConstEval, 1 - 4 arguments
    <ul>
    <li>RVL::UnaryLocalFunctionObjectScalarRedn</li>
    <li>RVL::BinaryLocalFunctionObjectScalarRedn</li>
    <li>RVL::TernaryLocalFunctionObjectScalarRedn</li>
    <li>RVL::QuaternaryLocalFunctionObjectScalarRedn</li>
    </ul>
    </li>
    </ul>
    </li>


    <li>Rn package:
    <ul>
    <li>RVL::RnArray - simple array child class of RVL::LocalDataContainter. Like GSL's Block, but with constructors, destructor, and access control</li>
    <li>RVL::RnSpace - completely implemented RVL::LocalSpace subclass, whose vectors own RVL::RnArray data structs</li>
    <li>RVL::GenMat - a simple RVL::LinearOp subclass implementing matrix multiplication, with a symmetric specialization (RVL::SymMat). These store the matrix elements in an RVL::RnArray, and act on any vectors of appropriate dimensions whose data structures are RVL::LocalDataContainer~s
    </li>
    </ul>

<li>ContentPackage package:  
<ul>
<li>RVL::ContentPackage - a <i>content package</i> is a pair consisting of (1) a <i>metadata</i> ("data about data") object, and (2) an array of <i>data</i> objects. This construct is naturally recursive: the metadata and/or data objects can themselves be content packages. This is an extremely flexible data structure which accommodates many types of scientific data. For example, an RnArray could have been (but was not) implemented as a content package with int metadata (the dimension) and scalar data. A matrix (or multivector, in Trilinos-speak) is naturally a content package with int metadata and RnArray data. A seismic trace in SEGY format is a content package with trace header metadata and float data. And so on. This concrete class template implements all such data structures and many more.</li>
<li>RVL::PackageContainer - essentially, a content package server. Subclass of RVL::DataContainer. Implements eval methods via RVL::PackageContainer::get, RVL::PackageContainer::put, and RVL::PackageContainer::reset methods, which provide sequential access to the RVL::ContentPackage array virtualized by this type. Designed for distributed or out-of-core data manipulation, for which only one (or a buffer's worth) of RVL::ContentPackage objects may be held in core at one time.</li>
<li>RVL::PackageContainerSpace - RVL::StdSpace instance whose buildDataContainer method returns a RVL::PackageContainer</li>
</ul>
</li>

<li>Functions package:
<ul>
<li>Implemented functions collection: local function object templates that carry out various standard array calculations. Some attention paid to efficiency - these may be used in reasonably large scale applications.
<ul>
<li>RVL::RVLCopy (float and double template specializations use memcpy)</li>
<li>RVL::RVLScale</li>
<li>RVL::RVLMax</li>
<li>RVL::RVLMin</li>
<li>RVL::RVLL2innerProd (specialization for std::complex scalar types implements Hermitian form)</li>
<li>RVL::RVLAddAccumulate</li>
<li>RVL::RVLAssignConst</li>
<li>RVL::RVLRandomize (with complex and int specializations)</li> 
<li>RVL::ASCIIReader, RVL::ASCIIWriter</li>
<li>RVL::BinaryReader, RVL::BinaryWriter</li>
<li>RVL::RVLBoxMaxStep (max step to boundary of box, used for bound constraints in optimization)</li>
<li>RVL::ElementwiseMultiply (dot-star)</li>
<li>RVL::ElementwiseDivision (dot-slash with safeguard)</li>
<li>RVL::ElementwiseSqrtAbs</li>
</ul>
</li>
<li>Functions-to-function-objects: wrap C functions to produce function objects which act by evaluating the function on each component of the inputs - one to six arguments
<ul>
<li>RVL::ScalarFO1</li>
<li>RVL::ScalarFO2</li>
<li>RVL::ScalarFO3</li>
<li>RVL::ScalarFO4</li>
<li>RVL::ScalarFO5</li>
<li>RVL::ScalarFO6</li>
</ul>
</li>
</ul>
</li>
    </ul>
<hr>
<a href="../../../doc/html/index.html">RVL Home Page</a>


*/
