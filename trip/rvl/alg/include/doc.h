/** \mainpage Algorithm

Invented by Tony Padula as part of his PhD thesis project (<a
 href="http://www.caam.rice.edu/tech_reports/2005_abstracts.html#TR05-11">Software
design for simulation-driven optimization</a>, TR 05-11, CAAM, Rice
University 2005). 

Abstracts iterative algorithm structure, provides a boolean framework for termination criteria.

Two main types:
<ul>
<li>RVLAlg::Algorithm - an absurdly simple abstraction. Public interface consists of one pure virtual method, run, which takes no arguments and returns no value - just "push the button".</li>
<li>RVLAlg::Terminator - almost as simple. Public interface also has one pure virtual method, query, which takes no arguments and returns a boolean, the answer to the question: "should this algorithm stop?"</li>
</ul>

Algorithms and Terminators can be combined to form a variety of other Algorithms and Terminators, which are subclasses of the abstract bases - lists, loops, conditional executions, etc. The class documentation describes exactly how this happens; here is a UML diagram of the principal Algorithm subclasses. Supplied subclasses of Terminator provide a boolean algebra, also various types of output and even interactive input ("would you like this algorithm to stop?").

\image html alg.png

<hr>
<a href="../../../doc/html/index.html">RVL Home Page</a>
*/
