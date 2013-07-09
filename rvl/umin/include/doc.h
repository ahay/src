/** \mainpage Unconstrained Optimization Package

Algorithms currently implemented:

<ul>
<li>RVLUmin::CGAlg - conjugate gradient iteration for SPD linear systems with trust region truncation. Combines RVLUmin::CGStep with RVLAlg::CountingThresholdIterationTable and RVLAlg::BallProjTerminator terminators in RVLAlg::LoopAlg. The RVLAlg::BallProjTerminator contribution makes this a component of a trust region algorithm.</li>
<li>RVLUmin::CGNEAlg - conjugate gradient iteration for solution of least squares problems, with trust region truncation. Avoids explicit formation of normal operator. Combines RVLUmin::CGNEStep with RVLAlg::VectorCountingThresholdIterationTable (displays iteration count, residual norm, and normal residual norm) and RVLAlg::BallProjTerminator terminators in RVLAlg::LoopAlg. The RVLAlg::BallProjTerminator contribution makes this a component of a trust region algorithm.</li>
<li>RVLUmin::LBFGSBT - limited memory BFGS implementation. Combines limited memory BFGS search direction comptation (RVLUmin::LBFGSDir) with geometric RVLUmin::BacktrackingLineSearchAlg to create a line search step (RVLUmin::UMinStepLS), and combines that with a RVLAlg::CountingIterationTable terminator in a RVLAlg::LoopAlg</li>
<li>RVLUmin::TRGNAlg - implements Steihaug-Toint truncated Gauss-Newton-Krylov algorithm. Combines RVLUmin::TRGNStep with RVLAlg::VectorCountingThresholdIterationTable (displays iteration number, residual, normal residual, predicted reduction, actual reduction, and trust radius) in a RVLAlg::LoopAlg.
<li>RVLUmin::PowerMethod - power iteration to estimate largest singular value (operator norm) of a linear operator. Uses a computable estimate for singular value relative error to stop. Combines RVLUmin::PowStep with RVLAlg::VectorCountingThresholdTable (displays iteration number, relative residual, and singular value estimate) in a RVLAlg::LoopAlg.
</ul>

See class documentation files for more information.

<hr>
<a href="../../../doc/html/index.html">RVL Home Page</a>

*/

