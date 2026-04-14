# splitEE_curvilinear

This is a companion software to the manuscript

[M. Caliari and F. Cassini. A tensor-based exponential integrator for
diffusion--reaction equations in common curvilinear coordinates,
arXiv preprint, arXiv:2604.11558, 2026](https://doi.org/10.48550/arXiv.2604.11558)

In particular, it contains all the functions and scripts needed to reproduce
the numerical experiments of the manuscript.

It is fully compatible with MathWorks(R) MATLAB(R) and GNU Octave.

We require some external functions for the execution of some scripts.
In particular:
* For the complete execution of the code in ```test_order_comparison.m```
  we need [KIOPS](https://gitlab.com/stephane.gaudreault/kiops).
* For the tensor-matrix operations in ```BS_schnakenberg_sphere.m``` and
  ```BS_DIB_cylinder.m``` we need the relevant functions in
  [KronPACK](https://github.com/caliarim/KronPACK).

