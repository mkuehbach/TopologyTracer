**References**
==============

The method
^^^^^^^^^^

**The method** and the key to its implementation.
 | Kühbach, M.
 | *On the Significance of the Long-Range Environment and Capillary Contributions*
 | *for Nucleating Abnormal Grain Growth and Recrystallization*
 | submitted to Acta Materialia, 2018

 | Kühbach, M., Mießen, C., Barrales-Mora, L. A., Gottstein, G.
 | *Simulation and data-analytics of sub-grain growth with consideration* 
 | *of stored elastic energy and anisotropic grain boundary properties*
 | Proceedings of the 6th International Conference on Recrystallization and Grain Growth
 | 2016 (ReX & GG 2016) in the Omni William Penn Hotel in Pittsburgh, PA, U.S.
 | *ISBN 978-1-119-32835-3*
 
 | Kühbach M., Barrales-Mora L.A., Mießen C., Gottstein G.: 
 | *Ultrafast analysis of individual grain behavior during grain growth by parallel computing*
 | Proceedings of the 36th Riso International Symposium on Materials Science 
 | http://dx.doi.org/10.1088/1757-899X/89/1/012031

**Further details** are provided within my dissertation entitled
 | Kühbach, M.
 | *Efficient Recrystallization Microstructure Modeling by Utilizing Parallel Computation*
 | *PhD thesis RWTH Aachen University*
 | Successfully defended and accepted for publication in August, 2017
 | Finally open source published in January, 2018
  

Contact
^^^^^^^
Feel free to contact me and allow me to identify whether the TopologyTracer supplies functionalities which can cater also your analyses needs! Interested? Please contact me_

 .. _me: m.kuehbach@mpie.de 
  
Specifically about GraGLeS
^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, the GraGLeS model is one of the very few grain coarsening simulation packages
available that is capable of stressing the TopologyTracer. 

 | Mießen, C. Velinov, N., Gottstein, G., Barrales-Mora, L. A.
 | *A highly efficient 3D level-set grain growth algorithm tailored for ccNUMA architecture*
 | Modelling and Simulation in Materials Science and Engineering
 | http://dx.doi.org/10.1088/1361-651x/aa8676
 
 | Mießen C., Liesenjohann M., Barrales-Mora L.A., Shvindlerman L.S., Gottstein G.
 | *An advanced level set approach to grain growth – Accounting for grain boundary anisotropy and finite triple junction mobility*
 | Acta Materialia, 2015, 99
 | http://dx.doi.org/10.1016/j.actamat.2015.07.040
 
 | Mießen, C.
 | *A massive parallel simulation approach to 2d and 3d grain growth*
 | *PhD thesis RWTH Aachen University*
 | Successfully defended in November, 2017
 

Third-party contributions
^^^^^^^^^^^^^^^^^^^^^^^^^
The TopologyTracer makes use of third-party contributions for some of its functionalities:

**RapidXML:** for the XML-based processing of control files and parameterization
 | Kalicinski, M.
 | http://rapidxml.sourceforge.net

**Poly2Tri:** for the triangularization of non-self-intersecting polygons
 | Liang, W.
 | http://sites-final.uclouvain.be/mema/Poly2Tri/

**Robust predicates** are required for numerical stable computational geometry. For these tasks Poly2Tri relies on 
 | Shewchuk, J. R.
 | *Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates*
 | Computational Geometry, 1997, 18, p305
 | http://dx.doi.org/10.1007/PL00009321

**Boost C++ library** for probing the project folder structure to obtain meta data automatically
 | http://www.boost.org/
 
**Eigen** open-source library for performing various geometrical operations and the fitting of circles to point clouds for principal curvature radius estimation.
 | http://eigen.tuxfamily.org/

The **detection of the intersection area** between arbitrary triangles and a circle was inspired by the algorithm idea described in
 | http://stackoverflow.com/questions/540014/compute-the-area-of-intersection-between-a-circle-and-a-triangle/

Furthermore, I would like to acknowledge the work of S. Strobl et. al. who provided an open source algorithm 
to the numerical robust computation of the **intersection volume of an arbitrary tetrahedron and with a sphere**:

 | Strobl, S., Formella, A., Pöschel, T.
 | *Exact calculation of the overlap volume of spheres and mesh elements*
 | Journal of Computational Physics, Vol 311, 2016, p158
 | http://dx.doi.org/10.1016/j.jcp.2016.02.003
 
Additionally, I acknowledge the tetrahedralization project **TetGen** by Hang Si, which enabled me at least to probe in a substantiated
manner the potential extension of the above-mentioned intersection problem in 3d:

 | Hang Si
 | TetGen A Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator
 | http://wias-berlin.de/software/tetgen/
 
Its source code is not part of the TopologyTracer due to licencing issues but can be linked to the program, for which
the interested user is referred to the source code.

