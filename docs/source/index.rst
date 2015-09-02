.. Documentation to the TopologyTracer
.. ===================================

1. Getting started
^^^^^^^^^^^^^^^^^^

The compilation of TopologyTracer utilizes standard Linux tools.

.. toctree::
   :maxdepth: 2
     
   setup
   

2. Creating input
^^^^^^^^^^^^^^^^^

It is necessary to have collected results in the form of a sequence of microstructure snapshots. In its first version data from the 2D anisotropic grain growth simulation packages are utilized.

.. toctree::
   :maxdepth: 2
   
   input


3. Tracing grains
^^^^^^^^^^^^^^^^^

Tracing, i.e. the process of postprocessing time-resolved data with TopologyTracer is invoked with the following command line call::

   mpiexec -n <proc> topologytracer <1> <2>
   
By default it is assumed that the executable is placed in the same folder as the raw data.

.. toctree::
   :maxdepth: 2
   
   tracing
   
   
4. Visualizing results
^^^^^^^^^^^^^^^^^^^^^^

In its current state the TopologyTracer does not contain a GUI yet. Instead, plain MPI I/O raw files are being generated and postprocessed with the MATLAB scripts available in the *scripts* folder.


References
^^^^^^^^^^

| Kuhbach M., Barrales-Mora L.A., Miessen C., Gottstein G.: 
| **Ultrafast analysis of individual grain behavior during grain growth by parallel computing** 
| Proceedings of the 36th Riso International Symposium on Materials Science 
| doi:10.1088/1757-899X/89/1/012031
   
   
Licence
^^^^^^^

The project is licenced under the GNU v2.0


.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


Questions, contributions? Please check contact markus.kuehbach@rwth-aachen.de