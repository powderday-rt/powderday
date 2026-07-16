Model Tests
**********

Here, we describe some tests of the model.  This growing suite of
tests are packaged with the repository so that developers can examine
the impact of any code updates.



SKIRT tests
============

 `SKIRT <http://www.skirt.ugent.be/root/index.html>`_ is an
 outstanding dust radiative transfer code developed by the Ghent group
 with a number of features that distinguish it from `powderday
 <https://github.com/powderday-rt/powderday.git>`_.  Here, we describe how to
 run model tests between the two codes.


 A few notes:

:Relevant Hash:
   
   The `powderday <https://github.com/powderday-rt/powderday.git>`_ hash for
   the most recent set of tests is `2fdf5e49074289382648a3be6a9478fcb9a0e84c <https://github.com/powderday-rt/powderday/commit/2fdf5e49074289382648a3be6a9478fcb9a0e84c>`_
   and the relevant `SKIRT
   <http://www.skirt.ugent.be/root/index.html>`_ hash is
   96e8be9761c3021498b847854bf84c3da3129555.
    

:Simulation Codes:

  The tests done here have been done for Gizmo simulations,
  Gasoline/Changa simulations with a tipsy output format, and Arepo.
  Note, that the Arepo test that ship with the code are currently
  configured to be run with the radiative transfer on the Voronoi mesh
  in `powderday <https://github.com/powderday-rt/powderday.git>`_, but
  on an octree grid in `SKIRT
  <http://www.skirt.ugent.be/root/index.html>`_.  Pull requests very
  welcome from community experts on setting up the appropriate ski
  file for `SKIRT <http://www.skirt.ugent.be/root/index.html>`_ with
  Voronoi mesh simulations.

:Downloads:

 The relevant
 parameter files (the ``.ski`` input file for `SKIRT
 <http://www.skirt.ugent.be/root/index.html>`_ and the
 ``parameters*.py`` files for `powderday
 <https://github.com/powderday-rt/powderday.git>`_ are packaged with the
 repository.
 

 The hydrodynamic simulations for the test runs that ship with
 `powderday <https://github.com/powderday-rt/powderday.git>`_ aren't included
 in the repository due to their size.  The gizmo simulation is the
 same mufasa zoom simulation as in the `Gadget/Gizmo quickstart
 <https://powderday.readthedocs.io/en/latest/quickstart.html#gadget-gizmo>`_
 section, and is available `here
 <https://www.dropbox.com/s/g6d47z3pm8l18p7/snapshot_134.hdf5?dl=0>`_.
 The Changa/Gasoline/Tipsy output is galaxy_00300 from the public `yt datasets page <https://yt-project.org/data/>`_
 
Actually Running the Tests: Powderday:
--------------

In what follows, we lead the user through the tests for a Gizmo type
simulation, though the process for testing the tipsy outputs is
similar.

To run the `powderday
<https://github.com/powderday-rt/powderday.git>`_ test code, run::

  pd_front_end.py tests/SKIRT/gizmo_mw_zoom parameters_master_gizmo parameters_model_gizmo

After having edited the ``parameters_master`` and ``parameters_model``
file for your specific paths.  

Actually Running the Tests: SKIRT:
--------------

`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ needs to know
the stellar and gas particle information from the hydrodynamic
simulation.  In ```analytics.py``, there's a function
``SKIRT_data_dump()`` that dumps these files to your  `powderday <https://github.com/powderday-rt/powderday.git>`_ output directory.  Here, it should have dumped files with a path like::
  tests/SKIRT/gizmo_mw_zoom/SKIRT.134.gas.particles.txt

`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ also requires an
input parameter file which can be created via a terminal input
following their `tutorials
<http://www.skirt.ugent.be/tutorials/index.html>`_.  We have created the ``.ski`` input file::
  tests/SKIRT/pd_test.dust.ski
  
`SKIRT <http://www.skirt.ugent.be/root/index.html>`_ requires paths
relative to the run directory.  After compiling `SKIRT
<http://www.skirt.ugent.be/root/index.html>`_, one possible setup might be::
  cd SKIRT/run
  mkdir pd_test
  cp <path_to_powderday>/tests/SKIRT/pd_test.dust.ski pd_test
  ../release/SKIRT/main/skirt pd_test/pd_test.dust.ski

In the ``run`` directory, this should create a file like:
``/home/desika.narayanan/SKIRT/run/test_i90_sed.dat``.  This has
chosen one particular viewing angle (which, in fact may not be the
same viewing angle as the `powderday
<https://github.com/powderday-rt/powderday.git>`_ simulation - in fact it was
arbitrarily chosen).  Then you can run::
  cd <path_to_powderday>
  python tests/SKIRT/skirt_sed_plot.py

and it should produce an image like the top left one in:

.. image :: images/powderday_skirt_comparison.png 
    :align: center

Note, there are analagous tests for the arepo and gasoline/changa
front ends that ship with `powderday
<https://github.com/powderday-rt/powderday.git>`_ as well that will produce the other panels in the aforementioned code comparison figure.
	  
Persistent Differences in the Models
--------------

While we have attempted to ensure as much of an apples-to-apples
comparison between codes, some differences are persistent that
manifest themselves in the emergent SEDs.

The input SEDs are different, as is evident in the UV portion of the comparisons. The former employ interpolated
BC03 stellar models at a lower resolution than the input MILES
spectral libraries for the `fsps
<https://code.google.com/p/fsps/source/checkout>`_ models that
`powderday <https://github.com/powderday-rt/powderday.git>`_ employs.



  


