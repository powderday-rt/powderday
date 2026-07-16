Getting Started
**********

Overview of Requirements
============

* **python>=3.7**

  * numpy (<1.24)
  * scipy
  * astropy (3.2.3)
  * unyt
  * h5py
  * scikit-learn
  * six
  * p_tqdm


* **compilers**
  
  * gcc
  * gfortran
  

* **Additional Packages (with Instructions Below)**
  
  * git  <http://git-scm.com/>
  * powderday <https://github.com/powderday-rt/powderday.git>
  * yt <http://yt-project.org>
  * python-fsps <https://dfm.io/python-fsps/current/>
  * Hyperion <http://www.hyperion-rt.org/>
  * Hyperion Dust Files <http://docs.hyperion-rt.org/en/stable/dust/dust.html>

Installation
============
    


Manual Installation
--------------

What follows is a self-contained installation manual, though for
problematic installs of any of the sub packages, it's definitely best
to visit the main docs on the main software site (which are always
linked below in each subsection).

.. _python:

python
--------------

`powderday <https://github.com/powderday-rt/powderday.git>`_ should work with python >=3.5 though is ideal with 3.8 (and some issues have been noted that may relate to python 3.7).
Please file an issue if you encounter one.

We very strongly recommend that the user set up a new python environment for the
`powderday <https://github.com/powderday-rt/powderday.git>`_ installation to avoid software conflicts.   This could look something like (assuming a ``conda`` installation of python)::

  >conda create --name pd_environment
  >source activate pd_environment

(And then when you want to exit the environment, you can type)::

  >source deactivate pd_environment

Then, whenever you're in the ``pd_environment``, everything you
install will remain contained to that particular installation of
python, and not conflict with any of your other installed packages.

.. _powderday:


powderday
--------------

Simply clone the latest and greatest from the repo::

  >git clone https://github.com/powderday-rt/powderday.git

To install, `cd` into the cloned repository and run the usual::

  >python setup.py install


.. _yt:



.. _Hyperion:

Hyperion
--------------

`Hyperion <http://www.hyperion-rt.org>`_ is the main work horse of
`powderday <https://github.com/powderday-rt/powderday.git>`_.  The full
directions for installation are well-described on the main
`Installation page for Hyperion
<http://docs.hyperion-rt.org/en/stable/installation/installation.html>`_.
Here, we summarize the installation which should get most users
through without any real difficulty.

There are two ways to install `Hyperion <http://www.hyperion-rt.org>`_.  The first is via ``conda``::

  >conda install -c conda-forge hyperion

Please note, though, that there is an issue with six no longer being
bundled with astropy that was fixed here:
https://github.com/hyperion-rt/hyperion/issues/219.  This said, at the
time of the last update of these docs (July 10th, 2020), this has not translated to the conda installation, meaning you will need to manually update all of the files listed here:

https://github.com/hyperion-rt/hyperion/issues/219#issuecomment-600036854  by replacing::

  >#from astropy.extern import six
  >import six

(for example, the files might be located in a location like:)::

  >home/desika.narayanan/miniconda3/envs/pd_test/lib/python3.6/site-packages/hyperion/filter/filter.py
  
The second and manual way to install `Hyperion
<http://www.hyperion-rt.org>`_ follows (note, for the manual installation you don't have to worry about the six replacement above):


#. First clone the main repository.::

     >git clone https://github.com/hyperion-rt/hyperion.git

#. Make sure that you have the correct modules loaded on your cluster.
   This will require a compiler, openmpi and HDF5.  For example, on
   the University of Florida HiPerGator supercomputing system, I would
   have::

   
  >module load intel/2020.0.166
  >module load git
  >module load hdf5/1.14.1
  >module load openmpi/4.1.5
  
#. Install the python module::

   >cd hyperion
   >python setup.py install


#. Ensure that if you type::

   >hyperion

it returns a sensible output.  It should return something along the lines of::

  >usage: hyperion [-h] [-f] [-m n_cores] input output
  >hyperion: error: too few arguments

If it can't find `Hyperion <http://www.hyperion-rt.org>`_, check the
the path that is near one of the last lines of the setup.py
installation (that is something associated with the number 755) and
make sure it's in your path.  It's most likely to be a python binaries
directory.

#. Install the submodules manually::

   >git submodule init
   >git submodule update

#. Install the Fortran binaries::

     > ./configure

or::

  >./configure --prefix=$HOME/local

or some such path if you aren't administrator on your computer.  Note
for this step you'll need your compilers, MPI and HDF5 installations
active (so, on a supercomputer you might need to load these modules
such as [for example, on the University of Florida HiPerGator
supercomputer])::

  
  
  >module load intel/2020.0.166
  >module load git
  >module load hdf5/1.14.1
  >module load openmpi/4.1.5
  

of course please be careful of mixing and matching compilers, and
ensuring that you have the same compilers loaded for all
installations.
  
#. Compile the code::

   > make
   > make install
   

Note this will take a while!  Make sure this works by typing at the command line::

  >hyperion_sph

which should return something like::

  >Usage: hyperion_sph [-f] input_file output_file


  .. _Hyperion_dust:

Hyperion Dust Files
--------------

Unless you've written your own dust files, you will likely want to use
the pre-compiled dust files developed by Tom Robitaille (though don't
ship with `Hyperion <http://www.hyperion-rt.org>`_ due to their size).


  >git clone https://github.com/hyperion-rt/hyperion-dust.git
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

If you want to use the PAH model in `powderday
<https://github.com/powderday-rt/powderday.git>`_, you'll additionally need
these files in the same dust directory.  To download, click on the link,
then click 'raw' on the right side of each page.

1. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/big.hdf5
2. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/vsg.hdf5
3. https://github.com/hyperion-rt/paper-galaxy-rt-model/blob/master/dust/usg.hdf5

Please note the caveat that the PAH files are generated using some
approxmations described in `Robitaille et
al. <http://www.aanda.org/articles/aa/abs/2012/09/aa19073-12/aa19073-12.html>`_,
and we encourage the user of these PAH files to read this paper,
especially section 3.4.2.


yt
--------------------

Next we need `yt <http://yt-project.org>`_ - to install this, clone the source and install::

  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .

Note, it is important to install this *after*  `Hyperion <http://www.hyperion-rt.org>`_.  This is because  if you used the conda installation of `Hyperion <http://www.hyperion-rt.org>`_ , then `yt <http://yt-project.org>`_ 3.x ships with it and auto-installs. However, powderday is no longer compatible with `yt <http://yt-project.org>`_ 3.x.






.. _python-fsps:

python-fsps
--------------

To install::

  >export SPS_HOME="/path/where/you/want/to/download/fsps"
  >git clone https://github.com/cconroy20/fsps.git $SPS_HOME
  >python -m pip install fsps

    
  
You can test the installation by opening python and typing::

>import fsps


Note, we no longer need to actually install `fsps
<https://code.google.com/p/fsps/source/checkout>`_ anymore.  This is
installed within `python-fsps <https://dfm.io/python-fsps/current/>`_
itself!  Please see the `python-fsps
<https://dfm.io/python-fsps/current/>`_ docs for how to change
underlying stellar physics when installing (i.e. the spectral
libraries or the stellar isochrones).





Troubleshooting your Installation
============

  .. numpy issues:

Numpy Issues
---------------

* np versions >=1.24 have deprecated float that causes (waves hands wildly) everything to break.   Roll back via::

  >conda install -c conda-forge numpy=1.23



  .. _fsps installation issues:

fsps Installation Issues
---------------
* One possibility can be that there are issues in compiling
   src/autosps.f90.  One solution is to replace RETURN with STOP in
   line 21.



  .. _python-fsps installation issues:

python-fsps installation issues
--------------
* With intel compilers (e.g., on the University of Florida HiPerGator system) you should try::
     
   >CC=icc F90=ifort python setup.py install

*  `python-fsps <https://dfm.io/python-fsps/current/>`_ can't find f2py
   
   f2py is a numpy package that is sometimes named f2py2.7 by numpy.
   At the same time, `python-fsps
   <https://dfm.io/python-fsps/current/>`_ expects it
   to be called f2py (as it sometimes is; for example in Anaconda).
   So, you might need to locate f2py (it ships with `yt
   <http://yt-project.org>`_, so if you for example use the `yt
   <http://yt-project.org>`_ python) you need to link the following
   files::

   >cd /Users/desika/yt-x86_64/bin
   >ln -s f2py2.7 f2py

   and::

   >cd /Users/desika/yt-x86_64/lib/python2.7/site-packages
   >ln -s numpy/f2py/ f2py

   This should hopefully fix it.


* Issues with 'f2py' in the  `python-fsps
   <https://dfm.io/python-fsps/current/>`_ installation:

   Numpy has made some changes to f2py in the 1.10.x version of numpy.
   The easiest fix is to use a non 1.10.* version of numpy (thanks to
   Ben Johnson for finding this).

*  `python-fsps <https://dfm.io/python-fsps/current/>`_ has mysterious
installation failures.  Often this has to do with a bad `FSPS
<https://github.com/cconroy20/fsps>`_ compilation. Even if it seems
like `FSPS <https://github.com/cconroy20/fsps>`_ has compiled, it may
not actually execute properly if the correct compilers aren't set in
the MakeFile.  Thanks to Ena Choi for pointing this one out.

  .. _hyperion installation issues:


Hyperion Installation Issues
---------------

  .. _yt installation issues:

   
yt Installation Issues
---------------

* If you have trouble with this installation, you may want to unload
your openmpi module that you previously had loaded for the `Hyperion
<http://www.hyperion-rt.org>`_ install.



* Another common trick to help the installation is to install with::

   >LDSHARED="icc -shared" CC=icc pip install -e .


* Finally, even if you're installing mostly everything else from
  source, there's no issue usually with installing yt via Conda.  This
  can often times work well with intel compilers, which yt can be a
  bit fussy about sometimes.::

    >conda install --channel conda-forge yt

   
System Specific Installation Notes
============

HiPerGator at the University of Florida
--------------

[1] The first set of instructions for the University of Florida
HiPerGator3.0 facility is to employ intel compilers, and to compile
everything manually.  This allows the greatest flexibility, as well as
the ability to use private forks of individual codes.

First, load up the compilers that we'll use throughout (though note: openmpi is not loaded until after yt is installed as yt will sometimes bork due to openmpi)::

  >module load intel/2020.0.166
  >module load git
  >module load hdf5/1.14.1


yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .

Note, if you have trouble, please see the troubleshooting below.

fsps and python-fsps:

The development version of python-fsps now includes the Fortran FSPS source code::

  >cd $HOME
  >git clone --recursive https://github.com/dfm/python-fsps.git

then in your .bashrc set the analog to::
  
  >export SPS_HOME=$HOME/python-fsps/src/fsps/libfsps
  
  >cd python-fsps
  >CC=icc F90=ifort python -m pip install .


Before going forward, please try::

  >python
  >import fsps

and ensure that it does not throw any errors.  If you get an error along the lines of::

  >ImportError: /blue/narayanan/desika.narayanan/conda/envs/test/lib/python3.8/site-packages/fsps/_fsps.cpython-38-x86_64-linux-gnu.so: undefined symbol: getenv_

then try to install via pip::

  >python -m pip install fsps
  
Next, before installing hyperion, lets make sure our openmpi is loaded::

    >module load openmpi/4.1.5


hyperion::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >pip install .
  >git submodule init
  >git submodule update

  >./configure --prefix=$HOME/local

  >make
  >make install

hyperion dust::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion-dust.git
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

  
powderday::

  >git clone https://github.com/powderday-rt/powderday.git
  >cd powderday
  >conda install numpy scipy cython h5py matplotlib psutil joblib six astropy scikit-learn ipython
  >python setup.py install

[2] The second set of instructions use gcc, but a manual installation of everything. Thanks to Prerak Garg for these.::

First, load up the compilers that we'll use throughout::

  >module load gcc/12 openmpi/4.1.5 hdf5/1.14.1 git libz

  
yt::

  >cd $HOME
  >git clone https://github.com/yt-project/yt
  >cd yt
  >pip install -e .


fsps and python-fsps

The development version of python-fsps now includes the Fortran FSPS source code::

  >cd $HOME
  >git clone --recursive https://github.com/dfm/python-fsps.git


then in your .bashrc set the analog to::
  
  >export SPS_HOME=$HOME/python-fsps/src/fsps/libfsps
  
  >cd python-fsps
  >CC=gcc F90=gfortran F77=gfortran python -m pip install .

Now load up openmpi::
  >ml  openmpi/4.1.1


hyperion::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion.git
  >cd hyperion
  >pip install .
  >git submodule init
  >git submodule update

  >./configure --prefix=$HOME/local

  >make
  >make install

hyperion dust::

  >cd $HOME
  >git clone https://github.com/hyperion-rt/hyperion-dust.git
  >cd hyperion-dust-0.1.0
  >python setup.py build_dust

  
powderday::

  >git clone https://github.com/powderday-rt/powderday.git
  >conda install numpy scipy cython h5py matplotlib psutil joblib six astropy scikit-learn ipython
  >cd powderday
  >python setup.py install




  
  

[3] The third set of instructions use gcc, and the conda installation
of `Hyperion <http://www.hyperion-rt.org>`_.  Thanks to Paul Torrey
for these.::

  >module load  gcc/12 openmpi/4.1.5 hdf5/1.14.1 git libz 
  >conda install -c conda-forge hyperion
  >python -c "import hyperion" (just to ensure no errors thrown)
  >hyperion (just to ensure command is found)
  >python -m pip install fsps
  >[set $SPS_HOME variable in .bashrc)
  >cd $HOME
  >git clone https://github.com/powderday-rt/powderday.git
  >cd powderday
  >python setup.py install

then fix import six line in the equivalent of all of these::

  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/util/validator.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/conf/conf_files.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/filter/filter.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/dust/dust_type.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/model_output.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/flared_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/alpha_disk.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/bipolar_cavity.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ulrich_envelope.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/power_law_envelope.py 
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/densities/ambient_medium.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/sed.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/model/image.py
  >vi /home/paul.torrey/.conda/envs/pd_gcc/lib/python3.8/site-packages/hyperion/grid/yt3_wrappers.py
