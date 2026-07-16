Detailed Description of Parameters
**********

There are two parameters files - parameters_master and
parameters_model.  The distinction between where a parameter goes is
somewhat arbitrary, but is based on the following.  Parameters_master
tends to contain parameters that will likely be shared amongst all
snapshots in a given galaxy run, while parameters_model are parameters
that might change from run to run (like the snapshot name or the
galaxy center).

The parameters files are (somewhat non-traditionally) written as
python code files because this enables you to embed little snippets in
them that may be particular to your galaxy run.  For example, for
gadget snapshot naming conventions, it can be useful to have a snippet along the lines of::

  snapshot_num = 20
  snapnum_str = '{:03d}'.format(snapshot_num)

  snapshot_name = 'snapshot_'+snapnum_str+'.hdf5'


That ensures a 3 digit snapshot number, common to many gadget-style
simulations.


parameters_master
============

Resolution Keywords
------------

:oref:

   Over Refinement of the Octree.  For particle-based codes, 1 means
   each data holding cell (a False) gets refined one additional time,
   even after the octree refinement criteria has stopped.  Very heavy
   on the memory but can enable higher pixel resolution for images.
   Default is 0.

:n_ref:
   
   Refinement criteria for octree refinement for particle-based codes.
   This is the threshold number of particles to refine over.  When
   nparticles > n_ref the octree refines further.  Default is 64.

:zoom_box_len:

   Side length to zoom in on.  Is +/- zoom_box_len from the center.
   Units are proper kpc.  So, a grid centered on [0,0,0] with
   zoom_box_len = 200 would extend from [-200,200] kpc in physical
   units at the redshift of the simulation.

:bbox_lim:

   Initial bounding box of grid for particle simulations
   (+/- bbox_lim).  Units are kpc.  This must encompass all of the
   particles in a simulation currently.  This just has to be a big
   number, but you want to be careful of making *too* large as
   precision limitations only allow for up to 20 levels of refinement.



Parallelization
------------

:n_processes:

   Number of pool processes to run the radiative transfer on.  Note,
   the stellar population synthesis will only run on as many
   processors as are on a core since its parallelization is pool.map
   (not MPI)

:n_MPI_processes:

    Number of MPI tasks to run. For TORQUE this is best set as the same 
    as n_processes, while for SLURM this may not be the case.


RT Information
------------

For all photon counts, a decent rule of thumb is 10-100x the number of
grid cells that you have, though of course you should check the
convergence properties of your simulation.

:n_photons_initial:

   Number of photons to use in main iterations (for the whole grid)
   for specific energy and temperature calculations.

:n_photons_imaging:

   Number of photons to use for the SED/image calculation

:n_photons_raytracing_sources:

   If raytracing is set (which is the default hard-coded into the
   code), number of raytracing photons to use for source emission.

:n_photons_raytracing_dust:

   Similar to n_photons_raytracing_sources but for dust emission.

:n_photons_DIG:
    
    Number of photons to use to calculate non-monochromatic SED to get the energy 
    dumped in each gas cell for DIG calculation
   
:FORCE_RANDOM_SEED:

    Boolean. True means the seed specified below will be used for random number
    generation in the Hyperion model. False means no seed will be set.

:seed:

    The seed with which to initialize random number generation in Hyperion. 
    Must be a negative integer.

Dust Information
------------

:dustdir:

   String. Path to where your dust files are.  String format -
   (e.g. '/home/desika/hyperion-dust-0.1.0/dust_files/')

:dustfile:
   
   String. Name of your main dust file.  String format -
   (e.g. 'd03_3.1_6.0_A.hdf5')

:PAH:

   Boolean - True means use model for PAHs, False means don't.


:dust_grid_type:

   String that is of one of the following values:
   ['dtm','rr','manual','li_bestfit'] that designates how we will
   determine the dust density in each octree grid cell.  The choices are:

   :dtm:
      A constant dust to metals ratio set by the keyword dusttometals_ratio
   :rr: Follows the best fit relation from Remy Ruyer
      2014A&A...563A..31R powerlaw relation between the the gas to
      dust ratio and metallicity (in which Xco is allowed to vary with
      the metallicity as well).
   :manual: Looks for a dust mass associated with the hydrodynamic
            simulation (i.e. if the hydrodynamic simulation carries
            dust information with it, as some flavors of gizmo and
            arepo do).  For gizmo and arepo currently assumes this
            information is stored in PartType0.
   :li_bestfit: Takes the best fit relation from the cosmological
                simulations of Li, Narayanan & Dave 2019 (MNRAS,
                490, 1425) that relates the dust to gas ratio to the
                metallicity for main sequence galaxies.  note - this
                applies galaxy integrated properties to individual cells for the radiative transfer simulations.

   :li_ml: Uses a the Extreme Randomized Trees learning algorithm from
           the aforementioned Li, Narayanan & Dave 2019 paper to
           predict the dust content based on the physical properties
           of the gas.  The algorithm here is modified from the
           original LND2019 paper in that it is trained on particles,
           as opposed to the aggregated properties of galaxies.

	   Note, for this algorithm to work, one needs to download the
	   training set datafile (a z=0 snapshot from the Simba
	   cosmological simulation: Dave et al. 2019, MNRAS, 486, 282)
	   and place it in the directory powderday/mlt.  This snapshot
	   is located here: https://www.astro.ufl.edu/~desika.narayanan/powderday_files

:dusttometals_ratio:

   Dust mass to metals mass ratio

:enforce_energy_range:

   Boolean. False ensures energy conservation.  But the emisivities
   may not be strictly correct if the energy in a cell is out of range
   of the emissivities.  True modifies the energy in the simulation,
   but ensures that the emissivities are consistent with the energy.
   See:
   <http://docs.hyperion-rt.org/en/latest/api/hyperion.model.Model.html?highlight=enforce#hyperion.model.Model.set_enforce_energy_range>

:SUBLIMATION:

    Boolean. If true, dust grains above the sublimation 
    temperature are automatically killed. Fast mode is to set this to False.

:SUBLIMATION_TEMPERATURE:

    The temperature in Kelvin above which dust grains are automatically killed. 
    Meaningless if SUBLIMATION == False.
   

Hydro Code Units
------------

Currently these are actually not used in `powderday
<https://github.com/powderday-rt/powderday.git>`_).  They remain in the
parameters file as a placeholder though as we may need them as an
over-ride if we find some HDF5 files don't contain this information.

:unit_mass:

   Mass code units for galaxy simulation.  Units: Msun/h

:unit_length:

   Length code unit for galaxy simulation.  Units: kpc/h

:unit_age:

   Stellar age units.  Units: Gyr/h

:unit_velocity:

   Velocity code unit for galaxy simulation.  Units: cm/s


Stellar SEDs Info
------------

:Force_BINNED:

   Boolean. True means force binning of the stellar SEDs (in bins of
   age and metallicity).  If False star particles below max_age_direct (next parameter)
   are individually represented in the stellar SED as opposed to broken up into bins.
   Note: Setting this to False can be slow to run, and hard on the memory.

:max_age_direct:
   
   Age (in Gyr) below which star particles will not be binned for getting stellar SEDs (works only if FORCE_BINNED is False) 
   To force all the star particles to be added without binning set this to an age greater than the maximum stellar age 
   of the galaxy (say 16 Gyr for example)


:imf_type:

   IMF parameter for stellar pops calculations.

   0. Salpeter
   1. Chabrier
   2. Kroupa
   3. Van Dokkum
   4. Dave

   Though note options 3 and 4 are currently not supported.

:imf1:
    
    Logarithmic slope of the IMF over the range 0.08 < M < 0.5. Only used if imf_type = 2. (Default: 1.3)    

:imf2:

    Logarithmic slope of the IMF over the range 0.5 < M < 1.0. Only used if imf_type = 2. (Default: 2.3)

:imf3:

    Logarithmic slope of the IMF over the range 1.0 < M < 120. Only used if imf_type = 2. (Default: 2.3)

:pagb:

    Weight given to post AGB stars.  1 is the default.

:add_agb_dust_model:
    
    Add circumstellar AGB dust model (100%); Villaume, Conroy & Jonson 2015

Nebular Emission Info
------------

:add_neb_emission:

    Boolean. If set to True, nebular line emission will be added. (under active development)

:use_cloudy_tables:
    
    If True, CLOUDY look up tables (dev. by Nell Byler) will be used to calculate nebular emission.
    Otherwise CLOUDY models are generated individually for each young star particle. 
    Note:  The lookup tables work only for stars particles below 10 Myr. (Default: True)

:use_cmdf:

    If True, star particles that have mass greater than cmdf_mas_mass (defined below) are broken down using a cluster mass
    distribution defined as dN/dM goes as M^(beta). This works irrespecitve of whether nebular emission is turned on or not. 
    The cmdf is set by the following parameters defined below: cmdf_min_mass, cmdf_max_mass, cmdf_bins and cmdf_beta.

:cmdf_min_mass:

    Minimum mass of the star clusters in units of log(Msun). Note: Results might be inconsistent if
    set lower than 3.5. (See Chandar et al.2014 for more info) (Default = 3.5)

:cmdf_max_mass:
    
    Maximum mass of the star clusters in units of log(Msun). (Default = 5.0). Note: Only star particles that
    have a mass greater than this parameter are broken down.

:cmdf_bins:

    The number of bins used for calulating the cluster mass distribution function (Default = 6.0)

:cdmf_beta:
    
    The power law exponent (beta) for calculating CMDF (dN/dM goes as M^(beta))

:use_age_distribution:
    
    If True, star particles with ages between age_dist_min and age_dist_max (next parameters) into divided into
    an ensemble of particles all of whom have the same properties except their age which is picked from a power law age 
    distribution of the form dN/dt is proportional to t^-0.65 (Imp: This can only be used if use_cmdf is also set to True).
    Note: The function has a bunch of tunable parameters that can be changed though we feel that their default values
    should be good enough for most cases. The function is located in cloudy_tools.py file under powderday/nebular_emission.

:age_dist_min: 
    
    Star particle above this age are sub-divided into an age distribution if use_age_distribution is set to True
    (Units: Gyr, Default = 3.e-3)
    
:age_dist_max:
    
    Star particles below this age are sub-divided into an age distribution if use_age_distribution is set to True
    (Units: Gyr, Default = 1.e-2)
    

:alpha_enhacement:

    If set, then the metallicity of star particles to [Fe/H] rather than the total metals.
    Since FSPS does not support non solar abundance ratios, this parameter can be used to mimic the
    hardening of the radiaiton field due to alpha-enhancement. (Default: False)

COMMON PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~

NOTE: These parmeters take three or four values as an input. 
They correspond to the value of the pararmeter for young_stars, Post-AGB stars , AGNs and Diffuse Ionized Gas (DIG) respectively.

:FORCE_gas_logu:
    
    If set, then we force the ionization parameter (gas_logu) to be gas_logu (next parameter) else, 
    it is taken to be variable and dependent on ionizing radiation from star particles. (Default: [False,False,False])

:gas_logu:
   
    Gas ionization parameter. This is only relevant if add_neb_emission is set to True and FORCE_gas_logu is set to True (Default: [-2.0,-2.0,-2.0])

:gas_logu_init:         

    Force the ionization parameter to increase/decrease by this value (Scale: log). Useful if you want to run tests (Default: [0.0,0.0,0.0])
    
:FORCE_gas_logz:
   
    If set, then we force the metallicity (gas_logz) to be gas_logz (next parameter) else, it is taken to be 
    the star particle's metallicity. (Default: [False,False,False])

:gas_logz:
   
    Gas metallicity in units of log(Z/Z_sun) only relevant if add_neb_emission = True 
    and FORCE_gas_logz = True (Default: [0.0,0.0,0.0])

:FORCE_logq:
    
    If set, then we force the number of ionizing photons to be source_logq (next parameter) 
    else, it is taken to be variable and dependent on ionizing radiation of the source. (Default: [False,False,False])

:source_logq:

    The number of ionizing photons emitted by the source in units of s^-1. Only relevant if add_neb_emission = True, 
    use_cloudy_tables = True and  FORCE_gas_logq = True (Default: [1.e47,1.e47,1.e47])
     
:FORCE_inner_radius:
    
    If set, then we force the inner radius of the cloud to be inner_radius (next parameter). This works only for young stars and Post-AGB stars. 
    IMP Note: For AGN we keep the inner radius fixed at whatever is set by inner_radius (next parameter) irrespective of what this parameter is set to. 
    (Default: [False,False,True])
    

:inner_radius:

    This sets the inner radius of the cloud in cm. This is used only when add_neb_emission = True,
    use_cloudy_tables = True and FORCE_inner_radius = True (Default: [1.e19,1.e19,2.777e+20])

:FORCE_N_O_Pilyugin:
    If set, then the Nitrogen abundances are set according to the N/O vs O/H relation from Pilyugin et al. 2012
    If FORCE_N_O ratio (next parameter) is set to True then this parameter is ignored. (Default: [False,False,False])

    
:FORCE_N_O_ratio: 
    
    If set, then we force the log of N/O ratio to be N_O_ratio (next parameter). 
    This can be used as a template fix adundance ratio of other elements (Default: False)

:N_O_ratio: 
    
    This sets the log of N/O ratio. This is used only when add_neb_emission = True,
    use_cloudy_tables = False, FORCE_N/O ratio = True and neb_abund = "direct" (Default: = -0.85)

:neb_abund:

    This sets the elemental abundances for generating CLOUDY models. 
    Available abundances are:
    
    :dopita:    
      Abundabces from Dopita (2001) with old solar abundances = 0.019 and ISM grains.
    
    :newdopita: 
      Abundances from Dopita (2013). Solar Abundances from Grevasse 2010 - z= 0.013
      includes smooth polynomial for N/O, C/O relationship functional form for He(z),
      new depletion and factors in ISM grains.
    
    :gutkin:    
      Abundances from Gutkin (2016) and PARSEC metallicity (Bressan+2012) based on 
      Grevesse+Sauvel (1998) and Caffau+2011 
    
    :direct:
      Abundances are taken directly from the simulation if possible. Defaults 
      to using "dopita" if there is an error. (Note: Works only for AGNs and star particles that are 
      added directly without binning. Make sure to set FORCE_BINNED to False)

   This is used only when add_neb_emission = True and use_cloudy_tables = True. (Default: ["dopita","dopita","dopita"])


Young Stars
~~~~~~~~~~~~~~~~~~~~~~

:add_young_stars:
    
    If set, the stars with age below HII_max_age (described below) are included when calculating nebular emission (Default: True)

:HII_Rinner_per_Rs:        		    

	Rinner for cloudy calculations is set to this value times the Stromgen Radius. For example, if set to 0.01 Rinner is taken 
	to be 1 % of Stromgren Radius. If FORCE_inner_radius (described above) is set to True for young stars then this is overridden 
	and the value set by the inner_radius (described above) is used. 
    
:HII_nh: 
   
    Gas hydrogen density for calcualting nebular emission in units if cm^-3. (Default = 1.e2)

:HII_min_age:

    Sets the minimum age limit for calculating nebular emission in units of Gyr. (Default = 1.e-3)

:HII_max_age:
   
    Sets the maximum age limit for calculating nebular emission in units of Gyr. (Default = 1.e-2)

:HII_escape_fraction:
    
    HII region escape fraction (Default = 0.0)

:HII_dust:
    
    If set, then dust grains are included in the CLOUDY model. We use grains orion command to add
    dust grains which specifies graphitic and silicate grains with a size distribution and abundance
    appropriate for those along the line of sight to the Trapezium stars in Orion (see CLOUDY documentation
    Hazy 1 for more info). (Default: False)
                                                                                                                            

Post-AGB stars
~~~~~~~~~~~~~~~~~~~~~~

:add_pagb_stars:

    If set, the Post-AGB stars are included when calculating nebular emission (Default: False)

:PAGB_N_enhancement: 

    Enhances the Nitrogen abundance Post-AGB stars by increasing the log(N/O) by this value. 
    Only relevant if add_neb_emission is set to True, use_cloudy_tables is set to False and 
    add_pagb_stars is set to True (Default = 0.4)  

:PAGB_C_enhancement:
    
    Enhances the Carbon abundance Post-AGB stars by increasing the log(C/O) by this value.
    Only relevant if add_neb_emission is set to True, use_cloudy_tables is set to False and 
    add_pagb_stars is set to True (Default = 0.4)
    
    
:PAGB_Rinner_per_Rs:        		    

	Rinner for cloudy calculations is set to this value times the Stromgen Radius. 
    For example, if set to 0.01 Rinner is taken to be 1 % of Stromgren Radius. 
    If FORCE_inner_radius (described above) is set to True for young stars then this 
    is overridden and the value set by the inner_radius (described above) is used. 
    
:PAGB_nh: 
   
    Gas hydrogen density for calcualting nebular emission in units if cm^-3. (Default = 1.e2)


:PAGB_min_age:

    Sets the minimum age limit for calculating nebular emission from post-AGB stars, in units of Gyr. 
    Only relevant if add_neb_emission is set to True, use_cloudy_tables is set to False and 
    add_pagb_stars is set to True (Default = 0.1)

:PAGB_max_age: 

    Sets the maximum age limit for calculating nebular emission from post-AGB stars, in units of Gyr.
    Only relevant if add_neb_emission is set to True, use_cloudy_tables is set to False and 
    add_pagb_stars is set to True (Default = 10)

:PAGB_escape_fraction:	    
	
	Fraction of H-ionizaing photons that escape the HII region. 
    This is used only when add_neb_emission = True and use_cloudy_tables = False and
    add_pagb_stars is set to True (Default = 0.0)

AGN
~~~~~~~~~~~~~~~~~~~~~~

:add_AGN_neb: 			            
	
	If set, AGNs are included when calculating nebular emission (Default: False)

:AGN_nh:					            
	
	Gas hydrogen density for calcualting nebular emission in units if cm^-3. (Default = 1.e3)

:AGN_num_gas:

	For CLOUDY calculations we use the distance weighted average metallicity of gas particles around the AGN. 
	The number of gas particles used for doing so is set by this parameter. 
	(Default: 32)


DIG
~~~~~~~~~~~~~~~~~~~~~~

:add_DIG_neb:
    
    If set, Contribution from DIG is included when calculating nebular emission (Default: False)

:DIG_nh:

    Gas hydrogen density for calcualting nebular emission in units of cm^-3. (Default: 10)

:DIG_min_logU:

    Only gas cells with ionization parameter greater than this are considered for DIG calculation. 
    This is done so as to speed up the calculation by ignoring the cells that do not have enough energy 
    to produce any substantial emission.  (Defualt: -6.0)

:use_black_sed:

    If set, Black (1987) ISRF is used as the input SED shape for DIG CLOUDY calculations 
    else, the input SED shape is calulated by by taking a distance weighted average of the CLOUDY 
    output spectrum of nearby young stars. The normalization of the SED is set by the total energy above the 
    lyman limit dumped in each cell. (Default: False)

:stars_max_dist:

    Only stars within this distance (Units: Kpc) are considered for getting the input spectrum shape. (Default = 1)

:max_stars_num:

     This sets the upper limit on the number of stars that are used for calculating the input spectrum shape.

DEBUGGING AND CLEAN UP
~~~~~~~~~~~~~~~~~~~~~~

:dump_emlines:        

    If True, The emission lines are saved in a file before going
    through the dust radiative transfer.  These are the cloudy
    computed emission line strengths, and are calculated for all lines
    cloudy calculates (i.e. not just those undergoing radiative
    transfer).  The format for the output is a wavelength array,
    followed by a (nlam+2) list for each nebular emission bearing particle.
    The +2 in the (nlam+2) list are the O/H ratio and the id of that particle.
    Where id = 0 , 1, 2 and 3 corresponds to young stars, PAGB stars, AGN 
    and DIG respectively.There is a convenience package in /convenience to help 
    read in this file.

    This can be used as a fast way getting emission lines for the
    purpose of debugging the code.  Naming convention:
    emlines.galaxy*.txt where * is the galaxy number. This works only
    when add_neb_emission is set to True (Default: False)

:cloudy_cleanup:
   
   If set to True, all the CLOUDY files will be deleted after the source addition is complete. 
   This is used only when add_neb_emission = True and use_cloudy_tables = True.
   Note that in case an error occurs, the files are not deleted even if this value is set to True.
   (Default: True)  
   
   
Birth Cloud Information
------------

:CF_on:

   Boolean.  If set to True, then enables the Charlot & Fall
   birthcloud models for all stars with age younger than
   birth_cloud_clearing_age.

:birth_cloud_clearing_age:

   Stars with age < birth_cloud_clearing_age have Charlot & Fall
   birthclouds (if CF_on == True).  Meaningless if CF_on == False.
   Units: Gyr.


Idealized Galaxy SED Parameters
-------------------------------

:Z_init:

   Forced metallicity increase in the newstar particles.  Useful for
   idealized galaxy simulations where the stars can form out of
   pristine gas.  Units are absolute (so 0.02 = Solar). Setting to 0
   (default) means that you use the stellar metallicities as they come
   in the simulation (i.e. for Cosmological simulations).

:disk_stars_age:

   Age in Gyr of disk stars for idealized simulations. Meaningless for
   cosmological simulations.  Note, if this is <=7, then these will
   live in Charlot & Fall birthclouds (if CF_on = True).

   Note, for Gadget simulations, stars are divided into newstars, disk
   stars and bulge stars.  For Tipsy outputs, the stars initalized
   with the simulation are auto-detected by their nonsensical ages,
   and assigned as disk stars.  So, if there are stars initalized with
   your Tipsy simulation, assign their ages (and metallicities below)
   as disk stars.

:bulge_stars_age:

   As disk_stars_age but for bulge stars.

:disk_stars_metals:

   Metallicity of disk stars in FSPS metallicity units.  See last page
   of FSPS manual for numbers.  (e.g. 20 = Solar for Padova + BaSeL
   tracks).  Meaningless for cosmological simulations.

:bulge_stars_metals:

   As disk_stars_metals but for bulge stars.


Stellar Ages and Metallicities
-------------------------------
   
:N_STELLAR_AGE_BINS:
   
   Number of bins to bin the stellar ages in (boundaries are the
   oldest and youngest star particles; linear bins in log(age)).

Black Holes
------------

:BH_SED:

    If true, `powderday <https://github.com/powderday-rt/powderday.git>`_ will 
    attempt to load black hole information from the snapshot.

:BH_eta:

    Used in calculating the black hole luminosity (bhluminosity = 
    BH_eta * mdot * c**2.)

:BH_model:

    BH model type, either Nenkova or other.

:BH_modelfile:

    The path to the Nenkova model file if BH_model is set to Nenkova. This file
    can be downloaded here and placed anywhere in the repository, as long as 
    the correct path is set in ``parameters_master``: 
    <https://www.clumpy.org/downloads/clumpy_models_201410_tvavg.hdf5>

:nenkova_params:

    Nenkova+ (2008) model parameters.


Images and SED Parameters
------------

:NTHETA:

   Number of polar angles to view galaxy at

:IMAGING:

   Must be set to ``True`` for 
   `powderday <https://github.com/powderday-rt/powderday.git>`_ to produce an image
   output file.

:filterdir:

   Directory where filter files are stored. They should be located in
   "/home/desika/powderday/filters/".

:filterfiles:

   A list of the names of all filters to be used. 
   `powderday <https://github.com/powderday-rt/powderday.git>`_ will run at each 
   wavelength in all specified filter files, and will produce a ``.hdf5`` file
   containing images convolved with each filter transmission function. Note 
   that this can be quite computationally intensive and scales with the number
   of wavelengths. Following the example in ``parameters_master``, additional 
   filters can be added to this list. In bash, ``cd`` into your ``filterdir`` 
   and use the following command to format the filenames for easy copying and 
   pasting into this list.
   
    .. code-block:: bash

       >>> shopt -s globstar; printf "#    '%s'\n" *.filter

:IMAGING_TRANSMISSION_FILTER:

   If enabled, filter convolution will be performed through 
   `Hyperion <http://www.hyperion-rt.org>`_ instead of through `powderday 
   <https://github.com/powderday-rt/powderday.git>`_. This is much faster, but is 
   still an experimental feature and does not seem to produce accurate 
   convolved images.


DEBUGGING
------------

You should probably never touch any of these.


parameters_model
============

:snapshot_name:

   String - currently the snapshot name of your galaxy run. (Naming
   will change as other front ends built).

:hydro_dir:

   Location of snapshots

:PD_output_dir:

   String - location of where `powderday
   <https://github.com/powderday-rt/powderday.git>`_ output files should go.

:Auto_TF_file:

   String - name of the TF logical file to be written (doesn't need a
   path - will go into PD_output_dir)

:Auto_dustdens_file:

   String - name of the dust density ascii file to be written (doesn't
   need a path - will go into PD_output_dir)

:inputfile:

   String - name of the input HDF5 (rtin) file for `powderday
   <https://github.com/powderday-rt/powderday.git>`_ to write before
   radiative transfer begins.

:outputfile:

   String - name of the output HDF5 (rtout) file after radiative transfer

:x_cent:

   Location in grid coordinates of the x-coordinate of the center of
   your galaxy.  Only pertinenet if MANUAL_CENTERING==True.  Otherwise
   ignored by `powderday <https://github.com/powderday-rt/powderday.git>`_.

:y_cent:

   As x_cent but for the y-coordinate

:z_cent:

   As x_cent but for the z-coordinate.
