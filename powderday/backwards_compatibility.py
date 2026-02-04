import powderday.config as cfg
import os.path
import numpy as np


def variable_set():
    try:
        cfg.par.FORCE_RANDOM_SEED
    except:
        cfg.par.FORCE_RANDOM_SEED  = None

    try:
        cfg.par.FORCE_BINNED
    except:
        cfg.par.FORCE_BINNED = True
        
    try:
        cfg.par.max_age_direct
    except:
        cfg.par.max_age_direct = 1.e-2
        
    try:
        cfg.par.imf1
    except:
        cfg.par.imf1 = 1.3

    try:
        cfg.par.imf2
    except:
        cfg.par.imf2 = 2.3

    try:
        cfg.par.imf3
    except:
        cfg.par.imf3 = 2.3
        
    try:
        cfg.par.use_cmdf
    except:
        cfg.par.use_cmdf = False
      
    try:
        cfg.par.use_cloudy_tables
    except:
        cfg.par.use_cloudy_tables = False
     
    try:
        cfg.par.cmdf_min_mass
    except:
        cfg.par.cmdf_min_mass = 3.5
    
    try:
        cfg.par.cmdf_max_mass
    except:
        cfg.par.cmdf_max_mass = 5.0

    try:
        cfg.par.cmdf_bins
    except:
        cfg.par.cmdf_bins = 6
    
    try:
        cfg.par.cmdf_beta
    except:
        cfg.par.cmdf_beta = -2.0
    
    try:
        cfg.par.use_age_distribution
    except:
        cfg.par.use_age_distribution = True

    try:
        cfg.par.age_dist_min
    except:
        cfg.par.age_dist_min = 3e-3

    try:
        cfg.par.age_dist_max
    except:
        cfg.par.age_dist_max = 1e-2

    try:
        if len(np.atleast_1d(cfg.par.FORCE_gas_logu)) == 1:
            cfg.par.FORCE_gas_logu = [cfg.par.FORCE_gas_logu,cfg.par.FORCE_gas_logu,cfg.par.FORCE_gas_logu]
    except:
        cfg.par.FORCE_gas_logu = [False, False, False]

    try:
        if len(np.atleast_1d(cfg.par.gas_logu)) == 1:
            cfg.par.gas_logu = [cfg.par.gas_logu,cfg.par.gas_logu,cfg.par.gas_logu]
    except:
        cfg.par.gas_logu = [-2.0,-2.0,-2.0]

    try:
        if len(np.atleast_1d(cfg.par.gas_logu_init)) == 1:
            cfg.par.gas_logu_init = [cfg.par.gas_logu_init,cfg.par.gas_logu_init,cfg.par.gas_logu_init]
    except:
        cfg.par.gas_logu_init = [0.0,0.0,0.0]
    
    try:
        if len(np.atleast_1d(cfg.par.FORCE_gas_logz)) == 1:
            cfg.par.FORCE_gas_logz = [cfg.par.FORCE_gas_logz,cfg.par.FORCE_gas_logz,cfg.par.FORCE_gas_logz]
    except:
        cfg.par.FORCE_gas_logz = [False,False,False]
    
    try:
        if len(np.atleast_1d(cfg.par.gas_logz)) == 1:
            cfg.par.gas_logz = [cfg.par.gas_logz,cfg.par.gas_logz,cfg.par.gas_logz]
    except:
        cfg.par.gas_logz = [0.0,0.0,0.0]

    try:
        if len(np.atleast_1d(cfg.par.FORCE_logz)) == 1:
            cfg.par.FORCE_logz = [cfg.par.FORCE_logz,cfg.par.FORCE_logz,cfg.par.FORCE_logz]
    except:
        cfg.par.FORCE_logq = [False,False,False]
    
    try:
        if len(np.atleast_1d(cfg.par.source_logq)) == 1:
            cfg.par.source_logq = [cfg.par.source_logq,cfg.par.source_logq,cfg.par.source_logq]

    except:
        cfg.par.source_logq = [1.e47,1.e47,1.e47]

    
    try:
        if len(np.atleast_1d(cfg.par.FORCE_inner_radius)) == 1:
            cfg.par.FORCE_inner_radius = [cfg.par.FORCE_inner_radius,cfg.par.FORCE_inner_radius,True]

    except:
        cfg.par.FORCE_inner_radius = [False, False, True]
    
    try:
        if len(np.atleast_1d(cfg.par.inner_radius)) == 1:
            cfg.par.inner_radius = [cfg.par.inner_radius,cfg.par.inner_radius,2.777e+20]
    except:
        cfg.par.inner_radius = [1.e19,1.e19,2.777e+20]

    try:
        if len(np.atleast_1d(cfg.par.FORCE_N_O_Pilyugin)) == 1:
            cfg.par.FORCE_N_O_Pilyugin = [cfg.par.FORCE_N_O_Pilyugin,cfg.par.FORCE_N_O_Pilyugin,cfg.par.FORCE_N_O_Pilyugin, cfg.par.FORCE_N_O_Pilyugin]
    except:
            cfg.par.FORCE_N_O_Pilyugin = [False, False, False, False]
    
    try:
        if len(np.atleast_1d(cfg.par.FORCE_N_O_ratio)) == 1:
            cfg.par.FORCE_N_O_ratio = [cfg.par.FORCE_N_O_ratio,cfg.par.FORCE_N_O_ratio,cfg.par.FORCE_N_O_ratio, cfg.par.FORCE_N_O_ratio]
    except:
        cfg.par.FORCE_N_O_ratio = [False, False, False, False]
        
    try:
        if len(np.atleast_1d(cfg.par.N_O_ratio)) == 1:
            cfg.par.N_O_ratio = [cfg.par.N_O_ratio,cfg.par.N_O_ratio,cfg.par.N_O_ratio,cfg.par.N_O_ratio]
    except:
        cfg.par.N_O_ratio = [-0.85,-0.85,-0.85, -0.85]


    try:
        if len(np.atleast_1d(cfg.par.neb_abund)) == 1:
            cfg.par.neb_abund = [cfg.par.neb_abund,cfg.par.neb_abund,cfg.par.neb_abund,cfg.par.neb_abund]
    except:
        cfg.par.neb_abund = ["dopita","dopita","dopita","dopita"]
    
    try:
        cfg.par.add_young_stars
    except:
        cfg.par.add_young_stars = True
    
    try:
        cfg.par.HII_Rinner_per_Rs
    except:
        try:
            cfg.par.HII_Rinner_per_Rs = cfg.par.Rinner_per_Rs
        except:
            cfg.par.HII_Rinner_per_Rs = 0.01

    try: 
        cfg.par.HII_nh
    except:
        cfg.par.HII_nh = 1.e2

    try:
        cfg.par.HII_min_age
    except:
        cfg.par.HII_min_age = 1.e-3
	
    try:
        cfg.par.HII_max_age
    except: 
        cfg.par.HII_max_age = 1.e-2
        
    
    try:
        cfg.par.HII_dust
    except:
        cfg.par.HII_dust = False
    
    try:
        cfg.par.HII_escape_fraction
    except:
        cfg.par.HII_escape_fraction = 0.0
        
    try:
        cfg.par.alpha_enhance
    except:
        try:
            cfg.par.HII_alpha_enhance
            cfg.par.alpha_enhance = cfg.par.HII_alpha_enhance
        except:
            cfg.par.alpha_enhance = False

    try:
        cfg.par.add_pagb_stars
    except:
        cfg.par.add_pagb_stars = False

    try:
        cfg.par.PAGB_min_age
    except:
        cfg.par.PAGB_min_age = 1.e-1

    try:
        cfg.par.PAGB_max_age
    except:
        cfg.par.PAGB_max_age = 14.

    try:
        cfg.par.PAGB_N_enhancement
    except:
        cfg.par.PAGB_N_enhancement = 0.4

    try:
        cfg.par.PAGB_C_enhancement
    except:
        cfg.par.PAGB_C_enhancement = 0.4
    
    try: 
        cfg.par.PAGB_Rinner_per_Rs
    except:
        try:
            cfg.par.PAGB_Rinner_per_Rs = cfg.par.Rinner_per_Rs
        except:
            cfg.par.PAGB_Rinner_per_Rs = 0.01

    try:
        cfg.par.PAGB_nh
    except:
        try:
            cfg.par.PAGB_nh = cfg.par.HII_nh
        except:
            cfg.par.PAGB_nh = 1.e2


    try:
        cfg.par.PAGB_escape_fraction
    except:
        try:
            cfg.par.PAGB_escape_fraction = cfg.par.HII_escape_fraction
        except:
            cfg.par.PAGB_escape_fraction = 0.0


    try:
        cfg.par.add_AGN_neb
    except:
        cfg.par.add_AGN_neb = False
        
    try:
        cfg.par.AGN_nh
    except:
        cfg.par.AGN_nh = 1.e3
        
    try:
        cfg.par.AGN_num_gas
    except:
        cfg.par.AGN_num_gas = 32    

    try:
        cfg.par.dump_emlines
    except:
        cfg.par.dump_emlines = False

    try:
        cfg.par.cloudy_cleanup
    except:
        cfg.par.cloudy_cleanup = True
        
    try:
        cfg.par.BH_SED
    except:
        cfg.par.BH_SED  = None

    try:
        cfg.par.IMAGING
    except:
        cfg.par.IMAGING  = None

    try:
        cfg.par.SED
    except:
        cfg.par.SED = True

    try:
        cfg.par.IMAGING_TRANSMISSION_FILTER
    except:
        cfg.par.IMAGING_TRANSMISSION_FILTER = False

    try:
        cfg.par.SED_MONOCHROMATIC
    except:
        cfg.par.SED_MONOCHROMATIC = False

    try:
        cfg.par.SKIP_RT
    except:
        cfg.par.SKIP_RT = False

    try:
        cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS 
    except:
        cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS = True

    try:
        cfg.par.n_MPI_processes
    except:
        cfg.par.n_MPI_processes = cfg.par.n_processes #default is to make the same as the number of pool processes

    try:
        cfg.par.SOURCES_RANDOM_POSITIONS
    except:
        cfg.par.SOURCES_RANDOM_POSITIONS = False

    try:
        cfg.par.SUBLIMATION
    except:
        cfg.par.SUBLIMATION = False

    try:
        cfg.par.SUBLIMATION_TEMPERATURE
    except:
        cfg.par.SUBLIMATION_TEMPERATURE = 1600

    try:
        cfg.model.TCMB
    except:
        cfg.model.TCMB = 2.73

    try:
        cfg.model.THETA
    except:
        cfg.model.THETA = 0

    try:
        cfg.model.PHI
    except:
        cfg.model.PHI = 0

    try:
        cfg.par.MANUAL_ORIENTATION
    except:
        cfg.par.MANUAL_ORIENTATION = False

    try:
        cfg.par.dust_grid_type
    except:
        cfg.par.dust_grid_type = 'dtm'

    try:
        cfg.par.BH_model
    except:
        cfg.par.BH_model = None

    try:
        cfg.par.BH_modelfile
    except:
        cfg.par.BH_modelfile = None

    try:
        cfg.par.BH_var
    except:
        cfg.par.BH_var = False

    try:
        cfg.par.FORCE_STELLAR_AGES
    except:
        cfg.par.FORCE_STELLAR_AGES = False

    try:
        cfg.par.FORCE_STELLAR_AGES_VALUE
    except:
        cfg.par.FORCE_STELLAR_AGES_VALUE = 1

    try:
        cfg.par.FORCE_STELLAR_METALLICITIES
    except:
        cfg.par.FORCE_STELLAR_METALLICITIES = False

    try:
        cfg.par.FORCE_STELLAR_METALLICITIES_VALUE
    except:
        cfg.par.FORCE_STELLAR_METALLICITIES_VALUE = 0.013
    
    try:
        cfg.par.NEB_DEBUG
    except:
        cfg.par.NEB_DEBUG = True

    try:
        cfg.par.filterdir
    except:
        try:
            cfg.par.filterdir = os.path.dirname(cfg.par.filter_file)+'/'
        except:
            cfg.par.filterdir = ''

    try:
        cfg.par.filterfiles
    except:
        cfg.par.filterfiles = ['pdfilters.dat']

    try:
        cfg.par.PAH_frac
    except:
        cfg.par.PAH_frac = {'usg': 0.0586, 'vsg': 0.1351, 'big': 0.8063}

    try:
        cfg.par.otf_extinction
    except:
        cfg.par.otf_extinction = False

    try:
        cfg.par.explicit_pah
    except:
        cfg.par.explicit_pah = False

    try:
        cfg.par.draine21_pah_model
    except:
        cfg.par.draine21_pah_model = False

    try:
        cfg.par.draine21_pah_model
    except:
        cfg.par.draine21_pah_model = False

    try:
        cfg.par.dust_density
    except:
        cfg.par.dust_density = 2.4
        
    try:
        cfg.par.add_DIG_neb 
    except:
        cfg.par.add_DIG_neb = False

    try:
        cfg.par.DIG_nh
    except:
        cfg.par.DIG_nh = 1e1

    try:
        cfg.par.DIG_min_logU
    except:
        cfg.par.DIG_min_logU = -6.0

    try:
        cfg.par.stars_max_dist
    except:
        cfg.par.stars_max_dist = 1.0


    try:
        cfg.par.max_stars_num
    except:
        cfg.par.max_stars_num = 20

    try:
        cfg.par.use_black_sed
    except:
        cfg.par.use_black_sed = False

    try:
        cfg.par.n_photons_DIG
    except:
        cfg.par.n_photons_DIG = 1.e8

    try:
        cfg.par.SKIRT_DATA_DUMP
    except:
        cfg.par.SKIRT_DATA_DUMP = True

    try:
        cfg.par.SAVE_NEB_SEDS
    except:
        cfg.par.SAVE_NEB_SEDS = False

    try:
        cfg.par.REMOVE_INPUT_SEDS
    except:
        cfg.par.REMOVE_INPUT_SEDS = False

    try:
        cfg.par.OTF_EXTINCTION_MRN_FORCE
    except:
        OTF_EXTINCTION_MRN_FORCE = False
        


    try:
        cfg.par.separate_into_dust_species
    except:
        cfg.par.separate_into_dust_species = True


    try:
        cfg.par.OTF_EXTINCTION_MRN_FORCE
    except:
        cfg.par.OTF_EXTINCTION_MRN_FORCE = False

        
    return cfg.par.FORCE_RANDOM_SEED, cfg.par.FORCE_BINNED, cfg.par.max_age_direct, cfg.par.imf1, cfg.par.imf2, cfg.par.imf3, cfg.par.use_cmdf, cfg.par.use_cloudy_tables, cfg.par.cmdf_min_mass, cfg.par.cmdf_max_mass, cfg.par.cmdf_bins, cfg.par.cmdf_beta, cfg.par.use_age_distribution, cfg.par.age_dist_min, cfg.par.age_dist_max, cfg.par.FORCE_gas_logu, cfg.par.gas_logu, cfg.par.gas_logu_init, cfg.par.FORCE_gas_logz, cfg.par.gas_logz, cfg.par.FORCE_logq, cfg.par.source_logq, cfg.par.FORCE_inner_radius, cfg.par.inner_radius, cfg.par.FORCE_N_O_Pilyugin, cfg.par.FORCE_N_O_ratio, cfg.par.N_O_ratio, cfg.par.neb_abund, cfg.par.add_young_stars, cfg.par.HII_Rinner_per_Rs, cfg.par.HII_nh, cfg.par.HII_min_age, cfg.par.HII_max_age, cfg.par.HII_dust, cfg.par.HII_escape_fraction, cfg.par.alpha_enhance, cfg.par.add_pagb_stars, cfg.par.PAGB_min_age, cfg.par.PAGB_max_age, cfg.par.PAGB_N_enhancement, cfg.par.PAGB_C_enhancement, cfg.par.PAGB_Rinner_per_Rs, cfg.par.PAGB_nh, cfg.par.PAGB_escape_fraction, cfg.par.add_AGN_neb, cfg.par.AGN_nh, cfg.par.AGN_num_gas, cfg.par.dump_emlines, cfg.par.cloudy_cleanup, cfg.par.BH_SED, cfg.par.IMAGING, cfg.par.SED, cfg.par.IMAGING_TRANSMISSION_FILTER, cfg.par.SED_MONOCHROMATIC, cfg.par.SKIP_RT, cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS, cfg.par.n_MPI_processes, cfg.par.SOURCES_RANDOM_POSITIONS, cfg.par.SUBLIMATION, cfg.par.SUBLIMATION_TEMPERATURE, cfg.model.TCMB, cfg.model.THETA, cfg.model.PHI, cfg.par.MANUAL_ORIENTATION, cfg.par.dust_grid_type, cfg.par.BH_model, cfg.par.BH_modelfile, cfg.par.BH_var, cfg.par.FORCE_STELLAR_AGES,cfg.par.FORCE_STELLAR_AGES_VALUE, cfg.par.FORCE_STELLAR_METALLICITIES, cfg.par.FORCE_STELLAR_METALLICITIES_VALUE, cfg.par.NEB_DEBUG, cfg.par.filterdir, cfg.par.filterfiles,  cfg.par.PAH_frac,cfg.par.otf_extinction, cfg.par.explicit_pah, cfg.par.draine21_pah_model, cfg.par.dust_density, cfg.par.add_DIG_neb, cfg.par.DIG_nh, cfg.par.DIG_min_logU, cfg.par.stars_max_dist, cfg.par.max_stars_num, cfg.par.use_black_sed, cfg.par.n_photons_DIG, cfg.par.SKIRT_DATA_DUMP, cfg.par.SAVE_NEB_SEDS, cfg.par.REMOVE_INPUT_SEDS,cfg.par.OTF_EXTINCTION_MRN_FORCE, cfg.par.separate_into_dust_species, cfg.par.OTF_EXTINCTION_MRN_FORCE
