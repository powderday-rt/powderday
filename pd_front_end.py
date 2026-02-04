#!/usr/bin/env python
# coding: utf-8
#Code:  pd_front_end.py

# =========================================================
# IMPORT STATEMENTS
# =========================================================
from __future__ import print_function
from powderday.front_end_tools import make_SED, make_image, make_DIG_SED,compute_ISRF_SED
from powderday.source_creation import direct_add_stars, add_binned_seds, BH_source_add, DIG_source_add
from powderday.analytics import stellar_sed_write, dump_data, SKIRT_data_dump, logu_diagnostic,dump_emlines,dump_NEB_SEDs
from astropy import constants
import fsps
from powderday.image_processing import add_transmission_filters, convolve
from powderday.m_control_tools import m_control_sph, m_control_enzo, m_control_arepo
import powderday.backwards_compatibility as bc
import powderday.error_handling as eh
import powderday.SED_gen as sg
from powderday.front_ends.front_end_controller import stream
from astropy import units as u
import powderday.config as cfg
import h5py
import matplotlib as mpl
import copy
import numpy as np
import sys
import gc
import random
import os

gc.set_threshold(0)

if (len(sys.argv)) == 4:
    script, pardir, parfile, modelfile = sys.argv

elif (len(sys.argv)) == 5:
    script, pardir, parfile, neb_parfile, modelfile = sys.argv


mpl.use('Agg')


sys.path.insert(0, pardir)
par = __import__(parfile)
model = __import__(modelfile)


cfg.par = par  # re-write cfg.par for all modules that read this in now
cfg.model = model
neb_param_file = False

if cfg.par.add_neb_emission and not cfg.par.use_cloudy_tables: 
    # Checking to see if additional parameter file for nebular emission is provided or not when neublar emission is set to be 
    # calculated without using lookup tables. If it is not provided then a warning is thrown.
    if len(sys.argv) == 5:
        neb_par = __import__(neb_parfile)
        cfg.neb_par = neb_par
        neb_param_file = True
    else:
        print("#-----------------------------------------------------------------------------------------------------#")
        print("Warning Nebular Emission is turned on but parmeters master file for nebular emission is not provided.")
        print("Reverting to using default values for including nebular emission. See documentation for more info")
        print("#------------------------------------------------------------------------------------------------------#")
    
# =========================================================
# CHECK FOR THE EXISTENCE OF A FEW CRUCIAL FILES FIRST
# =========================================================
eh.file_exist(model.hydro_dir+model.snapshot_name)
eh.file_exist(par.dustdir+par.dustfile)

# =========================================================
# Enforce Backwards Compatibility for Non-Critical Variables
# =========================================================

cfg.par.FORCE_RANDOM_SEED, cfg.par.FORCE_BINNED, cfg.par.max_age_direct, cfg.par.imf1, cfg.par.imf2, cfg.par.imf3, cfg.par.use_cmdf, cfg.par.use_cloudy_tables, cfg.par.cmdf_min_mass, cfg.par.cmdf_max_mass, cfg.par.cmdf_bins, cfg.par.cmdf_beta, cfg.par.use_age_distribution, cfg.par.age_dist_min, cfg.par.age_dist_max, cfg.par.FORCE_gas_logu, cfg.par.gas_logu, cfg.par.gas_logu_init, cfg.par.FORCE_gas_logz, cfg.par.gas_logz, cfg.par.FORCE_logq, cfg.par.source_logq, cfg.par.FORCE_inner_radius, cfg.par.inner_radius, cfg.par.FORCE_N_O_Pilyugin, cfg.par.FORCE_N_O_ratio, cfg.par.N_O_ratio, cfg.par.neb_abund, cfg.par.add_young_stars, cfg.par.HII_Rinner_per_Rs, cfg.par.HII_nh, cfg.par.HII_min_age, cfg.par.HII_max_age, cfg.par.HII_dust, cfg.par.HII_escape_fraction, cfg.par.alpha_enhance, cfg.par.add_pagb_stars, cfg.par.PAGB_min_age, cfg.par.PAGB_max_age, cfg.par.PAGB_N_enhancement, cfg.par.PAGB_C_enhancement, cfg.par.PAGB_Rinner_per_Rs, cfg.par.PAGB_nh, cfg.par.PAGB_escape_fraction, cfg.par.add_AGN_neb, cfg.par.AGN_nh, cfg.par.AGN_num_gas, cfg.par.dump_emlines, cfg.par.cloudy_cleanup, cfg.par.BH_SED, cfg.par.IMAGING, cfg.par.SED, cfg.par.IMAGING_TRANSMISSION_FILTER, cfg.par.SED_MONOCHROMATIC, cfg.par.SKIP_RT, cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS, cfg.par.n_MPI_processes, cfg.par.SOURCES_RANDOM_POSITIONS, cfg.par.SUBLIMATION, cfg.par.SUBLIMATION_TEMPERATURE, cfg.model.TCMB, cfg.model.THETA, cfg.model.PHI, cfg.par.MANUAL_ORIENTATION, cfg.par.dust_grid_type, cfg.par.BH_model, cfg.par.BH_modelfile, cfg.par.BH_var, cfg.par.FORCE_STELLAR_AGES,cfg.par.FORCE_STELLAR_AGES_VALUE, cfg.par.FORCE_STELLAR_METALLICITIES, cfg.par.FORCE_STELLAR_METALLICITIES_VALUE, cfg.par.NEB_DEBUG, cfg.par.filterdir, cfg.par.filterfiles,  cfg.par.PAH_frac, cfg.par.otf_extinction, cfg.par.explicit_pah, cfg.par.draine21_pah_model, cfg.par.dust_density, cfg.par.add_DIG_neb, cfg.par.DIG_nh, cfg.par.DIG_min_logU, cfg.par.stars_max_dist, cfg.par.max_stars_num, cfg.par.use_black_sed, cfg.par.n_photons_DIG, cfg.par.SKIRT_DATA_DUMP, cfg.par.SAVE_NEB_SEDS, cfg.par.REMOVE_INPUT_SEDS, cfg.par.OTF_EXTINCTION_MRN_FORCE, cfg.par.separate_into_dust_species, cfg.par.OTF_EXTINCTION_MRN_FORCE = bc.variable_set()

# If a seperate parameter file is provided for nebular emission then overwrite the relevant variables based on that.
if neb_param_file:
    cfg.par.use_cmdf, cfg.par.cmdf_min_mass, cfg.par.cmdf_max_mass, cfg.par.cmdf_bins, cfg.par.cmdf_beta, cfg.par.use_age_distribution, cfg.par.age_dist_min, cfg.par.age_dist_max, cfg.par.FORCE_gas_logu, cfg.par.gas_logu, cfg.par.gas_logu_init, cfg.par.FORCE_gas_logz, cfg.par.gas_logz, cfg.par.FORCE_logq, cfg.par.source_logq, cfg.par.FORCE_inner_radius, cfg.par.inner_radius, cfg.par.FORCE_N_O_Pilyugin, cfg.par.FORCE_N_O_ratio, cfg.par.N_O_ratio, cfg.par.neb_abund, cfg.par.add_young_stars, cfg.par.HII_Rinner_per_Rs, cfg.par.HII_nh, cfg.par.HII_min_age, cfg.par.HII_max_age, cfg.par.HII_dust, cfg.par.HII_escape_fraction, cfg.par.add_pagb_stars, cfg.par.PAGB_min_age, cfg.par.PAGB_max_age, cfg.par.PAGB_N_enhancement, cfg.par.PAGB_C_enhancement, cfg.par.PAGB_Rinner_per_Rs, cfg.par.PAGB_nh, cfg.par.PAGB_escape_fraction, cfg.par.add_AGN_neb, cfg.par.AGN_nh, cfg.par.AGN_num_gas, cfg.par.dump_emlines, cfg.par.cloudy_cleanup, cfg.par.NEB_DEBUG, cfg.par.add_DIG_neb, cfg.par.DIG_nh, cfg.par.DIG_min_logU, cfg.par.stars_max_dist, cfg.par.max_stars_num, cfg.par.SAVE_NEB_SEDS = cfg.neb_par.use_cmdf, cfg.neb_par.cmdf_min_mass, cfg.neb_par.cmdf_max_mass, cfg.par.cmdf_bins, cfg.neb_par.cmdf_beta, cfg.neb_par.use_age_distribution, cfg.neb_par.age_dist_min, cfg.neb_par.age_dist_max, cfg.neb_par.FORCE_gas_logu, cfg.neb_par.gas_logu, cfg.neb_par.gas_logu_init, cfg.neb_par.FORCE_gas_logz, cfg.neb_par.gas_logz, cfg.neb_par.FORCE_logq, cfg.neb_par.source_logq, cfg.neb_par.FORCE_inner_radius, cfg.neb_par.inner_radius, cfg.neb_par.FORCE_N_O_Pilyugin, cfg.neb_par.FORCE_N_O_ratio, cfg.neb_par.N_O_ratio, cfg.neb_par.neb_abund, cfg.neb_par.add_young_stars, cfg.neb_par.HII_Rinner_per_Rs, cfg.neb_par.HII_nh, cfg.neb_par.HII_min_age, cfg.neb_par.HII_max_age, cfg.neb_par.HII_dust, cfg.neb_par.HII_escape_fraction, cfg.neb_par.add_pagb_stars, cfg.neb_par.PAGB_min_age, cfg.neb_par.PAGB_max_age, cfg.neb_par.PAGB_N_enhancement, cfg.neb_par.PAGB_C_enhancement, cfg.neb_par.PAGB_Rinner_per_Rs, cfg.neb_par.PAGB_nh, cfg.neb_par.PAGB_escape_fraction, cfg.neb_par.add_AGN_neb, cfg.neb_par.AGN_nh, cfg.neb_par.AGN_num_gas, cfg.neb_par.dump_emlines, cfg.neb_par.cloudy_cleanup, cfg.neb_par.NEB_DEBUG, cfg.neb_par.add_DIG_neb, cfg.neb_par.DIG_nh, cfg.neb_par.DIG_min_logU, cfg.neb_par.stars_max_dist, cfg.neb_par.max_stars_num, cfg.neb_par.SAVE_NEB_SEDS

# =========================================================
# CHECK FOR COMPATIBLE PARAMETERS
# =========================================================
eh.check_parameter_compatibility()

# =========================================================
# GRIDDING
# =========================================================

fname = cfg.model.hydro_dir+cfg.model.snapshot_name
field_add, ds = stream(fname)

# figure out which tributary we're going to
ds_type = ds.dataset_type

from powderday.pah.pah_source_create import pah_source_add


sp = fsps.StellarPopulation()

#setting solar metallicity value based on isochrone
print(f'\n----------------------------------------------\nSetting solar metallicity value')

isochrone = str(sp.libraries[0].decode())

# As of commit #329774874cf40c04368ae300677576e3cae2c369 (August 11, 2022, https://github.com/dfm/python-fsps) 
# isochrone solar metallicity can now be accessed through the 'solar_metallicity' attribute of the sp object.
try:
    print('isochrone = ', isochrone)
    cfg.par.solar = sp.solar_metallicity
    print(f'solar metallicity = {cfg.par.solar}')

except:
    print("\nWARNING: Please update your python-fsps to at least commit #329774874cf40c04368ae300677576e3cae2c369 (August 11, 2022)")
    print("Powderday will no longer work with the older versions of python-fsps in the near future\n")
    Ziso = {'mist': ('mist', 0.0142), 'bsti': ('basti', 0.020),
            'gnva': ('geneva', 0.020), 'prsc': ('parsec', 0.01524),
            'pdva': ('padova', 0.019), 'bpss': ('bpass', 0.20)}
    iso, Zsun = Ziso[isochrone]
    print('isochrone = '+ str(iso))
    cfg.par.solar = Zsun
    print(f'solar metallicity = {cfg.par.solar}')
print('----------------------------------------------')


# define the options dictionary
options = {'gadget_hdf5': m_control_sph,
           'tipsy': m_control_sph,
           'enzo_packed_3d': m_control_enzo,
           'arepo_hdf5': m_control_arepo}

m_gen = options[ds_type]()
m, xcent, ycent, zcent, dx, dy, dz, reg, ds, boost = m_gen(fname, field_add)


#save the dataset_type for future use in reg
try: reg.parameters['dataset_type'] = ds.dataset_type
except AttributeError:
    reg.parameters={}
    reg.parameters['dataset_type'] = ds.dataset_type


# Get dust wavelengths. This needs to preceed the generation of sources
# for hyperion since the wavelengths of the SEDs need to fit in the
# dust opacities.
df = h5py.File(cfg.par.dustdir+cfg.par.dustfile, 'r')
o = df['optical_properties']
df_nu = o['nu']
df_chi = o['chi']
df.close()


# add sources to hyperion
stars_list, diskstars_list, bulgestars_list, reg = sg.star_list_gen(boost, dx, dy, dz, reg, ds, sp, m)
nstars = len(stars_list)

# figure out N_METAL_BINS:
fsps_metals = np.array(sp.zlegend)
N_METAL_BINS = len(fsps_metals)


#initializing the nebular diagnostic file newly
if cfg.par.add_neb_emission and cfg.par.NEB_DEBUG: logu_diagnostic(None,None,None,None,None,None,None,append=False)
if cfg.par.add_neb_emission and cfg.par.dump_emlines: dump_emlines(None,append=False)
if cfg.par.add_neb_emission and (cfg.par.SAVE_NEB_SEDS or cfg.par.add_DIG_neb): dump_NEB_SEDs(None, None, None, append=False)

if cfg.par.BH_SED == True:
    BH_source_add(m, reg, df_nu, boost)


if cfg.par.FORCE_BINNED == False:
    m = direct_add_stars(df_nu, stars_list, diskstars_list, bulgestars_list, ds.cosmological_simulation, m, sp)

# note - the generation of the SEDs is called within
# add_binned_seds itself, unlike add_newstars, which requires
# that sg.allstars_sed_gen() be called first.
m = add_binned_seds(df_nu, stars_list, diskstars_list,bulgestars_list, ds.cosmological_simulation, m, sp)



#set the random seeds
if cfg.par.FORCE_RANDOM_SEED == False:
    m.set_seed(random.randrange(0,10000)*-1)
else:
    m.set_seed(cfg.par.seed)

# save SEDs
# stars and black holes can't both be in the sim and write stellar SEDs to a file becuase they have different wavelength sizes
if (par.STELLAR_SED_WRITE == True) and not (par.BH_SED) and not (par.draine21_pah_model):
    stellar_sed_write(m)


if ds_type in ['gadget_hdf5','tipsy','arepo_hdf5'] and cfg.par.SKIRT_DATA_DUMP:
    SKIRT_data_dump(reg, ds, m, stars_list, bulgestars_list, diskstars_list, ds_type, sp)

nstars = len(stars_list)
nstars_disk = len(diskstars_list)
nstars_bulge = len(bulgestars_list)


'''
#EXPERIMENTAL FEATURES
if par.SOURCES_IN_CENTER == True:
    for i in range(nstars):
        stars_list[i].positions[:] =  np.array([xcent,ycent,zcent])
    for i in range(nstars_bulge):
        bulgestars_list[i].positions[:] =  np.array([xcent,ycent,zcent])
    for i in range(nstars_disk):
        diskstars_list[i].positions[:] = np.array([xcent,ycent,zcent])
if par.SOURCES_RANDOM_POSITIONS == True:
    print "================================"
    print "SETTING SOURCES TO RANDOM POSITIONS"
    print "================================"
    for i in range(nstars):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        stars_list[i].positions[:] = np.array([xpos,ypos,zpos])
    for i in range(nstars_bulge):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        bulgestars_list[i].positions[:] = np.array([xpos,ypos,zpos])
    for i in range(nstars_disk):
        xpos,ypos,zpos = np.random.uniform(-dx,dx),np.random.uniform(-dy,dy),np.random.uniform(-dz,dz)
        diskstars_list[i].positions[:] = np.array([xpos,ypos,zpos])
'''





# set up the CMB field -- place holder to put in haardt/madau eventually
'''
cmb = m.add_external_box_source()
cmb.temperature = cfg.model.TCMB
cmb_box_len = ds.quan(cfg.par.zoom_box_len,'kpc').in_units('cm').value
cmb.bounds = [[-cmb_box_len,cmb_box_len],[-cmb_box_len,cmb_box_len],[-cmb_box_len,cmb_box_len]]
pdb.set_trace()
L_CMB = (constants.sigma_sb*(cfg.model.TCMB*u.K)**4.).to(u.erg/u.cm**2/u.s)*4*(cmb_box_len*u.cm)**2 #get_J_CMB()
cmb.luminosity = L_CMB.cgs.value
'''

'''
energy_density_absorbed=energy_density_absorbed_by_CMB()
m.add_density_grid(density, dust, specific_energy=energy_density_absorbed)
m.set_specific_energy_type('additional')
'''

print('Done adding Sources')


print('Setting up Model')
m_imaging = copy.deepcopy(m)
m.conf.output.output_specific_energy = 'last'


print("Dumping grid information")



#if ds_type in ['gadget_hdf5','tipsy','arepo_hdf5']:
#    dump_data(reg, model)


if cfg.par.add_neb_emission and cfg.par.add_DIG_neb:
    make_DIG_SED(m, par, model)
    DIG_source_add(m, reg, df_nu,boost)
    print ("Removing the DIG energy dumped input SED file")
    os.remove(cfg.model.inputfile + '_DIG_energy_dumped.sed')
    os.remove(cfg.model.outputfile + '_DIG_energy_dumped.sed')
    if not cfg.par.SAVE_NEB_SEDS: dump_NEB_SEDs(None, None, None, append=False, clean_up=True)


if cfg.par.otf_extinction and cfg.par.draine21_pah_model:
    m.compute_isrf(True)
    compute_ISRF_SED(m, par, model)
    pah_source_add(ds,reg,m,boost)


if ds_type in ['gadget_hdf5','tipsy','arepo_hdf5']:
    dump_data(reg, model)

if cfg.par.SED:
    make_SED(m, par, model)
    if cfg.par.REMOVE_INPUT_SEDS:
        print ("Removing the input SED file")
        os.remove(cfg.model.inputfile+'.sed')

if cfg.par.IMAGING:
    make_image(m_imaging, par, model, dx, dy, dz)
