from __future__ import print_function
import numpy as np
import powderday.config as cfg
import pdb

import astropy.units as u
import astropy.constants as constants
from astropy import cosmology as cosmo

import fsps 
from datetime import datetime
from powderday.grid_construction import stars_coordinate_boost

from multiprocessing import Pool
from functools import partial
from itertools import repeat
try: #scipy versioning check
    from scipy.integrate import simps
except:
    from scipy.integrate import simpson as simps

from powderday.nebular_emission.cloudy_tools import calc_LogQ, age_dist, cmdf, get_nearest,convert_metals
from powderday.analytics import logu_diagnostic,dump_emlines
from powderday.nebular_emission.cloudy_model import get_nebular
from p_tqdm import p_map


#this is required to keep the reg as a strong reference.  for some
#reason in the star_list.append in star_list_gen, reg otherwise gets
#garbage collected.
import gc
gc.set_threshold(0)

# Lazily initialize FSPS
sp = None

class Stars:
    def __init__(self,mass,metals,positions,age,sed_bin=[-1,-1,-1],lum=-1,fsps_zmet=20,all_metals=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]):
        self.mass = mass
        self.metals = metals
        self.positions = positions
        self.age = age
        self.sed_bin = sed_bin
        self.lum = lum
        self.fsps_zmet = fsps_zmet
        self.all_metals = all_metals

    def info(self):
        return(self.mass,self.metals,self.positions,self.age,self.sed_bin,self.lum,self.fsps_zmet)
        
def star_list_gen(boost,dx,dy,dz,reg,ds,sp,m):
    print ('[SED_gen/star_list_gen]: reading in stars particles for SPS calculation')
    mass = reg["star","masses"].value
    positions = reg["star","coordinates"].value
    age = reg["stellar","ages"].value
    nstars = len(mass)
    el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe' ]

    try:                                                                                                                                                                
        metals = np.zeros((nstars,11))-10.0
        for i in range(11):
            if i == 0:
                el_str = ""
            else:
                el_str = "_"+el[i-1]
            metals[:, i] = reg["star","metals"+el_str].value
    except:
        metals = reg["star","metals"].value
    print ('number of new stars =',nstars)
    
    #calculate the fsps interpolated metallicity

    #if the metallicity has many fields, and not just global
    #metallicity then just extract the global metallicity
    if metals.ndim > 1:
        metals_tot = metals[:,0]
    else:
        metals_tot = metals

    print ('[SED_gen/star_list_gen:] Manually increasing the newstar metallicities by: ',cfg.par.Z_init)
    metals_tot += cfg.par.Z_init
    
    #ADVANCED FEATURE - if force_stellar_metallcities or force_stellar_ages are set, then we set to those values
    if cfg.par.FORCE_STELLAR_AGES:
        print ("[SED_GEN/stars_list_gen:]  FORCE_STELLAR_AGES is set to True: setting all stars to age: %e Gyr"%cfg.par.FORCE_STELLAR_AGES_VALUE)
        age = np.repeat(cfg.par.FORCE_STELLAR_AGES_VALUE,nstars)

    if cfg.par.FORCE_STELLAR_METALLICITIES:
        print ("[SED_GEN/stars_list_gen:]  FORCE_STELLAR_METALLICITIES is set to True: setting all stars to metallicity: %e "%cfg.par.FORCE_STELLAR_METALLICITIES_VALUE)
        metals_tot = np.repeat(cfg.par.FORCE_STELLAR_METALLICITIES_VALUE,nstars)



    zmet = fsps_metallicity_interpolate(metals_tot, sp)
    #mwd(zmet,mass,'zmet_distribution.png')

    #print '[SED_gen/star_list_gen: ] fsps zmet codes:',zmet

    #create the stars_list full of Stars objects
    stars_list = []
    
    if metals.ndim > 1:
        for i in range(nstars):
            stars_list.append(Stars(mass[i],metals_tot[i],positions[i],age[i],fsps_zmet=zmet[i],all_metals = metals[i]))
    else:
        for i in range(nstars):
            stars_list.append(Stars(mass[i],metals_tot[i],positions[i],age[i],fsps_zmet=zmet[i]))
    
    #boost stellar positions to grid center
    print ('boosting new stars to coordinate center')
    stars_list = stars_coordinate_boost(stars_list,boost)
    
    #orig_stars_list_len = len(stars_list)
    
    #ASSIGN DISK AND BULGE STARS - note, if these don't exist, it will
    #just make empty lists

   

    
    bulgestars_list = []
    diskstars_list = []

    
    #in principle, we should just be able to do the following blocks
    #if the particle types exist. the issue is that different groups
    #use PartType2 and 3 as 'filler' particle types, so they may exist
    #even if they don't correspond to disk/bulge stars.


    if ds.cosmological_simulation == False:

        #Disk Stars

        if ("diskstar","coordinates") in ds.derived_field_list:
            
            disk_positions = reg[("diskstar","coordinates")].value
            disk_masses =  reg[("diskstar","masses")].value
            nstars_disk = len(disk_masses)
     
            #create the disk_list full of DiskStars objects
            for i in range(nstars_disk):
                diskstars_list.append(Stars(disk_masses[i],cfg.par.solar,disk_positions[i],cfg.par.disk_stars_age))

            print ('boosting disk stars to coordinate center')    
            diskstars_list = stars_coordinate_boost(diskstars_list,boost)

        #orig_disk_stars_list_len = nstars_disk
            
       
        #Bulge Stars


        if ("bulgestarcoordinates") in ds.derived_field_list:
            bulge_positions = reg[("bulgestar","coordinates")].value
            bulge_masses =  reg[("bulgestar","masses")].value
            nstars_bulge = len(bulge_masses)
            
            #create the bulge_list full of BulgeStars objects
            
            for i in range(nstars_bulge):
                bulgestars_list.append(Stars(bulge_masses[i],cfg.par.solar,bulge_positions[i],cfg.par.bulge_stars_age))
                

            print ('boosting bulge stars to coordinate center')
            bulgestars_list = stars_coordinate_boost(bulgestars_list,boost)

    
    #remove stars that are outside the grid 
    stars_list,bulgestars_list,diskstars_list = remove_stars_outside_grid(stars_list,bulgestars_list,diskstars_list,m)



    #EXPERIMENTAL FEATURES
    if cfg.par.SOURCES_IN_CENTER == True:
        for i in range(nstars):
            stars_list[i].positions[:] =  np.array([0,0,0])
        if ("bulgestar","coordinates") in ds.derived_field_list:
            for i in range(nstars_bulge):
                bulgestars_list[i].positions[:] =  np.array([0,0,0])
            for i in range(nstars_disk):
                diskstars_list[i].positions[:] = np.array([0,0,0])

    if cfg.par.SOURCES_RANDOM_POSITIONS == True:
        print ("================================")
        print ("SETTING SOURCES TO RANDOM POSITIONS")
        print ("================================")
        for i in range(nstars):
            xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
            stars_list[i].positions[:] = np.array([xpos,ypos,zpos])

        if ("bulgestar","coordinates") in ds.derived_field_list:
            for i in range(nstars_bulge):
                xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
                bulgestars_list[i].positions[:] = np.array([xpos,ypos,zpos])
            for i in range(nstars_disk):
                xpos,ypos,zpos = np.random.uniform(-0.9*dx/2.,0.9*dx/2.),np.random.uniform(-0.9*dy/2.,0.9*dy/2.),np.random.uniform(-0.9*dz/2.,0.9*dz/2.)
                diskstars_list[i].positions[:] = np.array([xpos,ypos,zpos])

    return stars_list,diskstars_list,bulgestars_list,reg



def allstars_sed_gen(stars_list,cosmoflag,sp):


    #NOTE this part is just for the gadget simulations - this will
    #eventually become obviated as it gets passed into a function to
    #populate the stars_list with objects as we start to feed in new
    #types of simulation results.

    nstars = len(stars_list)
    
    #get just the wavelength array
    sp.params["tage"] = stars_list[0].age
    sp.params["imf_type"] = cfg.par.imf_type
    sp.params["pagb"] = cfg.par.pagb
    sp.params["sfh"] = 0
    sp.params["zmet"] = stars_list[0].fsps_zmet
    sp.params["add_neb_emission"] = cfg.par.add_neb_emission
    sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
    sp.params['gas_logu'] = cfg.par.gas_logu

    if cfg.par.FORCE_gas_logz == False:
        sp.params['gas_logz'] = np.log10(stars_list[0].metals/cfg.par.solar)
    else:
        sp.params['gas_logz'] = cfg.par.gas_logz

        '''
    sp = fsps.StellarPopulation(tage=stars_list[0].age,imf_type=cfg.par.imf_type,pagb = cfg.par.pagb,sfh=0,zmet=stars_list[0].fsps_zmet,
                                add_neb_emission = cfg.par.add_neb_emission, add_agb_dust_model=cfg.par.add_agb_dust_model)
                                '''
    spec = sp.get_spectrum(tage=stars_list[0].age,zmet=stars_list[0].fsps_zmet)
    nu = 1.e8*constants.c.cgs.value/spec[0]
    nlam = len(nu)

    nprocesses = np.min([cfg.par.n_processes,len(stars_list)]) #the pool.map will barf if there are less star bins than process threads

    t1=datetime.now()
    stars_sed_gen = p_map(newstars_gen, stars_list, num_cpus=nprocesses)
    t2=datetime.now()

    print ('Execution time for SED generation in Pool.map multiprocessing = '+str(t2-t1))
    
    cloudy_nlam = len(np.genfromtxt(cfg.par.pd_source_dir + "/powderday/nebular_emission/data/refLines.dat", delimiter=','))
    stellar_fnu = np.zeros([nstars,nlam])
    mfrac = np.zeros(nstars)
    line_em = np.zeros([nstars,cloudy_nlam])

    # stars_sed_gen returns three things for each star/star bin: the spectrum,the associated surviving stellar mass fraction for that SSP
    # and the line luminosities from CLOUDY
    for i in range(nstars):
        stellar_fnu[i,:] = stars_sed_gen[i][0]
        mfrac[i] = stars_sed_gen[i][1]
        line_em[i, :] = stars_sed_gen[i][2]

    stellar_nu = nu


    if cosmoflag == False:

        #calculate the SED for disk stars; note, this gets calculated
        #whether or not disk stars actually exist.  if they don't exist,
        #bogus values for the disk age and metallicity are assigned based
        #on whatever par.disk_stars_age and metallicity are.  it's no big
        #deal since these SEDs don't end up getting added to the model in
        #source_creation.  
        
        #note, even if there are no disk/bulge stars, these are still
        #created since they're completely based on input parameters in
        #parameters_master.  they just won't get used at a later point
        #as there will be no disk/bulge star positions to add them to.

        #dust_tesc is an absolute value (not relative to min star age) as the ages of these stars are input by the user

        # Load in the metallicity legend
        fsps_metals = np.array(sp.zlegend)

        sp.params["tage"] = cfg.par.disk_stars_age
        sp.params["imf_type"] = cfg.par.imf_type
        sp.params["imf1"] = cfg.par.imf1
        sp.params["imf2"] = cfg.par.imf2
        sp.params["imf3"] = cfg.par.imf3
        sp.params["pagb"] = cfg.par.pagb
        sp.params["sfh"] = 0
        sp.params["zmet"] = cfg.par.disk_stars_metals
        sp.params["add_neb_emission"] = cfg.par.add_neb_emission
        sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
        sp.params['gas_logu'] = cfg.par.gas_logu
        if cfg.par.FORCE_gas_logz == False:
            sp.params['gas_logz'] = np.log10(fsps_metals[cfg.par.disk_stars_metals]/cfg.par.solar)
        else:
            sp.params['gas_logz'] = cfg.par.gas_logz

        spec = sp.get_spectrum(tage=cfg.par.disk_stars_age,zmet=cfg.par.disk_stars_metals)
        disk_fnu = spec[1]
        
        #calculate the SED for bulge stars
        sp.params["tage"] = cfg.par.bulge_stars_age
        sp.params["imf_type"] = cfg.par.imf_type
        sp.params["imf1"] = cfg.par.imf1
        sp.params["imf2"] = cfg.par.imf2
        sp.params["imf3"] = cfg.par.imf3
        sp.params["pagb"] = cfg.par.pagb
        sp.params["sfh"] = 0
        sp.params["zmet"] = cfg.par.bulge_stars_metals
        sp.params["add_neb_emission"] = cfg.par.add_neb_emission
        sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
        sp.params['gas_logu'] = cfg.par.gas_logu
        if cfg.par.FORCE_gas_logz == False:
            sp.params['gas_logz'] = np.log10(fsps_metals[cfg.par.bulge_stars_metals]/cfg.par.solar)
        else:
            sp.params['gas_logz'] = cfg.par.gas_logz


        spec = sp.get_spectrum(tage=cfg.par.bulge_stars_age,zmet=cfg.par.bulge_stars_metals)
        bulge_fnu = spec[1]
    

    else: #we have a cosmological simulation

        disk_fnu = []
        bulge_fnu = []
        

    total_lum_in_sed_gen = 0.
    for i in range(stellar_fnu.shape[0]):
        total_lum_in_sed_gen += np.absolute(np.trapz(stellar_fnu[i,:],x=nu))

    print ('[SED_gen: ] total_lum_in_sed_gen = ',total_lum_in_sed_gen)

    #return positions,disk_positions,bulge_positions,mass,stellar_nu,stellar_fnu,disk_masses,disk_fnu,bulge_masses,bulge_fnu
    return stellar_nu,stellar_fnu,disk_fnu,bulge_fnu, mfrac, line_em


def newstars_gen(star_object):
    global sp
    if sp is None:
        sp = fsps.StellarPopulation()
    
    #the newstars (particle type 4; so, for cosmological runs, this is all
    #stars) are calculated in a separate function with just one argument so that it is can be fed 
    #into pool.map for multithreading.
    #sp = fsps.StellarPopulation()
    sp.params["tage"] = star_object.age
    sp.params["imf_type"] = cfg.par.imf_type
    sp.params["pagb"] = cfg.par.pagb
    sp.params["sfh"] = 0
    sp.params["zmet"] = star_object.fsps_zmet
    sp.params["add_neb_emission"] = False
    sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model

    #first figure out how many wavelengths there are
    
    spec = sp.get_spectrum(tage=star_object.age,zmet=star_object.fsps_zmet)
    nu = 1.e8*constants.c.cgs.value/spec[0]

    nlam = len(nu)
  
    minage = star_object.age

    tesc_age = np.log10((minage+cfg.par.birth_cloud_clearing_age)*1.e9)
    # Get the number of ionizing photons from SED

    #calculate the SEDs for new stars
    sp.params["tage"] = star_object.age
    sp.params["imf_type"] = cfg.par.imf_type
    sp.params["imf1"] = cfg.par.imf1
    sp.params["imf2"] = cfg.par.imf2
    sp.params["imf3"] = cfg.par.imf3
    sp.params["pagb"] = cfg.par.pagb
    sp.params["sfh"] = 0
    sp.params["zmet"] = star_object.fsps_zmet
    sp.params["add_neb_emission"] = False
    sp.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model


    if cfg.par.CF_on == True:
        sp.params["dust_type"] = 0
        sp.params["dust1"] = 1
        sp.params["dust2"] = 0
        sp.params["dust_tesc"] = tesc_age
        
    if cfg.par.dust_screen == True:
            sp.params["dust_type"]=0
            sp.params["dust1"]=cfg.par.dust1
            sp.params["dust2"]=cfg.par.dust2

    if cfg.par.alpha_enhance: #Setting Zstar based on Fe/H
        spec_noneb, mfrac_neb = alpha_enhance(star_object.all_metals[-1], star_object.fsps_zmet, star_object.age, tesc_age)
    else:
        spec_noneb = sp.get_spectrum(tage=star_object.age, zmet=star_object.fsps_zmet)
    
    f = spec_noneb[1]

    #NOTE: FSPS SSP/CSP spectra are scaled by *formed* mass, not current mass. i.e., the SFHs of the SSP/CSP are normalized such that 1 solar mass
    #is formed over the history. This means that the stellar spectra are normalized by the integral of the SFH =/= current
    #(surviving, observed, etc.) stellar mass. In simulations, we only know the current star particle mass. To get the formed mass for an SSP,
    #we generate the surviving mass fraction (sp.stellar_mass) to extrapolate the initial mass from the current mass, metallicity, and age
    #this 'mfrac' is used to scale the FSPS SSP luminosities in source_creation
    mfrac = sp.stellar_mass
    
    cloudy_nlam = len(np.genfromtxt(cfg.par.pd_source_dir + "/powderday/nebular_emission/data/refLines.dat", delimiter=','))
    line_em = np.zeros(cloudy_nlam)

    pagb = cfg.par.add_pagb_stars and cfg.par.PAGB_min_age <= star_object.age <= cfg.par.PAGB_max_age
    young_star = cfg.par.add_young_stars and cfg.par.HII_min_age <= star_object.age <= cfg.par.HII_max_age

    if cfg.par.add_neb_emission and cfg.par.use_cloudy_tables and star_object.age <= 1.e-2:
        LogZ = np.log10(star_object.metals/cfg.par.solar)
        Rin = 1.e19   # Rinner is fixed at 1.e19 cm for lookup tables
        nh = 1e3 # Hydrogen Density is fixed at 1.e3 cm for lookup tables
        LogQ = calc_LogQ(1.e8*constants.c.cgs.value/spec_noneb[0], spec_noneb[1]*constants.L_sun.cgs.value, 
                            mstar=star_object.mass/constants.M_sun.cgs.value, mfrac=mfrac)
        LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*nh*constants.c.cgs.value))

        sp.params['gas_logu'] = LogU
        sp.params['gas_logz'] = LogZ
        sp.params["add_neb_emission"] = True
        lam_neb, spec_neb = sp.get_spectrum(tage=star_object.age, zmet=star_object.fsps_zmet)
        line_em = sp.emline_luminosity
        f = spec_neb
        
    elif (cfg.par.add_neb_emission or cfg.par.use_cmdf) and (young_star or pagb) and not cfg.par.use_cloudy_tables:
        # For each star particle we break it into a collection or cluster of star particles which have the same property as the parent star particle
        # but their masses and ages follow a power-law distribution if use_cmdf and use_age_distribution are set to True respectively and it falls
        # under the constraints set in parameters_master. To do so we create 1D arrays that store the masses of each cluster (cluster_mass), num of particles
        # in that cluster (num_cluster), and age of that cluster (age_cluster). So for example (assuming default values), a 1e6 solar mass particle
        # will first be broken down into 6 particles with different masses ranging from 10^3.5 Msun to 10^5 Msun. Each of these particles will be further
        # broken down into 5 particles with their ages are distributed as per the age distribution. Thus in total, this one particle will be broken
        # down into 30 particles and these arrays will store the properties of all the 30 particles. This allows us to consider these as 30 individual
        # particles rest of the calculation and their fluxes are combined in end to get the final result for this one particle.
            
        cluster_mass = [np.log10(star_object.mass / constants.M_sun.cgs.value)]
        num_clusters = [1]
        age_clusters = [star_object.age]

        if star_object.mass / constants.M_sun.cgs.value > 10 ** cfg.par.cmdf_max_mass and cfg.par.use_cmdf:
            cluster_mass, num_clusters = cmdf(star_object.mass / constants.M_sun.cgs.value,
                                                int(cfg.par.cmdf_bins), cfg.par.cmdf_min_mass,
                                                cfg.par.cmdf_max_mass, cfg.par.cmdf_beta)
            age_clusters = []
            for k in range(len(cluster_mass)):
                age_clusters.append(star_object.age)

            if cfg.par.use_age_distribution:
                num_clusters_cmdf = num_clusters
                cluster_mass_cmdf = cluster_mass
                num_clusters = []
                cluster_mass = []
                age_clusters = []
                for k in range(len(cluster_mass_cmdf)):
                    num, t = age_dist(num_clusters_cmdf[k], star_object.age)
                    rescale = np.sum(num_clusters_cmdf[k])/np.sum(num)
                        
                    for l in range(len(num)):
                        if num[l] == 0:
                            continue
                        num_clusters.append(num[l])
                        cluster_mass.append(np.log10((10**cluster_mass_cmdf[k])*rescale))
                        age_clusters.append(t[l])

            cluster_mass = np.array(cluster_mass)
            num_clusters = np.array(num_clusters)
            age_clusters = np.array(age_clusters)
            
        f = np.zeros(nlam)

        if cfg.par.FORCE_gas_logz == False:
            LogZ = np.log10(stars_list[i].metals/cfg.par.solar)
        else:
            LogZ = cfg.par.gas_logz
        
        if cfg.par.dust_screen == True:
            sp.params["dust_type"]=0
            sp.params["dust1"]=cfg.par.dust1
            sp.params["dust2"]=cfg.par.dust2


        if cfg.par.CF_on == True:
            sp.params["dust_type"] = 0
            sp.params["dust1"] = 1
            sp.params["dust2"] = 0
            sp.params["dust_tesc"] = tesc_age


        for j in range(len(cluster_mass)):
            num_HII_clusters = num_clusters[j]
            age = age_clusters[j]
            neb_file_output = cfg.par.NEB_DEBUG

            sp.params["add_neb_emission"] = False

            if cfg.par.alpha_enhance: #Setting Zstar based on Fe/H
                spec, mfrac_neb = alpha_enhance(star_object.all_metals[-1], star_object.fsps_zmet, age, tesc_age)
                    
            else:
                spec = sp.get_spectrum(tage=age,zmet=star_object.fsps_zmet)
                mfrac_neb = sp.stellar_mass

            if cfg.par.add_neb_emission:
                # id_val = 0, 1, 2 for young stars, Post-AGB star and AGNs respectively.
                if young_star:
                    id_val = 0
                    Rinner_per_Rs = cfg.par.HII_Rinner_per_Rs
                    nh = cfg.par.HII_nh
                    escape_fraction  = cfg.par.HII_escape_fraction
                    
                elif pagb:
                    id_val = 1
                    Rinner_per_Rs = cfg.par.PAGB_Rinner_per_Rs
                    nh = cfg.par.PAGB_nh
                    escape_fraction  = cfg.par.PAGB_escape_fraction

                alpha = 2.5e-13 # Recombination Rate (assuming T = 10^4 K)

                if cfg.par.FORCE_gas_logu[id_val]:
                    LogU = cfg.par.gas_logu[id_val]
                    LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*nh))
                    Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                    
                elif cfg.par.FORCE_logq[id_val]:
                    LogQ = cfg.par.source_logq[id_val]
                    Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                    LogU = np.log10((10**LogQ)/(4*np.pi*Rs*Rs*nh*constants.c.cgs.value))

                else:
                    LogQ = calc_LogQ(1.e8*constants.c.cgs.value/spec[0], spec[1]*constants.L_sun.cgs.value, mstar=10**cluster_mass[j], mfrac=mfrac_neb)
                    Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)
                    LogU = np.log10((10**LogQ)/(4*np.pi*Rs*Rs*nh*constants.c.cgs.value))+cfg.par.gas_logu_init[id_val]
                    LogQ = np.log10((10 ** (3*LogU))*(36*np.pi*(constants.c.cgs.value**3))/((alpha**2)*nh))
                    Rs = ((3*(10 ** LogQ))/(4*np.pi*(nh**2)*alpha))**(1./3.)

                if cfg.par.FORCE_inner_radius[id_val]:
                     Rin = cfg.par.inner_radius[id_val]

                else:
                    Rin = Rinner_per_Rs*Rs

                if cfg.par.FORCE_gas_logz[id_val]:
                    LogZ = cfg.par.gas_logz[id_val]
                    
                else:
                    LogZ = np.log10(star_object.metals/cfg.par.solar)

                if cfg.par.FORCE_inner_radius[id_val]:
                    Rin = cfg.par.inner_radius[id_val]

                if neb_file_output:        
                    logu_diagnostic(LogQ, LogU, LogZ, Rs, 10**cluster_mass[j], num_HII_clusters, age, append=True)
                    neb_file_output = False

                try:
                    # Calculating ionizing photons again but for 1 Msun in order to scale the output for FSPS
                    LogQ_1 = calc_LogQ(1.e8 * constants.c.cgs.value / spec[0], spec[1] * constants.L_sun.cgs.value)

                    spec_neb, wave_line, line_lum = get_nebular(spec[0], spec[1], nh, star_object.all_metals, logq = LogQ, radius = Rin,
                                                logu = LogU, logz = LogZ, logq_1 = LogQ_1, Dust=cfg.par.HII_dust, abund=cfg.par.neb_abund[id_val],
                                            clean_up = cfg.par.cloudy_cleanup, index=id_val, efrac=escape_fraction)
                except ValueError as err:
                    # If the CLOUDY run crashes we switch to using FSPS lookup tables.
                    print ("WARNING: Switching to using lookup tables pre-packed with FSPS to calculate nebular emission for this particle.")
                    print ("WARNING: The emission line fluxes repoted may not be accurate if the particle lies outside the range of the lookup table paramters.")
                    sp.params['gas_logu'] = LogU
                    sp.params['gas_logz'] = LogZ
                    sp.params["add_neb_emission"] = True
                    lam_neb, spec_neb = sp.get_spectrum(tage=age, zmet=star_object.fsps_zmet)
                    line_lum = sp.emline_luminosity

            else:
                spec_neb = spec[1]

            weight = num_HII_clusters*(10**cluster_mass[j])/(star_object.mass/constants.M_sun.cgs.value)
            f = f + spec_neb*weight
            if cfg.par.add_neb_emission and cfg.par.dump_emlines:
                line_em = line_em + line_lum*weight

    stellar_fnu = f

    return stellar_fnu, mfrac, line_em


def get_gas_metals(ngas, reg):
    # This function outputs the metallicity (total as well as all the 10 elements tracked by the simulation)
    # for all the gas particles
    el = ['He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
    metals = np.zeros((ngas,11))-10.0
    try:
        for i in range(11):
            if i == 0:
                el_str = ""
            else:
                el_str = "_"+el[i-1]
            metals[:, i] = reg["gas","metals"+el_str].value
    except:
        metals[:,0] = reg["gas","metals"].value
    
    return metals


def get_nearest_gas_metals(all_gas_coordinates, particle_coordinates, all_gas_metals):
    # Getting N nearest gas particles to the AGN where N is defined by cfg.par.AGN_num_gas
    nearest_gas_dist, nearest_gas_id = get_nearest(all_gas_coordinates, particle_coordinates, num=cfg.par.AGN_num_gas)
    metals_avg = []
    
    # We take the distance weighted avearge of the metallicity of nearest N gas particles.
    # This is used as input to the CLOUDY model.
    for q in range(11):
        nearest_gas_metals = np.array(all_gas_metals[:,q][nearest_gas_id])
        metals_avg.append(np.sum(nearest_gas_metals*(1/nearest_gas_dist))/np.sum(1/nearest_gas_dist))

    return metals_avg
                                                                                                                

def get_agn_seds(agn_ids, reg):
    print ('Starting AGN SED generation')
    
    if cfg.par.dump_emlines:
        dump_emlines(None, None)

    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,len(agn_ids)])

    # Pre-calculating average metallicity and fluxes for all the BH particles
    nu = []
    fnu_in = []
    metals_avg = []
    for agn_id in agn_ids:
        fnu_in.append(reg["bh","sed"][agn_id,:].in_units("Lsun").value/reg["bh","nu"].in_units("Hz").value)
        nu.append(reg["bh","nu"].in_units("Hz").value)
        all_gas_coordinates = reg["gas","coordinates"].in_units('kpc').value
        agn_coordinates = reg["bh","coordinates"][agn_id].in_units('kpc').value
        all_gas_metals = get_gas_metals(len(all_gas_coordinates), reg)
        metals_avg.append(get_nearest_gas_metals(all_gas_coordinates, agn_coordinates, all_gas_metals))

    fnu_out = p_map(agn_sed, agn_ids, nu, fnu_in, metals_avg, num_cpus=nprocesses)
    fnu_out = np.atleast_2d(fnu_out)
    t2 = datetime.now()

    print ('Execution time for AGN SED generation = '+str(t2-t1))
    
    return fnu_out


def agn_sed(agn_id, nu, fnu, metals_avg):
    agn_id = int(agn_id)

    if cfg.par.add_neb_emission and cfg.par.add_AGN_neb:
        
        id_val = 2
            
        tot_metals = metals_avg[0]
        metals = metals_avg
        
        if cfg.par.FORCE_gas_logz[id_val]:
            LogZ = cfg.par.gas_logz[id_val]
        
        else:
            LogZ = tot_metals/cfg.par.solar

        Rin = cfg.par.inner_radius[id_val]

        if cfg.par.FORCE_gas_logu[id_val]:
            LogU = cfg.par.gas_logu[id_val]
            LogQ = np.log10((10**LogU)*(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))
        
        elif cfg.par.FORCE_logq[id_val]:
            LogQ = cfg.par.source_logq[id_val]
            LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))

        else:
            LogQ = calc_LogQ(nu, fnu*constants.L_sun.cgs.value)
            LogU = np.log10((10**LogQ)/(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))+cfg.par.gas_logu_init[id_val]
            LogQ = np.log10((10**LogU)*(4*np.pi*Rin*Rin*cfg.par.AGN_nh*constants.c.cgs.value))

        spec, wave_line, line_lum = get_nebular(1.e8 * constants.c.cgs.value / nu, fnu, cfg.par.AGN_nh, metals, logq = LogQ, radius = Rin, logu = LogU, 
                                        logz = LogZ, logq_1 = LogQ, Dust=False, abund=cfg.par.neb_abund[id_val],clean_up = cfg.par.cloudy_cleanup, index=id_val)

        if cfg.par.dump_emlines:
        # The stellar population returns the calculation in units of Lsun
            line_em = line_lum * 3.839e33  # Units: ergs/s
            # The last column in dump_emlines is reserved for age of the star particle. For AGN we just set it to -1 as a place holder.
            OH = metals[4]
            line_em = np.append(line_em, OH)
            line_em = np.append(line_em, 4)
            dump_emlines(line_em)

    else:
        spec = fnu

    return spec


def get_dig_seds(lam ,sed, logU, cell_widths, metals):

    print ('Starting DIG SED generation')
    
    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,len(cell_widths)])
    
    fnu_out = p_map(partial(dig_sed, spec_lam=lam), sed, logU, cell_widths, metals, num_cpus=nprocesses)
    fnu_out = np.atleast_2d(fnu_out)
    
    t2 = datetime.now()

    print ('Execution time for DIG SED generation = '+str(t2-t1))

    return fnu_out


def dig_sed(sspi, logU, cell_width, metal, spec_lam):    
    id_val = 3   

    spec, wave_line, line_lum = get_nebular(spec_lam, sspi, cfg.par.DIG_nh, metal, logu=logU, Cell_width=cell_width, Dust=False,
                                            abund=cfg.par.neb_abund[id_val], clean_up = cfg.par.cloudy_cleanup, index=id_val)

    if cfg.par.dump_emlines:
        # The stellar population returns the calculation in units of Lsun
        line_em = line_lum * 3.839e33  # Units: ergs/s
        OH = metal[4]
        line_em = np.append(line_em, OH)
        line_em = np.append(line_em, id_val)
        dump_emlines(line_em)
        
    return spec


def fsps_metallicity_interpolate(metals, sp):

    # takes a list of metallicities for star particles, and returns a
    # list of interpolated metallicities

    fsps_metals = np.array(sp.zlegend)

    nstars = len(metals)
    
    zmet = []

    for i in range(nstars):
        zmet.append(find_nearest_zmet(fsps_metals,metals[i]))
    
    return zmet
    
def find_nearest_zmet(array,value):
    # this is modified from the normal find_nearest in that it forces
    # the output to be 1 index higher than the true value since the
    # minimum zmet value fsps will take is 1 (not 0)

    idx = (np.abs(array-value)).argmin()
    
    return idx+1      



def remove_stars_outside_grid(stars_list,bulgestars_list,diskstars_list,m):
    
    #first get grid limits.  how we do this will depend on the type of
    #grid it is because hyperion stores this information differently
    #in the model object depending on the grid type.

    print("[SED_gen/remove_stars_outside_grid]: Searching to see if any stars are outside the dust grid, which will cause a crash in the radiative transfer modules")

    if m.__dict__['grid_type'] == 'vor':
        xmax = m.grid.xmax
        xmin = m.grid.xmin
        ymax = m.grid.ymax
        ymin = m.grid.ymin
        zmax = m.grid.zmax
        zmin = m.grid.zmin

    elif m.__dict__['grid_type'] == 'oct':
        xmax = m.grid.dx
        xmin = xmax*-1
        ymax = m.grid.dy
        ymin = ymax*-1
        zmax = m.grid.dz
        zmin = zmax*-1

    elif m.__dict__['grid_type'] == 'amr':
        xmax = m.grid.levels[0].grids[0].xmax
        xmin = m.grid.levels[0].grids[0].xmin
        ymax = m.grid.levels[0].grids[0].ymax
        ymin = m.grid.levels[0].grids[0].ymin
        zmax = m.grid.levels[0].grids[0].zmax
        zmin = m.grid.levels[0].grids[0].xmin
        for ilevel, level_ref in enumerate(m.grid.levels):
            for igrid, grid_ref in enumerate(level_ref.grids):
                if grid_ref.xmin<xmin:
                    xmin = grid_ref.xmin
                if grid_ref.xmax>xmax:
                    xmax = grid_ref.xmax
                if grid_ref.ymin<ymin:
                    ymin = grid_ref.ymin
                if grid_ref.ymax>ymax:
                    ymax = grid_ref.ymax
                if grid_ref.zmin<zmin:
                    zmin = grid_ref.zmin
                if grid_ref.zmax>zmax:
                    zmax = grid_ref.zmax

    star_idx_to_remove = []
    bulge_idx_to_remove = []
    disk_idx_to_remove = []

    total_mass = 0
    mass_removed = 0


    for i in range(len(stars_list)):
        if (stars_list[i].positions[0] > xmax) or \
           (stars_list[i].positions[0] < xmin) or \
           (stars_list[i].positions[1] > ymax) or \
           (stars_list[i].positions[1] < ymin) or \
           (stars_list[i].positions[2] > zmax) or \
           (stars_list[i].positions[2] < zmin):
           
            star_idx_to_remove.append(i)
            mass_removed += stars_list[i].mass


#        if cfg.par.dust_screen == True:
#            sp.params["dust_type"]=0
#            sp.params["dust1"]=cfg.par.dust1
#            sp.params["dust2"]=cfg.par.dust2


#        if cfg.par.CF_on == True:
#            sp.params["dust_type"] = 0
#            sp.params["dust1"] = 1
#            sp.params["dust2"] = 0
#            sp.params["dust_tesc"] = tesc_age


        if (len(bulgestars_list) > 0):
            if (bulgestars_list[i].positions[0] > xmax) or \
               (bulgestars_list[i].positions[0] < xmin) or \
               (bulgestars_list[i].positions[1] > ymax) or \
               (bulgestars_list[i].positions[1] < ymin) or \
               (bulgestars_list[i].positions[2] > zmax) or \
               (bulgestars_list[i].positions[2] < zmin):

                bulge_idx_to_remove.append(i)
                mass_removed += bulgestars_list[i].mass

        
        if (len(diskstars_list) > 0):
            if (diskstars_list[i].positions[0] > xmax) or \
               (diskstars_list[i].positions[0] < xmin) or \
               (diskstars_list[i].positions[1] > ymax) or \
               (diskstars_list[i].positions[1] < ymin) or \
               (diskstars_list[i].positions[2] > zmax) or \
               (diskstars_list[i].positions[2] < zmin):
            
                disk_idx_to_remove.append(i)
                mass_removed += diskstars_list[i].mass



        total_mass += stars_list[i].mass
    

    #now that we've figured out which stars to remove, actually remove
    #them from the lists.  remove them in reverse order so that we
    #don't destroy the indexing as soon as we remove a single star
    for idx in star_idx_to_remove[::-1]:
        stars_list.pop(idx)

    for idx in bulge_idx_to_remove[::-1]:
        bulgestars_list.pop(idx)

    for idx in disk_idx_to_remove[::-1]:
        diskstars_list.pop(idx)
    

    number_of_removed_stars = len(star_idx_to_remove) + len(bulge_idx_to_remove) + len(disk_idx_to_remove)
    mass_fraction_removed = mass_removed/total_mass
    print("[SED_gen/remove_stars_outside_grid:] removing %f stars because they are outside the dust grid" % number_of_removed_stars)
    print("[SED_gen/remove_stars_outside_grid:] this corresponds to %f of the total stellar mass in the volume " % mass_fraction_removed)

    return stars_list,bulgestars_list,diskstars_list
        
    
def alpha_enhance(star_Fe, star_zmet, age, tesc_age):
    Fe = star_Fe
    
    # Gizmo metallicity structure, photospheric abundances from Asplund et al. 2009:
    # Photospheric mass fraction of H = 0.7381
    # Photospheric mass fraction of Fe = 1.31e-3

    # Converting from mass fraction to atomic fraction of Fe
    # Taking atmoic mass of H = 1.008u
    # Taking atmoic mass of Fe = 55.845u
    FeH = (Fe/0.7381)*(1.008/55.845)

    # Solar atomic fraction of Fe. Calculated by substituting Fe = 1.31e-3 in the previous equation
    FeH_sol = 3.22580645e-5

    Logzsol = np.log10(FeH/FeH_sol)
                        
    sp1 = fsps.StellarPopulation(zcontinuous=1)
    sp1.params["tage"] = age
    sp1.params["imf_type"] = cfg.par.imf_type
    sp1.params["imf1"] = cfg.par.imf1
    sp1.params["imf2"] = cfg.par.imf2
    sp1.params["imf3"] = cfg.par.imf3
    sp1.params["pagb"] = cfg.par.pagb
    sp1.params["sfh"] = 0
    sp1.params["zmet"] = star_zmet
    sp1.params["add_neb_emission"] = False
    sp1.params["add_agb_dust_model"] = cfg.par.add_agb_dust_model
    sp1.params["logzsol"] = Logzsol

    if cfg.par.CF_on == True:
        sp1.params["dust_type"] = 0
        sp1.params["dust1"] = 1
        sp1.params["dust2"] = 0
        sp1.params["dust_tesc"] = tesc_age

        spec = sp1.get_spectrum(tage=age)
        mfrac_neb = sp1.stellar_mass

    return spec, mfrac_neb
