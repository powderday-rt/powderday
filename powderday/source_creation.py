from __future__ import print_function
#source_creation.py
#purpose: to add new star sources, disks and bulge stars
import powderday.config as cfg
import numpy as np
import powderday.SED_gen as sg
from datetime import datetime
import astropy.units as units
import astropy.constants as constants
from powderday.helpers import find_nearest
from powderday.analytics import dump_AGN_SEDs,dump_NEB_SEDs,dump_emlines
from hyperion.model import ModelOutput
from hyperion.grid.yt3_wrappers import find_order
from powderday.nebular_emission.cloudy_tools import get_DIG_sed_shape, get_DIG_logU
from tqdm import tqdm
from scipy import spatial

class Sed_Bins:
    def __init__(self,mass,metals,age,fsps_zmet,all_metals=[-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1]):
        self.mass = mass
        self.metals = metals
        self.age = age
        self.fsps_zmet = fsps_zmet
        self.all_metals = all_metals

def direct_add_stars(df_nu, stars_list, diskstars_list, bulgestars_list, cosmoflag, m, sp):

    print("--------------------------------\n")
    print("Adding unbinned stars to the grid\n")
    print("--------------------------------\n")

    totallum_newstars = 0.
    unbinned_stars_list = []

    for star in stars_list:
        if star.age <= cfg.par.max_age_direct:
            unbinned_stars_list.append(star)


    nstars = len(unbinned_stars_list)
    if (nstars == 0):
        print ("No unbinned stars to add")
        return m

    else:
        print ("Number of unbinned stars to added: ", nstars)
    
    stellar_nu, stellar_fnu, disk_fnu, bulge_fnu, mfrac, unbinned_line_em = sg.allstars_sed_gen(unbinned_stars_list, cosmoflag, sp)
    
    #SED_gen now returns an additional parameter, mfrac, to properly scale the FSPS SSP spectra, which are in units of formed mass, not current stellar mass
    pos_arr = []
    fnu_arr = []
    for i in range(nstars):
        nu = stellar_nu[:]
        fnu = stellar_fnu[i, :]
        pos = unbinned_stars_list[i].positions

        nu, fnu = wavelength_compress(nu, fnu, df_nu)
        # reverse the arrays for hyperion
        nu = nu[::-1]
        fnu = fnu[::-1]

        lum = np.absolute(np.trapz(fnu, x=nu))*unbinned_stars_list[i].mass/constants.M_sun.cgs.value/mfrac[i]
        lum *= constants.L_sun.cgs.value
        
        young_star = cfg.par.add_young_stars and cfg.par.HII_min_age <= unbinned_stars_list[i].age <= cfg.par.HII_max_age
        pagb = cfg.par.add_pagb_stars and cfg.par.PAGB_min_age <= unbinned_stars_list[i].age <= cfg.par.PAGB_max_age

        if cfg.par.add_neb_emission and (star or pagb):
            pos_arr.append(pos)
            fnu_arr.append(stellar_fnu[i, :])

            line_em = unbinned_line_em[i,:]

            # the stellar population returns the calculation in units of Lsun/1 Msun: 
            # https://github.com/dfm/python-fsps/issues/117#issuecomment-546513619
            line_em = line_em * (unbinned_stars_list[i].mass * units.g).to(units.Msun).value * 3.839e33 # Units: ergs/s
            OH = unbinned_stars_list[i].all_metals[4]

            line_em = np.append(line_em, OH)
            if young_star:
                line_em = np.append(line_em, 1)
            elif pagb:
                line_em = np.append(line_em, 2)
                
            dump_emlines(line_em)


        # add new stars
        totallum_newstars += lum
        #m.add_spherical_source(luminosity = lum,radius = 10.*const.rsun,spectrum = (nu,fnu),

        m.add_point_source(luminosity = lum,spectrum=(nu,fnu), position = pos)

    print('[source_creation/add_unbinned_newstars:] totallum_newstars = ', totallum_newstars)
    
    if cfg.par.add_neb_emission and (cfg.par.SAVE_NEB_SEDS or cfg.par.add_DIG_neb) and (len(pos_arr) != 0):
        dump_NEB_SEDs(stellar_nu, fnu_arr, pos_arr)

    if cosmoflag == False: add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,unbinned_stars_list,diskstars_list,bulgestars_list,m)
   
    m.set_sample_sources_evenly(True)
    
    return m


def add_bulge_disk_stars(df_nu,stellar_nu,stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m):
    print ('Non-Cosmological Simulation: Adding Disk and Bulge Stars:')
    

    nu = stellar_nu[:]
    nu,bulge_fnu = wavelength_compress(nu,bulge_fnu,df_nu)
    nu,disk_fnu = wavelength_compress(nu,disk_fnu,df_nu)


    #reverse the arrays for hyperion
    nu = nu[::-1]
    


    nstars_disk = len(diskstars_list)
    nstars_bulge = len(bulgestars_list)
    

    if nstars_disk >0: 

        fnu = disk_fnu[:]

        #reverse the arrays for hyperion
        fnu = fnu[::-1]


        print ('adding disk stars to the grid: adding as a point source collection')
        disksource = m.add_point_source_collection()
            
        
        disk_lum = np.absolute(np.trapz(fnu,x=nu))*diskstars_list[0].mass/constants.M_sun.cgs.value
        #since stellar masses are in cgs, and we need them to be in msun - we
        #multiply by mass to get the *total* luminosity of the stellar
        #cluster since int(nu,fnu) is just the luminosity of a 1 Msun single star
        disk_lum *= constants.L_sun.cgs.value
        disksource.luminosity = np.repeat(disk_lum,nstars_disk)
        
        disk_pos = np.zeros([len(diskstars_list),3])
 
        for i in range(len(diskstars_list)): disk_pos[i,:] = diskstars_list[i].positions
        disksource.position=disk_pos
        
        disksource.spectrum = (nu,fnu)
    
        print ('[source_creation/add_bulge_disk_stars:] totallum_disktars = ',disksource.luminosity[0])
    


    if nstars_bulge > 0:

        fnu = bulge_fnu[:]
        #reverse the arrays for hyperion
        fnu = fnu[::-1]



        print ('adding bulge stars to the grid: adding as a point source collection')
        bulgesource = m.add_point_source_collection()
        bulge_lum = np.absolute(np.trapz(fnu,x=nu))*bulgestars_list[0].mass/constants.M_sun.cgs.value
        bulge_lum *= constants.L_sun.cgs.value
        bulgesource.luminosity = np.repeat(bulge_lum,nstars_bulge)
        
        bulge_pos = np.zeros([len(bulgestars_list),3])
        for i in range(len(bulgestars_list)): bulge_pos[i,:] = bulgestars_list[i].positions
        bulgesource.position=bulge_pos
        bulgesource.spectrum = (nu,fnu)
        
        print ('[source_creation/add_bulge_disk_stars:] totallum_bulgetars = ',bulgesource.luminosity[0])

    

def add_binned_seds(df_nu,stars_list,diskstars_list,bulgestars_list,cosmoflag,m,sp):
    
    # calculate max and min ages
    minimum_age = 15 #Gyr - obviously too high of a number
    maximum_age = 0 #Gyr

    # calculate the minimum and maximum luminosity
    minimum_mass = 1e15*constants.M_sun.cgs.value #msun - some absurdly large value for a single stellar cluster
    maximum_mass = 0 #msun

    # calculate the minimum and maximum stellar metallicity
    minimum_metallicity = 1.e5 #some absurdly large metallicity
    maximum_metallicity = 0

    nstars = len(stars_list)
    for i in range(nstars):
        #if stars_list[i].metals[0] < minimum_metallicity: minimum_metallicity = stars_list[i].metals[0]
        #if stars_list[i].metals[0] > maximum_metallicity: maximum_metallicity = stars_list[i].metals[0]
        
        if stars_list[i].mass < minimum_mass: minimum_mass = stars_list[i].mass
        if stars_list[i].mass > maximum_mass: maximum_mass = stars_list[i].mass

        if stars_list[i].age < minimum_age: minimum_age = stars_list[i].age
        if stars_list[i].age > maximum_age: maximum_age = stars_list[i].age

    # If Flag is set we do not bin stars younger than the age set by max_age_unbinned_stars
    if not cfg.par.FORCE_BINNED and cfg.par.max_age_direct > minimum_age:
        minimum_age = cfg.par.max_age_direct + 0.001

    delta_age = (maximum_age-minimum_age)/cfg.par.N_STELLAR_AGE_BINS

    if delta_age <= 0: # If max age for direct adding stars is greater than the max age of stars in the galaxy then 
        return m       # exit the function since there are no stars left for binning

    # define the metallicity bins: we do this by saying that they are the number of metallicity bins in FSPS

    fsps_metals = np.array(sp.zlegend)
    N_METAL_BINS = len(fsps_metals)

    # note the bins are NOT metallicity, but rather the zmet keys in
    # fsps (i.e. the zmet column in Table 1 of the fsps manual)
    metal_bins = np.arange(N_METAL_BINS)+1

    # define the age bins in log space so that we maximise resolution around young stars
    age_bins = 10.**(np.linspace(np.log10(minimum_age),np.log10(maximum_age),cfg.par.N_STELLAR_AGE_BINS))

    #tack on the maximum age bin
    age_bins = np.append(age_bins,age_bins[-1]+delta_age)

   
    #define the mass bins (log)
    #note - for some codes, all star particles have the same mass.  in this case, we have to have a trap:
    if minimum_mass == maximum_mass or cfg.par.N_MASS_BINS == 0: 
        mass_bins = np.zeros(cfg.par.N_MASS_BINS+1)+minimum_mass
    else:
        delta_mass = (np.log10(maximum_mass)-np.log10(minimum_mass))/cfg.par.N_MASS_BINS
        mass_bins = np.arange(np.log10(minimum_mass),np.log10(maximum_mass),delta_mass)
        mass_bins = np.append(mass_bins,mass_bins[-1]+delta_mass)
        mass_bins = 10.**mass_bins
        
    print ('mass_bins = ',mass_bins)
    print ('metal_bins = ',metal_bins)
    print ('age_bins = ',age_bins)

    #has_stellar_mass is a 3D boolean array that's [wz,wa,wm] big and
    #says whether or not that bin is being used downstream for
    #creating a point source collection (i.e. that it actually has at
    #least one star cluster that falls into it)
    has_stellar_mass = np.zeros([N_METAL_BINS,cfg.par.N_STELLAR_AGE_BINS+1,cfg.par.N_MASS_BINS+1],dtype=bool)


    stars_in_bin = {} #this will be a dictionary that holds the list
    #of star particles that go in every [wz,wa,wm]
    #group.  The keys will be tuples that hold a
    #(wz,wa,wm) set that we will then use later to
    #speed up adding sources.
    
    for i in range(nstars):
        # If a particle was directly added using direct_add_stars() then it is skipped over.
        if not cfg.par.FORCE_BINNED and stars_list[i].age <= cfg.par.max_age_direct:
            continue
            
        wz = find_nearest(metal_bins,stars_list[i].fsps_zmet)
        wa = find_nearest(age_bins,stars_list[i].age)
        wm = find_nearest(mass_bins,stars_list[i].mass)
        
        stars_list[i].sed_bin = [wz,wa,wm]
        has_stellar_mass[wz,wa,wm] = True

        if (wz,wa,wm) in stars_in_bin:
            stars_in_bin[(wz,wa,wm)].append(i)
        else:
            stars_in_bin[(wz,wa,wm)] = [i]


    print ('assigning stars to SED bins')
    sed_bins_list=[]
    sed_bins_list_has_stellar_mass = []
       
    #we loop through age bins +1 because the max values were tacked
    #onto those bin lists. but for metal bins, this isn't the case, so
    #we don't loop the extra +1
    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                sed_bins_list.append(Sed_Bins(mass_bins[wm],fsps_metals[wz],age_bins[wa],metal_bins[wz]))
                if has_stellar_mass[wz,wa,wm] == True:
                    stars_metals = []
                    for star in stars_in_bin[(wz,wa,wm)]:
                        stars_metals.append(stars_list[star].all_metals)
                    stars_metals = np.array(stars_metals)
                    stars_metals = np.mean(stars_metals,axis=0)
                    sed_bins_list_has_stellar_mass.append(Sed_Bins(mass_bins[wm],fsps_metals[wz],age_bins[wa],metal_bins[wz],stars_metals))
   
    #sed_bins_list is a list of Sed_Bins objects that have the
    #information about what mass bin, metal bin and age bin they
    #correspond to.  It is unnecessary, and heavy computational work
    #to re-create the SED for each of these bins - rather, we can just
    #calculate the SED for the bins that have any actual stellar mass.
            
    print ('Running SPS for Binned SEDs')
    print ('calculating the SEDs for ',len(sed_bins_list_has_stellar_mass),' bins')
    
    binned_stellar_nu,binned_stellar_fnu_has_stellar_mass,disk_fnu,bulge_fnu,mfrac,binned_line_em_has_stellar_mass = sg.allstars_sed_gen(sed_bins_list_has_stellar_mass,cosmoflag,sp)
    

    #since the binned_stellar_fnu_has_stellar_mass is now
    #[len(sed_bins_list_has_stellar_mass),nlam)] big, we need to
    #transform it back to the a larger array.  this is an ugly loop
    #that could probably be prettier...but whatever.  this saves >an
    #order of magnitude in time in SED gen.  
    nlam = binned_stellar_nu.shape[0]
    cloudy_nlam = len(np.genfromtxt(cfg.par.pd_source_dir + "/powderday/nebular_emission/data/refLines.dat", delimiter=','))
    binned_stellar_fnu = np.zeros([len(sed_bins_list),nlam])
    binned_mfrac = np.zeros([len(sed_bins_list)])
    binned_line_em = np.zeros([len(sed_bins_list),cloudy_nlam])

    counter = 0
    counter_has_stellar_mass = 0
    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                if has_stellar_mass[wz,wa,wm] == True:
                    binned_mfrac[counter] = mfrac[counter_has_stellar_mass]
                    binned_stellar_fnu[counter,:] = binned_stellar_fnu_has_stellar_mass[counter_has_stellar_mass,:]
                    binned_line_em[counter,:] = binned_line_em_has_stellar_mass[counter_has_stellar_mass,:]
                    counter_has_stellar_mass += 1 
                counter+=1

    print(f'after selecting for ones with stellar mass: {np.shape(binned_stellar_fnu)}')


    
    '''
    #DEBUG trap for nans and infs
    if np.isinf(np.sum(binned_stellar_nu)):  pdb.set_trace()
    if np.isinf(np.sum(binned_stellar_fnu)): pdb.set_trace()
    if np.isnan(np.sum(binned_stellar_nu)): pdb.set_trace()
    if np.isnan(np.sum(binned_stellar_fnu)): pdb.set_trace()
    '''

    #now binned_stellar_nu and binned_stellar_fnu are the SEDs for the bins in order of wz, wa, wm 
    
    #create the point source collections: we loop through the bins and
    #see what star particles correspond to these.  if any do, then we
    #add them to a list, and create a point source collection out of
    #these


    print ('adding point source collections')
    t1=datetime.now()


    totallum = 0 
    totalmass = 0 
    counter=0
    fnu_arr = []
    pos_arr = []

    for wz in range(N_METAL_BINS):
        for wa in range(cfg.par.N_STELLAR_AGE_BINS+1):
            for wm in range(cfg.par.N_MASS_BINS+1):
                
                if has_stellar_mass[wz,wa,wm] == True:
                
                    source = m.add_point_source_collection()                    
                    
                    nu = binned_stellar_nu
                    fnu = binned_stellar_fnu[counter,:]

                    nu,fnu = wavelength_compress(nu,fnu,df_nu)

                    #reverse for hyperion
                    nu = nu[::-1]
                    fnu = fnu[::-1]

                    
                    #source luminosities
                    #here, each (wz, wa, wm) bin will have an associated mfrac that corresponds to the fnu generated for this bin
                    #while each star particle in the bin has a distinct mass, they all share mfrac as this value depends only on the age and Z of the star
                    #thus, there are 'counter' number of binned_mfrac values (to match the number of fnu arrays)
                    lum = np.array([stars_list[i].mass/constants.M_sun.cgs.value*constants.L_sun.cgs.value/binned_mfrac[counter] for i in stars_in_bin[(wz,wa,wm)]])
                    lum *= np.absolute(np.trapz(fnu,x=nu))
                    source.luminosity = lum

                    pagb = cfg.par.add_pagb_stars and cfg.par.PAGB_min_age <= age_bins[wa] <= cfg.par.PAGB_max_age
                    young_star = cfg.par.add_young_stars and cfg.par.HII_min_age <= age_bins[wa] <= cfg.par.HII_max_age
                    
                    for i in stars_in_bin[(wz,wa,wm)]:  totalmass += stars_list[i].mass
                    
                    #source positions
                    pos = np.zeros([len(stars_in_bin[(wz,wa,wm)]),3])

                    for i in range(len(stars_in_bin[(wz,wa,wm)])):
                        pos[i,:] = stars_list[stars_in_bin[(wz,wa,wm)][i]].positions
                        age_star = stars_list[stars_in_bin[(wz,wa,wm)][i]].age

                        if (cfg.par.add_neb_emission) and (young_star or pagb):
                            pos_arr.append(pos[i, :])
                            fnu_arr.append(binned_stellar_fnu[counter,:])
                            line_em = binned_line_em[counter,:]
                            # the stellar population returns the calculation in units of Lsun/1 Msun: 
                            # https://github.com/dfm/python-fsps/issues/117#issuecomment-546513619
                            line_em = line_em * (stars_list[stars_in_bin[(wz,wa,wm)][i]].mass * units.g).to(units.Msun).value * 3.839e33 # Units: ergs/s
                            OH = stars_list[stars_in_bin[(wz,wa,wm)][i]].all_metals[4]
                            line_em = np.append(line_em, OH)

                            if young_star:
                                line_em = np.append(line_em, 1)
                            elif pagb:
                                line_em = np.append(line_em, 2)

                            dump_emlines(line_em)

                            
                    source.position=pos
                    #source spectrum
                    source.spectrum = (nu,fnu)
                                    
                    totallum += np.sum(source.luminosity)

                    '''
                    if np.isnan(lum): 
                        print 'lum is a nan in point source collection addition. exiting now.'
                        sys.exit()
                    if np.isinf(lum): 
                        print 'lum is an inf in point source collection addition. exiting now.'
                        sys.exit()
                    '''
                counter+=1
    
    
    if cfg.par.add_neb_emission and (cfg.par.SAVE_NEB_SEDS or cfg.par.add_DIG_neb) and (len(pos_arr) != 0):
        dump_NEB_SEDs(binned_stellar_nu, fnu_arr, pos_arr)

    if cosmoflag == False: add_bulge_disk_stars(df_nu,binned_stellar_nu,binned_stellar_fnu,disk_fnu,bulge_fnu,stars_list,diskstars_list,bulgestars_list,m)

    m.set_sample_sources_evenly(True)

    t2=datetime.now()
    print ('[source_creation/add_binned_seds:] Execution time for point source collection adding = '+str(t2-t1))
    print ('[source_creation/add_binned_seds:] Total Luminosity of point source collection is: ',totallum)

    return m


def wavelength_compress(nu,fnu,df_nu):
    
    nu_inrange = np.logical_and(nu >= min(df_nu),nu <= max(df_nu))
    nu_inrange = np.where(nu_inrange == True)[0]
  
    compressed_nu = nu[nu_inrange]
    compressed_fnu = np.asarray(fnu)[nu_inrange]
    
    #get rid of all wavelengths below lyman limit
    dum_nu = compressed_nu*units.Hz
    dum_lam = constants.c.cgs/dum_nu
    dum_lam = dum_lam.to(units.angstrom)

    wll = np.where(dum_lam.value >= 912)[0] #where are lambda is above the lyman limit
    compressed_nu = compressed_nu[wll] 
    compressed_fnu = compressed_fnu[wll]
    
    return compressed_nu, compressed_fnu
    

def BH_source_add(m,reg,df_nu,boost):

    print("--------------------------------\n")
    print("Adding Black Holes to Source List in source_creation\n")
    print("--------------------------------\n")
 
    try:
        nholes = reg["bh","sed"].shape[0]

    except: 
        print('BH source creation failed. No BH found')
        nholes = 0


    if nholes == 0:
        print('BH source creation failed. No BH found')
    
    else:
        #temporary wavelength compress just to get the length of the
        #compressed nu for a master array
        dumnu,dumfnu = wavelength_compress(reg["bh","nu"].value,reg["bh","sed"][0,:].value,df_nu)
        master_bh_fnu = np.zeros([nholes,len(dumnu)])
         
        width = reg.right_edge-reg.left_edge

        agn_ids = []

        for i in range(nholes):  
            #since the BH was added in a front end, we don't know
            #if the hole is in the actual cut out region of the yt
            #dataset.  so we need to filter out any holes that
            #might not be in the simulation domain or have no luminosity.

            if ((reg["bh","coordinates"][i,0].in_units('kpc') <  (reg.center[0].in_units('kpc')+(0.5*width[0].in_units('kpc'))))
            and
            (reg["bh","coordinates"][i,0].in_units('kpc') >  (reg.center[0].in_units('kpc')-(0.5*width[0].in_units('kpc'))))
            and
            (reg["bh","coordinates"][i,1].in_units('kpc') <  (reg.center[1].in_units('kpc')+(0.5*width[1].in_units('kpc'))))
            and
            (reg["bh","coordinates"][i,1].in_units('kpc') >  (reg.center[1].in_units('kpc')-(0.5*width[1].in_units('kpc'))))
            and
            (reg["bh","coordinates"][i,2].in_units('kpc') <  (reg.center[2].in_units('kpc')+(0.5*width[2].in_units('kpc'))))
            and
            (reg["bh","coordinates"][i,2].in_units('kpc') >  (reg.center[2].in_units('kpc')-(0.5*width[2].in_units('kpc'))))
            and 
            (reg["bh","luminosity"][i].value > 0)):

                agn_ids.append(i)
        if len(agn_ids) == 0:
            return


        print ('Number AGNs in the cutout with non zero luminositites: ', len(agn_ids))

        fnu_arr = sg.get_agn_seds(agn_ids, reg)
        nu = reg["bh","nu"].value

        for j in range(len(agn_ids)):
            i = agn_ids[j]
            fnu = fnu_arr[j,:]
            nu, fnu = wavelength_compress(nu,fnu,df_nu)

            master_bh_fnu[i,:] = fnu

            print('Boosting BH Coordinates and adding BH #%d to the source list now'%i)
            # the tolist gets rid of the array brackets
            bh = m.add_point_source(luminosity = reg["bh","luminosity"][i].value.tolist(),
                                    spectrum = (nu,fnu),
                                    position = (reg["bh","coordinates"][i,:].in_units('cm').value-boost).tolist())

        dump_AGN_SEDs(nu,master_bh_fnu,reg["bh","luminosity"].value)


def DIG_source_add(m,reg,df_nu,boost):

    print("--------------------------------")
    print("Adding DIG to Source List in source_creation")
    print("--------------------------------")

    print ("Getting the specific energy dumped in each grid cell")

    try:
        rtout = cfg.model.outputfile + '_DIG_energy_dumped.sed'
        try: 
            grid_properties = np.load(cfg.model.PD_output_dir+"/grid_physical_properties."+cfg.model.snapnum_str+'_galaxy'+cfg.model.galaxy_num_str+".npz")
        except:
            grid_properties = np.load(cfg.model.PD_output_dir+"/grid_physical_properties."+cfg.model.snapnum_str+".npz")

        cell_info = np.load(cfg.model.PD_output_dir+"/cell_info."+cfg.model.snapnum_str+"_"+cfg.model.galaxy_num_str+".npz")
    except:
        print ("ERROR: Can't proceed with DIG nebular emission calculation. Code is unable to find the required files.") 
        print ("Make sure you have the rtout.sed, grid_physical_properties.npz and cell_info.npz for the corresponding galaxy.")

        return 

    m_out = ModelOutput(rtout)
    oct = m_out.get_quantities()
    grid = oct
    order = find_order(grid.refined)
    refined = grid.refined[order]
    quantities = {}
    for field in grid.quantities:
        quantities[('gas', field)] = grid.quantities[field][0][order][~refined]

    cell_width = cell_info["fw1"][:,0]
    mass = (quantities['gas','density']*units.g/units.cm**3).value * (cell_width**3)
    specific_energy = (quantities['gas','specific_energy']*units.erg/units.s/units.g).value
    specific_energy = (specific_energy * mass) # in ergs/s

    mask = np.where(mass != 0)[0] # Masking out all grid cells that have no gas mass
    pos = cell_info["fc1"][mask] - boost
    cell_width = cell_width[mask]
    met = grid_properties["grid_gas_metallicity"][:, mask]
    met = np.transpose(met)

    specific_energy = specific_energy[mask]

    mask = []
    logU_arr = []
    lam_arr = []
    fnu_arr = []
  
    sed_file_name = cfg.model.PD_output_dir+"neb_seds_galaxy_"+cfg.model.galaxy_num_str+".npz"
    data = np.load(sed_file_name)
    nu = data['nu']
    stars_fnu = data["fnu"]
    star_coordinates = data["positions"]
    
    print ("[Source creation (DIG)] Generating KD Tree")
    tree = spatial.KDTree(star_coordinates)
    
    print("----------------------------------------------------------------------------------")
    print("[Source creation (DIG)] Calculating ionizing spectrum for " + str(len(cell_width))+ " gas cells ")
    print("----------------------------------------------------------------------------------")

    for i in tqdm(range(len(cell_width))):
        
        # Getting shape of the incident spectrum by taking a distance weighted average of the CLOUDY output spectrum of nearby stars.
        # lam in Angstrom, fnu in Lsun/Hz
        lam, fnu = get_DIG_sed_shape(pos[i], cell_width[i], nu, stars_fnu, tree)

        # If the gas cell has no young stars within a specified distance (stars_max_dist) then skip it.
        if len(np.atleast_1d(fnu)) == 1:
            continue
        
        # Calulating the ionization parameter by extrapolating the specfic energy beyind the lyman limit 
        # using the SED shape calculated above.
        logU = get_DIG_logU(lam, fnu, specific_energy[i], cell_width[i])

        lam_arr.append(lam)
        fnu_arr.append(fnu)
        
        # Only cells with ionization parameter greater than the parameter DIG_min_logU are considered 
        # for nebular emission calculation. This is done so as to speed up the calculation by ignoring 
        # the cells that do not have enough energy to prduce any substantial emission
        
        if logU > cfg.par.DIG_min_logU:
            mask.append(i)
            fnu_arr.append(fnu)
            logU_arr.append(logU)
        
    cell_width = cell_width[mask]
    met = met[mask]
    fnu_arr = np.array(fnu_arr)
    logU = np.array(logU_arr)
    pos = pos[mask]


    if (len(mask)) == 0:
        print ("No gas particles fit the criteria for calculating DIG. Skipping DIG calculation")
        return

    print("----------------------------------------------------------------------------------")
    print("[Source creation (DIG)] Calculating nebular emission for " + str(len(mask)) + " gas cells that fit the criteria")
    print("----------------------------------------------------------------------------------")
    
    fnu_arr_neb = sg.get_dig_seds(lam, fnu_arr, logU, cell_width, met)
    for i in range(len(logU)):
        nu = 1.e8 * constants.c.cgs.value / lam
        fnu = fnu_arr_neb[i,:]
        nu, fnu = wavelength_compress(nu,fnu,df_nu)
        
        nu = nu[::-1]
        fnu = fnu[::-1]

        lum = np.trapz(fnu,x=nu)*constants.L_sun.cgs.value
        
        source = m.add_point_source()
        source.luminosity = lum # [ergs/s]
        source.spectrum = (nu,fnu)
        source.position = pos[i] # [cm]
