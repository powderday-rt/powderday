#FIELD NOMENCLATURE

#For the gadget front end, we have two types of fields: particle, and
#mesh.  the former is read in, defined, and then deposited onto the
#octree to create the mesh quantities.  The particle information is
#given the tuples (e.g.) ('gas','metallicity'), or ('dust','mass'),
#while the octree is given ('dust','smoothedmasses').  Note - there is
#one side rule to this: if we have PartType3 dust particles in the
#simulations, then those particles may be given the specific name
#('particle_dust','mass').

from __future__ import print_function
import numpy as np
import yt
#from yt.fields.particle_fields import add_volume_weighted_smoothed_field
import powderday.config as cfg
from powderday.mlt.dgr_extrarandomtree_part import dgr_ert

def gadget_field_add(fname, bounding_box=None, ds=None,add_smoothed_quantities=True):

    def _gassmoothinglength(field,data):
        return data[('PartType0', 'SmoothingLength')].in_units('pc')

    def _starmetals_00(field, data):
        el_dict = {'He': '01',
                   'C': '02',
                   'N': '03',
                   'O': '04',
                   'Ne': '05',
                   'Mg': '06',
                   'Si': '07',
                   'S': '08',
                   'Ca': '09',
                   'Fe': '10'}
        el_str = field.name[1]
        if '_' in el_str:
            el_name = field.name[1][field.name[1].find('_')+1:]
            el_num = el_dict[el_name]
        else:
            el_num = '00'
        return data[('PartType4', 'Metallicity_'+el_num)]

    def _starmetals(field, data):
        return data[('PartType4', 'Metallicity')]

    def _starcoordinates(field, data):
        return data[('PartType4', 'Coordinates')]

    def _starformationtime(field, data):
        return data[('PartType4', 'StellarFormationTime')]

    def _starmasses(field, data):
        return data[("PartType4", "Masses")]

    def _diskstarcoordinates(field, data):
        return data[('PartType2', 'Coordinates')]

    def _diskstarmasses(field, data):
        return data[("PartType2", "Masses")]

    def _bulgestarcoordinates(field, data):
        return data[('PartType3', 'Coordinates')]

    def _bulgestarmasses(field, data):
        return data[("PartType3", "Masses")]

    def _gasdensity(field, data):
        return data[('PartType0', 'Density')]

    def _gasmetals_00(field, data):
        el_dict = {'He': '01',
                   'C': '02',
                   'N': '03',
                   'O': '04',
                   'Ne': '05',
                   'Mg': '06',
                   'Si': '07',
                   'S': '08',
                   'Ca': '09',
                   'Fe': '10'}
        el_str = field.name[1]
        if '_' in el_str:
            el_name = field.name[1][field.name[1].find('_')+1:]
            el_num = el_dict[el_name]
        else:
            el_num = '00'
    
        return data[('PartType0', 'Metallicity_'+el_num)]

    def _gasmetals(field, data):
        return data[('PartType0', 'Metallicity')]

    def _gascoordinates(field, data):
        return data[('PartType0', 'Coordinates')]

    def _gasmasses(field, data):
        return data[('PartType0', 'Masses')]

    def _gasfh2(field, data):
        try: return data[('PartType0', 'FractionH2')]
        except: return data[('PartType0', 'metallicity')]*0.

    def _gassfr(field, data):
        return data[('PartType0', 'StarFormationRate')]
        
    def _gassmootheddensity(field, data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][('PartType0', 'density')]
        else:
            return data[("deposit","PartType0_smoothed_density")]
            
    def _gassmoothedmetals(field, data):
        if float(yt.__version__[0:3]) >= 4:
            try:
                el_str = field.name[1]
                if '_' in el_str:
                    el_name = field.name[1][field.name[1].find('_')+1:]+"_"
                else:
                    el_name = ""
                return data.ds.parameters['octree'][('PartType0', el_name+'metallicity')]
            except:
                return data.ds.parameters['octree'][('PartType0', 'metallicity')]
        else:
            try:
                el_str = field.name[1]
                if '_' in el_str:
                    el_name = field.name[1][field.name[1].find('_')+1:]+"_"
                else:
                    el_name = ""
                return data[("deposit","PartType0_smoothed_"+el_name+"metallicity")]
            except:
                return data[("deposit","PartType0_smoothed_metallicity")]

    def _gassmoothedmasses(field, data):
        if float(yt.__version__[0:3]) >= 4:
            return data.ds.parameters['octree'][('PartType0', 'Masses')]
        else:
            return data[('deposit', 'PartType0_mass')]

    def _metaldens_00(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity_00"])

    def _metaldens(field, data):
        return (data["PartType0", "Density"]*data["PartType0", "Metallicity"])

    def _metalmass_00(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity_00"].value))

    def _metalmass(field, data):
        return (data["PartType0", "Masses"]*(data["PartType0", "Metallicity"].value))

    def _metalsmoothedmasses(field, data):
        if float(yt.__version__[0:3]) >= 4:
            return (data.ds.parameters['octree'][('PartType0', 'Masses')]* data.ds.parameters['octree'][('PartType0','metallicity')])
        else:
            return (data[('deposit', 'PartType0_smoothed_metalmass')].value)

        
    def _dustmass_manual(field, data):
        if cfg.par.otf_extinction == True:
            return (data.ds.arr(data[("PartType3", "Masses")].value, 'code_mass'))
        else:
            return (data.ds.arr(data[("PartType0", "Dust_Masses")].value, 'code_mass'))

    def _dustmass_dtm(field,data):
        return (data["PartType0","metalmass"]*cfg.par.dusttometals_ratio)

    def _li_ml_dustmass(field,data):
        return (data.ds.arr(data.ds.parameters['li_ml_dustmass'].value,'code_mass'))
    
    def _dustmass_rr(field,data):
        #hard coded values from remy-ruyer table 1
        a = 2.21
        alpha = 2.02
        x_sun = 8.69
        
        x = 12.+np.log10(data["PartType0","Metallicity_00"]/cfg.par.solar * 10.**(x_sun-12.) )
        y = a + alpha*(x_sun-np.asarray(x))
        gas_to_dust_ratio = 10.**(y)
        dust_to_gas_ratio = 1./gas_to_dust_ratio
        return dust_to_gas_ratio * data["PartType0","Masses"]

    def _dustmass_li_bestfit(field,data):
        log_dust_to_gas_ratio = (2.445*np.log10(data["PartType0","Metallicity_00"]/cfg.par.solar))-(2.029)
        dust_to_gas_ratio = 10.**(log_dust_to_gas_ratio)
        return dust_to_gas_ratio * data["PartType0","Masses"]

    def _dustsmoothedmasses(field, data):
        if float(yt.__version__[0:3]) >= 4:
            if cfg.par.otf_extinction == True:
                dsm = ds.arr(data.ds.parameters['octree'][('PartType3','Masses')],'code_mass')
            else:
                dsm = ds.arr(data.ds.parameters['octree'][('PartType0','Dust_Masses')],'code_mass')
            #return (data.ds.parameters['octree'][('PartType0','Dust_Masses')])
        else:
            if cfg.par.otf_extinction == True:
                dsm = (data.ds.arr(data[("deposit", "PartType3_sum_Masses")].value,'code_mass'))
            else:
                dsm = (data.ds.arr(data[("deposit", "PartType0_sum_Dust_Masses")].value, 'code_mass'))


        return dsm



    def _li_ml_dustsmoothedmasses(field,data):
        if float(yt.__version__[0:3]) >= 4:
            return (data.ds.parameters['octree'][('PartType0', 'li_ml_dustmass')])
        else:
            return (data.ds.arr(data[("deposit", "PartType0_sum_li_ml_dustmass")].value, 'code_mass'))
        
    def _return_dust_mass(field,data):
        return data['dust','mass']

    def _stellarages(field, data):
        ad = data.ds.all_data()
        if data.ds.cosmological_simulation == False:
            simtime = data.ds.current_time.in_units('Gyr')
            simtime = simtime.value

            age = simtime-data.ds.arr(ad[('PartType4', 'StellarFormationTime')],'Gyr').value
            # make the minimum age 1 million years
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print(
                '[gadget2pd: ] Idealized Galaxy Simulation Assumed: Simulation time is (Gyr): ', simtime)
            print('--------------\n')
        else:
            yt_cosmo = yt.utilities.cosmology.Cosmology(hubble_constant=data.ds.hubble_constant,
                                                        omega_matter=data.ds.omega_matter,
                                                        omega_lambda=data.ds.omega_lambda)
            simtime = yt_cosmo.t_from_z(ds.current_redshift).in_units('Gyr').value # Current age of the universe
            scalefactor = data[('PartType4', 'StellarFormationTime')].value
            formation_z = (1./scalefactor)-1.
            formation_time = yt_cosmo.t_from_z(formation_z).in_units('Gyr').value
            age = simtime - formation_time
            # Minimum age is set to 1 Myr (FSPS doesn't work properly for ages below 1 Myr)
            age[np.where(age < 1.e-3)[0]] = 1.e-3

            print('\n--------------')
            print('[gadget2pd: ] Cosmological Galaxy Simulation Assumed: Current age of Universe is (Gyr): ', simtime)
            print('--------------\n')

        age = data.ds.arr(age, 'Gyr')
        return age
        
    def _starsmoothedmasses(field, data):
        return data[('deposit', 'PartType4_mass')]

    def _bhluminosity(field, data):
        ad = data.ds.all_data()
        mdot = ad[("PartType5", "BH_Mdot")]
        # give it a unit since usually these are dimensionless in yt

        mdot = data.ds.arr(mdot, "code_mass/code_time")


        c = yt.utilities.physical_constants.speed_of_light_cgs
        bhluminosity = (cfg.par.BH_eta * mdot * c**2.).in_units("erg/s")

        print("[front_ends/gadget2pd:] Generating the black hole luminosity")
        if cfg.par.BH_var:
            return bhluminosity * cfg.par.bhlfrac
        else:
            return bhluminosity

    def _bhcoordinates(field, data):
        return data["PartType5", "Coordinates"]

    def _bhsed_nu(field, data):
        bhluminosity = data["bh","luminosity"]
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu, bhlum = agn_spectrum(log_lum_lsun)
        # the last 4 numbers aren't part of the SED
        nu = nu[0:-4]
        nu = 10.**nu
        nu = yt.YTArray(nu, "Hz")

        return nu

    def _bhsed_sed(field, data):
        bhluminosity = data["bh","luminosity"]
        nholes = len(bhluminosity)

        # get len of nu just for the 0th hole so we know how long the vector is
        log_lum_lsun = np.log10(bhluminosity[0].in_units("Lsun"))
        nu, l_band_vec = agn_spectrum(log_lum_lsun)
        nu = nu[0:-4]
        n_nu = len(nu)

        bh_sed = np.zeros([nholes, n_nu])

        for i in range(nholes):

            log_lum_lsun = np.log10(bhluminosity[i].in_units("Lsun"))
            nu, l_band_vec = agn_spectrum(log_lum_lsun)

            l_band_vec = 10.**l_band_vec
            l_band_vec = l_band_vec[0:-4]
            for l in range(len(l_band_vec)):
                l_band_vec[l] = data.ds.quan(l_band_vec[l], "erg/s")

            bh_sed[i, :] = l_band_vec
        bh_sed = yt.YTArray(bh_sed, "erg/s")
        return bh_sed


    #functions used for OTF_Extinction
    def _dust_density(field,data):
        return data.ds.arr(data[('PartType3','Dust_Density')],'code_mass/code_length**3')

    def _size_with_units(field,data):
        return data.ds.parameters['size']
    
    def _particle_dust_numgrains(field,data):
        return data.ds.arr(data['PartType3','Dust_Size'])

    def _particle_dust_carbon_fraction(field,data):
        return (data['PartType3','Metallicity_02'])

    def _particle_dust_mass(field,data):
        #note - this is degenerate with _dustmass_manual.  we repeat
        #the field addition here in order to simplify PAH computation
        #downstream with a consistent nomenclature with other
        #(non-particle-based) front ends.
        return (data.ds.arr(data[("PartType3", "Masses")].value, 'code_mass'))

    def _particle_dust_coordinates(field,data):
        ad = data.ds.all_data()
        return (ad['PartType3','Coordinates'])


    # load the ds (but only if this is our first passthrough and we pass in fname)
    if fname != None:
        if float(yt.__version__[0:3]) >= 4:
            ds = yt.load(fname, bounding_box=bounding_box)
            
            #ds.sph_smoothing_style = "gather"
            ds.index
            ad = ds.all_data()


        else:
            ds = yt.load(fname,bounding_box=bounding_box,over_refine_factor=cfg.par.oref,n_ref=cfg.par.n_ref)
            ds.index
            ad = ds.all_data()
        
    #We're assuming that the particle type that the active dust is
    #in in PartType3 and adding it to the sph types so that it can
    #be deposited onto the octree
    if cfg.par.otf_extinction: ds._sph_ptypes = ('PartType0','PartType3')


    #if we're in the 4.x branch of yt, load up the octree for smoothing
    if float(yt.__version__[0:3]) >= 4:
        left = np.array([pos[0] for pos in bounding_box])
        right = np.array([pos[1] for pos in bounding_box])
        #octree = ds.octree(left, right, over_refine_factor=cfg.par.oref, n_ref=cfg.par.n_ref, force_build=True)
        octree = ds.octree(left,right,n_ref=cfg.par.n_ref)
        ds.parameters['octree'] = octree

    print ('BOUNDING BOX:', bounding_box, 'LEFT: ', left, 'RIGHT: ', right)

    # for the metal fields have a few options since gadget can have different nomenclatures
    ad = ds.all_data()
    if ('PartType4', 'Metallicity_00') in ds.derived_field_list:
        try:
            ds.add_field(('star','metals'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_He'), function=_starmetals_00, sampling_type='particle', units="code_metallicity")
            ds.add_field(('star','metals_C'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_N'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_O'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Ne'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Mg'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Si'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_S'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Ca'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('star','metals_Fe'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
        except:
            ds.add_field(('star','metals'), function=_starmetals_00, sampling_type='particle',units="code_metallicity")
    else:
        ds.add_field(('star','metals'), function=_starmetals, sampling_type='particle',units="code_metallicity")

    if ('PartType0', 'Metallicity_00') in ds.derived_field_list:
        try:
            ds.add_field(('gas','metals'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_He'), function=_gasmetals_00, sampling_type='particle', units="code_metallicity")
            ds.add_field(('gas','metals_C'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_N'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_O'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Ne'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Mg'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Si'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_S'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Ca'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
            ds.add_field(('gas','metals_Fe'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
        except:
            ds.add_field(('gas','metals'), function=_gasmetals_00, sampling_type='particle',units="code_metallicity")
        
        ds.add_field(('metal','dens'), function=_metaldens_00, sampling_type='particle',units="g/cm**3")
        # we add this as part type 0 (a non-general name) as it gets
        # smoothed immediately and that's all we end up using downstream
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass_00, sampling_type='particle',units="g")

    else:
        ds.add_field(('gas','metals'), function=_gasmetals, sampling_type='particle',units="code_metallicity")
        ds.add_field(('metal','dens'), function=_metaldens, sampling_type='particle',units="g/cm**3")
        ds.add_field(('PartType0', 'metalmass'), function=_metalmass, sampling_type='particle',units="g")

    #this line is deprecated and no longer used (and will throw an error in sufficiently new yt hashes)
    #metalmass_fn = add_volume_weighted_smoothed_field("PartType0", "Coordinates", "Masses",
    #                                             "SmoothingLength", "Density", "metalmass",
    #                                             ds.field_info)

    if add_smoothed_quantities == True: ds.add_field(('metal','smoothedmasses'), function=_metalsmoothedmasses, sampling_type='particle',units='code_metallicity')


    ds.add_field(('gas','masses'), function=_gasmasses, sampling_type='particle',units='g')
    ds.add_field(('gas','fh2'), function=_gasfh2, sampling_type='particle',units='dimensionless')
    ds.add_field(('gas','sfr'), function=_gassfr, sampling_type='particle',units='g/s')
    ds.add_field(('gas','smoothinglength'),function=_gassmoothinglength,sampling_type='particle',units='pc')


    # get the dust mass

    if cfg.par.dust_grid_type == 'dtm':
        ds.add_field(('dust','mass'), function=_dustmass_dtm,sampling_type='particle',units='code_mass')
    if cfg.par.dust_grid_type == 'manual':
        #if ('PartType0', 'Dust_Masses') in ds.derived_field_list:
        ds.add_field(('dust','mass'), function=_dustmass_manual, sampling_type='particle',units='code_mass')

        if cfg.par.otf_extinction:
            #we need to add this density field so that the masses can be projected onto the octree in _dustsmoothedmasses
            ds.add_field(('PartType3','density'),function=_dust_density,units='code_mass/code_length**3',sampling_type='particle')
            ds.add_deposited_particle_field(("PartType3", "Masses"), "sum")

            #just adding this so that we have access to it later for analytics
            ds.add_field(('particle_dust','numgrains'),function=_particle_dust_numgrains,units='dimensionless',sampling_type='particle')
            ds.add_field(('particle_dust','carbon_fraction'),function=_particle_dust_carbon_fraction,units='dimensionless',sampling_type='particle')
            ds.add_field(('particle_dust','mass'),function=_particle_dust_mass,units='code_mass',sampling_type='particle')
            ds.add_field(('particle_dust','coordinates'),function=_particle_dust_coordinates,units='code_length',sampling_type='particle')

        else:
            ds.add_deposited_particle_field(("PartType0", "Dust_Masses"), "sum")
            #this just saves (redundantly) for passive dust 'manual' models the dust mass in 'particle_dust','mass' tuple.
            ds.add_field(('particle_dust','mass'),function=_return_dust_mass,units='code_mass',sampling_type='particle')
    
        if add_smoothed_quantities == True: ds.add_field(('dust','smoothedmasses'), function=_dustsmoothedmasses, sampling_type='particle',units='code_mass')
            
    if cfg.par.dust_grid_type == 'rr':
        #ds.add_field(("dust','mass"),function=_dustmass_rr,sampling_type='particle',units='code_mass')
        ds.add_field(('dust','mass'), function=_dustmass_rr,sampling_type='particle',units='code_mass')
    if cfg.par.dust_grid_type == 'li_bestfit':
        ds.add_field(('dust','mass'),function=_dustmass_li_bestfit,sampling_type='particle',units='code_mass')

    #if we have the Li, Narayanan & Dave 2019 Extreme Randomized Trees
    #dust model in place, create a field for these so that
    #dust_grid_gen can use these dust masses
    if cfg.par.dust_grid_type == 'li_ml':
        #get the dust to gas ratio
        ad = ds.all_data()
        li_ml_dgr = dgr_ert(ad["PartType0","Metallicity_00"],ad["PartType0","StarFormationRate"],ad["PartType0","Masses"])
        li_ml_dustmass = ((10.**li_ml_dgr)*ad["PartType0","Masses"]).in_units('code_mass')
        #this is an icky way to pass this to the function for ds.add_field in the next line. but such is life.
        ds.parameters['li_ml_dustmass'] = li_ml_dustmass
        ds.add_field(('PartType0','li_ml_dustmass'),function=_li_ml_dustmass,sampling_type='particle',units='code_mass')
        #just adding a new field that is called 'dustmass' so that we
        #can save it later in analytics.  
        ds.add_field(("dust','mass"),function=_li_ml_dustmass,sampling_type='particle',units='code_mass')
        ds.add_deposited_particle_field(("PartType0","li_ml_dustmass"),"sum")
        if add_smoothed_quantities == True: 
            ds.add_field(("li_ml_dustsmoothedmasses"), function=_li_ml_dustsmoothedmasses, sampling_type='particle',units='code_mass')




    ds.add_field(('star','masses'), function=_starmasses, sampling_type='particle',units='g')
    ds.add_field(('star','coordinates'), function=_starcoordinates, sampling_type='particle',units='cm')
    ds.add_field(('star','formationtime'), function=_starformationtime, sampling_type='particle',units='dimensionless')
    
    ds.add_field(('stellar','ages'),function=_stellarages,sampling_type='particle',units='Gyr')
   
    if ('PartType2', 'Masses') in ds.derived_field_list:
        ds.add_field(('diskstar','masses'), function=_diskstarmasses, sampling_type='particle',units='g')
        ds.add_field(('diskstar','coordinates'), function=_diskstarcoordinates, sampling_type='particle',units='cm')

    if ('PartType3', 'Masses') in ds.derived_field_list:
        ds.add_field(('bulgestar','masses'), function=_bulgestarmasses, sampling_type='particle',units='g')
        ds.add_field(('bulgestar','coordinates'), function=_bulgestarcoordinates, sampling_type='particle',units='cm')

    if add_smoothed_quantities == True: ds.add_field(('star','smoothedmasses'), function=_starsmoothedmasses, sampling_type='particle',units='g')

    ds.add_field(('gas','density'), function=_gasdensity, sampling_type='particle',units='g/cm**3')
    # Gas Coordinates need to be in Comoving/h as they'll get converted later.
    ds.add_field(('gas','coordinates'), function=_gascoordinates, sampling_type='particle',units='cm')
    if add_smoothed_quantities == True:
        ds.add_field(('gas','smootheddensity'), function=_gassmootheddensity, sampling_type='particle',units='g/cm**3')
        ds.add_field(('gas','smoothedmasses'), function=_gassmoothedmasses, sampling_type='particle',units='g')
        #try:
        ds.add_field(('gas','smoothedmetals'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_He'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_C'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_N'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_O'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Ne'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Mg'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Si'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_S'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Ca'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        ds.add_field(('gas','smoothedmetals_Fe'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')
        #except:
        #    ds.add_field(('gassmoothedmetals'), function=_gassmoothedmetals, sampling_type='particle',units='code_metallicity')

    if cfg.par.BH_SED == True:
        if ('PartType5','BH_Mass') in ds.derived_field_list:
            nholes = len(ds.all_data()[('PartType5', 'BH_Mass')])
            print("The number of black holes is:",nholes)
            if nholes > 0:
                if cfg.par.BH_model == 'Nenkova':
                    from powderday.agn_models.nenkova import Nenkova2008
                    agn_spectrum = Nenkova2008(*cfg.par.nenkova_params).agn_spectrum
                else:
                    from powderday.agn_models.hopkins import agn_spectrum

                if cfg.par.BH_var:
                    from powderday.agn_models.hickox import vary_bhluminosity
                    cfg.par.bhlfrac = vary_bhluminosity(nholes)
                ds.add_field(("bh","luminosity"),function=_bhluminosity,sampling_type='particle',units='erg/s')
                ds.add_field(("bh","coordinates"),function=_bhcoordinates,sampling_type='particle',units="cm")
                ds.add_field(("bh","nu"),function=_bhsed_nu,sampling_type='particle',units='Hz')
                ds.add_field(("bh","sed"),function=_bhsed_sed,sampling_type='particle',units="erg/s")


            else:
                print('No black holes found (length of BH_Mass field is 0)')

                
    
    #OTF EXTINCTION
    if cfg.par.otf_extinction:

        if add_smoothed_quantities==True:

            print("==============================================\n")
            print("[front_ends/gadget2pd:] Entering OTF Field Addition\n")
            print("==============================================\n")
        

        
            #add the dust density field
            ds.add_field(('PartType3','density'),function=_dust_density,units='code_mass/code_length**3',sampling_type='particle')
            ad = ds.all_data()
            nsizes = ad['PartType3','Dust_Size'].shape[1] #number of dust size bins
            
            #now loop through th esize bins and project them each onto
            #their own octree.  we'll then, after projecting them into
            #individual grids, collate those grids back together into a
            #master size grid
            for isize in range(nsizes):
                
                ds.parameters['size'] = 0 #just to clear it out
                
                #we have to slice it before we add the field.  if you try to slice
                #in the field defintion (i.e., in _size_with_units), yt freaks
                ds.parameters['size'] = ad['PartType3','Dust_Size'][:,isize]
                
                #actually add the sliced field now.  we do this so that we can
                #prepare for depositing onto the octree. we call this a dummy size
                #since this is just there for a place holder to deposit into a dummy octree
                
                print('adding and depositing fields for dust size bin '+str(isize))
                ds.add_field(('PartType3','dummy_size_bin'+str(isize)),function=_size_with_units,sampling_type='particle',units='dimensionless',particle_type=True,force_override=True)

                #deposit onto the octree.
                if float(yt.__version__[0:3]) < 4:
                    ds.add_deposited_particle_field(('PartType3','dummy_size_bin'+str(isize)),"sum")
                    
                    print(np.sum(ad["PartType3","dummy_size_bin"+str(isize)]))
                else:
                    #just call a dummy octree so that its deposited
                    dum_size_octree = octree[('PartType3', 'dummy_size_bin'+str(isize))]
                    print(np.sum(dum_size_octree))


                #save in mater array that we lazily stuff into parameters.
                #Note this has to be done here; attempting to do this
                #downstream can cause the entire octree_of_sizes array to
                #take on the value of the -1 size bin.  this is likely due
                #to the order of operations of ds and octree construction
                #in zoom.py

                if isize == 0:
                    if float(yt.__version__[0:3]) < 4:
                        octree_of_sizes = np.zeros((ad[('deposit','PartType3_sum_dummy_size_bin'+str(isize))].shape[0],nsizes))
                    else:
                        octree_of_sizes = np.zeros((len(octree[('PartType3','dummy_size_bin'+str(isize))]),nsizes))

                if float(yt.__version__[0:3]) < 4:
                    octree_of_sizes[:,isize] = ad[('deposit','PartType3_sum_dummy_size_bin'+str(isize))]
                else:
                    octree_of_sizes[:,isize] = octree[('PartType3','dummy_size_bin'+str(isize))]



            #------------------------
            #DEBUG BLOCK

            #if cfg.par.OTF_EXTINCTION_MRN_FORCE == True:
            #    loga = np.linspace(-4,0,octree_of_sizes.shape[1])
            #    a = 10.**loga
            #    mrn_dn_da = a**(-3.5)
            
                
            #    da = [a[i+1]-a[i] for i in range(len(a)-1)]
            #    da.append(da[-1])
            #    mrn_dn = mrn_dn_da*da

            #    for i in range(octree_of_sizes.shape[0]):
            #        octree_of_sizes[i,:] = mrn_dn

            #------------------------

            octree_of_sizes[np.isnan(octree_of_sizes)] = 0 #just because the density can be zero in some extreme cases which screws up the particle deposition
            ds.parameters['octree_of_sizes']=octree_of_sizes
            ds.parameters['octree_carbon_fraction'] = octree['PartType3','Metallicity_02']

    return ds
