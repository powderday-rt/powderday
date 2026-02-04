import numpy as np
import os,h5py,pdb
import powderday.config as cfg
from astropy import units as u
from astropy import constants as const
from powderday.pah.pah_file_read import read_draine_file
from scipy.interpolate import interp1d,interp2d
from scipy.optimize import nnls
from tqdm import tqdm

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()

def get_Cabs(draine_directories,simulation_sizes,gsd):
    
    ncells = simulation_sizes.shape[0]
    n_simulation_sizes = simulation_sizes.shape[1]

    for file in os.listdir(draine_directories[0]):
        if file.startswith("iout_graD") and file.endswith("_0.00") and "ib" in file: Cabsfile_cation = draine_directories[0]+'/'+file
        if file.startswith("iout_graD") and file.endswith("_0.00") and "nb" in file: Cabsfile_neutral = draine_directories[0]+'/'+file
    PAH_list_cation = read_draine_file(Cabsfile_cation)    
    PAH_list_neutral = read_draine_file(Cabsfile_neutral)
 
   #i'm sure there's some awesome pythonic cool kids way to code this.
    #but i'm an IDL coder at heart.
    Cabs_cation = np.zeros([len(PAH_list_cation[0].size_list),len(PAH_list_cation[0].lam)])
    for i in range(len(PAH_list_cation)):
        #this works bc PAH_list is an n_sizes long list of PAH objects
        Cabs_cation[i,:] = PAH_list_cation[i].cabs

    Cabs_neutral = np.zeros([len(PAH_list_neutral[0].size_list),len(PAH_list_neutral[0].lam)])
    for i in range(len(PAH_list_neutral)):
        #this works bc PAH_list is an n_sizes long list of PAH objects
        Cabs_neutral[i,:] = PAH_list_neutral[i].cabs

        
    #We now have Cabs (for neutrals and cations) in terms of (size,nu),
    #and we want to get these in terms of just (nu). To do this, we need
    #to regrid (i.e. downsample) the first dimension so that the GSD is at
    #the same sizes as the simulation that was run.  Then we can just
    #multiply Cabs(size_regridded,nu) by dN, where dN is the number of
    #grains at a given size.
    
    # we can just take the very first size of the the cation PAH_list - they should all be the same.
    draine_sizes = PAH_list_cation[0].size_list*u.cm
    simulation_sizes = simulation_sizes.to(u.cm)
    draine_lam = PAH_list_cation[0].lam*u.micron

    '''
    #Cabs is in dimensions of draine_sizes,draine_lam.
    #resampling Cabs to a shape of simulation_sizes,draine_lam -- this
    #will allow us to deal with the fact that Cabs is given at individual
    #grain sizes (by multiplying by simulation GSD), while also keeping the draine_lam axis available for integration later.
    '''

    
    f_2d_interp_cation = interp2d(draine_sizes.value,draine_lam.value,Cabs_cation.T,kind = 'cubic')
    f_2d_interp_neutral = interp2d(draine_sizes.value,draine_lam.value,Cabs_neutral.T,kind = 'cubic')
    #Cabs_cation_regrid = f_2d_interp_cation(simulation_sizes.value,simulation_isrf_lam.value).T
    #Cabs_neutral_regrid = f_2d_interp_neutral(simulation_sizes.value,simulation_isrf_lam.value).T


    #Cabs regrid arrays are (n_sizes,nwavelengths,ncells) big
    Cabs_cation_regrid_sizes_lam_cells = np.empty([n_simulation_sizes,draine_lam.shape[0],ncells])
    Cabs_neutral_regrid_sizes_lam_cells = np.empty([n_simulation_sizes,draine_lam.shape[0],ncells])
    Cabs_cation_regrid_lam_cells = np.empty([draine_lam.shape[0],ncells])
    Cabs_neutral_regrid_lam_cells = np.empty([draine_lam.shape[0],ncells])
    gsd_normalized = np.empty([n_simulation_sizes,ncells])

    print("[pah/isrf_decompose]: resampling Cabs from the Draine size arrays to the simulation size arrays")

    for i in tqdm(range(ncells)):
        Cabs_cation_regrid_sizes_lam_cells[:,:,i] = f_2d_interp_cation(simulation_sizes[i,:],draine_lam.value).T
        Cabs_neutral_regrid_sizes_lam_cells[:,:,i] = f_2d_interp_neutral(simulation_sizes[i,:],draine_lam.value).T
    
        #we now have to get Cabs(nu) -- right now we have Cabs(size,nu).
        #To do this, we should ensure that the GSD is at the same sizes as the
        #draine files (and if not, regrid so that it is), and then multiply
        #Cabs(size,nu) * dN, where dN is the number of grains at a given size.
        
        #we normalize the grain size distribution by the integral over the grain sizes so that this is normalized for the dot product into Cabs
        

        gsd_normalized[:,i] = gsd[i,:]/np.trapz(gsd[i,:],simulation_sizes[i,:])


        Cabs_cation_regrid_lam_cells[:,i] = np.dot(gsd_normalized[:,i],Cabs_cation_regrid_sizes_lam_cells[:,:,i]) # now in terms of just wavelength,ncells (at simulation_isrf_lam wavelengths)
        Cabs_neutral_regrid_lam_cells[:,i] = np.dot(gsd_normalized[:,i],Cabs_neutral_regrid_sizes_lam_cells[:,:,i])
    
    #get the units right since they're not faithfully followed throughout

    Cabs_cation_regrid_lam_cells = (Cabs_cation_regrid_lam_cells)*u.cm**2
    Cabs_neutral_regrid_lam_cells = (Cabs_neutral_regrid_lam_cells)*u.cm**2


    #DEBUG DEBUG DEBUG - these 5 lines are garbage and remove them.
    #they're intended to only grab Cabs for the 0.1micron value since
    #that's whats in the draine paper.  we may wantto keep something
    #like this at the end, but not in it's current form.
    #Cabs_cation_regrid_lam_cells = Cabs_cation_regrid_lam_cells.value
    #Cabs_neutral_regrid_lam_cells = Cabs_neutral_regrid_lam_cells.value
    #for i in range(Cabs_cation_regrid_lam_cells.shape[1]):
    #    Cabs_cation_regrid_lam_cells[:,i] = Cabs_cation[19,:]
    #    Cabs_neutral_regrid_lam_cells[:,i] = Cabs_neutral[19,:]
        
    return Cabs_cation_regrid_lam_cells,Cabs_neutral_regrid_lam_cells



def get_beta_nnls(draine_directories, gsd, simulation_sizes, reg):

    #get the wavelengths of the simulation 
    f = h5py.File(cfg.model.outputfile + '_isrf.sed')

    #thiis gives us the list of iterations in the initial ISRF calculation.  
    iteration_list = [i for i in f.keys() if 'iteration_' in i]
    dset = f[iteration_list[-1]]
    simulation_isrf_nu = dset['ISRF_frequency_bins'][:] * u.Hz
    simulation_isrf_lam = (const.c/simulation_isrf_nu).to(u.micron)


    simulation_specific_energy_sum = dset['specific_energy_nu']*u.erg/u.s/u.g #is [n_nu, n_dust, n_cells] big  

    #get the simulation_isrf in units of erg/s
    grid_dust_masses = reg['dust','mass'].in_units('g').to_astropy() #getting the dust masses out of yt units and into astropy units
    simulation_specific_energy_sum *= grid_dust_masses.cgs #now in erg/s

    #clip values that are MC noise too high
    simulation_specific_energy_sum[simulation_specific_energy_sum.value > 1.e50] = np.median(simulation_specific_energy_sum)


    ncells = grid_dust_masses.shape[0]
        
    #convolve the simulation specific energy (ISRF) with the GSD to
    #get rid of the size dimension:
    simulation_specific_energy_gsd_convolved = np.zeros([simulation_specific_energy_sum.shape[0],simulation_specific_energy_sum.shape[2]])

    
    print("[isrf_decompose/get_beta_nnls]: Convolving the simulation specific energy grid with the dust types")
    for i in tqdm(range(ncells)):
        #x = simulation_specific_energy_sum[:,:,i]
        simulation_specific_energy_gsd_convolved[:,i] = np.dot(simulation_specific_energy_sum[:,:,i],gsd[i,:])
        simulation_specific_energy_gsd_convolved[:,i]/=np.sum(gsd[i,:])

    simulation_specific_energy_gsd_convolved *= u.erg/u.s #attach units back to it
    
    #we have read in the draine directories explicitly to ensure that the ordering of them is identical from pah_source_create
    isrf_files = []
    for directory in draine_directories:
        for file in os.listdir(directory):
            if file.startswith("isrf"):
                isrf_files.append(file)


    '''#note - this bit isn't formally needed, and even still, it reads
    in 2x iout files (compared to the isrf files) since there's a Cabs for
    the cation state of PAHs, and one for the neutral


    iout_U0_files = [] #just saving the U=0 files since we just need them to grab C_abs
    for directory in draine_directories:
    for file in os.listdir(directory):
        if file.startswith("iout_graD") and file.endswith("_0.00"):
            print(directory,file)
            iout_U0_files.append(file)
    
    '''
        

    #get the length of a basis ISRF vector
    data = np.loadtxt(draine_directories[0]+'/'+isrf_files[0],skiprows=7,usecols=(0,1))
    nlam = len(data[:,0])
    draine_lam = data[:,0]*u.micron
    basis_isrf_vectors = np.zeros([len(draine_directories),nlam])
 


    for counter, (directory,file) in enumerate(zip(draine_directories,isrf_files)):
        data = np.loadtxt(directory+'/'+file,skiprows=7,usecols=(0,1))
        basis_isrf_vectors[counter] = data[:,1]
        
    #add the units (as listed in the files)
    basis_isrf_vectors *= u.erg/u.cm**3

    #the draine vectors are in erg/cm**3 density.  we employ u_nu
    #(erg/cm^3) * c/4pi = I_nu (erg/s/cm*2/Hz) to get I_nu.  then we
    #multiply by an
    #arbitrary constant (1) to get rid of the cm^2.  the reason we can do
    #that is taht we only want the *relative* contributions of the basis
    #vectors to the local ISRF.  the normalization will get set later by
    #the grain size distribution anyways.

    basis_isrf_vectors *= const.c/(4.*np.pi)  #erg/s/cm**2/Hz
    basis_isrf_vectors = basis_isrf_vectors.to(u.erg/u.s/u.cm**2)
    basis_isrf_vectors *= 1*u.cm**2 #erg/s/Hz
    

    #4 now resample the hyperion ISRF to the wavelengths of the Draine
    #basis ISRFs so that we can NNLS
    f_1d_interp_lam = interp1d(simulation_isrf_lam.to(u.micron).value,simulation_specific_energy_gsd_convolved.cgs.T.value,kind='cubic')
    simulation_specific_energy_sum_regrid = f_1d_interp_lam(draine_lam.to(u.micron).value).T

    #the regridding can turn some wavelengths where there was 0 emission
    #(from too low photon count simulations) to negative, so we zero these
    #back out.
    simulation_specific_energy_sum_regrid[simulation_specific_energy_sum_regrid < 0] = 0

    #in the interpolation we lost our units, so lets get them back
    simulation_specific_energy_sum_regrid *= u.erg/u.s


    
    
    
    simulation_sizes = np.broadcast_to(simulation_sizes,(ncells,simulation_sizes.shape[0]))*u.micron
    gsd = gsd.value
    
    
    Cabs_cation_regrid,Cabs_neutral_regrid = get_Cabs(draine_directories,simulation_sizes,gsd)
     
    if cfg.par.SKIP_LOGU_CALC == False:
        logU_grid = get_logU(simulation_specific_energy_sum_regrid,Cabs_cation_regrid,Cabs_neutral_regrid,draine_lam,reg)
    else:
        print("[pah/isrf_decompose:] SKIP_LOGU_CALC is set to True: Assuming logU across the grid is 0")
        logU_grid = np.zeros(ncells)








    #5. then nnls!  with nnls, we can then sum the PAH components for each
    #cell accordingly.  note - because the ISRF computed from hyperion has
    #the infrared component saved, we need to cut off our ISRF for both
    #the hyperion model and draine basis functions at some wavelength
    #before thermal IR emission gets big, like maybe 10 micron.  also may
    #be useful to renormalize things so that we don't have 10s of orders
    #of mag difference bewteen the ISRF field and basis vectors.
    
    beta_nnls = np.zeros([basis_isrf_vectors.shape[0],simulation_specific_energy_sum_regrid.shape[1]])
    ncells = simulation_specific_energy_sum_regrid.shape[1]



    x = basis_isrf_vectors
    y = simulation_specific_energy_sum_regrid

    #cut off everything after 1 micron
    idx = (np.abs(draine_lam.to(u.micron).value - 1)).argmin()
    x = x[:,0:idx]
    y = y[0:idx,:]
    


    
    for i in tqdm(range(ncells)):
        beta_nnls[:,i] = nnls(x.T,y[:,i])[0]
        isrf_lum = np.trapz(simulation_specific_energy_sum_regrid[:,i]/draine_lam,draine_lam)
        nnls_lum = np.trapz(np.dot(x.T,beta_nnls[:,i])/draine_lam[0:idx],draine_lam[0:idx])
        beta_nnls[:,i]*=isrf_lum.value/nnls_lum.value
    

    return beta_nnls,logU_grid


def get_logU(cell_isrf,Cabs_cation,Cabs_neutral,draine_lam,reg):
    #get the cell sizes since ISRF is not per cm^2 (and it needs to be
    #so that we can convert to erg/cm^3)

    
    #for some reason the units attached in reg.parameters needs to be reattached for it to carry through in the upcoming arithmetic
    cell_sizes = reg.parameters['cell_size'].value*u.cm


    cell_isrf = (cell_isrf/(cell_sizes**2.))
    

    #get into erg/cm^3 
    cell_isrf /=const.c
    cell_isrf = cell_isrf.to(u.erg/u.cm**3)

    h_ref = 1.958e-12*u.erg/u.s
    draine_nu = (const.c/draine_lam).to(u.Hz)



    #DEBUG DEBUG DEBUG
    #cell_isrf = np.ones(cell_isrf.shape)*1.043e-12*u.erg/u.cm**3


    
    cell_isrf = (cell_isrf.T/draine_nu[::-1]).T
    
    #eq. 5 from Draine et al. 2021, 917, 3, ApJ
    #U = int(d_nu u * c* C_abs)/h_ref


    print('[pah/isrf_decompose/get_logU:] Computing logU for PAH calculation')

    y = cell_isrf * const.c*Cabs_neutral/h_ref
    U= (np.trapz(y,draine_nu[::-1],axis=0)).decompose()
    
    #DEBUG DEBUG DEBUG
    '''
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.hist(np.log10(U.value).flatten())
    ax.set_yscale('log')
    fig.savefig('hist_logu.png',dpi=300)
    '''
    
    #just to make any numerical isseues with U<~0 not impact our logU calc
    U[U<=0] = 1.e-10
    logU = np.log10(U)

    
    return logU

    
