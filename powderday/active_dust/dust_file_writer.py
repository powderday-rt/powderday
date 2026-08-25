#Based on routines originally written by Qi Li, University of Florida

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import powderday.config as cfg
import astropy.units as u
import astropy.constants as constants
from hyperion.dust import IsotropicDust
import os


def find_nearest(array,value):
    idx = (np.abs(np.array(array)-value)).argmin()
    return idx

def Q_unit_transform(Q):
	# wavelength [micron] -> 1/wavelength [1/micron]
	# linear -> log10 
	Q[:,0] = np.log10(1/Q[:,0])
	Q[:,1] = np.log10(Q[:,1])
	return Q

def Q_interp_wl(wlen,Q0):
	Q = interpolate.interp1d(Q0[:,0],Q0[:,1],fill_value='extrapolate')
	Q_wl = Q(wlen)
	return Q_wl

def Qext_tab_load(fin_xtab = cfg.par.pd_source_dir+'/powderday/active_dust/grain_size.txt',fin_gra=cfg.par.pd_source_dir+'/powderday/active_dust/Gra_Optical/Gra_LD93_',fin_sil=cfg.par.pd_source_dir+'/powderday/active_dust/Sil_Optical/Sil_LD93_'):
	# Return
	#  xtab, tuple Qtab = (Qabs_C, Qabs_Si, Qsca_C, Qsca_Si)
	xtab = np.log10(np.loadtxt(fin_xtab))
	ntab = len(xtab)
	Qabs_C = [[] for i in range(ntab)]
	Qsca_C = [[] for i in range(ntab)]
	Qabs_Si = [[] for i in range(ntab)]
	Qsca_Si = [[] for i in range(ntab)]

	# read Qext tables into arrays
	for i in range(ntab):
		num = str(i).zfill(2)
		Qabs_C[i] = np.loadtxt(fin_gra + num,usecols=(0,1))
		Qabs_C[i] = Q_unit_transform(Qabs_C[i])
		Qabs_Si[i] = np.loadtxt(fin_sil + num,usecols=(0,1))
		Qabs_Si[i] = Q_unit_transform(Qabs_Si[i])
		Qsca_C[i] = np.loadtxt(fin_gra + num,usecols=(0,2))
		Qsca_C[i] = Q_unit_transform(Qsca_C[i])
		Qsca_Si[i] = np.loadtxt(fin_sil + num,usecols=(0,2))
		Qsca_Si[i] = Q_unit_transform(Qsca_Si[i])

	return xtab,(Qabs_C,Qabs_Si,Qsca_C,Qsca_Si)


def Qext_get(x,wlen,cfrac,xtab,Qtab):
	# Input
	#  x: array of log10(grain_radii) of the simulation
	#  wlen: wavelengths of desired extinction curve
	#  cfrac: mass fraction of graphite
	#  xtab: array of log10(grain_radii) of the optical property table
    #  Qtab: tabulated extinction coefficient
	# Return
	#  tuple (Q_absorption, Q_scatter)

	ntab = len(xtab)
	dim = (len(x),len(wlen)) # dim: wavelength x grain size
	wlen = np.log10(1. / wlen)
	print(wlen)

 	# original tables
	Qabs_C,Qabs_Si,Qsca_C,Qsca_Si = Qtab[0],Qtab[1],Qtab[2],Qtab[3]

	# extended in wlen
	Qabs_C_wl = np.zeros((ntab,dim[1]))
	Qsca_C_wl = np.zeros((ntab,dim[1]))
	Qabs_Si_wl = np.zeros((ntab,dim[1]))
	Qsca_Si_wl = np.zeros((ntab,dim[1]))

	# fully extended table (additionally extended in grain size)
	Qabs_C_2D = np.zeros(dim)
	Qsca_C_2D = np.zeros(dim)
	Qabs_Si_2D = np.zeros(dim)
	Qsca_Si_2D = np.zeros(dim)


	# align arrays across grain sizes by linear interpolation
	for i in range(ntab):
		Qabs_C_wl[i,:] = Q_interp_wl(wlen,Qabs_C[i])
		Qabs_Si_wl[i,:] = Q_interp_wl(wlen,Qabs_Si[i])
		Qsca_C_wl[i,:] = Q_interp_wl(wlen,Qsca_C[i])
		Qsca_Si_wl[i,:] = Q_interp_wl(wlen,Qsca_Si[i])

	# bilinear interpolate to obtain the extended table
	for i in range(dim[1]):
		Q = interpolate.interp1d(xtab,Qabs_C_wl[:,i],fill_value='extrapolate')
		Qabs_C_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qabs_Si_wl[:,i],fill_value='extrapolate')
		Qabs_Si_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qsca_C_wl[:,i],fill_value='extrapolate')
		Qsca_C_2D[:,i] = Q(x)
		Q = interpolate.interp1d(xtab,Qsca_Si_wl[:,i],fill_value='extrapolate')
		Qsca_Si_2D[:,i] = Q(x)


	return (cfrac*10**Qabs_C_2D + (1.0-cfrac)*10**Qabs_Si_2D, cfrac*10**Qsca_C_2D + (1.0-cfrac)*10**Qsca_Si_2D)



def extinction_law(x,dsf,wlen,cfrac,t_Qext,t_Qext_V):
	# Input
	#  x: array of log10(grain_radii) of the simulation
	#  dsf: counts (not mass!) of dust grains within bins with centers x
	#  wlen: wavelengths of desired extinction curve
	#  cfrac: mass fraction of graphite
    #  t_Qext: tuple (Q_absorption, Q_scatter)
	# return
	#   A/Av, R (albedo)

	a = 10**x
	Qabs, Qsca = t_Qext[0], t_Qext[1]
	Qabs_V, Qsca_V = t_Qext_V[0].transpose()[0], t_Qext_V[0].transpose()[0]
	Qext = Qabs + Qsca
	Qext_V = Qabs_V + Qsca_V
		
	A = np.zeros(len(wlen))
	Asca = np.zeros(len(wlen))
	Av = 2.5 * np.log10(np.e) * np.pi * np.sum(Qext_V*a**2*dsf)
	for i in range(len(wlen)):
		Asca[i] = 2.5 * np.log10(np.e) * np.pi * np.sum(Qsca[:,i]*a**2*dsf)
		A[i] = 2.5 * np.log10(np.e) * np.pi * np.sum(Qext[:,i]*a**2*dsf)
	R = Asca / A
	A #/= Av
	print(Av)
	return A, R, Qext


#Li & Draine (2001) two-composition opacities. : the carbonaceous component is evaluated at
#every grain size from the same Draine et al. (2021) cross sections the
#PAH emission model.  PAH-like at small sizes, blending smoothly
#into graphite at large sizes, with the Hensley & Draine (2023) ionized
#fraction .  this allows the radiative transfer and the PAH emission to use 
#identical grains at every size. Silicate absorption and all scattering
#come from the Mie tables (Laor & Draine 1993 optical constants, the
#same Draine-family lineage). Outside the wavelength range of the
#carbonaceous cross-section data the mixed Mie values are retained.
CFRAC = 0.54                           # carbonaceous mass fraction (diagnostic co-add only)
N_SPECIES = 2                          # per size bin: graphite, silicate
SPECIES_ORDER = ('graphite', 'silicate')   # flat index = i_size*N_SPECIES + i_species
_CARB_XSEC_RANGE_UM = (0.0912, 1000.)  # validity of the carbonaceous data
_PS = []


def _f_ion(a_um):
    """Hensley & Draine (2023) standard ionized PAH fraction vs size."""
    a_h = 10.e-4   # 10 Angstrom in micron
    return 1. - 1. / (1. + a_um / a_h)


def _carbonaceous_kappa_abs(grain_sizes_this_bin, dN, wlen, rho):
    """dN-weighted carbonaceous absorption opacity per gram for one
    size bin, from the Draine et al. (2021) cross sections with the
    Hensley & Draine (2023) ionized fraction."""
    import pah_spec
    if not _PS:
        _PS.append(pah_spec.PahSpec())
    ps = _PS[0]
    lam = wlen.to(u.micron)
    kabs = np.zeros([len(grain_sizes_this_bin), len(lam)])
    for j, loga in enumerate(grain_sizes_this_bin):
        a = (10.**loga) * u.micron
        ion_cabs, neu_cabs = ps.calc_c_abs(lam, a.to(u.AA))
        fi = _f_ion(a.value)
        cabs = (1. - fi) * neu_cabs[0] + fi * ion_cabs[0]
        m_grain = (4. / 3.) * np.pi * a.to(u.cm)**3 * rho
        kabs[j, :] = (cabs / m_grain).to(u.cm**2 / u.g).value
    return np.sum(kabs * dN[:, np.newaxis], axis=0) / np.sum(dN)


def dust_file_writer(nsizes):

#if __name__ == "__main__":

    #grain size range that we're modeling.  we set it up as a linspace
    #array so that it can have two purposes: (a) to serve as an input
    #array for a test run where we compare against a known size
    #distribution function read in from DNSF_example.txt, and then (b)
    #to set the limits for the grain size bins that we're going to
    #create for powderday
    x = np.linspace(cfg.par.otf_extinction_log_min_size,cfg.par.otf_extinction_log_max_size,41)
    #x=np.linspace(-4,0,41)
        
    #avelength array that we're modeling: 0.1-1000 micron
    wlen = 1. / np.logspace(-4,3,201)*u.micron
    ASSUMED_DENSITY_OF_DUST = 2.4*u.g/u.cm**3

    nu = (constants.c/wlen).to(u.Hz)
    
    #load an example dust size function for testing against
    dsf = np.loadtxt(cfg.par.pd_source_dir+'powderday/active_dust/DNSF_example.txt')

    #assumed graphite fraction
    cfrac = 0.54

    #compute quantities for exmaple dsf above
    xtab, Qtab = Qext_tab_load()
    t_Qext = Qext_get(x,wlen.value,cfrac,xtab,Qtab)
    t_Qext_V = Qext_get(x,np.array([0.551]),cfrac,xtab,Qtab)
    Aext, R, dumQext = extinction_law(x,dsf,wlen.value,cfrac,t_Qext,t_Qext_V)
        
    #nowbreak up the size distribution into bins and then
    #scale wiht the loaded up DSF to see if their co-added
    #extinction laws look reasonable or not.
    
    nbins = nsizes
    
    #array that holds the left edge, right edge, and edges of bins in between.  we use this to set the left and right edge arrays

    #the hydro simulation reports grain counts at nsizes discrete
    #sizes, so each dust file has to describe a bin *centred* on one of
    #those sizes.  taking the edges as the midpoints between neighbouring
    #simulation sizes gives bin i the range [x_i - dx/2, x_i + dx/2],
    #which both centres the within-bin opacity average on the size the
    #bin represents and makes the simulation-size to dust-file mapping
    #one-to-one across the whole grain size range.
    sim_sizes = np.linspace(cfg.par.otf_extinction_log_min_size,
                            cfg.par.otf_extinction_log_max_size,nsizes)
    if nsizes > 1:
        dx = sim_sizes[1]-sim_sizes[0]
    else:
        dx = 1.
    edges = np.linspace(sim_sizes[0]-dx/2.,sim_sizes[-1]+dx/2.,nbins+1)
    grain_size_left_edge_array = edges[0:-1]
    grain_size_right_edge_array = edges[1::]

    #grain_size_left_edge_array = np.linspace(np.min(x),np.max(x),nbins)
    #grain_size_right_edge_array = []
    outfile_filenames = []
    species_of_file = []

    Aext_array = np.zeros([len(wlen.value),nbins])

    #set up an array for storing the weighted fractional contribution
    #each size bin contributes to the final extinction curve. note -
    #this is only done for testing here by comparing these co-added
    #extinction laws back to the original one computed from dsf
    frac = np.zeros(nbins)
    
    #loop through the bins (but not the right most one, hence the -1)
    for counter,i in enumerate(range(nbins)):
        
        grain_sizes_this_bin = np.linspace(grain_size_left_edge_array[i],grain_size_right_edge_array[i],41)#this 41 is an arbitrary choice
        
        #save the right edge of the bin
        #grain_size_right_edge_array.append(grain_size_left_edge_array[i+1])
        

        #used to figure out what bin size of the sample dsf the current size bin we're on corresponds to
        idx = find_nearest(x,grain_size_left_edge_array[i])
        
        #assuming a flat distribution of sizes within each bin
        temp_dsf = np.repeat(1.e59,len(grain_sizes_this_bin))
        
        #set the weighted fraction to scale the binned final dust size distribution by  for testing
        frac[i] = dsf[idx]

        #calculate Qext for the bins
        t_Qext = Qext_get(grain_sizes_this_bin,wlen.value,cfrac,xtab,Qtab)
        t_Qext_V = Qext_get(grain_sizes_this_bin,np.array([0.551]),cfrac,xtab,Qtab)
        
        #compute the stuff!
        temp_Aext,albedo,temp_Qext = extinction_law(grain_sizes_this_bin,temp_dsf,wlen.value,cfrac,t_Qext,t_Qext_V)
        
        #save the extinction laws to an array, just for testing
        Aext_array[:,i] = temp_Aext
        
        
        #calculate kappa
        #kappa (a) = 3 * Qext (a) / ( 4 * a * rho_grain)

        kappa_a_lambda = np.zeros(temp_Qext.shape)
        for i in range(kappa_a_lambda.shape[1]):
            kappa_a_lambda[:,i] = 3.*temp_Qext[:,i]/(4.* ((10.**grain_sizes_this_bin)*u.micron).to(u.cm)*ASSUMED_DENSITY_OF_DUST)
        #kappa_lambda = np.sum(kappa_a_lambda,axis=0)

        #to go from kappa_a_lambda-->kappa_lambda, we need to
        #integrate over a size distribution (as we want the
        #distribution-weighted average over the size range we're
        #modeling here).  this represents a minor inconsistency in
        #these methods as it {\it a priori} assumes a size
        #distribution.  we will assume MRN (dn/da = a^-3.5, so dn = a^-2.5)

        dN = 1400.*(10.**grain_sizes_this_bin)**(-2.5)
        #dN-weighted mean of kappa(a) over the sizes sampled within the
        #bin: sum(kappa * dN) / sum(dN). The weights and the opacities
        #are evaluated at the same size samples, so this is the opacity
        #per gram of an MRN-distributed population within the bin.
        kappa_lambda = np.sum(kappa_a_lambda*dN[:,None],axis=0)/np.sum(dN)

        #per-species opacities.  here we write a pure-graphite and a
        #pure silicate dust file for each size bin.  the per-cell,
        #per-size carbonaceous fraction is then applied downstream through
        #the density grid (tributary_dust_add), so a silicate-rich cell
        #carries no bump and a carbon-rich cell a strong one, taken from
        #the simulation's own grain populations.  carbonaceous absorption
        #uses the Draine (2021)/Hensley-Draine (2023) cross sections where
        #valid and the Mie graphite tables outside that range; silicate
        #and all scattering come from the Mie tables.
        a_cm_bin = ((10.**grain_sizes_this_bin)*u.micron).to(u.cm)
        lamv = wlen.to(u.micron).value
        in_range = ((lamv >= _CARB_XSEC_RANGE_UM[0])
                    & (lamv <= _CARB_XSEC_RANGE_UM[1]))
        t_Q_sil = Qext_get(grain_sizes_this_bin, wlen.value, 0.0, xtab, Qtab)
        t_Q_gra = Qext_get(grain_sizes_this_bin, wlen.value, 1.0, xtab, Qtab)
        kabs_sil_a = np.zeros(t_Q_sil[0].shape)
        ksca_sil_a = np.zeros(t_Q_sil[1].shape)
        kabs_gra_a = np.zeros(t_Q_gra[0].shape)
        ksca_gra_a = np.zeros(t_Q_gra[1].shape)
        for _i in range(kabs_sil_a.shape[1]):
            kabs_sil_a[:,_i] = 3.*t_Q_sil[0][:,_i]/(4.*a_cm_bin*ASSUMED_DENSITY_OF_DUST)
            ksca_sil_a[:,_i] = 3.*t_Q_sil[1][:,_i]/(4.*a_cm_bin*ASSUMED_DENSITY_OF_DUST)
            kabs_gra_a[:,_i] = 3.*t_Q_gra[0][:,_i]/(4.*a_cm_bin*ASSUMED_DENSITY_OF_DUST)
            ksca_gra_a[:,_i] = 3.*t_Q_gra[1][:,_i]/(4.*a_cm_bin*ASSUMED_DENSITY_OF_DUST)
        kabs_sil = np.sum(kabs_sil_a*dN[:,None],axis=0)/np.sum(dN)
        ksca_sil = np.sum(ksca_sil_a*dN[:,None],axis=0)/np.sum(dN)
        kabs_gra_mie = np.sum(kabs_gra_a*dN[:,None],axis=0)/np.sum(dN)
        ksca_gra = np.sum(ksca_gra_a*dN[:,None],axis=0)/np.sum(dN)
        kabs_car = _carbonaceous_kappa_abs(grain_sizes_this_bin, dN, wlen,
                                           ASSUMED_DENSITY_OF_DUST)

        #assemble each pure species: chi = kappa_abs + kappa_sca, albedo
        #= kappa_sca/chi.
        kabs_gra = np.where(in_range, kabs_car, kabs_gra_mie)
        chi_gra = kabs_gra + ksca_gra
        albedo_gra = np.where(chi_gra > 0, ksca_gra/chi_gra, 0.)
        chi_sil = kabs_sil + ksca_sil
        albedo_sil = np.where(chi_sil > 0, ksca_sil/chi_sil, 0.)

        #(Aext_array for the diagnostic co-add plot is already set above,
        #before the inner loop clobbers the loop variable.)

        #----------------------------------
        #write the two per-species HDF5 files for this size bin.
        #flat index = i_size*N_SPECIES + i_species, species order
        #(graphite, silicate) -- a C-order flatten of (n_sizes, n_species),
        #matched exactly by frac_grid in tributary_dust_add and by the
        #kappa read-back in isrf_decompose.
        #----------------------------------
        if not os.path.exists(cfg.model.PD_output_dir+'/dust_files/'):
            os.makedirs(cfg.model.PD_output_dir+'/dust_files/')
        for i_species, (chi_sp, alb_sp) in enumerate(
                ((chi_gra, albedo_gra), (chi_sil, albedo_sil))):
            flat = counter*N_SPECIES + i_species
            d = IsotropicDust(nu.value, alb_sp, chi_sp)
            filename = (cfg.model.PD_output_dir
                        + '/dust_files/binned_dust_sizes.'+str(flat)+'.hdf5')
            outfile_filenames.append(filename)
            species_of_file.append(SPECIES_ORDER[i_species])
            d.write(filename)



    #----------------------------------
    #save the metadata
    #----------------------------------
    grain_size_right_edge_array = np.asarray(grain_size_right_edge_array)

    x = grain_size_left_edge_array[0:-1]
    y = grain_size_right_edge_array
    z = np.asarray(outfile_filenames)
    #np.savetxt('dust_files/binned_dust_sizes.key',np.transpose([grain_size_left_edge_array[0:-1],grain_size_right_edge_array,np.asarray(outfile_filenames)]))

    #outfile_filenames is now n_sizes*N_SPECIES long, flattened size-major
    #(species order graphite, silicate).  n_sizes / n_species / species are
    #stored so every reader reshapes the flat dust-type axis from one
    #source of truth rather than assuming a count.
    np.savez(cfg.model.PD_output_dir+'/dust_files/binned_dust_sizes.npz',
             grain_size_left_edge_array = grain_size_left_edge_array,
             grain_size_right_edge_array = grain_size_right_edge_array,
             outfile_filenames = outfile_filenames,
             species_of_file = np.asarray(species_of_file),
             n_sizes = nbins, n_species = N_SPECIES,
             species_order = np.asarray(SPECIES_ORDER))



    #----------------------------------
    #Making some plots just for testing
    #----------------------------------
    final_Aext = np.average(Aext_array,axis=1,weights=frac)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(wlen,final_Aext,label='coadded')
    ax.loglog(wlen,Aext,label='original')

    plt.legend(loc=2)

    fig.savefig('extinction.png',dpi=300)
        
