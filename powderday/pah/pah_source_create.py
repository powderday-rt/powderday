#DEBUG NEED TO NOT JUST USE 0TH 3 INDICES AS BELOW, BUT NEED TO MAKE THIS DYNAMIC BOTH HERE AND IN PAH_SPEC



import numpy as np
from powderday.pah.pah_file_read import read_draine_file
from powderday.helpers import find_nearest
from astropy import units as u
from astropy import constants as constants
import powderday.config as cfg
import pdb
from tqdm import tqdm
from powderday.pah.isrf_decompose import get_beta_nnls,get_isrf
import os,glob
import multiprocessing as mp
from functools import partial
from datetime import datetime
from unyt import unyt_quantity,unyt_array
import pah_spec  # import Helena Richie's pah_spec model for SPA calculations
from powderday.pah.isrf_decompose import get_Cabs,get_logU
from scipy.interpolate import interp1d
import tqdm
#This global variable will hold the model inside each worker process
_worker_ps_model = None



import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as constants
# Make sure pah_spec is imported as ps_module or similar if 'pah_spec' variable name is used for the module
import pah_spec 


import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as constants
import pah_spec

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as constants
# Ensure pah_spec is imported
import pah_spec
    

def debug_compare_inputs(ps, sim_lam, sim_u_lambda, sim_gsd_neu, sim_gsd_ion, sim_sizes):
    """
    Plots and prints a comparison between the Simulation inputs for a single cell
    and the pah_spec reference defaults.
    """
    # ==========================================
    # 1. SETUP REFERENCE DATA (PAH_SPEC DEFAULTS)
    # ==========================================
    
    # Reference ISRF: strictly 10 * ps.u_lambda_arr as requested
    ref_u_lambda = 10 * ps.u_lambda_arr
    
    ref_lam = ps.wavelength_u_arr
    if not hasattr(ref_lam, 'unit'): 
        ref_lam *= u.micron
    
    # Reference GSD
    ref_gsd_neu = ps.size_dist_neu
    ref_gsd_ion = ps.size_dist_ion
    
    # Reference Sizes: Use module constant explicitly to avoid AttributeError
    ref_sizes = pah_spec.GRAIN_SIZES
    if not hasattr(ref_sizes, 'unit'): 
        ref_sizes *= u.angstrom # pah_spec defaults are typically Angstroms here
    else:
        ref_sizes = ref_sizes.to(u.angstrom)

    # ==========================================
    # 2. CALCULATE DIAGNOSTICS
    # ==========================================
    
    # Helper to integrate energy density
    def get_integrated_energy(lam, u_lam):
        # Sort by wavelength for integration
        sort_idx = np.argsort(lam)
        # Trapz integration
        return np.trapz(u_lam[sort_idx], lam[sort_idx])

    # Ensure units for integration
    sim_lam_cm = sim_lam.to(u.cm)
    sim_u_lam_cgs = sim_u_lambda.to(u.erg / u.cm**4)
    
    ref_lam_cm = ref_lam.to(u.cm)
    ref_u_lam_cgs = ref_u_lambda.to(u.erg / u.cm**4)

    E_dens_sim = get_integrated_energy(sim_lam_cm, sim_u_lam_cgs)
    E_dens_ref = get_integrated_energy(ref_lam_cm, ref_u_lam_cgs)
    
    ratio_U = (E_dens_sim / E_dens_ref).decompose().value

    # ==========================================
    # 3. PLOTTING
    # ==========================================
    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    
    # User Preference: Font sizes
    plt.rcParams.update({'axes.labelsize': 16, 'xtick.labelsize': 16, 
                         'ytick.labelsize': 16, 'legend.fontsize': 14})

    # --- Plot 1: Radiation Fields (ISRF) ---
    # Plot Reference
    ax[0].loglog(ref_lam.to(u.micron), ref_u_lambda.to(u.erg/u.cm**4), 
                 label=r'Reference ($10 \times$ Default)', color='black', linestyle='--')
    
    # Plot Simulation
    # Ensure sim_lam is sorted for plotting if it isn't already
    sort_idx_sim = np.argsort(sim_lam)
    ax[0].loglog(sim_lam[sort_idx_sim].to(u.micron), sim_u_lambda[sort_idx_sim].to(u.erg/u.cm**4), 
                 label='Simulation Cell (ergcm4)', color='red')
    
    ax[0].set_xlabel(r'Wavelength [$\mu$m]')
    ax[0].set_ylabel(r'$u_\lambda$ [erg cm$^{-4}$]')
    ax[0].legend()
    # User Preference: No title

    # --- Plot 2: Grain Size Distributions (Neutrals) ---
    # Plot Reference
    # Note: pah_spec defaults are often per H. 
    ax[1].loglog(ref_sizes.to(u.angstrom), ref_gsd_neu, 
                 label='Reference GSD (Neutral)', color='black', linestyle='--')
    
    # Plot Simulation
    sim_sizes_ang = sim_sizes.to(u.angstrom)
    
    # Check if shapes match for plotting. If sim_gsd is a single value per bin (expected), plot directly.
    ax[1].loglog(sim_sizes_ang, sim_gsd_neu, 
                 label='Simulation GSD (Neutral)', color='red', marker='o')

    ax[1].set_xlabel(r'Grain Size $a$ [$\AA$]')
    ax[1].set_ylabel(r'Size Dist ($dn/da/n_H$)')
    ax[1].legend()
    # User Preference: No title

    plt.tight_layout()
    # User Preference: dpi=300
    plt.savefig('debug_pah_comparison.png', dpi=300)
    # plt.show() # Uncomment if you have a display

    # ==========================================
    # 4. PRINT DIAGNOSTICS
    # ==========================================
    print("\n" + "="*40)
    print(f"DEBUG DIAGNOSTICS: Sim vs Reference")
    print("="*40)
    print(f"Total Integrated Energy Density (approx U):")
    print(f"  Reference (10*default): {E_dens_ref:.3e}")
    print(f"  Simulation Cell:        {E_dens_sim:.3e}")
    print(f"  Ratio (Sim / Ref):      {ratio_U:.5f}")
    
    if ratio_U < 1e-4:
         print("  --> [CRITICAL] Simulation ISRF is extremely faint!")
    
    print("-" * 20)
    print(f"Grain Size Distribution (Neutral) Max Value:")
    print(f"  Reference Max: {np.max(ref_gsd_neu):.3e}")
    print(f"  Simulation Max: {np.max(sim_gsd_neu):.3e}")
    
    if np.max(sim_gsd_neu) > 1e10:
        print("  --> [CRITICAL] Sim GSD is huge (~1e50?). This looks like Total Number.")
        print("                 pah_spec expects abundance per H (~1e-7).")
        print("                 Check division by (n_H * Volume).")
    elif np.max(sim_gsd_neu) < 1e-20:
        print("  --> [CRITICAL] Sim GSD is extremely small. Check units (cm^-3 vs code units).")
    
    print("="*40 + "\n")
    


def get_whole_ceil(n,near):
    nn = np.divide(n,np.linspace(1,np.ceil(n/near),int(np.ceil(n/near))))
    return(nn[nn%1==0][-1])




#SERIES OF FUNCTIONS FOR PARALLEL PROCESSING THE SPA CALCULATIONS

#Fuunction designed to run once per CPU core when the pool starts to call pah_spec
import numpy as np
import multiprocessing as mp
from astropy import constants
from astropy import units as u
import pah_spec

# Global variable for the worker process
_worker_ps_model = None

def _init_worker():
    """
    Initialize the pah_spec model once per core to avoid reloading 
    it for every single cell.
    """
    global _worker_ps_model
    import pah_spec
    _worker_ps_model = pah_spec.PahSpec()


# ==========================================
# PARALLEL HELPER FUNCTIONS (Part 2/3)
# ==========================================

def _process_cell_task(args):
    """
    Worker process that replicates the exact physics and unit logic 
    of the serial implementation.
    """
    # Unpack arguments
    # Note: These inputs are expected to be Astropy Quantities (with units)
    (wav_input,          # [micron] Wavelengths (reversed)
     u_lambda_input,     # [erg/cm^4] Specific Energy Density (reversed)
     gsd_slice,          # [dimensionless] Grain counts for the PAH bins
     n_H,                # [cm^-3] Gas density
     cell_vol,           # [cm^3] Cell volume
     sim_sizes,          # [cm] Simulation grain sizes (sliced to PAH bins)
     f_ion               # [dimensionless] Ionization fraction (sliced)
     ) = args     
    

    norm_factor = (n_H.value * cell_vol.value)
    if norm_factor <= 0:
        # Return zero arrays with correct units if cell is empty
        zero_arr = np.zeros(len(wav_input)) * u.erg / (u.cm * u.s)
        return zero_arr, zero_arr
    
    # Calculate distributions
    #cell_size_dist_neu = gsd_slice.value / norm_factor / sim_sizes.value * (1. - f_ion.value)
    #cell_size_dist_ion = gsd_slice.value / norm_factor / sim_sizes.value * f_ion.value

    cell_size_dist_neu = (gsd_slice.value / norm_factor) * (1. - f_ion.value)
    cell_size_dist_ion = (gsd_slice.value / norm_factor) * f_ion.value
    
    # ------------------------------------------------------------------
    # GENERATE SPECTRUM
    # ------------------------------------------------------------------
    # Inputs:
    #   wavelength_arr: [micron] Quantity
    #   u_lambda_arr:   [erg/cm^4] Quantity
    #   size_dist_neu:  [float array] Abundance (per H)
    #   size_dist_ion:  [float array] Abundance (per H)
    
    spec_neu, spec_ion = _worker_ps_model.generate_spectrum(
        wavelength_arr=wav_input, 
        u_lambda_arr=u_lambda_input, 
        size_dist_neu=cell_size_dist_neu, 
        size_dist_ion=cell_size_dist_ion
    )
    
    # ------------------------------------------------------------------
    # SCALE BACK TO LUMINOSITY
    # ------------------------------------------------------------------
    # Serial Code Reference:
    # neutral_grid_PAH_luminosity[...] = spectrum_neu * (n_H.value * vol.value)
    
    # Result is Total Luminosity per wavelength [erg / (cm * s)]
    lum_neu = spec_neu * norm_factor
    lum_ion = spec_ion * norm_factor
    
    return lum_neu, lum_ion


# ==========================================
# PARALLEL DRIVER FUNCTION (Part 3/3)
# ==========================================

def compute_grid_PAH_luminosity_SPA_parallel(cell_list, gsd, reg, simulation_sizes, ds, draine_directories, f_ion):

    #The number of PAH sizes considered by pah_spec
    n_pah_sizes = len(pah_spec.GRAIN_SIZES)

    # ---------------------------------------------------------
    # 1. PREPARE ISRF
    # ---------------------------------------------------------
    simulation_specific_energy_gsd_convolved, simulation_isrf_nu, simulation_isrf_lam = get_isrf(gsd, reg)
    cell_isrf = simulation_specific_energy_gsd_convolved.cgs.value.T * u.erg / u.Hz

    #convert yt-->astropy units
    cell_sizes = reg.parameters['cell_size'].in_units('cm').value * u.cm

    #Convert E_nu [erg/Hz] -> Energy Density u_nu [erg/cm^3/Hz]
    cell_vol = cell_sizes**3
    u_nu = cell_isrf.T / cell_vol
    
    # Convert u_nu [per Hz] -> u_lambda [per cm]
    lam = simulation_isrf_lam
    jacobian = constants.c / (lam**2)
    cell_isrf_ergcm4 = u_nu.T * jacobian
    # Verify strict unit compliance 
    cell_isrf_ergcm4 = cell_isrf_ergcm4.to(u.erg / u.cm**4)

    # ---------------------------------------------------------
    # 2. PREPARE DENSITY & CONSTANTS
    # ---------------------------------------------------------
    n_H = ds.arr(reg['PartType3', 'Dust_GasDensity'], 'code_mass/code_length**3').in_units('g/cm**3').value * u.g / u.cm**3
    n_H /= constants.m_p.cgs
    n_H = n_H.to(u.cm**-3)

    n_cells = len(cell_list)
    
    # Instantiate temp model just to get wavelength array dimensions
    temp_ps = pah_spec.PahSpec() 
    n_out_wav = len(temp_ps.emission_wavelengths)
    
    # Flip wavelengths (ascending order for pah_spec)
    wav_input = simulation_isrf_lam.to(u.micron)[::-1]
    
    # Slice simulation inputs constant across cells
    sim_sizes_sliced = simulation_sizes[0:n_pah_sizes]
    f_ion_sliced = f_ion[0:n_pah_sizes]

    # ---------------------------------------------------------
    # 3. BUILD TASK LIST
    # ---------------------------------------------------------
    print(f"[SPA Parallel] Preparing tasks for {n_cells} cells...")
    tasks = []
    
    for i in range(n_cells):
        
        # Specific inputs for this cell
        # ISRF (flipped to ascending wavelength)
        u_lambda_this = cell_isrf_ergcm4[i, :][::-1]
        
        # GSD counts for the PAH bins
        gsd_slice = gsd[i, :][0:n_pah_sizes]

        # Geometry
        n_H_this = n_H[i]
        vol_this = cell_sizes[i]**3 # Quantity [cm^3]

        task_tuple = (
            wav_input,          # Quantity [micron]
            u_lambda_this,      # Quantity [erg/cm^4]
            gsd_slice,          # Array (counts)
            n_H_this,           # Quantity [cm^-3]
            vol_this,           # Quantity [cm^3]
            sim_sizes_sliced,   # Quantity [cm]
            f_ion_sliced        # Array (fraction)
        )
        tasks.append(task_tuple)

    # ---------------------------------------------------------
    # 4. RUN PARALLEL POOL
    # ---------------------------------------------------------
    n_procs = cfg.par.n_processes
    print(f"[SPA Parallel] Launching pool with {n_procs} workers...")

    with mp.Pool(processes=n_procs, initializer=_init_worker) as pool:
        results = list(tqdm.tqdm(pool.imap(_process_cell_task, tasks), total=n_cells))

    # ---------------------------------------------------------
    # 5. UNPACK RESULTS
    # ---------------------------------------------------------
    neutral_grid_PAH_luminosity = np.zeros((n_cells, n_out_wav)) * u.erg / (u.cm * u.s)
    ion_grid_PAH_luminosity = np.zeros((n_cells, n_out_wav)) * u.erg / (u.cm * u.s)

    for i, (res_neu, res_ion) in enumerate(results):
        neutral_grid_PAH_luminosity[i, :] = res_neu
        ion_grid_PAH_luminosity[i, :] = res_ion

    grid_PAH_luminosity = neutral_grid_PAH_luminosity + ion_grid_PAH_luminosity

    return grid_PAH_luminosity, neutral_grid_PAH_luminosity, ion_grid_PAH_luminosity



def compute_grid_PAH_luminosity_SPA_serial(cell_list, gsd, reg, simulation_sizes, ds, draine_directories, f_ion):
    
    #IN PRACTICE THIS SERIAL CODE ISN'T INTENDED TO BE USED IN
    #PRODUCTION. IT'S JUST HERE BECAUSE IT'S 1000X EASIER TO DEVELOP
    #AND DEBUG SERIAL CODE, AND THEN, ONCE VERIFIED, CONVERT TO
    #PARALLEL.


    n_pah_sizes = len(pah_spec.GRAIN_SIZES)

    
    # Get the ISRF for all cells
    simulation_specific_energy_gsd_convolved, simulation_isrf_nu, simulation_isrf_lam = get_isrf(gsd, reg)

    cell_isrf = simulation_specific_energy_gsd_convolved.cgs.T
    cell_sizes = reg.parameters['cell_size'].value * u.cm

    
    # 1. Convert E_nu [erg/Hz] -> Energy Density u_nu [erg/cm^3/Hz]
    cell_vol = cell_sizes**3
    cell_isrf = cell_isrf.value*u.erg/u.Hz
    u_nu = cell_isrf.T / cell_vol
    
    u_nu = u_nu.to(u.erg / u.cm**3 / u.Hz)
    
    # 2. Convert u_nu -> u_lambda 
    lam_cm = simulation_isrf_lam.to(u.cm)
    jacobian = constants.c / (lam_cm**2)
    cell_isrf_ergcm4 = u_nu.T * jacobian

    cell_isrf_ergcm4 = cell_isrf_ergcm4.to(u.erg / u.cm**4)


    # Initialize pah_spec
    ps = pah_spec.PahSpec()
    
    n_cells = len(cell_list)
    n_output_wav = len(ps.emission_wavelengths)
    
    # Pre-allocate output grids
    neutral_grid_PAH_luminosity = np.zeros((n_cells, n_output_wav)) * u.erg / (u.cm * u.s)
    ion_grid_PAH_luminosity = np.zeros((n_cells, n_output_wav)) * u.erg / (u.cm * u.s)

    # Calculate n_H
    n_H = ds.arr(reg['PartType3', 'Dust_GasDensity'], 'code_mass/code_length**3').in_units('g/cm**3').value * u.g / u.cm**3
    n_H /= constants.m_p.cgs
    n_H = n_H.to(u.cm**-3)


    for counter, cell in enumerate(cell_list):
        print(counter)

        # Skip empty cells to prevent divide-by-zero instability
        if n_H[counter].value < 1e-5: 
            continue

        # Get the ISRF for this cell (ascending wavelength order for pah_spec)
        isrf_lam = simulation_isrf_lam.to(u.micron)[::-1]
        cell_isrf_this = cell_isrf_ergcm4[counter, :][::-1]

        
        # DIAGNOSTIC 1: Check inputs for instability
        print(f"  [Cell {counter}] n_H: {n_H[counter]:.3e}")
        print(f"  [Cell {counter}] ISRF max: {np.max(cell_isrf_this.value):.3e}")

        
        # Original normalization approach - divide by (n_H * volume * size)
        # This converts grain counts to something like dn/da per H per volume
        #cell_size_dist_neu = gsd[counter, :][0:n_pah_sizes].value / (n_H[counter].value * cell_sizes[counter].value**3) / simulation_sizes[0:n_pah_sizes].value * (1. - f_ion[0:n_pah_sizes].value)
        #cell_size_dist_ion = gsd[counter, :][0:n_pah_sizes].value / (n_H[counter].value * cell_sizes[counter].value**3) / simulation_sizes[0:n_pah_sizes].value * f_ion[0:n_pah_sizes].value

        cell_size_dist_neu = gsd[counter, :][0:n_pah_sizes].value / (n_H[counter].value * cell_sizes[counter].value**3)  * (1. - f_ion[0:n_pah_sizes].value)
        cell_size_dist_ion = gsd[counter, :][0:n_pah_sizes].value / (n_H[counter].value * cell_sizes[counter].value**3)  * f_ion[0:n_pah_sizes].value


        # DIAGNOSTIC 3: Check normalization passed to spectrum generator
        print(f"  [Cell {counter}] Max Grain/H input: {np.max(cell_size_dist_neu):.3e}")
        
        # DIAGNOSTIC 4: Check Physical Cell Size
        # (Add this near the other print statements)
        vol_cm3 = cell_sizes[counter].value**3
        size_pc = cell_sizes[counter].to(u.pc).value
        
        print(f"--- CELL {counter} GEOMETRY CHECK ---")
        print(f"  Size:   {size_pc:.4f} pc")
        print(f"  Volume: {vol_cm3:.3e} cm^3")
        print(f"  n_H:    {n_H[counter].value:.3e} cm^-3")
        print(f"  Total H atoms (n_H * Vol): {(n_H[counter].value * vol_cm3):.3e}")

        
        # Generate the spectrum using pah_spec
        spectrum_neu, spectrum_ion = ps.generate_spectrum(
            wavelength_arr=isrf_lam,
            u_lambda_arr=cell_isrf_this,
            size_dist_neu=cell_size_dist_neu,
            size_dist_ion=cell_size_dist_ion
        )
        

        # Scale back by (n_H * volume) to get total luminosity per wavelength
        neutral_grid_PAH_luminosity[counter, :] = spectrum_neu * (n_H[counter].value * cell_sizes[counter].value**3)
        ion_grid_PAH_luminosity[counter, :] = spectrum_ion * (n_H[counter].value * cell_sizes[counter].value**3)


        
        # --- DEBUG CALL HERE (Only for the first cell) ---
        debug_compare_inputs(
            ps,
            isrf_lam,
            cell_isrf_this,
            cell_size_dist_neu,
	    cell_size_dist_ion,
            simulation_sizes[0:n_pah_sizes]
        )


        # DIAGNOSTIC 2: Compare raw sum vs integrated sum
        
        # Current method (Sum of L_lambda)
        raw_sum = np.sum(neutral_grid_PAH_luminosity[counter, :].value + 
                         ion_grid_PAH_luminosity[counter, :].value)
        
        # Integration approximation (L_lambda * lambda) -> erg/s
        # (This is rough integration, but gets the order of magnitude right)
        wavelengths_cm = ps.emission_wavelengths.to(u.cm).value
        integrated_sum = np.sum(
            (neutral_grid_PAH_luminosity[counter, :].value + 
             ion_grid_PAH_luminosity[counter, :].value) * wavelengths_cm
        )
        
        print(f"  [Cell {counter}] Raw Sum (erg/s/cm): {raw_sum:.3e}")
        print(f"  [Cell {counter}] Integrated (erg/s):   {integrated_sum:.3e}")

                
	#ENDDEBUG

        
    grid_PAH_luminosity = neutral_grid_PAH_luminosity + ion_grid_PAH_luminosity

    #returns are in erg/cm/s
    return grid_PAH_luminosity, neutral_grid_PAH_luminosity, ion_grid_PAH_luminosity



def compute_grid_PAH_luminosity(cell_list,beta_nnls,grid_of_sizes,numgrains,draine_sizes,draine_lam,f_ion,neutral_PAH_reference_objects,ion_PAH_reference_objects,
                                logU,basis_logU_values, draine_bins_idx):

    #these are re-defined for each pool thread.  when they're
    #returned, they'll be packed into a master list with one extra
    #dimension for the thread number.
    neutral_grid_PAH_luminosity = np.zeros([len(cell_list),len(draine_lam)])
    ion_grid_PAH_luminosity = np.zeros([len(cell_list),len(draine_lam)])


    for counter,cell in enumerate(cell_list):
        print(cell)
        beta_cell = beta_nnls[:,cell]
        beta_cell = beta_cell/np.max(beta_cell)
        
        #need to make a temporary (for this cell) PAH_list that is
        #just n_draine_sizes long that is convolved with beta_nnls
        neutral_pah_grid = np.zeros([len(draine_sizes),len(draine_lam)])
        ion_pah_grid = np.zeros([len(draine_sizes),len(draine_lam)])

        logU_cell = logU[cell]
                
        for j in np.flatnonzero(beta_cell): #for j in range(len(Beta_cell))
        #for j in range(len(beta_cell)):

            #find the logU that is closest to what draine has computed in their basis computations
            nearest_logU_idx = find_nearest(logU_cell,np.asarray(basis_logU_values))
            
            #here, identify the correct reference object that
            #corresponds to both the beta_cell (j), but also which
            #logU value it is.

            neutral_PAH_list = neutral_PAH_reference_objects[j,nearest_logU_idx]
            temp_neutral_pah_grid = np.array([x.lum for x in neutral_PAH_list])
            temp_neutral_pah_grid *= beta_cell[j]
            neutral_pah_grid += temp_neutral_pah_grid #this is the running summation of the (n_sizes,n_lam) pah grid for the i-th cell
            
            ion_PAH_list = ion_PAH_reference_objects[j,nearest_logU_idx]
            temp_ion_pah_grid = np.array([x.lum for x in ion_PAH_list])
            temp_ion_pah_grid *= beta_cell[j]
            ion_pah_grid += temp_ion_pah_grid #this is the running summation of the (n_sizes,n_lam) pah grid for the i-th cell
            
        #set the PAH luminosity of the cell to be the dot product of
        #the Draine luminosities (i.e., pah_grid[draine_bins_idx,:] which has
        #dimensions (simulation_sizes,wavelengths)) with the actual
        #grain size distribution in that cell (i.e.,
        #grid_of_sizes[i_cell,:]). note, we take the transpose of
        #grid_of_sizes to get the dimensions to match up correctly for the dot product

        #note - we're also folding in the ionized fraction 
        neutral_grid_PAH_luminosity[counter,:] = np.dot(neutral_pah_grid[draine_bins_idx,:].T*(1.-f_ion), grid_of_sizes.T[:,cell])
        ion_grid_PAH_luminosity[counter,:] = np.dot(ion_pah_grid[draine_bins_idx,:].T*f_ion, grid_of_sizes.T[:,cell])
        
        
    #particle_PAH_luminosity = np.dot(pah_grid[draine_bins_idx,:].T,numgrains.T[:,cell])


    grid_PAH_luminosity = neutral_grid_PAH_luminosity + ion_grid_PAH_luminosity
    

    return grid_PAH_luminosity,neutral_grid_PAH_luminosity,ion_grid_PAH_luminosity


def get_PAH_lum_cdf(nu_reverse,fnu,wpah_nu_reverse,grid_PAH_luminosity):
    lum_list = []

    for i in range(grid_PAH_luminosity.shape[0]): 
        fnu_reverse = fnu[i,:][::-1]
        lum = (np.absolute(np.trapz(nu_reverse[wpah_nu_reverse].cgs.value,fnu_reverse[wpah_nu_reverse])).item()*u.Lsun).to(u.erg/u.s).value
        lum_list.append(lum)

    lum_list = np.asarray(lum_list)
    lum_list[lum_list == 0] = np.min(lum_list[lum_list > 0])

    #omg this is so janky. but it works.
    percentile_list = []

    loglum_bins = np.linspace(np.min(np.log10(lum_list)),np.max(np.log10(lum_list)),100)
    for loglum in loglum_bins:
        percentile_list.append(np.sum(lum_list[lum_list > 10.**loglum])/np.sum(lum_list))
        
    percentile_idx = find_nearest(np.asarray(percentile_list),cfg.par.percentile_LPAH_to_include)
    lum_to_cut_below = 10.**(loglum_bins[percentile_idx])

    useful_idxs = np.where(lum_list >= lum_to_cut_below)[0]


    return useful_idxs


def pah_source_add(ds,reg,m,boost):
    
    LUM_FLOOR = 1.e20 #erg/s -- just some small value compared to the ~few Lsun we typically get in a cell

    
    #first - establish where we're working
    draine_directories = []
    print('powderday/pah/pah_source_create]: reading from the following Draine PAH directories')
    for it in os.scandir(cfg.par.draine_data_dir):
        if it.is_dir():
            print(it.path)
            draine_directories.append(it.path)


    #first establish the grain size distribution and sizes from the
    #hydro simulation
    grid_of_sizes = ds.parameters['reg_grid_of_sizes']

    #these are the size bins from the hydro sim
    simulation_sizes = (ds.parameters['grain_sizes_in_micron']*u.micron)


    
    #second, use the Hensley & Draine fitting formula to determine
    #f_ion as a function of size: #f_ion(a) = 1 - 1/(1 + a/10 A)
    

    f_ion = 1.- 1./(1+(simulation_sizes.to(u.angstrom)/(10.*u.angstrom)))
    

    #determine q_PAH for analysis and save it to parameters for
    #writing out 
    ad = ds.all_data()


    #compute the PAH mass
    #first set the dust density as 2.4 g/cm**3 (the assumed density) 
    dust_density = np.ones(grid_of_sizes.shape)*ds.quan(cfg.par.dust_density,'g/cm**3')
    mass_per_bin = dust_density* np.pi * 3./4 * unyt_array.from_astropy(simulation_sizes.to(u.cm))**3 * (ds.parameters['reg_grid_of_sizes_graphite']*ds.parameters['reg_grid_of_sizes_aromatic_fraction'])
    #ensure mass_per_bin is in g since we'll lose the unit later
    mass_per_bin = mass_per_bin.in_units('g')
    
    #pah_idx = np.where(simulation_sizes.to(u.angstrom).value <= 13) #this is the sizes q_pah usually measures; 1000 carbon atoms and eq. 22 of Hensley and Draine
    #mass_per_bin_only_pah_idx = np.zeros(mass_per_bin.shape)
    #mass_per_bin_only_pah_idx[:,pah_idx]=mass_per_bin[:,pah_idx]
    m_pah = ds.arr(np.sum(mass_per_bin,axis=1),'g')
    reg.parameters['m_pah'] = m_pah
    q_pah = m_pah.in_units('g')/reg['dust','mass'].in_units('g')
    #in case we have an errant cell with no dust information (super rare corner case)
    q_pah[np.isnan(q_pah)] = np.min(q_pah[~np.isnan(q_pah)])
    reg.parameters['q_pah'] = q_pah


    
    #dN_pah = np.sum(ds.parameters['reg_grid_of_sizes_graphite']*ds.parameters['reg_grid_of_sizes_aromatic_fraction'],axis=1)
    dN_pah = np.sum(ds.parameters['reg_grid_of_sizes_graphite'],axis=1)
    dN_total = np.sum(ds.parameters['reg_grid_of_sizes'],axis=1)

    #compute the mass weighted grain size distributions for comparison in analytics.py 
    #try: #for mesh based code
    #particle_mass_weighted_gsd = np.average(reg['dust','numgrains'],weights=reg['dust','mass'],axis=0)
    #grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','mass'],axis=0)
    #except:
    particle_mass_weighted_gsd = np.average(reg['particle_dust','numgrains'],weights=reg['dust','mass'],axis=0)
    try: #for octree    
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','smoothedmasses'],axis=0)
    except:
        grid_mass_weighted_gsd = np.average(grid_of_sizes,weights=reg['dust','mass'],axis=0)

    #second, read the information from the Draine files. We can do
    #this just for an arbitrary file in one of the draine_directories
    #since the size bins and wavelengths are all the same.
    temp_filename = glob.glob(draine_directories[0]+'/*iout_graD16*nb*_0.00')[0]
    temp_PAH_list = read_draine_file(temp_filename)
    draine_sizes = temp_PAH_list[0].size_list
    draine_lam = temp_PAH_list[0].lam*u.micron


    #get the logU and beta_nnls for the local ISRF
    beta_nnls,logU = get_beta_nnls(draine_directories,grid_of_sizes,simulation_sizes,reg)
    #just sto save it through analytics
    reg.parameters['logU'] = logU
    
    #read in a single draine directory, get the files, and from this
    #grab the logU values that the draine heating rates have been
    #computed for. note that this assumes that every single basis ISRF
    #sed shape has been computed for the exact same logU files.  pd
    #will throw some error around here if for some reason this isn't
    #true and it tries to read in a file that doesn't exist.
    temp_neutral_draine_files = np.sort(glob.glob(draine_directories[0]+'/*iout_graD16*nb*_*.*0'))

    basis_logU_values = []
    for file in temp_neutral_draine_files:
        basis_logU_values.append(float(file[-4::]))
    

    
    for draine_directory in draine_directories:
        neutral_logU_iout_files = np.sort(glob.glob(draine_directory+'/*iout_graD16*nb*_*.*0'))
        ion_logU_iout_files = np.sort(glob.glob(draine_directory+'/*iout_graD16*ib*_*.*0'))
        
    try:
        print("attempting to read in Draine reference objects from file: "+cfg.par.draine_data_dir+'/draine_reference_objects.npz')
        data = np.load(cfg.par.draine_data_dir+'/draine_reference_objects.npz',allow_pickle=True)
        neutral_PAH_reference_objects = data['neutral_PAH_reference_objects']
        ion_PAH_reference_objects = data['ion_PAH_reference_objects']


    except:
        print("[pah/pah_source_create:] Draine Reference Objects file does not exist: Creating now [may take 5 min]")

        
        #we will build the PAH reference objects as [n_draine_directories
        #[i.e., SED shapes], nlogU, len(temp_PAH_list)].  we can then
        #reference the correct SED shape and logU idx when actually
        #computing the SEd.  i'm sure there's some cool
        #pandas/pythonic/object oriented way to do this but whatever. this
        #country was built on endlessly cumbersome multidimensional
        #arrays.
              
        neutral_PAH_reference_objects = np.zeros([len(draine_directories),len(basis_logU_values),len(temp_PAH_list)],dtype=object)
        ion_PAH_reference_objects = np.zeros([len(draine_directories),len(basis_logU_values),len(temp_PAH_list)],dtype=object)
        
        print("[pah/pah_source_create:] building the reference PAH list for neutrals and ions")

        for draine_directory_idx,draine_directory in enumerate(draine_directories):

            neutral_logU_iout_files = np.sort(glob.glob(draine_directory+'/*iout_graD16*nb*_*.*0'))
            ion_logU_iout_files = np.sort(glob.glob(draine_directory+'/*iout_graD16*ib*_*.*0'))
        
            for logU_idx in range(len(basis_logU_values)):
                neutral_file = neutral_logU_iout_files[logU_idx]
                ion_file = ion_logU_iout_files[logU_idx]
                
                print('[pah/pah_source_create: ] processing PAH file: '+neutral_file)
                neutral_PAH_reference_objects[draine_directory_idx,logU_idx,:] = np.asarray(read_draine_file(neutral_file))
                
                print('[pah/pah_source_create: ] processing PAH file: '+ion_file)
                ion_PAH_reference_objects[draine_directory_idx,logU_idx,:] = np.asarray(read_draine_file(ion_file))

        np.savez(cfg.par.draine_data_dir+'/draine_reference_objects.npz',neutral_PAH_reference_objects=neutral_PAH_reference_objects,ion_PAH_reference_objects=ion_PAH_reference_objects)
    

        
    
    #print("[pah/pah_source_create:] building the reference PAH list for neutrals")
    #for counter,neutral_file in tqdm(enumerate(neutral_logU_iout_files)):
    #    neutral_PAH_reference_objects[counter,:] = np.asarray(read_draine_file(neutral_file))

    #print("[pah/pah_source_create:] building the reference PAH list for ions")
    #for counter,ion_file in tqdm(enumerate(ion_logU_iout_files)):
    #    ion_PAH_reference_objects[counter,:] = np.asarray(read_draine_file(ion_file))
            
              #third, on a cell-by-cell basis, interpolate the luminosity for
              #each grain size bin, and multiply by the number of grains in that
    #bin
    ncells = grid_of_sizes.shape[0]

    total_PAH_luminosity = np.zeros(len(temp_PAH_list[0].lam))
    
    
    #in regions where the radiation field has been poorly sampled (due
    #to low photon count) we can have beta_nnls for the whole cell is
    #0.  then, due to the normalization of beta_nnls in get_beta_nnls,
    #this means NaNs.  so we take those cells and assume equipartition
    #in the draine basis functions.
    beta_nnls[np.isnan(beta_nnls)] = 1./beta_nnls.shape[0]




    #find the indices of the Draine sizes that best match those that are in the simulation
    Draine_simulation_idx_left_edge_array = []
    for size in simulation_sizes.to(u.cm).value:
        idx0 = find_nearest(draine_sizes,size)
        #if draine_sizes[idx0] > size: idx0 -=1
        
        #this is really the nearest point in the Draine sizes to the
        #simulation_sizes. 
        Draine_simulation_idx_left_edge_array.append(idx0)



    #get the indices for where the Draine size bins match ours 
    size_arange = np.arange(len(simulation_sizes))
    draine_bins_idx = np.asarray(Draine_simulation_idx_left_edge_array)[size_arange]

    #pah_grid = np.array([x.lum for x in temp_PAH_list])


    #initialize the process pool and build the chunks
    t1 = datetime.now()
    nprocesses = np.min([cfg.par.n_processes,ncells]) #pool.map will barf in the corner case that we have less cells than cores


    cell_list = np.arange(ncells)

    #chunking to speed up multiprocessing: since the processes are so
    #quick, we can lose a factor of 50% time in just spawning new
    #threads.  it saves a ton of time to chunk up the work and send it all off once.

    #set the number of chunks to be divisble evenly by the number of
    #cells: this will make the concatenation below work for the
    #grid_PAH_luminosities.  this will force a small slowdown if
    #nchunks>nprocessors, but it's not a big penalty.

    nchunks=nprocesses
    print("nchunks = ",nchunks)
    nchunks = int(get_whole_ceil(len(cell_list),nchunks))
    print("modified nchunks = ",nchunks)

    chunk_start_indices = []
    chunk_start_indices.append(0) #the start index is obviously 0
    #this should just be int(ncells/nchunks) but in case ncells < nchunks, we need to ensure that this is at least  1
    delta_chunk_indices = np.max([int(len(cell_list) / nchunks),1])
    print ('delta_chunk_indices = ',delta_chunk_indices)

    for n in range(1,nchunks):
        chunk_start_indices.append(chunk_start_indices[n-1]+delta_chunk_indices)

    list_of_chunks = []
    for n in range(nchunks):
        cells_list_chunk = cell_list[chunk_start_indices[n]:chunk_start_indices[n]+delta_chunk_indices]
        #if we're on the last chunk, we might not have the full list included, so need to make sure that we have that here
        if n == nchunks-1:
            cells_list_chunk = cell_list[chunk_start_indices[n]::]
        list_of_chunks.append(cells_list_chunk)


    
    print("Computing the PAH luminosities for every cell given its grain size distribution and logU. Entering Pool.map multiprocessing.")
    p = mp.Pool(processes = nprocesses)
    dum_numgrains = reg['particle_dust','numgrains'].value 


    




    pah_grid_of_sizes = ds.parameters['reg_grid_of_sizes_graphite']*ds.parameters['reg_grid_of_sizes_aromatic_fraction']

    #ignore the calculation for any sizes bigger than PAH traditional sizes to avoid LIR double counting issue (issue #241 on github)
    sizes = 10.**(np.linspace(cfg.par.otf_extinction_log_min_size,cfg.par.otf_extinction_log_max_size,pah_grid_of_sizes.shape[1])) #micron
    max_pah_size = 2.e-3 #micron -- using Eq 21 from Narayanan et al. 2023 and assuming PAHs <= 1000 carbon atoms
    not_pah_bin_indices = np.where(sizes >= max_pah_size)[0]
    pah_grid_of_sizes[:,not_pah_bin_indices]=1.e-30

    
    #pah_grid_of_sizes = ds.parameters['reg_grid_of_sizes_graphite']


    if cfg.par.PAH_SPA == False:
        dum  = compute_grid_PAH_luminosity(cell_list,
                                           beta_nnls = beta_nnls,
                                           grid_of_sizes = pah_grid_of_sizes.value,
                                           numgrains = dum_numgrains,
                                           draine_sizes = draine_sizes,
                                           draine_lam = draine_lam.value,
                                           f_ion=f_ion,
                                           neutral_PAH_reference_objects = neutral_PAH_reference_objects,
                                           ion_PAH_reference_objects = ion_PAH_reference_objects,
                                           logU = logU,
                                           basis_logU_values = basis_logU_values,
                                           draine_bins_idx = draine_bins_idx)
        

        grid_PAH_luminosity = dum[0]*u.erg/u.s
        grid_neutral_PAH_luminosity = dum[1]*u.erg/u.s
        grid_ion_PAH_luminosity = dum[2]*u.erg/u.s
        pah_lam = draine_lam

    else:

        
        grid_PAH_luminosity, grid_neutral_PAH_luminosity, grid_ion_PAH_luminosity = compute_grid_PAH_luminosity_SPA_parallel(
            cell_list, 
            grid_of_sizes, 
            reg, 
            simulation_sizes, 
            ds, 
            draine_directories, 
            f_ion
        )
        
        #grid_PAH_luminosity, grid_neutral_PAH_luminosity,grid_ion_PAH_luminosity = compute_grid_PAH_luminosity_SPA_serial(cell_list,grid_of_sizes,reg,simulation_sizes,ds,draine_directories,f_ion)
        #get the units of wavelength back out - get the emission wavelengths
        temp_ps = pah_spec.PahSpec()
        SPA_emission_wavelengths = temp_ps.emission_wavelengths
        grid_PAH_luminosity *= SPA_emission_wavelengths.to(u.cm)
        grid_neutral_PAH_luminosity *= SPA_emission_wavelengths.to(u.cm)
        grid_ion_PAH_luminosity *= SPA_emission_wavelengths.to(u.cm)

        pah_lam = SPA_emission_wavelengths.to(u.micron)
        
    '''

#THIS COMMENTED BLOCK IS FOR DOING THE PAH LUMINOSITY COMPUTATIONS IN
PARALLEL WITH PARTIAL.  HONESTLY, IT DOESN'T SEEM A LOT FASTER, AND
WITH THE OVERHEAD IT CAN ACTUALLY CAUSE SLOWDOWNS AT TIME.  WITH SOME
TESTING IT MAY BE WORTH RE-INTRODUCING.

    dum  = p.map(partial(compute_grid_PAH_luminosity,
                         beta_nnls = beta_nnls,
                         grid_of_sizes = pah_grid_of_sizes.value,
                         numgrains = dum_numgrains,
                         draine_sizes = draine_sizes,
                         draine_lam = draine_lam.value,
                         f_ion=f_ion,
                         neutral_PAH_reference_objects = neutral_PAH_reference_objects,
                         ion_PAH_reference_objects = ion_PAH_reference_objects,
                         logU = logU,
                         basis_logU_values = basis_logU_values,
                         draine_bins_idx = draine_bins_idx),[arg for arg in list_of_chunks])
    

    
    #this is some crazy business here, so let me explain.  dum returns
    #a tuple that is nprocesses big.  each element of this tuple is 3
    #elements long, each of which is (ncells/nprocessors,n_draine_lam)
    #long.  the 3 corresponds to [total, neutral, ions]. for exmaple,
    #the [0][0] element of dum is the first PAH emission spectrum
    #chunk (ncells/nprocessors , n_wavelengths) for the total PAH
    #spectrum.  the [0][1] corresponds to the neutrals for the first
    #chunk, and [0][2] the ions.  hence, the following :does a list
    #comprehension on all of the chunks [ so that we have a list of
    #lists, where each sublist is a chunk], nparray's it, and then
    #concatenates on the 0th axis to make one grand array.  

    grid_PAH_luminosity = np.concatenate( np.asarray([dum[i][0] for i in range(len(dum))] ),axis=0)
    grid_neutral_PAH_luminosity = np.concatenate( np.asarray([dum[i][1] for i in range(len(dum))] ),axis=0)
    grid_ion_PAH_luminosity = np.concatenate( np.asarray([dum[i][2] for i in range(len(dum))] ),axis=0)

    '''

    t2 = datetime.now()
    print ('Execution time for PAH dot producting [is that a word?] across the grid = '+str(t2-t1))


    grid_PAH_luminosity[np.isnan(grid_PAH_luminosity)] = 0
    grid_neutral_PAH_luminosity[np.isnan(grid_neutral_PAH_luminosity)] = 0
    grid_ion_PAH_luminosity[np.isnan(grid_ion_PAH_luminosity)] = 0

    nu = (constants.c/pah_lam).to(u.Hz)
    #the units here are Lsun/Hz - this is to be consistent with our
    #stellar fnu addition later. The SEDs of individiual sources all
    #end up getting renormalized by the luminosity, so the exact units
    #don't matter as long as they're consistent across all the sources
    #(and types of sources) being added to the grid.


    fnu = np.divide((grid_PAH_luminosity).to(u.Lsun).value,nu.to(u.Hz).value)

    #Because the Draine templates include re-emission, but we want to
    #add the PAHs as sources only, we restrict to the PAH range.
    nu_reverse = nu[::-1]

    #this is in here to set up the testing/debugging infrastructure
    #for a scenario where we only want to include certain wavelengths.
    #We set it to [0.1,1e3] as a default to catch all the emission.
    nu_pah2 = (constants.c/(0.1*u.micron)).to(u.Hz) #start of the pah range
    nu_pah1 = (constants.c/(1.e3*u.micron)).to(u.Hz) #end of pah range
    wpah_nu_reverse = np.where( (nu_reverse.value < nu_pah2.value) & (nu_reverse.value > nu_pah1.value))[0]




    #set a fnu floor since the 0's can propagate to NaNs
    fnu_floor = np.min(fnu[fnu>0])
    fnu[fnu==0]=fnu_floor

    #to reduce memory requirements, we can't really add a PAH source
    #for every single cell.  so we just do so for cells that are at a
    #luminosity such that the CDF (sum(L>L_threshold) > 99% of the
    #total luminosity).  
 
    only_important_PAH_idx = get_PAH_lum_cdf(nu_reverse,fnu,wpah_nu_reverse,grid_PAH_luminosity)

    
    for i in only_important_PAH_idx:#range(grid_PAH_luminosity.shape[0]): #np.arange(2500)

        fnu_reverse = fnu[i,:][::-1]



        lum = (np.absolute(np.trapz(nu_reverse[wpah_nu_reverse].cgs.value,fnu_reverse[wpah_nu_reverse])).item()*u.Lsun).to(u.erg/u.s).value
        print(lum)

        if lum <= LUM_FLOOR: lum = LUM_FLOOR #just a jamky variable
                                            #defined at the top of
                                            #this function to define a
                                            #lowest luminosity so that
                                            #we don't add PAH cells
                                            #with 0 luminosity
        #reversing arrays to make nu increasing, and therefore correct for hyperion addition

        print(lum)
        
        m.add_point_source(luminosity=lum,spectrum=(nu_reverse[wpah_nu_reverse].value,fnu_reverse[wpah_nu_reverse]),position=reg.parameters['cell_position'][i,:].in_units('cm').value-boost)



    if cfg.par.draine21_pah_grid_write: #else, the try/except in analytics.py will get caught and will just write a single -1 to the output npzfile
        reg.parameters['grid_PAH_luminosity'] = grid_PAH_luminosity
        reg.parameters['grid_neutral_PAH_luminosity'] = grid_neutral_PAH_luminosity
        reg.parameters['grid_ion_PAH_luminosity'] = grid_ion_PAH_luminosity

    reg.parameters['PAH_lam'] = pah_lam.value

    total_PAH_luminosity =np.sum(grid_PAH_luminosity,axis=0)
    total_neutral_PAH_luminosity = np.sum(grid_neutral_PAH_luminosity,axis=0)
    total_ion_PAH_luminosity = np.sum(grid_ion_PAH_luminosity,axis=0)

    reg.parameters['total_PAH_luminosity'] = total_PAH_luminosity
    reg.parameters['total_neutral_PAH_luminosity'] = total_neutral_PAH_luminosity
    reg.parameters['total_ion_PAH_luminosity'] = total_ion_PAH_luminosity

    reg.parameters['only_important_PAH_idx'] = only_important_PAH_idx



    grid_PAH_L_lam = grid_PAH_luminosity/pah_lam
    integrated_grid_PAH_luminosity = np.trapz((grid_PAH_luminosity/pah_lam.value),pah_lam.value,axis=1)
    integrated_grid_neutral_PAH_luminosity = np.trapz((grid_neutral_PAH_luminosity/pah_lam.value),pah_lam.value,axis=1)
    integrated_grid_ion_PAH_luminosity = np.trapz((grid_ion_PAH_luminosity/pah_lam.value),pah_lam.value,axis=1)

    reg.parameters['integrated_grid_PAH_luminosity'] = integrated_grid_PAH_luminosity
    reg.parameters['integrated_grid_neutral_PAH_luminosity'] = integrated_grid_neutral_PAH_luminosity
    reg.parameters['integrated_grid_ion_PAH_luminosity'] = integrated_grid_ion_PAH_luminosity


    #save some information for dumping into analytics
    reg.parameters['q_pah'] = q_pah
    reg.parameters['particle_mass_weighted_gsd'] = particle_mass_weighted_gsd
    reg.parameters['grid_mass_weighted_gsd'] = grid_mass_weighted_gsd
    reg.parameters['simulation_sizes'] = simulation_sizes


    #just for funzies save the beta 
    for i in range(beta_nnls.shape[1]): beta_nnls[:,i]/=np.max(beta_nnls[:,i])
    reg.parameters['beta_nnls'] = beta_nnls
