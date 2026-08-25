#This file holds the function that adds the active dust grids.  This
#file is called from the tributaries, and contains all the common
#lines of code for adding dust grids once the grid_of_sizes has been
#established.

from __future__ import print_function
import numpy as np
import powderday.config as cfg
from hyperion.dust import SphericalDust
import pdb
from powderday.helpers import find_nearest
from powderday.active_dust.dust_file_writer import *

def active_dust_add(ds,m,grid_of_sizes,nsizes,dustdens,specific_energy,refined=[False],grid_of_sizes_graphite = [-1], grid_of_sizes_silicates = [-1], grid_of_sizes_aromatic_fraction = [-1]):

        #go ahead and call the active dust writer to write dust
        #extinction files at the exact sizes of the hydro sim.  this
        #will help later convolution if needed in the PAH modules
        dust_file_writer(nsizes)
        

        #first, save the grid_of_sizes to the ds.paramteters so we can carry it around
        ds.parameters['reg_grid_of_sizes'] = grid_of_sizes #named 'reg_grid_of_sizes' 
        ds.parameters['reg_grid_of_sizes_graphite'] = grid_of_sizes_graphite
        ds.parameters['reg_grid_of_sizes_silicate'] = grid_of_sizes_silicates
        ds.parameters['reg_grid_of_sizes_aromatic_fraction'] = grid_of_sizes_aromatic_fraction

        #empty size bins are intentionally left at zero here.  cells
        #that carry dust mass but have no resolved grain-size
        #distribution are handled below at the fraction level, where
        #they fall back to the grid-integrated size distribution so
        #their dust mass is still partitioned across bins and conserved.


        #now load the mapping between grain bin and filename for the lookup table
        data = np.load(cfg.model.PD_output_dir+'/dust_files/binned_dust_sizes.npz')
        grain_size_left_edge_array = data['grain_size_left_edge_array']
        grain_size_right_edge_array  = data['grain_size_right_edge_array']
        dust_filenames = data['outfile_filenames']

        nbins = len(grain_size_left_edge_array)

        #dust_file_writer may split each size bin into per-species dust
        #files (graphite, silicate); n_species tells us how many, and the
        #flat dust-type ordering is size-major (i_size*n_species+i_species).
        #older dust files predate this and are single-species.
        n_species = int(data['n_species']) if 'n_species' in data.files else 1




        #find which sizes in the hydro simulation correspond to the
        #pre-binned extinction law sizes from dust_file_writer.py

        #each dust file describes a bin centred on one of the simulation
        #grain sizes, so the size a bin represents is its centre, not its
        #left edge.  matching on the centre keeps every simulation size
        #bin mapped to exactly one dust file.
        dust_file_to_grain_size_mapping_idx = []
        x=np.linspace(cfg.par.otf_extinction_log_min_size,cfg.par.otf_extinction_log_max_size,nsizes)
        grain_size_center_array = 0.5*(grain_size_left_edge_array+grain_size_right_edge_array)
        for i in range(nbins):
                dust_file_to_grain_size_mapping_idx.append(find_nearest(x,grain_size_center_array[i]))


        #set up the frac array that is nbins big.  this is the
        #fractional contribution of each dust file bin which is based
        #on the total number of grains in the grid in that bin.

        #frac =np.zeros([dustdens.shape[0],nbins])

        dsf_grid = np.zeros([dustdens.shape[0],nbins])
        frac_grid = np.zeros([dustdens.shape[0],nbins])
        debug_nearest_extinction_curve = np.zeros([nbins])


        if cfg.par.OTF_EXTINCTION_MRN_FORCE == True:
                if n_species > 1:
                        raise NotImplementedError(
                            "[tributary_dust_add:] OTF_EXTINCTION_MRN_FORCE "
                            "carries no silicate/graphite information and is "
                            "incompatible with per-species dust files "
                            "(n_species=%d). Regenerate single-species dust "
                            "files or disable MRN_FORCE." % n_species)
                grid_sum = np.zeros(nbins)

                #how DNSF was set up.  not needed other than for testing
                #x=np.linspace(-4,0,41)
                #load an example dust size function for testing against
                dsf = np.loadtxt(cfg.par.pd_source_dir+'/powderday/active_dust/mrn_dn.txt')#DNSF_example.txt')
                if dsf.shape[0] != nsizes:
                        raise Exception("[tributary_dust_add:] You have enabled the experimental feature OTF_EXTINCTION_MRN_FORCE. Here, the MRN distribution that we are assuming does not have the same shape as the grid size distribution which can cause trouble if the Draine PAH model is enabled.  Therefore, please re-run your MRN grid generator in [active_dust/mrn_test_writer] with the following number of grid sizes: ",nsizes)

                        #raise Exception("[tributary_dust_add:] You have enabled the experimental feature OTF_EXTINCTION_MRN_FORCE.")
#  Here, the MRN distribution that we are assuming does not have the same shape as the grid size distribution which can cause trouble if the Draine PAH model is enabled.  Therefore, please re-run your MRN grid generator in [active_dust/] with the following number of grid sizes")



                #nbins = len(grain_size_left_edge_array)


                for i in range(nbins):
                        #find the index bounds in x that we want to interpolate between
                        idx0 = find_nearest(x,grain_size_left_edge_array[i])
                        if x[idx0] > grain_size_left_edge_array[i]: idx0 -= 1
                        idx1 = idx0+1
                
                        dsf_interp = np.interp(grain_size_left_edge_array[i],[x[idx0],x[idx1]],[dsf[idx0],dsf[idx1]])
                
                        #this sets the fraction of each bin size we need (for the
                        #entire grid!)
                        dsf_grid[:,i] = dsf_interp
                        grid_sum[i] = np.sum(dsf_grid[:,i])
                        debug_nearest_extinction_curve[i] = dsf_interp


                #the forced MRN distribution is a grain-number
                #distribution (dn per size bin) and is spatially
                #uniform, so convert it to a dust-mass fraction per
                #size bin (mass ~ dn * a^3) and give every cell that
                #same size partition.
                a_bin3 = (10.**grain_size_left_edge_array)**3.
                bin_mass = dsf_grid[0,:] * a_bin3
                frac = bin_mass/np.sum(bin_mass)

                for i in range(nbins):
                        frac_grid[:,i] = frac[i]

                #in an octree the refined (parent) cells carry no
                #density, so make sure they are never assigned any dust.
                if np.sum(refined) > 0:
                        frac_grid[np.asarray(refined) != 0,:] = 0.

                '''
                import matplotlib.pyplot as plt
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.plot(x,dsf,label='dsf')
                ax.plot(grain_size_left_edge_array,frac_grid[0,:],label='frac_grid')
                ax.plot(grain_size_left_edge_array,grid_sum,label='grid_sum')
                ax.plot(grain_size_left_edge_array,debug_nearest_extinction_curve,label='d_n_e_c')
                ax.set_yscale('log')
                plt.legend()
                fig.savefig('junk.png',dpi=300)
                
                import pdb
                pdb.set_trace()
                '''

                #------------------------    
        
        else:


                #each hyperion dust type corresponds to a single grain-size
                #bin, so a cell's total dust density has to be divided
                #among the bins according to how much of the cell's dust
                #mass lives at each size.  the grid stores grain *counts*,
                #and the mass in a bin is proportional to N * a^3, so the
                #per-cell size distribution is the grain count weighted by
                #a^3.  a single grain material density is assumed and
                #cancels in the per-cell normalisation below.
                a_micron = 10.**x
                bin_size_idx = np.asarray(dust_file_to_grain_size_mapping_idx)
                a_bin3 = (a_micron[bin_size_idx])**3.

                if n_species > 1:
                        #graphite and silicate are separate hyperion dust
                        #types (one file each per size bin), so the per-cell
                        #mass is partitioned over BOTH size and species from
                        #the simulation's own grain populations.
                        
                        #species order (graphite, silicate) matches
                        #dust_file_writer; the third axis is flattened
                        #size-major at the hyperion boundary below.
                        gra = np.asarray(grid_of_sizes_graphite)[:,bin_size_idx]
                        sil = np.asarray(grid_of_sizes_silicates)[:,bin_size_idx]
                        cell_mass = np.stack([gra,sil],axis=2) * a_bin3[None,:,None]
                        cell_mass_total = np.sum(cell_mass,axis=(1,2))
                        nonempty = cell_mass_total > 0
                        frac = np.sum(cell_mass,axis=0) / np.sum(cell_mass)   #(nbins,n_species) fallback
                        fg = np.zeros_like(cell_mass)
                        fg[nonempty] = cell_mass[nonempty]/cell_mass_total[nonempty,None,None]
                        fg[~nonempty] = frac
                        if np.sum(refined) > 0:
                                fg[np.asarray(refined) != 0] = 0.
                        #C-order flatten of (nbins, n_species) -> size-major,
                        #matching the dust-file flat index i_size*n_species+i_species
                        frac_grid = fg.reshape(fg.shape[0], nbins*n_species)
                else:
                        grain_counts = grid_of_sizes[:,bin_size_idx]
                        cell_mass = grain_counts * a_bin3[None,:]

                        #normalise within each cell so the size fractions sum
                        #to one and the cell's total dust mass is conserved.
                        cell_mass_total = np.sum(cell_mass,axis=1)
                        nonempty = cell_mass_total > 0

                        #grid-integrated mass fraction per size bin, used both
                        #as a reference and as the fallback size distribution
                        #for cells that carry dust mass but have no resolved
                        #grain sizes (so their mass is still partitioned).
                        frac = np.sum(cell_mass,axis=0) / np.sum(cell_mass)

                        frac_grid[nonempty,:]  = cell_mass[nonempty,:] / cell_mass_total[nonempty,None]
                        frac_grid[~nonempty,:] = frac

                        #in an octree the refined (parent) cells carry no
                        #density, so make sure they are never assigned any dust.
                        if np.sum(refined) > 0:
                                frac_grid[np.asarray(refined) != 0,:] = 0.

        #now add the dust grids to hyperion. frac_grid has one column per
        #hyperion dust type (nbins*n_species, flattened size-major) and its
        #width must equal the number of dust files, in the same order.
        assert frac_grid.shape[1] == len(dust_filenames), \
            "[tributary_dust_add:] frac_grid width %d != n dust files %d" % (
                frac_grid.shape[1], len(dust_filenames))
        for bin in range(len(dust_filenames)):
                file = dust_filenames[bin]

                d = SphericalDust(file)


                m.add_density_grid(dustdens*frac_grid[:,bin],d,specific_energy=specific_energy)
                        #m.add_density_grid(dustdens*frac[bin],d,specific_energy=specific_energy)

                        
        #finally, re-save the grid_of_sizes and grain sizes to the ds.paramteters so we can carry it around
        ds.parameters['reg_grid_of_sizes'] = grid_of_sizes #named 'reg_grid_of_sizes'
        ds.parameters['grain_sizes_in_micron'] = 10.**(x)

