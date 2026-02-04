from __future__ import print_function
from astropy import constants
from powderday.image_processing import add_transmission_filters, convolve
from astropy import units as u
import powderday.config as cfg
import numpy as np

def make_SED(m, par, model):
    # set up the SEDs and images
    
    if cfg.par.SED_MONOCHROMATIC == True:

        # since all sources have the same spectrum just take the nu
        # from the input SED from the first source

        monochromatic_nu = m.sources[0].spectrum['nu'] * u.Hz
        monochromatic_lam = (constants.c / monochromatic_nu).to(u.micron).value[::-1]

        if cfg.par.FIX_SED_MONOCHROMATIC_WAVELENGTHS == True:
            # idx = np.round(np.linspace(np.min(np.where(monochromatic_lam > cfg.par.SED_MONOCHROMATIC_min_lam)[0]),\
            ##                           np.max(np.where(monochromatic_lam < cfg.par.SED_MONOCHROMATIC_max_lam)[0]),\
            #                           cfg.par.SED_MONOCHROMATIC_nlam))

            idx = np.where((monochromatic_lam > cfg.par.SED_MONOCHROMATIC_min_lam) & (
                        monochromatic_lam < cfg.par.SED_MONOCHROMATIC_max_lam))[0]
            monochromatic_lam = np.take(monochromatic_lam, list(idx))
        m.set_monochromatic(True, wavelengths=monochromatic_lam)
        m.set_raytracing(True)
        m.set_n_photons(initial=par.n_photons_initial,
                        imaging_sources=par.n_photons_imaging,
                        imaging_dust=par.n_photons_imaging,
                        raytracing_sources=par.n_photons_raytracing_sources,
                        raytracing_dust=par.n_photons_raytracing_dust)

        m.set_n_initial_iterations(3)
        m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)
        m.set_copy_input(False)
        sed = m.add_peeled_images(sed=True, image=False)

        if cfg.par.MANUAL_ORIENTATION == True:
            sed.set_viewing_angles(np.array(cfg.par.THETA), np.array(cfg.par.PHI))

        else:
            sed.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist() * par.NPHI,
                                   np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))
        sed.set_track_origin('basic')

        if cfg.par.SKIP_RT == False:
            m.write(model.inputfile + '.sed', overwrite=True)
            if par.n_MPI_processes > 1:
                m.run(model.outputfile + '.sed', mpi=True, n_processes=par.n_MPI_processes, overwrite=True)
            else:
                m.run(model.outputfile + '.sed', mpi=False, overwrite=True)

        print(
            '[pd_front_end]: Beginning RT Stage: Calculating SED using a monochromatic spectrum equal to the input SED')

    else:

        m.set_raytracing(True)
        m.set_n_photons(initial=par.n_photons_initial, imaging=par.n_photons_imaging,
                        raytracing_sources=par.n_photons_raytracing_sources,
                        raytracing_dust=par.n_photons_raytracing_dust)
        m.set_n_initial_iterations(7)
        m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

        sed = m.add_peeled_images(sed=True, image=False)
        sed.set_wavelength_range(2500, 0.001, 1000.)

        if cfg.par.MANUAL_ORIENTATION == True:
            sed.set_viewing_angles(np.array(cfg.par.THETA), np.array(cfg.par.PHI))
        else:
            sed.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist(
            ) * par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))
        sed.set_track_origin('basic')

        print('[pd_front_end]: Beginning RT Stage: Calculating SED using a binned spectrum')

        # Run the Model
        if cfg.par.SKIP_RT == False:
            m.write(model.inputfile + '.sed', overwrite=True)
            if par.n_MPI_processes > 1:
                m.run(model.outputfile + '.sed', mpi=True,
                      n_processes=par.n_MPI_processes, overwrite=True)
            else:
                m.run(model.outputfile + '.sed', mpi=False, overwrite=True)



def make_DIG_SED(m, par, model):
    m.set_raytracing(True)
    m.set_n_photons(initial=par.n_photons_DIG, imaging=par.n_photons_DIG,
                    raytracing_sources=par.n_photons_DIG,
                    raytracing_dust=par.n_photons_DIG)
    m.set_n_initial_iterations(7)
    m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_wavelength_range(2500, 0.001, 1000.)

    sed.set_viewing_angles(np.linspace(0, 90, 1).tolist(), np.repeat(np.linspace(0, 90, 1), 1))
    sed.set_track_origin('basic')

    print('[pd_front_end]: Beginning RT Stage: For DIG calculation')

    # Run the Model
    m.write(model.inputfile + '_DIG_energy_dumped.sed', overwrite=True)
    if par.n_MPI_processes > 1:
        m.run(model.outputfile + '_DIG_energy_dumped.sed', mpi=True, n_processes=par.n_MPI_processes, overwrite=True)
    else:
        m.run(model.outputfile + '_DIG_energy_dumped.sed', mpi=False, overwrite=True)
    
    print('[pd_front_end]: RT Stage For DIG calculation has ended')

def compute_ISRF_SED(m, par, model):
    m.set_raytracing(True)
    m.set_n_photons(initial=par.n_photons_initial, imaging=par.n_photons_imaging,
                    raytracing_sources=par.n_photons_raytracing_sources,
                    raytracing_dust=par.n_photons_raytracing_dust)
    m.set_n_initial_iterations(7)
    m.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

    sed = m.add_peeled_images(sed=True, image=False)
    sed.set_wavelength_range(2500, 0.001, 1000.)
    #sed.set_viewing_angles(np.linspace(0, 90, 1).tolist(), np.repeat(np.linspace(0, 90, 1), 1))
    sed.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist(
            ) * par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))
    sed.set_track_origin('basic')

    print('[pd_front_end]: Beginning RT Stage: For ISRF calculation')

    # Run the Model
    m.write(model.inputfile + '_isrf.sed', overwrite=True)
    if par.n_MPI_processes > 1:
        m.run(model.outputfile + '_isrf.sed', mpi=True,n_processes=par.n_MPI_processes, overwrite=True)
    else:
        m.run(model.outputfile + '_isrf.sed', mpi=False, overwrite=True)

    
def make_image(m_imaging, par, model,dx,dy,dz):
    print("Beginning Monochromatic Imaging RT")
    
    if cfg.par.IMAGING_TRANSMISSION_FILTER == False:
        
        # read in the filters file
        try:
            filter_data = [np.loadtxt(par.filterdir+f) for f in par.filterfiles]
        except:
            raise ValueError("Filters not found. You may be running above changeset 'f1f16eb' with an outdated parameters_master file. Please update to the most recent parameters_master format or ensure that the'filterdir' and 'filterfiles' parameters are set properly.")
            
        # Extract and flatten all wavelengths in the filter files

        wavs = []
        for single_filter in par.filterfiles:
            single_filter_data = np.loadtxt(par.filterdir+'/'+single_filter)
            wavs.append(single_filter_data[:,0])


        wavs = np.unique(np.asarray(wavs)) #remove duplicates for efficiency

        #wavs = [wav[0] for single_filter in filter_data for wav in single_filter]
        #wavs = list(set(wavs))      # Remove duplicates, if they exist

        m_imaging.set_monochromatic(True, wavelengths=wavs)
        m_imaging.set_raytracing(True)
        m_imaging.set_n_photons(initial=par.n_photons_initial,
                                imaging_sources=par.n_photons_imaging,
                                imaging_dust=par.n_photons_imaging,
                                raytracing_sources=par.n_photons_raytracing_sources,
                                raytracing_dust=par.n_photons_raytracing_dust)
    else:
        m_imaging.set_n_photons(initial=par.n_photons_initial, imaging=par.n_photons_imaging)

    m_imaging.set_n_initial_iterations(7)
    m_imaging.set_convergence(True, percentile=99., absolute=1.01, relative=1.01)

    image = m_imaging.add_peeled_images(sed=True, image=True)

    if cfg.par.IMAGING_TRANSMISSION_FILTER == True:
        add_transmission_filters(image)

    if cfg.par.MANUAL_ORIENTATION == True:
        image.set_viewing_angles(np.array(cfg.par.THETA), np.array(cfg.par.PHI))
    else:
        image.set_viewing_angles(np.linspace(0, 90, par.NTHETA).tolist()*par.NPHI, np.repeat(np.linspace(0, 90, par.NPHI), par.NPHI))

    image.set_track_origin('basic')
    image.set_image_size(cfg.par.npix_x, cfg.par.npix_y)
    image.set_image_limits(-dx/2., dx/2., -dy/2., dy/2.)

    m_imaging.write(model.inputfile+'.image', overwrite=True)
    if par.n_MPI_processes > 1:
        m_imaging.run(model.outputfile+'.image', mpi=True, n_processes=par.n_MPI_processes, overwrite=True)
    else:
        m_imaging.run(model.outputfile+'.image', mpi=False, overwrite=True)
    
    convolve(model.outputfile+'.image', par.filterfiles, filter_data)
