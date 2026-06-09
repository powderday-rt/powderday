import powderday.config as cfg
import numpy as np
from astropy import units as u
from hyperion.model import ModelOutput
import h5py


def add_transmission_filters(image):

    for i in range(len(cfg.par.filterfiles)):
        lam, throughput = np.loadtxt(cfg.par.filterfiles[i], unpack=True)
        f = image.add_filter()
        f.name = 'dum'
        lam /= (1.+cfg.par.TRANSMISSION_FILTER_REDSHIFT)
        f.spectral_coord = lam * u.micron
        f.transmission = throughput * u.percent
        f.detector_type = 'energy'
        f.alpha = 1
        f.central_spectral_coord = np.mean(lam)*u.micron

    return None


def convolve(image_file, filterfilenames, filter_data):
    # Load the model output object
    m = ModelOutput(image_file)
    # Get the image
    image = m.get_image(units='ergs/s')
    if(image.val.shape[0]>1):
        for inc in range(image.val.shape[0]):
            # Get image bounds for correct scaling
            w = image.x_max * u.cm
            w = w.to(u.kpc)
        
            # This is where the convolved images will go
            image_data = []
        
            # List the filters that shouldn't be used in convolution
            skip_conv = ['arbitrary.filter', 'pdfilters.dat']
        
            # Loop through the filters and match wavelengths to those in the image
            for i in range(len(filterfilenames)):
        
                # Skip "arbitrary.filter" if it is selected
                if filterfilenames[i] in skip_conv:
                    print(" Skipping convolution of default filter")
                    continue
        
                print("\n Convolving filter {}...".format(filterfilenames[i]))
                wavs = filter_data[i][:, 0]
        
                # Figure out which indices of the image wavelengths correspond to
                # this filter
                indices = []
                for wav in wavs:
                    diffs = np.abs(image.wav - wav)
        
                    # Make sure the closest wavelength is *really* close --- there
                    # could be rounding errors, but we don't want to accidentally grab
                    # the wrong wavelength
                    if min(diffs) <= 1e-10:
                        indices.append(diffs.argmin())
        
                if len(indices) != len(wavs):
                    raise ValueError(
                        "Filter wavelength mismatch with available image wavelengths")
        
                # Get the monochromatic images at each wavelength in the filter
                images = [image.val[inc, :, :, j] for j in indices]
                print(' Found {} monochromatic images'.format(len(images)))
        
                # Show wavelengths and weights from filter file
                wavelengths = [image.wav[j] for j in indices]
                weights = filter_data[i][:, 1]
        
                print('\n Wavelength              Weight')
                print(' ----------              ------')
                for k in range(len(wavelengths)):
                    print('  {:.2E}              {:.2E}'.format(wavelengths[k],
                                                                weights[k]))
        
                # Apply appropriate transmissivities from filter file
                image_data.append(np.average(images, axis=0, weights=weights))        
            # Save the image data and filter information as an .hdf5 file
            f = h5py.File(cfg.model.PD_output_dir+f"convolved.inc.{inc}" +
                          cfg.model.snapnum_str+".hdf5", "w")
            f.create_dataset("image_data", data=image_data)
            f['image_data'].attrs['width'] = w.value
            f['image_data'].attrs['width_unit'] = np.bytes_('kpc')
        
            # Don't add the names of filters that were skipped
            trimmed_names = list(set(filterfilenames) - set(skip_conv))
            # Encode in utf8 so hdf5 file can process and append names (doesn't accept strings)
            names = [n.encode('utf8') for n in trimmed_names]
            f.create_dataset("filter_names", data=names)
        
            for i in range(len(filterfilenames)):
                f.create_dataset(filterfilenames[i], data=filter_data[i])
            f.close()
    else:
        # Get image bounds for correct scaling
        w = image.x_max * u.cm
        w = w.to(u.kpc)
    
        # This is where the convolved images will go
        image_data = []
    
        # List the filters that shouldn't be used in convolution
        skip_conv = ['arbitrary.filter', 'pdfilters.dat']
    
        # Loop through the filters and match wavelengths to those in the image
        for i in range(len(filterfilenames)):
    
            # Skip "arbitrary.filter" if it is selected
            if filterfilenames[i] in skip_conv:
                print(" Skipping convolution of default filter")
                continue
    
            print("\n Convolving filter {}...".format(filterfilenames[i]))
            wavs = filter_data[i][:, 0]
    
            # Figure out which indices of the image wavelengths correspond to
            # this filter
            indices = []
            for wav in wavs:
                diffs = np.abs(image.wav - wav)
    
                # Make sure the closest wavelength is *really* close --- there
                # could be rounding errors, but we don't want to accidentally grab
                # the wrong wavelength
                if min(diffs) <= 1e-10:
                    indices.append(diffs.argmin())
    
            if len(indices) != len(wavs):
                raise ValueError(
                    "Filter wavelength mismatch with available image wavelengths")
    
            # Get the monochromatic images at each wavelength in the filter
            images = [image.val[0, :, :, j] for j in indices]
            print(' Found {} monochromatic images'.format(len(images)))
    
            # Show wavelengths and weights from filter file
            wavelengths = [image.wav[j] for j in indices]
            weights = filter_data[i][:, 1]
    
            print('\n Wavelength              Weight')
            print(' ----------              ------')
            for k in range(len(wavelengths)):
                print('  {:.2E}              {:.2E}'.format(wavelengths[k],
                                                            weights[k]))
    
            # Apply appropriate transmissivities from filter file
            image_data.append(np.average(images, axis=0, weights=weights))    
        # Save the image data and filter information as an .hdf5 file
        f = h5py.File(cfg.model.PD_output_dir+"convolved." +
                      cfg.model.snapnum_str+".hdf5", "w")
        f.create_dataset("image_data", data=image_data)
        f['image_data'].attrs['width'] = w.value
        f['image_data'].attrs['width_unit'] = np.bytes_('kpc')
    
        # Don't add the names of filters that were skipped
        trimmed_names = list(set(filterfilenames) - set(skip_conv))
        # Encode in utf8 so hdf5 file can process and append names (doesn't accept strings)
        names = [n.encode('utf8') for n in trimmed_names]
        f.create_dataset("filter_names", data=names)
    
        for i in range(len(filterfilenames)):
            f.create_dataset(filterfilenames[i], data=filter_data[i])
        f.close()
