import astropy.io.fits as ts
from TheCannon import dataset,apogee
from astropy.table import Table
import numpy as np
import os
from astropy.io import fits
import pickle
import AnniesLasso_2 as tc
import time
from os import listdir
from os.path import isfile, join

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.time import Time




## attention!!
# The observer locations are different for different telescopes

"""
APOGEE-Mexicon

    lon = -105.820417
    lat = 32.780361

('-0.0077d', '51.4826')

"""


def MJD2BJD(mjd, target_coord, site_location=(-105.820417,32.780361)):
    """do the conversion
mjd -- input mjd, scale is utc
target_coord -- the coordinate of the target, in astropy.coord format,
to caculate the light travel time
site_location -- location of the telescope, to make accurate calcualtion of
light travel time and tdb conversion. Not very important here. The default value
is for Greenwich Observatory.
"""
    t = Time(mjd, format='mjd', scale='utc', location=site_location)
    # calculate light travel time
    ltt = t.light_travel_time(target_coord)
    # print(t, ltt)
    # convert t to tdb, and add light travel time

    # You can use JD or MJD. It's different
    # Maybe MJD?---BMJD

    t_out = (t.tdb + ltt).jd

    return t_out



## training the cannon

pkl_file = open('wl.pkl', 'rb')
wl = pickle.load(pkl_file)
pkl_file.close()



training_set_path = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/reference_labels.csv"
training_set_spectrum_dir = "/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/Cannon Experiment python 3.5/Data/"


training_set = Table.read("reference_labels.csv")



def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


N = len(training_set)
keep = np.ones(N, dtype=bool)

training_set_flux = []
training_set_ivar = []
training_set_error = []

for i, row in enumerate(training_set):
    image_path = os.path.join(training_set_spectrum_dir, row["ID"])
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False
        continue
    print("{}/{}: {}".format(i + 1, N, image_path))
    image = fits.open(image_path)
    flux = image[1].data
    flux_err = image[2].data
    badpix = get_pixmask(flux, flux_err)
    ivar = 1.0/flux_err**2
    error = flux_err
    # badpix is a array and the length is 8575
    flux[badpix] = np.median(flux)
    ivar[badpix] = 0.0
    training_set_flux.append(flux)
    training_set_ivar.append(ivar)
    training_set_error.append(error)

training_set_flux = np.array(training_set_flux)
training_set_ivar = np.array(training_set_ivar)
training_set_error = np.array(training_set_error)


assert all(keep)

# Construct model.
model = tc.L1RegularizedCannonModel(
    training_set, training_set_flux, training_set_ivar,threads=8)
model.s2 = 0
model.regularization = 0
model.vectorizer = tc.vectorizer.NormalizedPolynomialVectorizer(training_set,
    tc.vectorizer.polynomial.terminator(("Teff_{corr}", "logg_{corr}", "[M/H]_{corr}"), 2))

model.train()



# function
def get_pixmask(flux, err):
    bad_flux = ~np.isfinite(flux)
    bad_err = (~np.isfinite(err)) | (err <= 0)
    bad_pixels = bad_err | bad_flux
    return bad_pixels


##########
#Run them on our stars and count time
##########


def fitting_ve(name):

    image_path = name
    if not os.path.exists(image_path):
        print("{}/{} could not be found: {}".format(i + 1, N, image_path))
        keep[i] = False

    # We only store flux,ivar,inf_flux,parameters,parameters_new,parameters_sim,ve(n*3)(include ve, ve_new,ve_sim)
    try:
        image = fits.open(image_path, ignore_missing_end=True)
        dat = Table.read(image_path)

        flux = image[1].data
        flux_err = image[2].data

        flux = np.atleast_2d(flux)
        flux_err = np.atleast_2d(flux_err)


    except IOError:

        print("opts. This one fail")
        em =0

    else:

        em =1

        badpix = get_pixmask(flux, flux_err)
        ivar = 1.0 / flux_err ** 2
        error = flux_err
        # badpix is a array and the length is 8575
        flux = np.array(flux, dtype=np.float64)
        ivar = np.array(ivar, dtype=np.float64)

        flux[badpix] = np.median(flux)
        ivar[badpix] = 0.0

        flux = np.array(flux)
        ivar = np.array(ivar)

        # normalize flux:
        # value

        tr_ID = image_path

        test_labels_all_i = np.array([5000, 1, 1])

        ds = dataset.Dataset(wl, tr_ID, flux, ivar,
                             test_labels_all_i, tr_ID, flux, ivar)

        ds.ranges = [[371, 3192], [3697, 5997], [6461, 8255]]

        # set sudo-continuous spectrum
        pseudo_tr_flux, pseudo_tr_ivar = ds.continuum_normalize_training_q \
            (q=0.90, delta_lambda=50)

        # set mask
        contmask = ds.make_contmask(pseudo_tr_flux, pseudo_tr_ivar, frac=0.07)

        # get continuous mask

        ds.set_continuum(contmask)

        # fit the normalized-spectrum in the continuous region

        cont = ds.fit_continuum(3, "sinusoid")

        # Obtain the normalized flux
        norm_tr_flux, norm_tr_ivar, norm_test_flux, norm_test_ivar = \
            ds.continuum_normalize(cont)

        norm_tr_flux = np.atleast_2d(norm_tr_flux)

        if len(norm_tr_flux[:,0])<3:
            em=0
        else:
            nothing=1

        # infer labels


        # inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)


        # Use inferred labels from the combined spectra:


        inf_labels = model.fit(norm_tr_flux, norm_tr_ivar)
        # only use the inf labels from the combined spectra

        com = len(inf_labels[:, 0])

        inf_labels_com = inf_labels[0, :]

        inf_labels = []
        for z in range(0, com):
            inf_labels.append(inf_labels_com)

        inf_labels = np.array(inf_labels)

        v = model.vectorizer.get_label_vector(inf_labels)
        inf_flux = np.dot(v, model.theta.T)
        opt_flux, parameters = model.fitting_spectrum_parameters_single \
            (norm_tr_flux, norm_tr_ivar, inf_flux)


        # calculate chi-squared!

        chi_inf = (norm_tr_flux-inf_flux)**2*norm_tr_ivar
        chi_inf = np.sum(chi_inf,axis=1)

        chi_mix = (norm_tr_flux-opt_flux)**2*norm_tr_ivar
        chi_mix = np.sum(chi_mix,axis=1)



        ve = (parameters[:, 2] - parameters[:, 0]) / (parameters[:, 0] + parameters[:, 1] + parameters[:, 2]) * 4144.68

        ve_un = model.uncertainty

        # old
        a0 = parameters
        a1 = ve
        a2 = ve_un

        # covariance matrix for abc
        a3 = model.un_cov

        # spectra

        a4 = norm_tr_flux
        a5 = norm_tr_ivar
        a6 = inf_flux
        a7 = opt_flux

        # inf_labels are from the
        a8 = inf_labels

        a9 = chi_inf

        a10 = chi_mix

        # VHELIO
        a11 = np.array(dat[0]["VHELIO"])

        # Fiber

        a12 = np.array(dat[0]["FIBER"])

        # Files

        # BJD

        RA = image[0].header["RA"]

        DEC = image[0].header["DEC"]

        SNR = image[0].header["SNR"]

        MJD = dat[0]["MJD"]

        c = SkyCoord(RA, DEC, frame='icrs', unit='deg')

        BJD = MJD2BJD(MJD, c)

        a13 = np.array(BJD)

        # calculate chi-squared:


        try:
            # save them

            # pay attention to the fits file saving

            path_fits_i = image_path.replace("/Volumes/Data_2TB/Data/DR13_rc/apStar-r6-",
                                             "/Users/caojunzhi/Desktop/Data/dr13_red_clump/")

            print("saving files" + path_fits_i)

            hdu = fits.PrimaryHDU(data=a0)
            hdu.header[
                'COMMENT'] = "Simple orange juice"

            # add header info

            hdu.header['SNR'] = SNR
            hdu.header['RA'] = RA
            hdu.header['DEC'] = DEC

            hdu.writeto(path_fits_i, clobber=True)

            ts.append(path_fits_i, a1)
            ts.append(path_fits_i, a2)
            ts.append(path_fits_i, a3)
            ts.append(path_fits_i, a4)
            ts.append(path_fits_i, a5)
            ts.append(path_fits_i, a6)
            ts.append(path_fits_i, a7)
            ts.append(path_fits_i, a8)

            ts.append(path_fits_i, a9)
            ts.append(path_fits_i, a10)
            ts.append(path_fits_i, a11)
            ts.append(path_fits_i, a12)
            ts.append(path_fits_i, a13)

        except OSError:
            print("fail")
            em=0

    return em



start_time = time.time()


pkl_file = open('Red_clump_dr13_ori.pkl', 'rb')
rc_path = pickle.load(pkl_file)
pkl_file.close()

length = len(rc_path)
print(length)
print(rc_path)

rc_path = rc_path[1:100]

path_fits_dr13 = []
path_ori_dr13 = []
save_path = "/Users/caojunzhi/Desktop/Data/dr13_red_clump/"

for i in range(0,len(rc_path)):

    original_i = rc_path[i]

    # opt
    # input ori
    print("doing star %d"%i)
    em = fitting_ve(name=original_i)


stop_time = time.time()

print("The time we use %.2f s"%(stop_time-start_time))















