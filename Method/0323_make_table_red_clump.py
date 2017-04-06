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



fail = 0

class make_table():

    def read_data(self):


        pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_rc_fits.pkl', 'rb')
        path_fits = pickle.load(pkl_file)
        pkl_file.close()


        pkl_file = open('/Users/caojunzhi/Desktop/NYU/Laboratory/task 2016.8.1-12.23/My codes/AnniesLasso_v2/dr13_rc_ori.pkl', 'rb')
        path_ori = pickle.load(pkl_file)
        pkl_file.close()


        N = len(path_fits)
        print(N)


        star_name = []
        star_visit = []

        BJD = []
        inf_labels = np.array([0,0,0])
        chi_inf = []
        chi_mix = []

        parameters= np.array([0,1,0])

        un_cov = np.zeros((3, 3))

        fiber= []
        VBARY = []
        SHIFT = []

        SNR = []



        for i in range(0,N):

            try:

                print("reading star %d  %.2f %%" % (i + 1, (i + 1) / N * 100))

                star = fits.open(path_fits[i])

                dat = Table.read(path_ori[i])

                ni = len(star[1].data)

                star_name_i = path_fits[i].replace(".fits", "")
                star_name_i = star_name_i.replace("/Users/caojunzhi/Desktop/Data/dr13_red_clump/", "")

                SNR_i = []
                star_name_ii = []
                for j in range(0, ni - 2):
                    star_name_ii = np.append(star_name_ii,star_name_i)
                    SNR_i = np.append(SNR_i, star[0].header["SNR"])



                star_visit_i = np.array(dat[0]["FILE"])


                chi_inf_i = star[9].data[2:ni]


                chi_mix_i = star[10].data[2:ni]


                BJD_i = star[13].data


                fiber_i = star[12].data


                VBARY_i = star[11].data


                SHIFT_i = star[1].data[2:ni]



                inf_labels_i = star[8].data[2:ni, :]


                parameters_i = star[0].data[2:ni, :]


                un_cov_i = star[3].data[:, :, 2:ni]



            except IndexError:
                print("This one fails")



            except OSError:
                print("This one fails")




            else:
                SNR = np.append(SNR,SNR_i)

                star_visit = np.append(star_visit, star_visit_i)

                star_name = np.append(star_name,star_name_ii)

                chi_inf = np.append(chi_inf, chi_inf_i)

                chi_mix = np.append(chi_mix, chi_mix_i)

                BJD = np.append(BJD, BJD_i)

                fiber = np.append(fiber, fiber_i)

                VBARY = np.append(VBARY, VBARY_i)

                SHIFT = np.append(SHIFT, SHIFT_i)

                inf_labels = np.vstack((inf_labels, inf_labels_i))

                parameters = np.vstack((parameters, parameters_i))

                un_cov = np.dstack((un_cov, un_cov_i))




        n_visit = len(un_cov[0,0,:])



        star_name = np.array(star_name)
        self.star_name = star_name

        SNR = np.array(SNR)
        self.SNR =SNR

        star_visit = np.array(star_visit)
        star_visit = star_visit.ravel()
        self.star_visit = star_visit

        chi_inf = np.array(chi_inf)
        chi_inf = chi_inf.ravel()
        self.chi_inf = chi_inf

        chi_mix = np.array(chi_mix)
        chi_mix = chi_mix.ravel()
        self.chi_mix = chi_mix

        BJD = np.array(BJD).ravel()
        self.BJD = BJD

        fiber = np.array(fiber).ravel()
        self.fiber = fiber

        VBARY = np.array(VBARY).ravel()
        self.VBARY = VBARY

        SHIFT = np.array(SHIFT).ravel()
        self.SHIFT =SHIFT

        inf_labels = inf_labels[1:n_visit,:]
        self.inf_labels = inf_labels

        parameters = parameters[1:n_visit,:]
        self.parameters = parameters


        un_cov = un_cov[:,:,1:n_visit]
        self.un_cov = un_cov

        print("check shape")

        print(n_visit - 1)

        print(star_name.shape,SNR.shape,star_visit.shape,chi_inf.shape,chi_mix.shape)
        print(BJD.shape,fiber.shape,inf_labels.shape,parameters.shape,un_cov.shape,VBARY.shape,SHIFT.shape)

    # Input the path of the table:

    def write_table(self):

        # save them in the header


        path = "/Users/caojunzhi/Downloads/upload_20170322/red_clump_dr13.fits"

        self.table_path = path

        prihdr = fits.Header()
        prihdr['COMMENT'] = "id MJD HJD VELIO VELIOUN Shift ShiftUN Shiftabs ShiftabsUN"

        prihdu = fits.PrimaryHDU(data=self.un_cov,header=prihdr)

        # Table list


        col1 = fits.Column(name='APOGEEID', format="25A", array=self.star_name)
        col2 = fits.Column(name='VISIT', format="25A", array=self.star_visit)
        col3 = fits.Column(name='BJD', format="E", array= self.BJD)

        col4 = fits.Column(name='TEFF', format='E', array=self.inf_labels[:,0])
        col5 = fits.Column(name='LOGG', format="E", array=self.inf_labels[:,1])
        col6 = fits.Column(name='FEH', format='E', array=self.inf_labels[:,2])

        col7 = fits.Column(name='A', format="E", array=self.parameters[:,0])
        col8 = fits.Column(name='B', format="E", array=self.parameters[:,1])
        col9 = fits.Column(name='C', format="E", array=self.parameters[:,2])

        col10 = fits.Column(name='CHIINF', format="E", array=self.chi_inf)
        col11 = fits.Column(name='CHIMIX', format="E", array=self.chi_mix)

        col12 = fits.Column(name='VBARY', format="E", array=self.VBARY)
        col13 = fits.Column(name='VSHIFT', format="E", array=self.SHIFT/1000)

        col14 = fits.Column(name='FIBER', format="E", array=self.fiber)
        col15 = fits.Column(name='SNR', format="E", array=self.SNR)



        cols = fits.ColDefs(
            [col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13, col14, col15])

        tbhdu = fits.BinTableHDU.from_columns(cols)

        thdulist = fits.HDUList([prihdu, tbhdu])

        print("saving table")

        thdulist.writeto(path, clobber=True)
        thdulist.close()

    def check(self,path):

        star = fits.open(path)
        table = Table.read(path)
        self.star = star
        self.table = table

        print(table[0])
        print(star[0].data.shape)




# construct table:

start_time = time.time()

model = make_table()
model.read_data()
model.write_table()

path = "/Users/caojunzhi/Downloads/upload_20170322/red_clump_dr13.fits"
#model.check(path)

print("The number of fails %d"%fail)

stop_time = time.time()

print("The time we use is %.2f"%(stop_time-start_time))

