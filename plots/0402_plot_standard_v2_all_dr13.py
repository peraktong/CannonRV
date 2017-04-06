import numpy as np
from astropy.table import Table
from astropy.io import fits

import matplotlib.pyplot as plt
import matplotlib
import pickle
from matplotlib import cm
from numpy.random import randn

# table path

path = "/Users/caojunzhi/Downloads/upload_20170330_all/dr13.fits"


star = fits.open(path)
table = Table.read(path)

"""
There are 13 columns in the table:

1. 'APOGEEID' -- The name of the star
2. 'VISIT' -- The name of the visit file
3. BJD -- Barycentric JD

Inferred labels are from the Cannon. The spectra we use are from the first combined spectra
(There are two combined spectra for each star, which are obtained by two different methods)

: (1) global weighting, where each visit spectrum is weighted by its (S/N)2, and
(2) pixel-by-pixel weighting, where each pixel is weighted by its (S/N)2.

4. TEFF
5. LOGG
6. FEH

The abc parameters for each visit:

7. A -- parameter a
8. B -- parameter b
9. C -- parameter c

10. CHIINF -- chi-squared for the inferred flux from the cannon (a=0,b=1,c=0)
11. CHIMIX -- chi-squared for the mixed flux from the abc fit.

12. VBARY -- The barycentric Velocity(km/s) from the APOGEE team.
13. VSHIFT -- The velocity shift from the abc fit(km/s)

14. FIBER -- Fiber ID
15. SNR  -- SNR of the visit
####
The covariance matrix of the abc fit is in HDU0 data, which is
a 3*3*N 3-d matrix. N is the number of visits.
###
"""

# read covariance matrix from the abc fit:

un_cov = star[0].data[:,:,0]

#print(un_cov)


# read the velocity shift from the abc fit
v_shift = table["VSHIFT"]
#print(v_shift.shape)

########################
#Read table and plot to check.

class plot():

    def read_table(self):
        path = "/Users/caojunzhi/Downloads/upload_20170330_all/dr13.fits"

        star = fits.open(path)
        table = Table.read(path)

        print("reading data")

        # read it:

        un_cov = star[0].data
        self.un_cov = un_cov

        a = table["A"]
        b = table["B"]
        c = table["C"]
        self.a = a
        self.b = b
        self.c = c

        mask = 2*b>a+c
        self.mask = mask


        name = table["APOGEEID"]
        self.name = name



        SHIFT = table["VSHIFT"]
        self.shift = SHIFT


        VBARY = table["VBARY"]
        self.VBARY = VBARY


        teff = table["TEFF"]
        self.teff = teff

        logg = table["LOGG"]
        self.logg = logg


        feh = table["FEH"]
        self.feh = feh

        self.chi_inf = table["CHIINF"]
        self.chi_mix = table["CHIMIX"]

        self.BJD = table["BJD"]

        self.fiber = table["FIBER"]

        self.SNR =table["SNR"]



    def plot_teff_logg(self):

        # only show visits with 2b>a+c
        mask = self.mask

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]

        # shift is in km/s

        shift = self.shift[mask]

        a = self.a
        b = self.b
        c = self.c

        bac = (2*b-a-c)[mask]



        font = {'family': 'normal',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)

        alpha = 0.3
        #ax1

        mask2 = (shift>-1)&(shift<1)

        logg = logg[mask2]
        teff = teff[mask2]
        feh = feh[mask2]

        shift = shift[mask2]

        mask2 = (logg>-1)&(logg<6)&(teff>3000)&(teff<6000)&(feh>-3)&(feh<1)


        logg = logg[mask2]
        teff = teff[mask2]
        shift = shift[mask2]



        ax1.set_ylabel("$Teff \quad (K)$", fontsize=20)
        ax1.set_xlabel("$logg \quad (dex)$", fontsize=20)

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(logg,teff, marker='o', c=shift,
                    vmin=np.min(shift), vmax=np.max(shift), alpha=alpha, cmap="RdBu")


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])

        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$RV\quad shifts \quad (Km/s)$", fontsize=20)


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True



        f.suptitle(r'\textbf{Teff vs Logg for the stars in DR13}', fontsize=30,weight="bold")

        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "Teff_logg" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


    def plot_teff_feh(self):

        # only show visits with 2b>a+c
        mask = self.mask

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]

        shift = self.shift[mask]


        a = self.a
        b = self.b
        c = self.c

        bac = (2*b-a-c)[mask]





        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)

        alpha = 0.3

        mask2 = (shift > -1) & (shift < 1)

        logg = logg[mask2]
        teff = teff[mask2]
        feh = feh[mask2]
        shift = shift[mask2]

        mask2 = (logg > -1) & (logg < 6) & (teff > 3000) & (teff < 6000) & (feh > -3) & (feh < 1)


        feh = feh[mask2]
        teff = teff[mask2]
        shift = shift[mask2]

        #ax1

        ax1.set_ylabel("$Teff \quad (K)$", fontsize=20,weight="bold")
        ax1.set_xlabel("$FeH \quad (dex)$", fontsize=20,weight="bold")

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(feh,teff, marker='o', c=shift,
                    vmin=np.min(shift), vmax=np.max(shift), alpha=alpha, cmap="RdBu")


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$RV \quad shifts \quad (Km/s)$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{Teff vs FeH for the stars in DR13}', fontsize=30)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "Teff_feh" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()



    def plot_teff_logg_bac(self):

        # only show visits with 2b>a+c
        mask = self.mask

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]

        shift = self.shift[mask]

        a = self.a
        b = self.b
        c = self.c

        bac = (2*b-a-c)[mask]



        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3
        #ax1

        low = -1
        up = 3

        mask2 = ((bac>-1) & (bac<3))

        logg = logg[mask2]
        teff = teff[mask2]
        feh = feh[mask2]
        bac = bac[mask2]

        mask2 = (logg > -1) & (logg < 6) & (teff > 3000) & (teff < 6000) & (feh > -3) & (feh < 1)

        logg = logg[mask2]
        teff = teff[mask2]
        bac = bac[mask2]

        # plt.cm.viridis


        ax1.set_ylabel('$Teff \quad (K)$', fontsize=20,weight="bold")
        ax1.set_xlabel('$Logg \quad (dex)$', fontsize=20,weight="bold")

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(logg,teff, marker='o', c=bac,
                         vmin=low, vmax=up, alpha=alpha, cmap="viridis")


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$2b-a-c$ ", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{Teff vs Logg for the stars in DR13}', fontsize=30)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "Teff_logg_2bac" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


    def plot_teff_feh_bac(self):

        # only show visits with 2b>a+c
        mask = self.mask

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]

        shift = self.shift[mask]


        a = self.a
        b = self.b
        c = self.c

        bac = (2*b-a-c)[mask]


        low = -1
        up = 3

        mask2 = ((bac > -1) & (bac < 3))

        logg = logg[mask2]
        teff = teff[mask2]
        feh = feh[mask2]
        bac = bac[mask2]

        mask2 = (logg > -1) & (logg < 6) & (teff > 3000) & (teff < 6000) & (feh > -3) & (feh < 1)


        feh = feh[mask2]
        teff = teff[mask2]
        bac = bac[mask2]

        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3
        #ax1

        ax1.set_ylabel('$Teff \quad (K)$', fontsize=20,weight="bold")
        ax1.set_xlabel('$FeH \quad (dex)$', fontsize=20,weight="bold")

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(feh,teff, marker='o', c=bac,
                         vmin=low, vmax=up, alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$2b-a-c$", fontsize=20,weight="bold")



        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{Teff vs FeH for the stars in DR13}', fontsize=30)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "Teff_feh_2bac" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()






    def plot_shift_bjd(self):
        mask = self.mask

        shift =self.shift[mask]
        BJD = self.BJD[mask]
        feh = self.feh[mask]

        mask = (shift>-4)&(shift<4)

        shift = shift[mask]
        BJD =BJD[mask]
        feh = feh[mask]



        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3
        #ax1

        ax1.set_xlabel('$BJD$', fontsize=20,weight="bold")
        ax1.set_ylabel('$RV shift \quad (Km/s)$ ', fontsize=20,weight="bold")

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(BJD,shift, marker='o', c=feh,
                    vmin=np.min(feh), vmax=np.max(feh), alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$Fe/H \quad (dex)$", fontsize=20,weight="bold")



        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{RV shift vs BJD for the stars in DR13}', fontsize=30)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "RV_shift_vs_BJD" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


    def plot_rv_fiber(self):
        mask = self.mask

        a = self.a[mask]
        b = self.b[mask]
        c = self.c[mask]
        fiber = self.fiber[mask]
        SNR = self.SNR[mask]

        portion = (c+a)/(a+b+c)

        RV = (c - a) / (a + b + c) * 4.14468


        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3
        #ax1

        mask2 = ((SNR>30) & (SNR<300)&(RV>-4)&(RV<4))

        SNR = SNR[mask2]
        fiber = fiber[mask2]
        RV = RV[mask2]



        ax1.set_xlabel('$FiberID$', fontsize=20,weight="bold")
        ax1.set_ylabel('$RV shift \quad (Km/s)$', fontsize=20,weight="bold")

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(fiber,RV, marker='o', c=SNR,
                    vmin=30, vmax=300, alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$SNR$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{RV shifts vs FiberID for the stars in DR13}', fontsize=30)

        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "RV_shift_vs_Fiber" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()





    def plot_ac_fiber(self):
        mask = self.mask

        a = self.a[mask]
        b = self.b[mask]
        c = self.c[mask]
        fiber = self.fiber[mask]

        portion = (c+a)/(a+b+c)

        RV = (c - a) / (a + b + c) * 4.14468


        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        mask2 = (RV>-1)&(RV<1)

        fiber = fiber[mask2]
        portion = portion[mask2]
        RV = RV[mask2]



        alpha = 0.3
        #ax1

        ax1.set_xlabel('$FiberID$', fontsize=20,weight="bold")
        ax1.set_ylabel('$(c+a)/(a+b+c)$ ', fontsize=20,weight="bold")

        axes = plt.gca()

        axes.set_ylim([-1,1])

        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(fiber,portion, marker='o', c=RV,
                    vmin=np.min(RV), vmax=np.max(RV), alpha=alpha, cmap='RdBu')


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$RV \quad shifts \quad (Km/s)$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{(c+a)/(a+b+c) vs FiberID for the stars in DR13}', fontsize=30)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "ac_vs_Fiber" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()



    def histogram_shift_abc(self):


        a = self.a
        b = self.b
        c = self.c

        RV = (c-a)/(a+b+c)*4.14468

        # add a mask: only show results with 2b>a+c

        mask = 2*b>a+c

        a = a[mask]
        b = b[mask]
        c = c[mask]
        RV = RV[mask]



        font = {'weight': 'bold', 'size': 15}
        matplotlib.rc('font', **font)


        f, ((ax1, ax2), (ax3, ax4)) = \
            plt.subplots(2, 2)

        colors = ["cyan",'b', 'g', 'r']
        name = ["RV","a", "b", "c"]

        # histogram of rv
        #ax1

        rms_RV = (np.nansum(RV*RV)/len(RV))**0.5
        rms_a = (np.nansum(a * a) / len(a)) ** 0.5
        rms_b = (np.nansum(b*b) / len(b)) ** 0.5
        rms_c = (np.nansum(c * c) / len(c)) ** 0.5

        ax1.hist(RV, bins=40, color=colors[0], label="%s RMS = %.2f $Km/s$"%(name[0],rms_RV))


        #ax1.set_title('Histogram of Radial velocity shifts', fontsize=30)
        ax1.set_xlabel('$values \quad of\quad radial\quad velocity\quad shifts \quad (Km/s)$', fontsize=15,weight="bold")
        ax1.set_ylabel('$Number$', fontsize=15,weight="bold")
        ax1.legend(prop={'size': 15})

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)


        # histogram of a
        #ax2

        ax2.hist(a, bins=40, color=colors[1], label="$%s RMS = %.2f$"%(name[1],rms_a))


        #ax2.set_title('Histogram of parameter a', fontsize=30)
        ax2.set_xlabel('$values\quad of\quad parameter\quad a$', fontsize=15,weight="bold")
        ax2.set_ylabel('$Number$', fontsize=15,weight="bold")
        ax2.legend(prop={'size': 15})

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)


        # histogram of b
        #ax3

        ax3.hist(b, bins=40, color=colors[2], label="$%s RMS = %.2f$"%(name[2],rms_b))
        ax3.legend(prop={'size': 15})


        #ax3.set_title('Histogram of paramete b', fontsize=30)
        ax3.set_xlabel("$values\quad of\quad parameter\quad b$", fontsize=15,weight="bold")
        ax3.set_ylabel('$Number$', fontsize=15,weight="bold")

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)



        # histogram of c
        #ax4

        ax4.hist(c, bins=40, color=colors[3], label="$%s RMS = %.2f$"%(name[3],rms_c))
        ax4.legend(prop={'size': 15})


        #ax4.set_title('Histogram of parameter c', fontsize=30)
        ax4.set_xlabel("$values\quad of\quad parameter\quad c$", fontsize=15,weight="bold")
        ax4.set_ylabel('$Number$', fontsize=15,weight="bold")

        # add vertical grey line
        # ax1.plot((wl[index], wl[index]), (0.5, 1 + 0.5 * N), 'k-', linewidth=1.5)


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{Histogram of RV shifts, a, b and c for the stars in DR13}', fontsize=25)

        f.legends
        #f.suptitle("Histogram of RV shifts, a, b and c by using the absorption lines")


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 8.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "histogram_rv_shift" + ".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


    # RV before after


    def plot_RV_rms_before_after_teff(self):

        mask = self.mask

        shift =self.shift[mask]
        VBARY = self.VBARY[mask]

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]


        # From the average (c+a)/(a+b+c)
        # Do put a mask here
        mask = self.mask

        # add points with the same fiberid together
        name = self.name[mask]
        target = list(set(name))

        VBARY = self.VBARY[mask]
        shift =self.shift[mask]

        #SNR = self.SNR[mask]

        fusion_new = []

        # name+rms_old and rms_new + Teff logg feh

        for i in range(0,len(target)):

            print("Doing %.2f %%"%(i/len(target)*100))

            index = np.where(name == target[i])
            index = np.array(index)
            index = index.ravel()

            rms_old_i = np.std(VBARY[index])

            rms_new_i = np.std(VBARY[index]+shift[index])

            teff_i = np.nanmedian(teff[index])
            logg_i = np.nanmedian(logg[index])
            feh_i = np.nanmedian(feh[index])

            fusion_new.append([target[i],rms_old_i,rms_new_i,teff_i,logg_i,feh_i])

        fusion_new = np.array(fusion_new)



        self.fusion_new = fusion_new



        # portion+fiber+rv

        # name = fusion_new[:, 0]
        rms_old = np.array(fusion_new[:,1],dtype=np.float32).ravel()
        rms_new = np.array(fusion_new[:,2],dtype=np.float32).ravel()

        # use int
        teff = np.array(fusion_new[:,3],dtype=np.float16).ravel()


        ## Calculate good/bad rate:

        good = np.sum(np.array(rms_new<rms_old,dtype=int))
        bad = np.sum(np.array(rms_new > rms_old, dtype=int))
        tot = good+bad

        good = good/tot
        bad = bad/tot

        self.good = good

        self.bad = bad




        mask2 = (rms_old < 1)&(rms_new<1)

        rms_old = rms_old[mask2]
        rms_new = rms_new[mask2]
        teff = teff[mask2]




        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3

        low = 4500
        up =5300

        mask2 = (teff>4500)&(teff<5300)



        teff = teff[mask2]
        rms_old = rms_old[mask2]
        rms_new = rms_new[mask2]


        #ax1

        ax1.plot(rms_old,rms_old,"k",alpha=0.3,linewidth=0.3)

        ax1.set_xlabel('$RMS\quad of\quad RVs\quad before\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")
        ax1.set_ylabel('$RMS\quad of\quad RVs\quad after\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")

        ax1.plot([],[],"ro",label=" $%.2f\quad percent\quad of\quad the\quad stars\quad have\quad smaller\quad RV\quad RMS\quad after\quad the\quad correction$"%(self.good*100))

        ax1.legend()
        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(rms_old,rms_new, marker='o', c=teff,
                    vmin=low, vmax=up, alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$Teff \quad (K)$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{RMS of RVs before vs after the correction for the stars in DR13}', fontsize=25)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 11.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "RV_rms_before_after_teff" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()



    def plot_RV_rms_before_after_logg(self):

        mask = self.mask

        shift =self.shift[mask]
        VBARY = self.VBARY[mask]

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]


        fusion_new =self.fusion_new


        # name = fusion_new[:, 0]
        rms_old = np.array(fusion_new[:,1],dtype=np.float32).ravel()
        rms_new = np.array(fusion_new[:,2],dtype=np.float32).ravel()
        logg = np.array(fusion_new[:,4],dtype=np.float16).ravel()




        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)

        mask2 = ((logg>2) & (logg<3.25))

        logg = logg[mask2]
        rms_old = rms_old[mask2]
        rms_new = rms_new[mask2]

        mask2 = (rms_old < 1) & (rms_new < 1)

        rms_old = rms_old[mask2]
        rms_new = rms_new[mask2]
        logg = logg[mask2]

        alpha = 0.3
        #ax1

        ax1.plot(rms_old,rms_old, "k", alpha=0.3, linewidth=0.3)

        ax1.set_xlabel('$RMS\quad of\quad RVs\quad before\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")
        ax1.set_ylabel('$RMS\quad of\quad RVs\quad after\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")


        ax1.plot([],[],"ro",label=" $%.2f\quad percent\quad of\quad the\quad stars\quad have\quad smaller\quad RV\quad RMS\quad after\quad the\quad correction$"%(self.good*100))

        ax1.legend()
        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(rms_old,rms_new, marker='o', c=logg,
                    vmin=np.min(logg), vmax=np.max(logg), alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$logg \quad (dex)$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{RMS of RVs before vs after the correction for the stars in DR13}', fontsize=25)


        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 11.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "RV_rms_before_after_logg" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()


    def plot_RV_rms_before_after_feh(self):

        mask = self.mask

        shift =self.shift[mask]
        VBARY = self.VBARY[mask]

        teff = self.teff[mask]
        logg = self.logg[mask]
        feh = self.feh[mask]


        fusion_new =self.fusion_new


        # name = fusion_new[:, 0]
        rms_old = np.array(fusion_new[:,1],dtype=np.float32).ravel()
        rms_new = np.array(fusion_new[:,2],dtype=np.float32).ravel()
        feh = np.array(fusion_new[:,5],dtype=np.float16).ravel()

        # only show rms<1 km/s

        mask2 = (rms_old < 1) & (rms_new < 1)

        rms_old = rms_old[mask2]
        rms_new = rms_new[mask2]
        feh = feh[mask2]



        font = {'family': 'normal',
                'weight': 'bold',
                'size': 14}

        matplotlib.rc('font', **font)


        f, ax1 = plt.subplots(1,1)


        alpha = 0.3

        #ax1

        ax1.plot(rms_old,rms_old, "k", alpha=0.3, linewidth=0.3)

        ax1.set_xlabel('$RMS\quad of\quad RVs\quad before\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")
        ax1.set_ylabel('$RMS\quad of\quad RVs\quad after\quad the\quad correction \quad (km/s)$', fontsize=20,weight="bold")


        ax1.plot([],[],"ro",label=" $%.2f\quad percent\quad of\quad the\quad stars\quad have\quad smaller\quad RV\quad RMS\quad after\quad the\quad correction$"%(self.good*100))

        ax1.legend()


        f.subplots_adjust(right=0.8)


        pl = ax1.scatter(rms_old,rms_new, marker='o', c=feh,
                    vmin=np.min(feh), vmax=np.max(feh), alpha=alpha, cmap=plt.cm.viridis)


        cbar_ax = f.add_axes([0.85, 0.15, 0.02, 0.7])
        cb = f.colorbar(pl, cax=cbar_ax)

        cb.solids.set_edgecolor("face")


        cb.set_label("$FeH \quad (dex)$", fontsize=20,weight="bold")


        matplotlib.rcParams['text.usetex'] = True
        matplotlib.rcParams['text.latex.unicode'] = True

        f.suptitle(r'\textbf{RMS of RVs before vs after the correction for the stars in DR13}', fontsize=25)

        # save them:

        fig = matplotlib.pyplot.gcf()

        # adjust the size based on the number of visit

        fig.set_size_inches(14.5, 11.5)

        save_path = "/Users/caojunzhi/Downloads/upload_20170330_all/" + "RV_rms_before_after_feh" +".png"
        fig.savefig(save_path, dpi=500)

        plt.close()



model = plot()
model.read_table()



# Teff logg feh
model.plot_teff_logg()

model.plot_teff_feh()
model.plot_teff_logg_bac()
model.plot_teff_feh_bac()



# vs fiber+bjd



model.plot_shift_bjd()

model.plot_rv_fiber()
model.plot_ac_fiber()




# histogram
model.histogram_shift_abc()



#VBARY vs



model.plot_RV_rms_before_after_teff()
model.plot_RV_rms_before_after_logg()
model.plot_RV_rms_before_after_feh()




