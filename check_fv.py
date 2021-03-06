import numpy
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
from pylab import *

#################################
#Adjust physical properties here
#################################

boxsize = 0.1   # Mpc/h
N_fast = 10000  # particles
N_slow = 100    # particles
sigma_fast = 10.0 # km/s
sigma_slow = 2.0 # km/s
OmegaM = 1.0
#################################

Mpc2m = 3.08567758e22  #m
m2Mpc = 1./Mpc2m
Msun2kg = 1.98855e30 #kg
kg2Msun = 1./Msun2kg
m2km = 0.001
h_mass = 1.6737237e-27 #kg
G = 6.674e-11   # SI
H0 = 100.0        # km/s / (Mpc/h)
Msun2Gadget = 1e-10
G = G * (m2km**2.) * (m2Mpc) / (kg2Msun * Msun2Gadget) # (Mpc/h) (km/s)^2 / (1e10Msun/h)
print "G (internal unit) = ",G

k = 1.38064852e-23 #SI

rho_crit_0 = 3.* H0**2 / (8.*pi*G)  # (1e10 Msun/h)/(Mpc/h)^3
mass_p = rho_crit_0*OmegaM*boxsize**3/(N_fast+N_slow)
soft_r = 0.02*(mass_p/rho_crit_0)**(1./3.)


header_struct = numpy.dtype([
        ('npart', numpy.int32, 6),                #number of particles of each type in this file */
        ('mass', numpy.float64, 6), # mass of particles of each type. If 0, then the masses are explicitly
                                    #   stored in the mass-block of the snapshot file, otherwise they are omitted */
        ('time', numpy.float64, 1),                  # time of snapshot file */
        ('redshift', numpy.float64, 1), # redshift of snapshot file */
        ('flag_sfr', numpy.int32, 1),                 # flags whether the simulation was including star formation */
        ('flag_feedback', numpy.int32,1),           # flags whether feedback was included (obsolete) */
        ('npartTotal', numpy.uint32, 6),   # total number of particles of each type in this snapshot. This can be
                                #   different from npart if one is dealing with a multi-file snapshot. */
        ('flag_cooling', numpy.int32, 1),             # flags whether cooling was included  */
        ('num_files', numpy.int32, 1),                # number of files in multi-file snapshot */
        ('BoxSize', numpy.float64, 1),              # box-size of simulation in case periodic boundaries were used */
        ('Omega0', numpy.float64, 1),                # matter density in units of critical density */
        ('OmegaLambda', numpy.float64, 1),           # cosmological constant parameter */
        ('HubbleParam', numpy.float64, 1),           # Hubble parameter in units of 100 km/sec/Mpc */
        ('flag_stellarage', numpy.int32, 1),          # flags whether the file contains formation times of star particles */
        ('flag_metals', numpy.int32, 1),              # flags whether the file contains metallicity values for gas and star
                                #   particles */
        ('npartTotalHighWord', numpy.uint32, 6),   # High word of the total number of particles of each type */
        ('flag_entropy_instead_u', numpy.int32, 1),   # flags that IC-file contains entropy instead of u */
        ('flag_doubleprecision', numpy.int32, 1),     # flags that snapshot contains double-precision instead of single precision */
        ('flag_lpt_ics', numpy.int32, 1),             # flag to signal that IC file contains 2lpt initial conditions */
        ('lpt_scalingfactor', numpy.float32, 1),      # scaling factor for 2lpt initial conditions */
        ('flag_tracer_field', numpy.int32, 1),        # flags presence of a tracer field */
        ('composition_vector_length', numpy.int32,1),        # specifies the length of the composition vector (0 if not present)  */
        ('fill',numpy.string_,40)]) # fills to 256 Bytes */


def KrookWuModel(minv,maxv,stepv,tau,beta):
    v = numpy.arange(minv,maxv,stepv,dtype=numpy.float32)
    K = 1.0 - numpy.exp(tau/-6.0)
    f = numpy.exp(-0.5*v*v/(K*beta*beta)) * ((5.0*K-3.0)/K + (1.0-K)/(K*K)*v*v/(beta*beta))/(2.0*(2.0*numpy.pi*K*beta*beta)*numpy.sqrt(2.0*numpy.pi*K*beta*beta))
    return (v,f)

def read(filename):
    print "Opening ", filename
    fp = open(filename, "rb")
    # begin header
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    #print dummy
    header = numpy.fromfile(fp,dtype=header_struct,count=1)
    #print header
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    # end header

    # begin position
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    pos = numpy.fromfile(fp,dtype=numpy.float32,count = 3*header[0]['npart'][1]).reshape((header[0]['npart'][1],3))
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    #end position

    #begin velocity
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    vel = numpy.fromfile(fp,dtype=numpy.float32,count = 3*header[0]['npart'][1]).reshape((header[0]['npart'][1],3))
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    #end velocity

    vv = numpy.sqrt(vel[:,0]*vel[:,0]+vel[:,1]*vel[:,1]+vel[:,2]*vel[:,2])
    histogram = numpy.histogram(vv,bins=100)
    x = numpy.empty(100,dtype=numpy.float64)
    for i in range(100):
        x[i] = 0.5*(histogram[1][i]+histogram[1][i+1])
    y = numpy.array(histogram[0]).astype(numpy.float64)
    fig = figure()
    ax = fig.add_subplot(111)
    ax.plot(x,y*(histogram[1][1]-histogram[1][0]))
    ax.set_title(filename)
    fig.savefig(filename+".pdf")   # save the figure to file
    close(fig)    # close the figure

    #begin PID
    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    
    ID = numpy.fromfile(fp,dtype=numpy.int32,count = header[0]['npart'][1])

    dummy = numpy.fromfile(fp,dtype=numpy.int32,count=1)
    #end PID

    # check = 8060
    # index = numpy.where(ID == check)[0]
    # print ID[index],pos[index],vel[index,:]
    
    # #begin mass
    # dummy.tofile(fp)
    # data = mass_p*numpy.ones(N_slow+N_fast,dtype=numpy.float32)
    # data.tofile(fp)
    # dummy.tofile(fp)
    # #end mass
    
    fp.close()

def main():
    for i in range(20):
        filename = sys.argv[1].strip()+"%03d"%(i)
        if not os.path.isfile(filename):
            filename = sys.argv[1].strip()
        read(filename)
    #print KrookWuModel(0,100.0,0.01,6.0,7.0)
    return 0

if __name__ == "__main__":
    main()
