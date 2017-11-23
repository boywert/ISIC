import numpy
from pylab import *

#################################
#Adjust physical properties here
#################################

boxsize = 100.0   # Mpc/h
N_fast = 1000000  # particles
N_slow = 10000    # particles
sigma_fast = 10.0 # km/s
sigma_slow = 2.0 # km/s
Mhalo = 1.0       # 1e10 Msun/h

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
Mhalo_added = Mhalo*N_slow/N_fast
R_fast = 0.5*G*Mhalo/(sigma_fast*sigma_fast)
R_slow = 0.5*G*Mhalo_added/(sigma_slow*sigma_slow)
mass_p = Mhalo/N_fast
rho_crit_0 = 3.* H0**2 / (8.*pi*G)  # (1e10 Msun/h)/(Mpc/h)^3
soft_r = 0.02*(mass_p/rho_crit_0)**(1./3.)
print "R_fast = ", R_fast
print "R_slow = ", R_slow
print "mass_p = ", mass_p
print "r_soft = ", soft_r

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

header = numpy.zeros(1,dtype=header_struct)
pointer_hdr = header[0]
pointer_hdr['npart'][1] = N_fast+N_slow
pointer_hdr['mass'][1] = mass_p
pointer_hdr['time'] = 0.0
pointer_hdr['redshift'] = 10000000.0
pointer_hdr['npartTotal'][1] = N_fast+N_slow
pointer_hdr['num_files'] = 1
pointer_hdr['BoxSize'] = boxsize
pointer_hdr['Omega0'] = 1.0
pointer_hdr['OmegaLambda'] = 0.0
pointer_hdr['HubbleParam'] = 0.0

dummy = numpy.int32(256)
fp = open("test.bin", "wb")

# begin header
dummy.tofile(fp)
header[0].tofile(fp)
dummy.tofile(fp)
# end header

# begin position
dummy.tofile(fp)
phi = numpy.random.uniform(0,2*numpy.pi,N_fast)
costheta = numpy.random.uniform(-1.0,1.0,N_fast)
u = numpy.random.uniform(0,1,N_fast)
theta = numpy.arccos( costheta )
r = R_fast * u

data = (r * sin( theta) * cos( phi ) + boxsize/2).astype(numpy.float32)
#print data
data.tofile(fp)
fp.seek(N_slow*4,1)
data = (r * sin( theta) * sin( phi ) + boxsize/2).astype(numpy.float32)
#print data
data.tofile(fp)
fp.seek(N_slow*4,1)
data = (r * cos( theta ) + boxsize/2).astype(numpy.float32)
#print data
data.tofile(fp)

phi = numpy.random.uniform(0,2*numpy.pi,N_slow)
costheta = numpy.random.uniform(-1.0,1.0,N_slow)
u = numpy.random.uniform(0,1,N_slow)
theta = numpy.arccos( costheta )
r = R_slow * u

data = (r * cos( theta ) + boxsize/2).astype(numpy.float32)
data.tofile(fp)
fp.seek(-4*(N_fast+2*N_slow),1)
data = (r * sin( theta) * sin( phi ) + boxsize/2).astype(numpy.float32)
#print data
data.tofile(fp)
fp.seek(-4*(N_fast+2*N_slow),1)
data = (r * sin( theta) * cos( phi ) + boxsize/2).astype(numpy.float32)
#print data
data.tofile(fp)
fp.seek(4*2*(N_fast+N_slow),1)
#print data
dummy.tofile(fp)
del(phi,costheta,u,theta,data)

#end position

#begin velocity
dummy.tofile(fp)

data = numpy.random.normal(0., sigma_fast, N_fast).astype(numpy.float32)
vx = data
print data
data.tofile(fp)
data = numpy.random.normal(0., sigma_slow, N_slow).astype(numpy.float32)
print data
data.tofile(fp)
data = numpy.random.normal(0., sigma_fast, N_fast).astype(numpy.float32)
vy = data
print data
data.tofile(fp)
data = numpy.random.normal(0., sigma_slow, N_slow).astype(numpy.float32)
print data
data.tofile(fp)
data = numpy.random.normal(0., sigma_fast, N_fast).astype(numpy.float32)
vz = data
print data
data.tofile(fp)
data = numpy.random.normal(0., sigma_slow, N_slow).astype(numpy.float32)
print data
data.tofile(fp)
dummy.tofile(fp)
#end velocity
#vv = numpy.sqrt(vx*vx+vy*vy+vz*vz)
#print vv
#begin PID
dummy.tofile(fp)
data = numpy.arange(N_slow+N_fast,dtype=numpy.int32)
data.tofile(fp)
dummy.tofile(fp)
#end PID

fp.close()

