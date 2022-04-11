# -----------------------------------------------------------------------------
#
# CASE: 2D URANS equations - Flat Plate 
#       -- with kernels: rhs.py genRhs.py
#
# -----------------------------------------------------------------------------

# ======================================================================== LOAD
from tkinter.tix import Tree
import dnami as dn # dNami kernel

from dnami import np, sys # re-import non dnami modules
import os

# ===================================================================  FOLDERS

# -- restarts 
try :
    os.mkdir('./restarts/')
except FileExistsError:
    pass

# ========================================================================= ASK

# Parameters for the case ...
alpha = dn.cst(2.5)
Ma    = dn.cst(0.2)
Re    = dn.cst(5000000)
Pr    = dn.cst(0.71)
gamma = dn.cst(1.4)
T0    = dn.cst(300)
P0    = dn.cst(1.)
R     = dn.cst(287.05)
Rho0  = dn.cst((1+0.2*Ma**2)**(-2.5))
dimPcr = False

if dimPcr:
	Cv = dn.cst(1.0/(gamma-1.0))
else:
	Cv = dn.cst(1.0/(Ma**2*gamma*(gamma-1.0)))

filtr_amp = dn.cst(0.1)    # filter amplitude

# ... in time ...
with_dt = dn.cst(1e-3)
#nitmax  = 50000 # for actual run
nitmax  = 1000  # for test case 

# ... in space ...
#L = dn.cst(2.*np.pi) 
with_length = [2.3333333333,1]      # domain length in each direction
#with_grid   = [nx-2*hlo,64]   # number of points in each direction

# ... as fast as possible!
with_proc     = [2,2] # mpi proc. topology

# ===================================================================== PREPARE

dtree = dn.create_tree()

# .. assign user-defined values
dtree['eqns']['coeff'][0][1] = dn.cst(0.0000186)
dtree['eqns']['coeff'][1][1] = dn.cst(1.0/( (gamma-1.0)*Ma**2*Re*Pr ))
dtree['eqns']['coeff'][2][1] = dn.cst(gamma-1.)
dtree['eqns']['coeff'][3][1] = dn.cst(1.0/Cv)
dtree['eqns']['coeff'][4][1] = dn.cst(Ma*(gamma*R*T0)**0.5)
dtree['eqns']['coeff'][5][1] = dn.cst(0.1355)                 #Cb1
dtree['eqns']['coeff'][6][1] = dn.cst(0.622)                  #Cb2
dtree['eqns']['coeff'][7][1] = dn.cst(2.0/3)                  #sigma
dtree['eqns']['coeff'][8][1] = dn.cst(0.41)                   #k
dtree['eqns']['coeff'][9][1] = dn.cst(0.1355/0.41**2+(1+0.622)/2.0/3) #Cw1
dtree['eqns']['coeff'][10][1] = dn.cst(0.3)                   #Cw2
dtree['eqns']['coeff'][11][1] = dn.cst(2)                     #Cw3
dtree['eqns']['coeff'][12][1] = dn.cst(7.1)                   #Cv1
dtree['eqns']['coeff'][13][1] = dn.cst(1)                     #Ct1
dtree['eqns']['coeff'][14][1] = dn.cst(2)                     #Ct2
dtree['eqns']['coeff'][15][1] = dn.cst(1.1)                   #Ct3
dtree['eqns']['coeff'][16][1] = dn.cst(2)                     #Ct4
dtree['eqns']['coeff'][17][1] = dn.cst(1.0/2.0/3)             #sigmaI



# .. shortcut key handles
numerics = dtree['num']
grid     = dtree['grid']['size']
geom     = dtree['grid']['geom']
mpi      = dtree['mpi']['split']


#grid['nxgb'] = with_grid[0] 
#grid['nygb'] = with_grid[1]
#grid['nzgb'] = 1

geom['Lx'] = with_length[0] 
geom['Ly'] = with_length[1] 
geom['Lz'] = 0.

mpi['nxpr'] = with_proc[0] 
mpi['nypr'] = with_proc[1] 
mpi['nzpr'] = 1

# .. start the message passing interface
dtree = dn.start_mpi(dtree)
dMpi = dtree['mpi']['dMpi']

nx  = dMpi.nx
ny  = dMpi.ny

hlo  = numerics['hlo']
with_grid = [545-2*hlo,385-2*hlo]

grid['nxgb'] = with_grid[0] 
grid['nygb'] = with_grid[1]
grid['nzgb'] = 1

# .. create the *hlocomputational grid and write to file
dtree = dn.create_grid(dtree)
dn.dnami_io.hello_world(dtree)

# define useful aliases
xloc, yloc = geom['xloc'], geom['yloc']
Lx  , Ly   = geom['Lx']  , geom['Ly']
dx  , dy   = geom['dx']  , geom['dy']


numerics['tint']['tstep'] = with_dt 
dt = numerics['tint']['tstep']
numerics['filtr']['eps'] = filtr_amp 

# .. allocate tree
large = 10000 #no cache blocking in this example
dtree['libs']['cache blocking'] = [large,large,large]
dtree = dn.allocate(dtree)

# - Primitive variables
rh = dtree['eqns']['qvec']['views']['rho']
ux = dtree['eqns']['qvec']['views']['u']
uy = dtree['eqns']['qvec']['views']['v']
et = dtree['eqns']['qvec']['views']['et']

q  = dtree['eqns']['qvec']['views']['q'] 

# - Metrics
d  = dtree['eqns']['qvec']['views']['d']
ksi=dtree['eqns']['qvec']['views']['ksi']
eta=dtree['eqns']['qvec']['views']['eta']

with open('x_coord.dat','r') as f:
    dat=f.readlines()
    dat=np.array(dat)
    dat=dat.T
    for j in np.arange(0,ny+2*hlo):
        ksi[j,:]=dat
with open('y_coor.dat','r') as f:
    dat=f.readlines()
    dat=np.array(dat).T
    for i in np.arange(0,nx+2*hlo):
        eta[:,i]=dat
#u_wall=dtree['eqns']['qvec']['views']['d']

#deltay = dtree['eqns']['qvec']['views']['deltay']
deltaxI = dtree['eqns']['qvec']['views']['deltaxI']
deltayI = dtree['eqns']['qvec']['views']['deltayI']

# - Store variables aliases if any
if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	qstored  = dtree['eqns']['qvec']['views']['qstored'] 

# ================================================================== FUNCTIONS 

def sound_speed():
    e = et - .5*(ux*ux+uy*uy)
    T = (1./alpha)*( e*Ma*Ma )
    c = np.sqrt( T*(1.+1./alpha) )/Ma
    return c	

# ================================================================== INITIALISE

# initial clock
ti = dn.cst(0.0)
ni = 1
trstart = dn.cst(0.)

#init thermo
T0   = dn.cst(1.0)
P0   = dn.cst(1.0)/(Ma**2*gamma)
Rho0 = dn.cst(1.0)#P0/T0*Ma**2*gamma

#numpy slice refering to the core of the domain
dom = np.s_[hlo:nx+hlo,hlo:ny+hlo]

rh[dom] = Rho0
ux[dom] = (np.sin(xloc[:, np.newaxis, np.newaxis])
			*np.cos(yloc[np.newaxis, :, np.newaxis]))

# -- Swap 
dMpi.swap(q,hlo,dtree)

if 'qstored' in dtree['eqns']['qvec']['views'].keys():
	dn.dnamiF.stored(intparam,fltparam,data)	
	dMpi.swap(qstored,hlo,dtree)

# -- Write the first restart
dn.dnami_io.write_restart(0,ti,0,dtree)

# ========================================================================= RUN

intparam,fltparam,data = (dtree['libs']['fort']['integers'],
						  dtree['libs']['fort']['floats'],
						  dtree['libs']['fort']['data'])

mod_filter = 1
mod_output = 1000
mod_info   = 100


for n in range(1,nitmax+1):
    ti = ti + dt

    # - RK loop
    for nrk in range(1,4):
        intparam[7] = nrk
        dMpi.swap(	q,hlo,dtree)

        if 'qstored' in dtree['eqns']['qvec']['views'].keys():
        	dn.dnamiF.stored(intparam,fltparam,data)
        	dMpi.swap(	qstored,hlo,dtree)

        dn.dnamiF.time_march(intparam,fltparam,data)	

    # - Filter
    if np.mod(n,mod_filter) == 0:

        dMpi.swapXc(q,hlo,dtree)
        dn.dnamiF.filter(1,intparam,fltparam,data)
        dMpi.swapYc(q,hlo,dtree)
        dn.dnamiF.filter(2,intparam,fltparam,data)
        dMpi.swapZc(q,hlo,dtree)
        dn.dnamiF.filter(3,intparam,fltparam,data)

    # - Output restarts 
    if np.mod(n,mod_output) == 0:
            dn.dnami_io.write_restart(n,ti,0,dtree)

    # - Output information
    if np.mod(n,mod_info) == 0:

        if dMpi.ioproc:
            print('____________________________________________________________')
            print('iteration',n,' with time t =',ti)
            sys.stdout.flush()
        dn.dnami_io.globalMinMax(dtree,rh,'r')
        dn.dnami_io.globalMinMax(dtree,ux,'u')
        dn.dnami_io.globalMinMax(dtree,uy,'v')
        dn.dnami_io.globalMinMax(dtree,uz,'w')
        dn.dnami_io.globalMinMax(dtree,et,'et')

# ----------------------------------------------------------------------------

# -- Grab the max value of rho-rho0 at end of run

if dMpi.iMpi:        
    maxval = np.amax(rh[:]-Rho0)
    MPI    = dMpi.MPIlib
    erra   = dMpi.comm_torus.reduce(maxval,op=MPI.MAX,root=0)
    if dMpi.ioproc:
        np.savetxt('out.dat',np.asarray([erra]))
else:
    erra = np.amax(rh[:]-Rho0)
    np.savetxt('out.dat',np.asarray([erra]))

