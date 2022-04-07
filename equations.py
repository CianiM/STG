# =============================================================================
# 3D navier stokes equations  
# =============================================================================
from re import L
import sympy as sym
# number of dimensions
dim = 2 

# coefficients ////////////////////////////////////////////////////////////////
coefficients = {'visc'       : 1, 
                'kappa'      : 2,  
                'gamma_m1'   : 3,
                'CvInv'      : 4,
                'u_0'        : 5,
                'Cb1'        : 6,
                'Cb2'        : 7,
                'sigma'      : 8,
                'k'          : 9,
                'Cw1'        : 10,
                'Cw2'        : 11,
                'Cw3'        : 12,
                'Cv1'        : 13,
                'Ct1'        : 14,
                'Ct2'        : 15,
                'Ct3'        : 16,
                'Ct4'        : 17,
                'sigmaI'     : 18,
                'esse'       : 19,
                'L_ref'      : 20
                }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                  'u' : 2,  
                  'v' : 3, 
                  'et': 4,
                  'nut': 5}
                  
varsolved = ['rho','u','v','et','nut']

# of which
consvar      = [2,3,4,5] # are conservative variables         

# derived local variables /////////////////////////////////////////////////////

# -- Two divergence forms for inside and outside first derivative
divops  = ' deltaxI*( [u]_1x ) + deltayI*( [v]_1y )  '
ddivops = ' deltaxI*( {u}_1x ) + deltayI*( {v}_1y )  '

# -- Srain tensor
S =   { 'uu' : ' deltaxI*( [u]_1x )',
        'uv' : ' 0.5_wp*( deltayI*( [u]_1y ) + deltaxI*( [v]_1x ) ) ',
        'vu' : ' 0.5_wp*( deltaxI*( [v]_1x ) + deltayI*( [u]_1y ) ) ',
        'vv' : ' deltayI*( [v]_1y ) ' }
#s=''

#for i in S.keys():
#    s = s + '('+ S[i] +')**2 + '
s= '(( '+ '('+S['uu'] +')**2 +' + '('+ S ['vv']+')**2 +' +'('+S['uv'] + ')**2+ ' +'('+S['vu'] +')**2'+')*2)**0.5'
print(s)

varloc       = {'p': 'gamma_m1*rho*(e)',
                'e': '(et-0.5_wp*(u*u+v*v))',
                'chi': '( nut/visc ) ',
                'fv1': '( chi**3/( chi**3 + Cv1**3) )',
                'nu' : '( nut*fv1 ) ',
                'visc_t' : '( visc + rho*nu )',
                'fv2': '( 1-chi/( 1 + chi*fv1) )',
                'ft2': '( Ct3*exp(-Ct4*chi**2) )',
                'fw' : '( gg*( 1+Cw3**6 )/( gg**6+Cw3**6) ) ',
                'gg' : '( rr + Cw2*( rr**6-rr ) )',
                'rr' : '( nut/(SS*k**2*eta**2) )',
                'SS' : '( stemp+nut/(k**2*eta**2) )',
                #'stemp' : s, 
                'T': 'CvInv*(e)',
                #'symm':'( ( sign( 1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) ) )' ,
                'wall':' dabs( 1-symm )'}

varstored   = { 'd'  : {'symb' : "d" ,
                                 'ind': 1 ,
                                 'static' : True}, # normal distance
                'eta' : {'symb' : "eta" ,
                                  'ind' : 2 ,
                                  'static' : True},
                'ksi' : {'symb' : "ksi" ,
                                  'ind' :3,
                                  'static' : True},
                'stemp' : {'symb' : "stemp",
                                    'ind' : 4,
                                    'static' : True},
                'deltaxI' : {'symb' : '1.0_wp /( [ ksi ]_1x )' ,
                                      'ind' : 5,
                                      'static' : True},
                'deltayI' : {'symb' : '1.0_wp /(  [ eta ]_1y )' ,
                                      'ind' : 6 ,
                                      'static' : True},
                'symm'    : {'symb' : '( ( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) )',
                                      'ind': 7,
                                      'static': True}} 

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt'  
          ,  'u' : 'd(rho u)/dt',
             'v' : 'd(rho v)/dt', 
             'et': 'd(rho et)/dt',
             'nut':'d(rho nut)/dt'}
       

# .. name tags to use for intermediate variables created by the constructor
locname_dif = {'rho': 'dif_rho',
               'u'  : 'dif_rhou',
               'v'  : 'dif_rhov',
               'et' : 'dif_et  ',
               'nut': 'dif_nut'}

locname_conv = {'rho': 'conv_rho',
                'u'  : 'conv_rhou',
                'v'  : 'conv_rhov',
                'et' : 'conv_et  ',
                'nut': 'conv_nut '}

locname_bc  = {'rho': 'bc_rho',
               'u'  : 'bc_u',
               'v'  : 'bc_v',
               'et' : 'bc_et  ',
               'nut': 'bc_nut '}   
locname_bc_conv= { 'rho' : 'conv_rho',
                   'u'   : 'conv_u' ,
                   'v'   : 'conv_v',
                   'et'  : 'conv_et',
                   'nut' : 'conv_nut' }

locname_bc_dif = { 'rho' : 'dif_rho',
                   'u'   : 'dif_u',
                   'v'   : 'dif_v',
                   'et'  : 'dif_et',
                   'nut' : 'dif_nut' }             

# RHS terms ///////////////////////////////////////////////////////////////////

# Euler 

Fx = {'rho' : 'rho*u         ',
      'u'   : 'rho*u*u  + p  ',
      'v'   : 'rho*u*v       ',
      'et'  : '(rho*et + p)*u ',
      'nut' : 'rho*u*nut-sigmaI*(visc+rho*nut)*( { nut }_1x )*deltaxI ' }

Fy = {'rho' : 'rho*v         ',
      'u'   : 'rho*v*u       ', 
      'v'   : 'rho*v*v  + p  ', 
      'et'  : '(rho*et + p)*v ',
      'nut' : 'rho*v*nut-sigmaI*(visc+rho*nut)*deltayI*( { nut }_1y ) ' } 

Src_conv={}

for key in Fx.keys():
    Src_conv[key]= 'deltaxI*( [ ' +Fx[key] + ' ]_1x )' + ' + ' + 'deltayI*( [ '+Fy[key]+ ' ]_1y ) '

 
######################################
#                                    #
# Navier-Stokes Diffusive terms only # 
#                                    #
######################################

Fx = {'u'   : ' - visc_t *( 2.0_wp * deltaxI*( {u}_1x ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )',
      'v'   : ' - visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x ) )', 
      
      'et'  : ' - kappa*deltaxI*( {T}_1x ) '
              ' - u*(visc_t *( 2.0_wp *deltaxI*( {u}_1x ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'
              ' - v*(visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x )))'}

Fy = {'u'   : ' - visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x ))  ',
      'v'   : ' - visc_t *( 2.0_wp * deltayI*( {v}_1y ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )', 
      'et'  : ' - kappa*deltayI*( {T}_1y )'
              ' - u*(visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x )))'
              ' - v*(visc_t *( 2.0_wp * deltayI*( {v}_1y ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'}
       
# -- Divergence formulation

Src_dif  = {}

for key in Fx.keys():
    Src_dif[key]= 'deltaxI*( [ ' + Fx[key] +' ]_1x )' + ' + ' + 'deltayI *( [ '+ Fy[key]  +' ]_1y ) '

Src_dif['nut'] = ' -Cb2*sigmaI*( (deltaxI)**2*( [ rho*nut ]_1x )*( [ nut ]_1x )+ (deltayI)**2*( [ rho*nut ]_1y )*( [ nut ]_1y ) ) \
                  - Cb1*(1-ft2)*SS*rho*nut + (Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2 '


########################################################################################################## 
#----------------------------------------BOUNDARY CONDITIONS---------------------------------------------#
##########################################################################################################


######################################
#                                    #
#----------Symmetry, j1--------------#  
#                                    #
######################################

Src_BC_Symm_conv ={}

Fx = {'rho' : ' rho*u ',
       'u'  : ' rho*u*u + p ',
       'v'  : ' rho*u*v ',
       'et' : ' (rho*et+ p)*u ',
       'nut': ' rho*u*nut-sigmaI*(visc+rho*nut)*( { nut }_1x )*deltaxI'}


for key in Fx.keys():
      Src_BC_Symm_conv[key] = 'deltaxI*( [ ' + Fx[key] + ' ]_1x )'

Src_BC_Symm_dif={}

Fx = { 'u'  : ' -4.0_wp/3.0_wp*visc_t*( {u}_1x )*deltaxI',
       'v'  : ' -visc_t*( {v}_1x )*deltaxI ',
       'et' : ' -4.0_wp/3.0_wp*visc_t*( {u}_1x )*deltaxI*u -visc_t*( {v}_1x )*deltaxI*v -kappa*( {T}_1x )*deltaxI'
      }
for key in Fx.keys():
      Src_BC_Symm_dif[key] = 'deltaxI*( [ ' + Fx[key] +' ]_1x )'

Src_BC_Symm_dif['nut']=' -Cb2*sigmaI*( (deltaxI)**2*( [ rho*nut ]_1x )*( [ nut ]_1x ) ) - Cb1*(1-ft2)*SS*rho*nut + (Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2 '

######################################
#                                    #
#-------------Wall, j1---------------#  
#                                    #
######################################
Src_BC_Wall={}

Src_BC_Wall_conv = { 'rho' : ' rho*( [v]_1y )*deltayI',
                     'et'  : ' ( rho*et+p )* ( [v]_1y )*deltayI'}

Src_BC_Wall_dif = { 'et'  : ' (-visc_t*( [u]_1y )*( [u]_1y ) - 4.0_wp/3.0_wp*visc_t*( [v]_1y )*( [v]_1y ) -kappa*[T]_2y )*deltayI**2 -kappa*[T]_2x'}


##--Building boundary conditions--##
Src_BC_conv={}
Src_BC_dif={}
##--Boundary conditions for j1--
Src_BC_conv['j1'] = {}
Src_BC_dif['j1'] = {}

for key in Src_BC_Symm_conv.keys():
      Src_BC_conv['j1'][key]='( '+ Src_BC_Symm_conv[key] + ' )*symm'

for key in Src_BC_Symm_dif.keys():
      Src_BC_dif['j1'][key]= '( '+ Src_BC_Symm_dif[key] + ' )*symm'

for key in Src_BC_Wall_conv.keys():
      Src_BC_conv['j1'][key] = Src_BC_conv['j1'][key] + ' + ( '+ Src_BC_Wall_conv[key]+' )*wall'

Src_BC_dif['j1']['et']= Src_BC_dif['j1']['et'] + ' + ( '+Src_BC_Wall_dif['et']+' )*wall'

##--Physical BC
Src_phyBC = {}
Src_phyBC['j1'] = { 'u' : ' symm*u ',
                    'v' : '0.0_wp',
                    'nut':' summ*nut' }

##--Symmetric boundary conditions for jmax--##

#...
#Src_BC_conv['jmax']={}
#Src_BC_dif['jmax']={}
#for key in Src_BC_Symm_conv.keys():
#      Src_BC_conv['jmax'][key]= Src_BC_Symm_conv[key]
#for key in Src_BC_Symm_dif.keys():
#      Src_BC_dif['jmax'][key]= Src_BC_Symm_dif[key]
#...




######################################
#                                    #
#----Inlet i1, Subsonic inflow-------#  
#                                    #
######################################


##-- Charactertistics --##
from CharsForConsLaw import characteristics
from genNSBC import sympy2dNami

Char={}
Char=characteristics('Euler')

x,y,t=sym.symbols(['x','y','t'],Real=True)
rho,u,v,et,p=sym.symbols(['rho','u','v','et','p'],Real=True)
rhou ,rhov ,rhow ,rhoet = sym.symbols(['rhou' ,'rhov' ,'rhow' ,'rhoet'],Real=True)
gamma=sym.Symbol('gamma')
rho = sym.Function('rho')(x,y,t)
u = sym.Function('u')(x,y,t)
v = sym.Function('v')(x,y,t)
p = sym.Function('p')(x,y,t)

p=(gamma-1)*(rho*et-(rho*u)**2/(2*rho))

Q=sym.Matrix([[rho],
             [u],
             [v],
             [et]])

Q_CS=sym.Matrix([[rho],
                [rho*u],
                [rho*v],
                [rho*et]])
M=Q_CS.jacobian(Q)

Li_BC_i1_in = Char['xi'][0].copy()
Li_BC_i1=[]
for i in Li_BC_i1_in:
      Li_BC_i1.append(sympy2dNami(i))
#--> Li_BC_imax[0] = L3 <--#
#--> Li_BC_imax[1] = L2 <--#
#--> Li_BC_imax[2] = L1 <--#
#--> Li_BC_imax[3] = L5 <--#
#--------------------------#
#Li_BC_i1[2] = Li_BC_i1[3]        #--->Costant U-velocity at inlet
Li_BC_i1[3] = '-'+Li_BC_i1[2]        #--->Constant Pressure at inlet  --> L5=-L1
Li_BC_i1[0] = ' 0.0_wp '          #--->Costant V-velocity at inlet
Li_BC_i1[1] = ' 0.0_wp '          #--->Costant T          at inlet
#Li_BC_i1[1] = ' 0.0_wp '         #--->Costant entropy at inlet

##--Physical boundary conditions--##

Src_phyBC['i1'] = { 'rho' : 'Rho0',
                    'v' : '0.0_wp',
                    #'T' : 'T0',  
                    'nut': '0.0_wp'}

D = { 1 : '(1.0_wp/c**2*('+ Li_BC_i1[1]+'0.5_wp*('+Li_BC_i1[3] + ' + ' + Li_BC_i1[2]+'))', #0.0_wp
      2 : '(0.5_wp*( '+Li_BC_i1[3] + ' + ' + Li_BC_i1[2]+'))',  #0.0_wp
      3 : '(1.0_wp/(2*rho*c)*( '+ Li_BC_i1[3] + ' - ' + Li_BC_i1[2]+'))',
      4 : Li_BC_i1[0]}

Src_BC_conv_i1 = {   'u'   : 'u*'+D[1]+'+ rho*'+D[3]+'rho*u*( [ v ]_1y )',
                     'et'  : '0.5_wp*(u**2)*' + D[1]+'+1.0_wp/gamma_m1*'+D[2]+'+rho*u*'+D[3]+'+rho*v*'+D[4]+' [ (rho*et+p)*v ]_1y'}


Src_BC_dif_i1 = {}

Fx = {'u'   : ' - visc_t *( 2.0_wp * deltaxI*( {u}_1x ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  ) )',
      'et'  : ' - kappa*deltaxI*( {T}_1x ) '
              ' - u*(visc_t *( 2.0_wp *deltaxI*( {u}_1x ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'
              ' - v*(visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x )))'}

Fy = {'u'   : ' - visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x ))  ',
      'et'  : ' - kappa*deltayI*( {T}_1y )'
              ' - u*(visc_t *( deltayI*( {u}_1y ) + deltaxI*( {v}_1x )))'
              ' - v*(visc_t *( 2.0_wp * deltayI*( {v}_1y ) - 2.0_wp/3.0_wp * ( '+ ddivops +'  )))'}

for key in Fx.keys():
    Src_BC_dif_i1[key]= 'deltaxI*( [ ' + Fx[key] +' ]_1x )' + ' + ' + 'deltayI *( [ '+ Fy[key]  +' ]_1y ) '

######################################
#                                    #
#---Outflow imax, Costant pressure---#  
#                                    #
######################################

Li_BC_imax_out     = Char['xi'][0].copy()
velchar_BC_imaxout  = Char['xi'][1].copy()

Li_BC_imax=[]
for i in Li_BC_imax_out:
      Li_BC_imax.append(sympy2dNami(i))
velchar_BC_imax= []
for i in velchar_BC_imaxout:
      velchar_BC_imax.append(sympy2dNami(i))

#--> Li_BC_imax[0] = L3 <--#
#--> Li_BC_imax[1] = L2 <--#
#--> Li_BC_imax[2] = L1 <--#
#--> Li_BC_imax[3] = L5 <--#

Li_BC_imax[2]= '-'+Li_BC_imax[3]  #--->Costant Pressure at the outlet

D = { 1 : '(1.0_wp/c**2*('+ Li_BC_imax[1]+'0.5_wp*('+Li_BC_imax[3] + ' + ' + Li_BC_imax[2]+'))',
      2 : '(0.5_wp*( '+Li_BC_imax[3] + ' + ' + Li_BC_imax[2]+'))',  #0.0_wp
      3 : '(1.0_wp/(2*rho*c)*( '+ Li_BC_imax[3] + ' - ' + Li_BC_imax[2]+'))',
      4 : Li_BC_imax[0]}

Src_BC_conv_imax = { 'rho' : Src_conv['rho']+'+'+ D[1],
                     'u'   : 'u*'+D[1]+'+ rho*'+D[3]+'+ [ rho*u*v ]_1y',
                     'v'   : 'v*'+D[1]+'+ rho*'+D[4]+'+ [ rho*v*v ]_1y',
                     'et'  : '0.5_wp*(u**2+v**2)*'+D[1]+'+1.0_wp/gamma_m1*'+D[2]+'+rho*u*'+D[3]+'+rho*v*'+D[4]+' [ (rho*et+p)*v ]_1y'}

Src_BC_dif_imax = {}
for key in Src_dif.keys():
      Src_BC_dif_imax[key]= Src_dif[key]

######################################
#                                    #
#----Outflow jmax, Non reflective----#  
#                                    #
######################################

Mi_BC_jmax_sym     = Char['eta'][0].copy()
velchar_BC_jmaxsym = Char['eta'][1].copy()

Mi_BC_jmax=[]
for i in Mi_BC_jmax_sym:
      Mi_BC_jmax.append(sympy2dNami(i))
velchar_BC_jmax= []
for i in velchar_BC_jmaxsym:
      velchar_BC_jmax.append(sympy2dNami(i))
#Mi_BC_jmax[2]= -Mi_BC_jmax[3]  #--->Costant Pressur at the outlet
Mi_BC_jmax[2] ='esse*c*(1-M_jmax*M_jmax)/L_ref*( p - P0)'  #---> Non reflective condition

D = { 1 : '(1.0_wp/c**2*('+ Mi_BC_jmax[1]+'0.5_wp*('+Mi_BC_jmax[3] + ' + ' + Mi_BC_jmax[2]+'))',
      2 : '(0.5_wp*( '+Mi_BC_jmax[3] + ' + ' + Mi_BC_jmax[2]+'))',
      3 : '(1.0_wp/(2*rho*c)*( '+ Mi_BC_jmax[3] + ' - ' + Mi_BC_jmax[2]+'))',
      4 : Mi_BC_jmax[0]}

Src_BC_conv_jmax = { 'rho' : Src_conv['rho']+'+'+ D[1],
                     'u'   : 'u*'+D[1]+'+ rho*'+D[3]+'[ rho*u*v ]_1y',
                     'v'   : 'v*'+D[1]+'+ rho*'+D[4]+'[ rho*v*v ]_1y',
                     'et'  : '0.5_wp*(u**2+v**2)*'+D[1]+'+1.0_wp/gamma_m1*'+D[2]+'+rho*u*'+D[3]+'+rho*v*'+D[4]+' [ (rho*et+p)*v ]_1y'}

Src_BC_dif_jmax = {}
for key in Src_dif.keys():
      Src_BC_dif_jmax[key]= Src_dif[key]


varbc = { 'u_wall' : {'symb'  : ' u_wall ', 
                      'ind'   : 1 ,
                      'static': False,
                      'face'  : 'j1'},
          'v_wall' : {'symb'  : 'v_wall ',
                      'ind'   : 2 ,
                      'static': False,
                      'face'  : 'j1'},
          'nut_wall':{'symb'  : 'nut_wall ',
                      'ind'   : 3 ,
                      'static': False,
                      'face'  : 'j1'}}           