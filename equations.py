# =============================================================================
# 3D navier stokes equations  
# =============================================================================

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
                'sigmaI'     : 18
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
                'symm':' ( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) ) ' ,
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
                'symm'    : {'symb' : '( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) )',
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

 

# -- Conservative formulation

#Src_conv = {'rho' : '[ '+Fx['rho']+' ]_1x' + ' + ' + '[ '+Fy['rho']+' ]_1y' + ' + ' + '[ '+Fz['rho']+' ]_1z ',
#            'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y' + ' + ' + '[ '+Fz['u']  +' ]_1z ',
#            'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y' + ' + ' + '[ '+Fz['v']  +' ]_1z ',
#            'w'   : '[ '+Fx['w']  +' ]_1x' + ' + ' + '[ '+Fy['w']  +' ]_1y' + ' + ' + '[ '+Fz['w']  +' ]_1z ',
#            'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y' + ' + ' + '[ '+Fz['et'] +' ]_1z ',
#            'nut': ' [ '+Fx['nut']+' ]_1x' + ' + ' + '[ '+Fy['nut']+' ]1_y' + ' + ' + '[ '+Fz['nut']+' ]1_z '}

# -- Skew symmetric formulation

#Src_skew = {'rho' : '0.5_wp*( [rho*u]_1x + [rho*v]_1y + [rho*w]_1z '
#                   '+   u*[rho]_1x + v*[rho]_1y + w*[rho]_1z '
#                   '+ rho*( '+divops+' ) )',
#        'u'   : '0.5_wp*( [rho*u*u]_1x + [rho*u*v]_1y + [rho*u*w]_1z '
#                   '+   u*[rho*u]_1x + v*[rho*u]_1y + w*[rho*u]_1z '
#                   '+rho*u*( '+divops+' ) ) + [p]_1x ',
#        'v'   : '0.5_wp*( [rho*u*v]_1x + [rho*v*v]_1y + [rho*v*w]_1z '
#        '+   u*[rho*v]_1x + v*[rho*v]_1y + w*[rho*v]_1z '
#                   '+rho*v*( '+divops+' ) ) + [p]_1y ',
#        'w'   : '0.5_wp*( [rho*w*u]_1x + [rho*w*v]_1y + [rho*w*w]_1z '
#                   '+   u*[rho*w]_1x + v*[rho*w]_1y + w*[rho*w]_1z '
#                   '+rho*w*( '+divops+' ) ) + [p]_1z ',
#        'et'  : '0.5_wp*( [rho*et*u]_1x + [rho*et*v]_1y + [rho*et*w]_1z '
#                   '+   u*[rho*et]_1x + v*[rho*et]_1y + w*[rho*et]_1z '
#                   '+rho*et*( '+divops+') ) + [p*u]_1x + [p*v]_1y + [p*w]_1z '}

######################################
##Navier-Stokes Diffusive terms only## 
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
    Src_dif[key]= 'deltaxI*( [ ' + Fx[key] +' ]_1x )' + ' + ' + 'deltayI *( [ '+ Fy[key]  +' ]_1y )'

#Src_dif  = {'u'   : '[ '+Fx['u']  +' ]_1x' + ' + ' + '[ '+Fy['u']  +' ]_1y' + ' + ' + '[ '+Fz['u']  +' ]_1z ',
#            'v'   : '[ '+Fx['v']  +' ]_1x' + ' + ' + '[ '+Fy['v']  +' ]_1y' + ' + ' + '[ '+Fz['v']  +' ]_1z ',
#            'et' : ' [ '+Fx['et'] +' ]_1x' + ' + ' + '[ '+Fy['et'] +' ]_1y' + ' + ' + '[ '+Fz['et'] +' ]_1z '}



Src_dif['nut'] = ' -Cb2*sigmaI*( (deltaxI)**2*( [ rho*nut ]_1x )*( [ nut ]_1x ) + (deltayI)**2*( [ rho*nut ]_1y )*( [ nut ]_1y ) ) - Cb1*(1-ft2)*SS*rho*nut + (Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2 '
#Src_SA = {'nut':' -Cb2*sigmaI*( ( [ rho*nut ]_1x )*( [ nut ]_1x )+ deltayI**2*( [ rho*nut ]_1y )*( [ nut ]_1y ) ) - Cb1*(1-ft2)*( SS )*rho*nut + (Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/d**2 '}


########################################################################################################## 
#########################################BOUNDARY CONDITIONS##############################################
##########################################################################################################

##Characteristic equations##

from CharsForConsLaw import characteristics

Char={}



##--Symm--##

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

##--Wall--##
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
      Src_BC_conv['j1'][key]= '( '+Src_BC_Symm_conv[key]+' )*symm '+ '( ' + Src_BC_Wall_conv[key]+' )*wall'
      Src_BC_dif['j1'][key] = '( '+Src_BC_Symm_dif[key] +' )*symm '+ '( ' + Src_BC_Wall_dif[key] +' )*wall'

##--Boundary conditions for jmax--##


######
Src_BC_conv['jmax']={}
Src_BC_dif['jmax']={}
for key in Src_BC_Symm_conv.keys():
      Src_BC_conv['jmax'][key]= Src_BC_Symm_conv[key]
      Src_BC_dif['jmax'][key]= Src_BC_Symm_dif[key]

Src_BC ={}
Src_BC['j1'] = { 'u' : ' symm*u ',
                 'v' : '0.0_wp',
                 'nut':' summ*nut' }
##--Boundary conditions for i1--##




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