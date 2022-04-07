import sympy as sym
import sys

sym.init_printing(use_latex=True,wrap_line=False)

def characteristics(eqn):

    #################################################################################
    #
    #    Some automative characteristic derivations for systems of conservation laws
    #
    #################################################################################
    
    #####################################
    #  Deriving the quasi-linear form
    #####################################
    
    # Example of the Euler's equations in cartesian coordinates
    
    # General variables for the Euler equations:
    
    x,y,z,t     = sym.symbols(['x','y','z','t']   ,Real=True)
    xi,eta,zeta  = sym.symbols(['xi','eta','zeta'],Real=True)    

    rho,u  ,v  ,w  ,et,Bx,By,Bz     = sym.symbols(['rho','u'  ,'v'  ,'w'  ,'et' , 'Bx','B_y','B_z'],Real=True)
    mu0 = sym.Symbol('mu_0',Real=True)
    mu0 = 1.0

    vel = sym.Matrix([u,v,w])
    B = sym.Matrix([Bx,By,Bz])
    Bnorm = B.dot(B)
    sym.pprint(mu0)

    rhou ,rhov ,rhow ,rhoet = sym.symbols(['rhou' ,'rhov' ,'rhow' ,'rhoet'],Real=True)
    x_xi,y_xi,z_xi          = sym.symbols(['x_xi','y_xi','z_xi'],Real=True)    
    x_eta,y_eta,z_eta       = sym.symbols(['x_eta','y_eta','z_eta'],Real=True)  
    x_zeta,y_zeta,z_zeta    = sym.symbols(['x_zeta','y_zeta','z_zeta'],Real=True)  
    J = sym.Symbol('J',real=True)

    # Equation of states:
    #IG EOS
    
    gamma = sym.Symbol('gamma')
    p     = (gamma-1)*(rho*et-(rhou**2+rhov**2)/(2*rho))
    
    #Arbitrary EOS:
    # p   = sym.Function('p')(rho*et-0.5*(rhou**2+rhov**2)/rho,rho)
    # This general definition will introduce the Grüneisen index in the following derivations.
    
    # Defining the flux vector "Flx":
    #                  ---  
    #   ∂   ->    -->  ---     ->
    #   ──( f ) + Div( Flx ) = 0
    #   ∂t         
    # and the quasi-linear form and the "A matrix, A = [Ax,Ay,Az]""
    #
    #             -----------      
    #   ∂   ->    -----------  ->    ->
    #   ──( f ) + (A. \nabla)( f ) = 0
    #   ∂t         
    #
    

    # Flx_x for cartesian coordinates
    if eqn == 'Euler':
        Flx_x   = sym.Matrix([[rho*u],
                             [rho*u*u+p],
                             [rho*v*u],
                             [(rho*et+p)*u]])
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
        Flx_y   = sym.Matrix([[rho*v],
                             [rho*u*v],
                             [rho*v*v+p],
                             [(rho*et+p)*v]])
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])    
        # Flx_x for curvilinear coordinates    
   
        # Flx_x   = sym.Matrix([[rho*        ( u*x_xi + v*x_eta  )/J         ],
        #                       [(rho*u*     ( u*x_xi + v*x_eta  )+p*x_xi)/J ],
        #                       [(rho*v*     ( u*x_xi + v*x_eta  )+p*x_eta)/J],
        #                       [((rho*et+p)*( u*x_xi + v*x_eta  ))/J       ]])    
        #                      # [rho*w],
        #                      # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])    
        
        # Flx_y   = sym.Matrix([[rho*        ( u*y_xi+v*y_eta  )/J          ],
        #                       [(rho*u*     ( u*y_xi+v*y_eta  ) +p*y_xi )/J],
        #                       [(rho*v*     ( u*y_xi+v*y_eta  ) +p*y_eta)/J],
        #                       [((rho*et+p)*( u*y_xi+v*y_eta  ))/J        ]])    
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]]) 
    
     #   Flx_z   = sym.Matrix([[rho* ( u*x_zeta + v*y_zeta  )],
     #                        [rho*u*( u*x_zeta + v*y_zeta  )+p*x_zeta],
     #                        [rho*v*( u*x_zeta + v*y_zeta  )+p*y_zeta],
     #                        [(rho*et+p)*( u*x_zeta + v*y_zeta  )]])    
     #                        # [rho*w],
     #                        # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])    

        Q    = sym.Matrix([[rho],
                       [rhou],
                       [rhov],
                       # [rhow],
                       [rhoet]])

        p = sym.Symbol('p')
    
        Q_CS = sym.Matrix([[rho],
                       [rho*u],
                       [rho*v],
                       # [rho*w],
                       [0.5*rho*(u*u+v*v)+1./(gamma-1.)*p]])
                       # [rho*w],
                       # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
    
        Q2    = sym.Matrix([[rho],
                       [u],
                       [v],
                       # [w],
                       [p]])

    elif eqn == 'MHD':
        Flx_x   = sym.Matrix([[rho*u],
                             [rho*u*u+p + Bnorm/(2*mu0)-Bx**2/mu0],
                             [rho*v*u - Bx*By/mu0],
                             [rho*w*u - Bx*Bz/mu0],
                             [By*u-Bx*v],
                             [Bz*u-Bx*w],
                             [(rho*et+p + Bnorm/(2*mu0))*u + Bx*B.dot(vel)]])
        sym.pprint(Flx_x)
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
        # Flx_y   = sym.Matrix([[rho*v],
        #                      [rho*u*v - By*Bx/mu0],
        #                      [rho*v*v+p + Bnorm/(2*mu0)-By**2/mu0],
        #                      [rho*w*v - By*Bz/mu0],
        #                      [By*u-Bx*v],
        #                      [Bz*u-Bx*w],
        #                      [(rho*et+p + Bnorm/(2*mu0) )*v + Bx*B.dot(vel) ]])
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])  
        Q    = sym.Matrix([[rho],
                       [rhou],
                       [rhov],
                       [rhow],
                       [By],
                       [Bz],
                       [rhoet]])

        p = sym.Symbol('p')
    
        Q_CS = sym.Matrix([[rho],
                       [rho*u],
                       [rho*v],
                       [rho*w],
                       [By],
                       [Bz],
                       [0.5*rho*(u*u+v*v+w*w)+1./(gamma-1.)*p]])
                       # [rho*w],
                       # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
    
        Q2    = sym.Matrix([[rho],
                       [u],
                       [v],
                       [w],
                       [By],
                       [Bz],
                       [p]])
         
    else:
        print('Eqn not implemented yet.')
        import sys;sys.exit()                                          
        
    prim = [rho*u ,rho*v,rho*w,u       ,v       ,w,rho*et]
    cons = [rhou  ,rhov ,rhow, rhou/rho,rhov/rho,rhow/rho, rhoet] 
    switch = {}
    switchb= {}
    for p,c in zip(prim,cons):
        switch[p] = c
        switchb[c]= p
    
    
    Flx_x = Flx_x.subs(switch,simultaneous=True,evaluate=False)
    Flx_y = Flx_y.subs(switch,simultaneous=True,evaluate=False)    

    
    # Ax, aY and Az are obtained from the Jacobian of the Flx_x, Flx_y and Flx_z compenant 
    Ax = Flx_x.jacobian(Q)
    Ay = Flx_y.jacobian(Q)     



    Ax = Ax.subs(switchb,simultaneous=True,evaluate=False) # to switch back to "velocities" expressions
    Ay = Ay.subs(switchb,simultaneous=True,evaluate=False) # to switch back to "velocities" expressions




    
    
    M   = Q_CS.jacobian(Q2)

    Mm1 = M.inv() 

    Atild = Mm1*Ax*M
    Btild = Mm1*Ay*M    

    c = sym.Symbol('c')
    p     = (gamma-1)*(rho*et-(u**2+v**2)*rho/2)
    
    cexp = sym.simplify(gamma*p/rho)

    Atild = sym.simplify(Atild)
    Btild = sym.simplify(Btild)

    eq = c**2/gamma/p*rho
    
    etsub = sym.solve(eq-1,et)
    
    Atild = Atild.subs(et,etsub[0])
    Atild = sym.simplify(Atild)

    Btild = Btild.subs(et,etsub[0])
    Btild = sym.simplify(Btild)
    
    rho = sym.Function('rho')(x,y,z,t)
    u = sym.Function('u')(x,y,z,t)
    v = sym.Function('v')(x,y,z,t)
    w = sym.Function('w')(x,y,z,t)
    By = sym.Function('By')(x,y,z,t)
    Bz = sym.Function('Bz')(x,y,z,t)
    p = sym.Function('p')(x,y,z,t)
    
    if eqn == 'Euler':
        Phi_x = sym.Matrix([rho,u,v,p]).diff(x)
        Phi_y = sym.Matrix([rho,u,v,p]).diff(y)
    elif eqn == 'MHD':    
        Phi_x = sym.Matrix([rho,u,v,w,By,Bz,p]).diff(x)
        Phi_y = sym.Matrix([rho,u,v,w,By,Bz,p]).diff(y)

    Li,charSpeeds_x, Pleft_x   = charspeed_Lis(Atild,Phi_x)
    Mi,charSpeeds_y, Pleft_y = charspeed_Lis(Btild,Phi_y)

    # rho = sym.Function('rho')(xi,eta,zeta,t)
    # u = sym.Function('u')    (xi,eta,zeta,t)
    # v = sym.Function('v')    (xi,eta,zeta,t)
    # w = sym.Function('w')    (xi,eta,zeta,t)
    # p = sym.Function('p')    (xi,eta,zeta,t)    
    
    # Phi_xi = sym.Matrix([rho,u,v,p]).diff(xi)
    # Phi_eta = sym.Matrix([rho,u,v,p]).diff(eta)
    
    # Li,charSpeeds_xi, Pleft_xi   = charspeed_Lis(Atild,Phi_xi)
    # Mi,charSpeeds_eta, Pleft_eta = charspeed_Lis(Btild,Phi_eta)
 
    # sym.pprint(Li)
    
    return {'xi' :[Li,charSpeeds_x , Pleft_x],
            'eta':[Mi,charSpeeds_y, Pleft_y]}


#
#  Get characteristic speeds and eigenvectors
# 
def charspeed_Lis(Atild,Phi):


    (Sm1,D) = Atild.diagonalize()
    S = Sm1.inv()

    LeftEV = [ v[2] for v in Atild.left_eigenvects()]
    
    EV     = []
    for v in Atild.left_eigenvects():
      for i in range(0,v[1]):
        EV.append(v[0])
    
    
    Pleft = []
    for row in LeftEV:
      for a in row: 
        Pleft.append(a)             
    
    Pleft = sym.Matrix(Pleft)
    
    LiSym = Pleft*Phi
    
    LiFinal = sym.Matrix([e*l for e,l in zip(EV,LiSym)])

    return LiFinal,EV,Pleft

if __name__ == '__main__':
    # name,fbeg,fend,nstep = sys.argv
    sym.pprint(characteristics('Euler')['xi'])   # Pour tester, j'affiche les Li (pour la direction xi)
   
