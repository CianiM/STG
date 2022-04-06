import sympy as sym
sym.init_printing(use_latex=True, wrap_line=False)
rho,u,et=sym.symbols(['rho','u','et'],Real=True)
u=rho*u/rho
Flx_x=sym.Matrix([[rho],[(rho*u)**2/rho],[rho*et]])

Q=sym.Matrix([[rho],[rho*u],[rho*et]])

Ax=Flx_x.jacobian(Q)
sym.pprint(Ax)
