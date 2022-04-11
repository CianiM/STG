import sys
import re
import numpy as np
from genKer import  genrk3, genrk3update, genFilter, genBC, append_Rhs
import os


wp = 'float64'
  
from  equations import *

# Set install path (needed by genKer)
instpath = os.environ['INSTALLPATH']
incPATH  = instpath+'/src_for/includes/gen/'

hlo_glob = 4

def main():
      
    from genKer import rhs_info    
    #rhs = rhs_info(dim,wp,hlo_glob,incPATH,varsolved,varname, 
    #               consvar=consvar,varstored=varstored,varloc=varloc,varbc=varbc,
    #               coefficients=coefficients)
    rhs=rhs_info()

# Generate LHS:
    genrk3(len(varsolved)      ,rhs=rhs) 
    genrk3update(len(varsolved),rhs=rhs)

# Generate RHS:
#    append_Rhs(Src_conv, 5, 4, rhsname, locname_conv, update=False,rhs=rhs)
#    append_Rhs(Src_dif , 3 ,2 , rhsname, locname_dif , update=True ,rhs=rhs,stored=True)  
    append_Rhs(Src_rhs, 5, 4, rhsname, locname_rhs, update=False,rhs=rhs)              
# Generate Filters (if required):      
#    genFilter(5,4, len(varsolved),rhs=rhs)     #No filter required

# Generate BCs:
#    genBC(Src_conv ,5,4,rhsname , locname_conv, update=False,rhs=rhs,stored=False)
#    genBC(Src_dif  ,3 ,2 ,rhsname , locname_dif , update=True ,rhs=rhs,stored=True)
    genBC(Src_rhs ,5,4,rhsname , locname_rhs, update=False,rhs=rhs,stored=False)
    #--i1
    genBC(Src_BC_phy_i1,3,2, rhsname, locname_bc, setbc=[True,{'Inlet_surf':{'i1':['q']}}],update=False, rhs=rhs)
    genBC(Src_BC_rhs_i1  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Inlet_surf':{'i1':['rhs']}}] , update=False, rhs=rhs )
#    genBC(Src_BC_conv['i1'] , 3,2,rhsname, locname_bc_conv, setbc=[True, {'Subsonic_inflow':{'i1':['rhs']}}],  update=False, rhs=rhs )
#    genBC(Src_BC_dif['i1'], 3,2, rhsname, locname_bc_dif, setbc=[True, {'Subsonic_inflow':{'i1':['rhs']}}], update=True, rhs=rhs )   
    #--j1
    genBC(Src_BC_phy_j1  ,3,2,rhsname , locname_bc, setbc=[True,{'Low_surf':{'j1':['q']}}]  , update=False,rhs=rhs)
    genBC(Src_BC_rhs_j1  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Low_surf':{'j1':['rhs']}}] , update=False, rhs=rhs )
#    genBC(Src_BC_conv['j1'] , 3,2,rhsname, locname_bc_conv, setbc=[True, {'Symm_Wall':{'j1':['rhs']}}],  update=False, rhs=rhs )
#    genBC(Src_BC_dif['j1'], 3,2, rhsname, locname_bc_dif, setbc=[True, {'Symm_Wall':{'j1':['rhs']}}], update=True, rhs=rhs )
    #--imax
    genBC(Src_BC_rhs_imax  ,3,2,rhsname , locname_bc_rhs, setbc=[True,{'Out_surf':{'imax':['rhs']}}] , update=False, rhs=rhs )
#    genBC(Src_BC_conv['imax'] , 3,2,rhsname, locname_bc_conv, setbc=[True, {'Costant_pressure':{'imax':['rhs']}}],  update=False, rhs=rhs )
#    genBC(Src_BC_dif['imax'], 3,2, rhsname, locname_bc_dif, setbc=[True, {'Costant_pressure':{'imax':['rhs']}}], update=True, rhs=rhs )
    #--jmax
    genBC(Src_BC_rhs_jmax  ,3,2,rhsname , locname_rhs, setbc=[True,{'Top_surf':{'jmax':['rhs']}}] , update=False, rhs=rhs )
#    genBC(Src_BC_conv['jmax'], 3,2,rhsname, locname_bc_conv, setbc=[True,{'Non_reflective':{'jmax':['rhs']}}], update=False, rhs=rhs )
#    genBC(Src_BC_dif['jmax'], 3,2,rhsname, locname_bc_dif, setbc=[True,{'Non_reflective':{'jmax':['rhs']}}], update=False, rhs=rhs )
  
    
# Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
