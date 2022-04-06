import sys
import re
import numpy as np
from genKer import  genrk3, genrk3update, genFilter, genBC, append_Rhs
import os

from generate.equations import Src_BC_dif 

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
    append_Rhs(Src_conv, 5, 4, rhsname, locname_conv, update=False,rhs=rhs)
    append_Rhs(Src_dif , 3 ,2 , rhsname, locname_dif , update=True ,rhs=rhs,stored=True)  
    #append_Rhs(Src_SA , 3 ,2 , rhsname, locname_dif , update=True ,rhs=rhs,stored=True)                        
    #append_Rhs(divF, 9,8, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 7,6, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 5,4, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           
    #append_Rhs(divF, 3,2, rhsname, vnamesrc_divFx, update=False,rhs=rhs,stored=True)                           

# Generate Filters (if required):      
#    genFilter(5,4, len(varsolved),rhs=rhs)     #No filter required

# Generate BCs:
    genBC(Src_conv ,5,4,rhsname , locname_conv, update=False,rhs=rhs,stored=False)
    genBC(Src_dif  ,3 ,2 ,rhsname , locname_dif , update=True ,rhs=rhs,stored=True)
    
    #--j1
    genBC(Src_BC['j1']  ,3,2,rhsname , locname_bc, setbc=[True,{'Wall_BC':{'j1':['q']}}]  , update=False,rhs=rhs)
    genBC(Src_BC_conv['j'] , 3,2,rhsname, locname_bc_conv, setbc=[True, {'Symm_Wall':{'j1':['rhs']}}],  update=False, rhs=rhs )
    genBC(Src_BC_dif['j1'], 3,2, rhsname, locname_bc_dif, setbc=[True, {'Symm_Wall':{'j1':['rhs']}}]], update=False, rhs=rhs )


    #--jmax
    genBC(Src_BC_)
    # genBC(Src_phybc_rhs ,5,4,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'jmax':['rhs']}}], update=False,rhs=rhs)
    # genBC(Src_phybc_qmax,5,4,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'jmax':['q'  ]}}], update=False,rhs=rhs)
    
    
#    # j1
#    genBC(Save_eqns['Src_conv'] ,11,10,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'j1':['rhs']}}]  , update=False,rhs=rhs)
#    # genBC(Src_phybc_q1  ,5,4,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'j1':['q'  ]}}]  , update=False,rhs=rhs)
#
#    # jmax
#    genBC(Save_eqns['Src_conv'] ,11,10,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'jmax':['rhs']}}], update=False,rhs=rhs)
#    # genBC(Src_phybc_qmax,5,4,rhsname , vnamesrc_bc, setbc=[True,{'wall':{'jmax':['q'  ]}}], update=False,rhs=rhs)

    # i1
    # genBC(Src_phybc_rhs ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'i1':['rhs']}}]  , update=False,rhs=rhs)
    # genBC(Src_phybc_q1  ,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'i1':['q'  ]}}]  , update=False,rhs=rhs)
    # imax
    # genBC(Src_phybc_rhs ,5,,4,rhsname , locname_bc, setbc=[True,{'wall':{'imax':['rhs']}}], update=False,rhs=rhs)      
    # genBC(Src_phybc_qmax,5,4,rhsname , locname_bc, setbc=[True,{'wall':{'imax':['q'  ]}}], update=False,rhs=rhs)      



# Extract RHS info:
    rhs.export()

if __name__ == '__main__':
    main()
