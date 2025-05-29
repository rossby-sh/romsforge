# -*- coding: utf-8 -*-
"""
Created on Sat Jan 21 15:47:04 2023

@author: birostris
@email : birostris36@gmail.com

Name : 
Reference :
Description :
"""

import numpy as np
from scipy.interpolate import griddata


def stretching(Vstretching,theta_s,theta_b,Layer_N,kgrid=1):
    
    # Vstretching=MyVar['Vstretching']
    # theta_s=MyVar['Theta_s']
    # theta_b=MyVar['Theta_b']
    # N=MyVar['Layer_N']    
    Np=Layer_N+1
    if Vstretching==1:
        ds=1/Layer_N
        if kgrid==1:
            nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-.5
            s=(lev-Layer_N)*ds
        
        if theta_s>0:
            Ptheta=np.sinh(theta_s*s)/np.sinh(theta_s)
            Rtheta=np.tanh(theta_s*(s+0.5))/(2.0*np.tanh(0.5*theta_s))-0.5
            C=(1.0-theta_b)*Ptheta+theta_b*Rtheta
        else:
            C=s 
        
    elif Vstretching==2:
        
        alfa=1.0
        beta=1.0
        ds=1.0/Layer_N
        if kgrid==1:
            Nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-.5
            s=(lev-Layer_N)*ds
        
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
            if theta_b>0:
                Cbot=-1.0+np.sinh(theta_b*(s+1.0))/np.sinh(theta_b)
                weigth=(s+1.0)**alfa*(1.0+(alfa/beta)*(1.0-(s+1)**beta))
                C=weigth*Csur+(1.0-weigth)*Cbot
            else:
                C=Csur
    elif Vstretching==4:
        ds=1.0/Layer_N
        if kgrid==1:
            Nlev=Np
            lev=np.arange(Layer_N+1)
            s=(lev-Layer_N)*ds
        else:
            nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-0.5
            s=(lev-Layer_N)*ds
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=-s**2
        
        if theta_b>0:
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
        
        
    elif Vstretching==5:
        if kgrid==1:
            nlev=Np
            lev=np.arange(Layer_N+1)
            s=-(lev*lev-2.0*lev*Layer_N+lev+Layer_N*Layer_N-Layer_N)/(Layer_N*Layer_N-Layer_N)-\
                0.01*(lev*lev-lev*Layer_N)/(1.0-Layer_N)
            s[0]=-1.0
        else:
            Nlev=Layer_N
            lev=np.arange(1,Layer_N+1)-0.5
            s=-(lev*lev-2.0*lev*Layer_N+lev+Layer_N*Layer_N-Layer_N)/(Layer_N*Layer_N-Layer_N)-\
                0.01*(lev*lev-lev*Layer_N)/(1.0-Layer_N)
        if theta_s>0:
            Csur=(1.0-np.cosh(theta_s*s))/(np.cosh(theta_s)-1.0)
        else:
            Csur=-s**2
        if theta_b>0:
            Cbot=(np.exp(theta_b*Csur)-1.0)/(1.0-np.exp(-theta_b))
            C=Cbot
        else:
            C=Csur
    return s,C
    
            
def ztosigma(var,z,depth):
    Ns,Mp,Lp=z.shape
    Nz=len(depth)
    vnew=np.zeros([Ns,Mp,Lp])
    for ks in range(Ns):
        sigmalev=np.squeeze(z[ks,:,:])
        thezlevs=0*sigmalev
        for kz in range(Nz):
            thezlevs[sigmalev>depth[kz]]=thezlevs[sigmalev>depth[kz]]+1
        if np.max(thezlevs)>=Nz or np.min(thezlevs)<=0:
            print("min sigma level = "+str(np.min(z))+' - min z level = '+\
                  str(np.min(depth)))
            print("max sigma level = "+str(np.max(z))+' - max z level = '+\
                  str(np.max(depth)))            
        thezlevs=thezlevs.astype('int32')
        imat,jmat=np.meshgrid(np.arange(1,Lp+1),np.arange(1,Mp+1))
        pos=Nz*Mp*(imat-1)+Nz*(jmat-1)+thezlevs
        z1,z2=depth[thezlevs-1],depth[thezlevs]
        tmp_var=var.transpose().flatten()
        v1=tmp_var[pos-1].reshape(Mp,Lp)
        v2=tmp_var[pos].reshape(Mp,Lp)
        vnew[ks,:,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vnew

def ztosigma_1d(var,z,depth):
    Ns,Lp=z.shape
    Nz=len(depth)
    vnew=np.zeros([Ns,Lp])
    for ks in range(Ns):
        sigmalev=np.squeeze(z[ks,:])
        thezlevs=0*sigmalev
        for kz in range(Nz):
            thezlevs[sigmalev>depth[kz]]=thezlevs[sigmalev>depth[kz]]+1
        if np.max(thezlevs)>=Nz or np.min(thezlevs)<=0:
            print("min sigma level = "+str(np.min(z))+' - min z level = '+\
                  str(np.min(depth)))
            print("max sigma level = "+str(np.max(z))+' - max z level = '+\
                  str(np.max(depth))) 
                
        thezlevs=thezlevs.astype('int32')
        
        jmat= np.arange(1,Lp+1)
        pos=Nz*(jmat-1)+thezlevs
        
        z1,z2=depth[thezlevs-1],depth[thezlevs]
        tmp_var=var.transpose().flatten()
        v1=tmp_var[pos-1].reshape(Lp)
        v2=tmp_var[pos].reshape(Lp)
        vnew[ks,:]=(((v1-v2)*sigmalev+v2*z1-v1*z2)/(z1-z2))
    return vnew














        
def zlevs(Vtransform, Vstretching,theta_s, theta_b, hc, N,igrid, h, zeta):
    from copy import deepcopy
    
    # for get section
    if len(h.shape)!=2:
        h=h.reshape([1,len(h)])
        zeta=zeta.reshape([1,len(zeta)])
    Np=N+1;
    Lp,Mp=h.shape
    L=Lp-1;
    M=Mp-1;
    
    hmin=np.min(h);
    hmax=np.max(h);

    # Compute vertical stretching function, C(k):

    if igrid==5:
        kgrid=1
        s,C=stretching(Vstretching, theta_s, theta_b, N, kgrid)
    else:
        kgrid=0
        s,C=stretching(Vstretching, theta_s, theta_b, N, kgrid)

    #  Average bathymetry and free-surface at requested C-grid type.

    if igrid==1:
        hr=deepcopy(h)
        zetar=deepcopy(zeta)
    elif igrid==2:
        hp=0.25*(h[:L,:M]+h[1:Lp,:M]+h[:L,1:Mp]+h[1:Lp,1:Mp])
        zetap=0.25*(zeta[:L,:M]+zeta[1:Lp,:M]+zeta[:L,1:Mp]+zeta[1:Lp,1:Mp])
    elif igrid==3:
        hu=0.5*(h[:L,:Mp]+h[1:Lp,:Mp])
        zetau=0.5*(zeta[:L,:Mp]+zeta[1:Lp,:Mp])
    elif igrid==4:
        hv=0.5*(hp[:Lp,:M]+h[:Lp,1:Mp])
        zetav=0.5*(zeta[:Lp,:M]+zeta[:Lp,1:Mp])        
    elif igrid==5:
        hr=deepcopy(h)
        zetar=deepcopy(zeta)


    # Compute depths (m) at requested C-grid location.
    if Vtransform==1:
        if igrid==1:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hr
                z[:,:,k]=z0+zetar*(1+z0/hr)
        elif igrid==2:
            z=np.zeros([zetap.shape[0],zetap.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hp
                z[:,:,k]=z0+zetap*(1+z0/hp)
        elif igrid==3:
            z=np.zeros([zetau.shape[0],zetau.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hu
                z[:,:,k]=z0+zetau*(1+z0/hu)
        elif igrid==4:
            z=np.zeros([zetav.shape[0],zetav.shape[-1],N])
            for k in range(N):
                z0=(s[k]-C[k])*hc+C[k]*hv
                z[:,:,k]=z0+zetav*(1+z0/hv)
        elif igrid==5:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],Np])
            z[:,:,0]=-hr
            for k in range(1,Np):
                z0=(s[k]-C[k])*hc+C[k]*hr
                z[:,:,k]=z0+zetar*(1+z0/hr)
    
    elif Vtransform==2:
        if igrid==1:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hr)/(hc+hr)
                z[:,:,k]=zetar+(zeta+hr)*z0
        elif igrid==2:
            z=np.zeros([zetap.shape[0],zetap.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hp)/(hc+hp)
                z[:,:,k]=zetap+(zetap+hp)*z0
        elif igrid==3:
            z=np.zeros([zetau.shape[0],zetau.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hu)/(hc+hu)
                z[:,:,k]=zetau+(zetau+hu)*z0
        elif igrid==4:
            z=np.zeros([zetav.shape[0],zetav.shape[-1],N])
            for k in range(N):
                z0=(hc*s[k]+C[k]*hv)/(hc+hv)
                z[:,:,k]=zetav+(zetav+hv)*z0
        elif igrid==5:
            z=np.zeros([zetar.shape[0],zetar.shape[-1],Np])
            for k in range(0,Np):
               z0=(hc*s[k]+C[k]*hr)/(hc+hr)
               z[:,:,k]=zetar+(zetar+hr)*z0
    z=np.squeeze(np.transpose(z,[2,0,1]))
    
    return z



def get_section(ncG,ncD,vname,lon_rng,lat_rng,tindx=0):
        
    if vname =='u':
        mask=ncG['mask_u']
        igrid=3
    elif vname =='v':
        mask=ncG['mask_v']
        igrid=4
    else:
        mask=ncG['mask_rho']
        igrid=1
        
    lat_rho,lon_rho=ncG['lat_rho'][:].data,ncG['lon_rho'][:].data
    M,L=lon_rho.shape
    dl=1.5*max(np.max(np.abs(lon_rho[1:M,:]-lon_rho[:M-1,:])),
                   np.max(np.abs(lon_rho[:,1:L]-lon_rho[:,:L-1])),
                   np.max(np.abs(lat_rho[1:M,:]-lat_rho[:M-1,:])),
                   np.max(np.abs(lat_rho[:,1:L]-lat_rho[:,:L-1])))
            
    if dl>10:
        dl=1.5        
    # ncG.variables.keys()
    
    ncVar={'Layer_N':len(ncD['Cs_r'][:]),'Vtransform':ncD['Vtransform'][:].data,\
           'Vstretching':ncD['Vstretching'][:].data,'Theta_s':ncD['theta_s'][:].data,\
               'Theta_b':ncD['theta_b'][:].data,'Tcline':ncD['Tcline'][:].data}
    
    lon_rng_re=[lon_rng[0]-dl,lon_rng[-1]+dl]
    print(lat_rng)

    lat_rng_re=[lat_rng[0]-dl,lat_rng[-1]+dl]
    sub= (lon_rho>=lon_rng_re[0]) & (lon_rho<=lon_rng_re[-1]) & (lat_rho>=lat_rng_re[0]) & (lat_rho<=lat_rng_re[-1])
    ival,jval=np.sum(sub,axis=0),np.sum(sub,axis=1)
    imin,imax=np.min(np.where(ival!=0)),np.max(np.where(ival!=0))
    jmin,jmax=np.min(np.where(jval!=0)),np.max(np.where(jval!=0))

    # Get subgrids
    lon_rho=lon_rho[jmin:jmax+1,imin:imax+1];
    lat_rho=lat_rho[jmin:jmax+1,imin:imax+1];
    sub=sub[jmin:jmax+1,imin:imax+1];
    mask=mask[jmin:jmax+1,imin:imax+1];

    if len(lon_rng)==1:
        print('N-S section at longitude: '+str(lon_rng[0]))
        if len(lat_rng)==1:
            raise
        elif len(lat_rng)==2:
            latsec=np.arange(lat_rng[0],lat_rng[-1],dl)
        lonsec=np.zeros_like(latsec)+lon_rng[0]
    elif len(lat_rng)==1:
        print('E-W section at latitude: '+str(lat_rng[0]))
        if len(lon_rng)==2:
            lonsec=np.arange(lon_rng[0],lon_rng[-1],dl)
        latsec=np.zeros_like(lonsec)+lat_rng[0]
    elif len(lon_rng)==2 and len(lat_rng)==2:
        Npts= np.max(np.array(abs(lon_rng[-1]-lon_rng[0])/dl,\
                       abs(lat_rng[-1]-lat_rng[0])/dl))
        if lon_rng[0]==lon_rng[-1]:
            lonsec=lon_rho[0]+np.zeros([1,Npts+1])
        else:
            lonsec=np.arange(lon_rng[0],lon_rng[-1],(lon_rng[-1]-lon_rng[0])/Npts)
        if lat_rng[0]==lat_rng[-1]:
            latsec=lat_rng[0]+np.zeros([1,Npts+1])
        else:
            latsec=np.arange(lat_rng[0],lat_rng[-1],(lat_rng[-1]-lat_rng[0])/Npts)
    elif len(lon_rho)!=len(lat_rho):
        raise

    Npts=len(lonsec)

    # Get the subgrid
    sub=np.zeros_like(lon_rho)
    for i in range(Npts):
        sub[(lon_rho>lonsec[i]-dl)&(lon_rho<lonsec[i]+dl)&(lat_rho>latsec[i]-dl)&\
            (lat_rho<latsec[i]+dl)]=1

    # get the coefficients of the objective analysis
    londata=lon_rho[sub==1]
    latdata=lat_rho[sub==1]

    # coef=oacoef(londata,latdata,lonsec,latsec,100e3);
    # Get the mask
    mask=mask[sub==1]
    m1=griddata((londata,latdata),mask,(lonsec,latsec),'nearest');

    londata=londata[mask==1]
    latdata=latdata[mask==1]

    # Get the vertical levels
    topo=ncG['h'][:]
    topo=topo[jmin:jmax+1,imin:imax+1]
    topo=topo[sub==1]
    topo=topo[mask==1]
    topo=griddata((londata,latdata),topo,(lonsec,latsec),'linear');

    zeta=ncD['zeta'][tindx,jmin:jmax+1,imin:imax+1]
    zeta=zeta[sub==1]
    zeta=zeta[mask==1]
    zeta=griddata((londata,latdata),zeta,(lonsec,latsec),'linear');

    Z=zlevs(ncVar['Vtransform'], ncVar['Vstretching'],ncVar['Theta_s'],\
            ncVar['Theta_b'], ncVar['Tcline'], ncVar['Layer_N'],igrid, topo, zeta)
    N,M=Z.shape
    # Loop on the vertical levels
    VAR=np.zeros_like(Z)
    for i in range(N):
        var=np.squeeze(ncD[vname][tindx,i,jmin:jmax+1,imin:imax+1])
        var=var[sub==1]
        var=var[mask==1]
        var=griddata((londata,latdata),var,(lonsec,latsec),'linear');
        VAR[i,:]=m1*var
    X=np.zeros([N,M])
    if len(lat_rng)==1 and len(lon_rng)==2:
        for i in range(N):
            X[i,:]=lonsec
    elif len(lat_rng)==2 and len(lon_rng)==1:
        for i in range(N):
            X[i,:]=latsec
    elif len(lat_rng)==2 and len(lon_rng)==2:
        for i in range(N):
            X[i,:]=lonsec
        print('!!! Dist --> In development !!!')
    return X,Z,VAR

            
            
            
    