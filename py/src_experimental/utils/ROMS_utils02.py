# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:53:22 2024

@author: shjo
"""

def rho2u_3d(var):
    N,Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:,:L]+var[:,:,1:Lp])
    return var_u
def rho2v_3d(var):
    N,Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:,:M,:]+var[:,1:Mp,:])
    return var_v
def rho2u_2d(var):
    Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:L]+var[:,1:Lp])
    return var_u
def rho2v_2d(var):
    Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:M,:]+var[1:Mp,:])
    return var_v
    
def rho2u_2d(var):
    N,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:L]+var[:,1:Lp])
    return var_u
def rho2u_3d(var):
    N,Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*(var[:,:,:L]+var[:,:,1:Lp])
    return var_u
def rho2u_4d(var):
    T,N,Mp,Lp=var.shape
    L=Lp-1
    var_u=0.5*( var[:,:,:,:L]+var[:,:,:,1:Lp] )
    return var_u
def rho2v_3d(var):
    T,Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:,:M,:]+var[:,1:Mp,:])
    return var_v
def rho2v_4d(var):
    T,N,Mp,Lp=var.shape
    M=Mp-1
    var_v=0.5*(var[:,:,:M,:]+var[:,:,1:Mp,:])
    return var_v