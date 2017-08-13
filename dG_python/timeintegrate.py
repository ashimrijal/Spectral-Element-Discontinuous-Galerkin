def ADER_Wave_dG(u,v,D,NP,num_element,dx,w,x,t,r0,rn,dt,rho,mu):

    import rate
    import numpy as np
        
    wave_dG_return = rate.wave_dG(u,v,D,NP,num_element,dx,w,r0,rn,rho,mu,x,t)
    Klu = wave_dG_return['Hu']
    Klv = wave_dG_return['Hv']
    
    fac = dt
    K3u = fac * Klu
    K3v = fac * Klv
    
    for k in range(1,NP+1): 
        fac = fac * (dt) /(k+1)
        wave_dG_return = rate.wave_dG(Klu,Klv,D,NP,num_element,dx,w,r0,rn,rho,mu,x,t)
        Klu = wave_dG_return['Hu']
        Klv = wave_dG_return['Hv']
        
        K3u = K3u + fac * Klu
        K3v = K3v + fac * Klv

    Hu = u + K3u
    Hv = v + K3v

    return {'Hu':Hu,'Hv':Hv}



def ADER_Wave_1D_GL(u_1,v_1,u_2,v_2,S,psi,D,NP,num_element1,num_element2,dx1,dx2,w,x,t,r0,rn,dt,rho1,mu1,rho2,mu2,alpha,Tau_0,fric_law):
    
    import rate
    import numpy as np
    
    wave_1D_GL_return = rate.wave_1D_Friction(u_1,v_1,u_2,v_2,S,psi,D,NP,num_element1,num_element2,dx1,dx2,w,r0,rn,rho1,mu1,rho2,mu2, x,t,alpha,Tau_0,fric_law)
    Klu_1 = wave_1D_GL_return['Hu_1']
    Klv_1 = wave_1D_GL_return['Hv_1']
    Klu_2 = wave_1D_GL_return['Hu_2']
    Klv_2 = wave_1D_GL_return['Hv_2']
    Kld   = wave_1D_GL_return['H_d']
    Klpsi   = wave_1D_GL_return['H_psi']

    fac = dt
    K3u_1 = fac * Klu_1
    K3v_1 = fac * Klv_1
    K3u_2 = fac * Klu_2
    K3v_2 = fac * Klv_2
    K3d =  fac * Kld
    K3psi =  fac * Klpsi
    
    for k in range(1,NP+1): 
        fac = fac * (dt) /(k+1)
        wave_1D_GL_return = rate.wave_1D_Friction(Klu_1,Klv_1,Klu_2,Klv_2,S,psi,D,NP,num_element1,num_element2,dx1,dx2,w,r0,rn,rho1,mu1,rho2,mu2,x,t,alpha,Tau_0,fric_law)
        Klu_1 = wave_1D_GL_return['Hu_1']
        Klv_1 = wave_1D_GL_return['Hv_1']
        Klu_2 = wave_1D_GL_return['Hu_2']
        Klv_2 = wave_1D_GL_return['Hv_2']
        Kld   = wave_1D_GL_return['H_d']
        Klpsi   = wave_1D_GL_return['H_psi']
    
        K3u_1 = K3u_1 + fac * Klu_1
        K3v_1 = K3v_1 + fac * Klv_1
        K3u_2 = K3u_2 + fac * Klu_2
        K3v_2 = K3v_2 + fac * Klv_2
        K3d= K3d + fac * Kld
        K3psi = K3psi + fac * Klpsi

    Hu_1 = u_1 + K3u_1
    Hv_1 = v_1 + K3v_1
    Hu_2 = u_2 + K3u_2
    Hv_2 = v_2 + K3v_2
    H_d = S + K3d
    H_psi = psi + K3psi
    
    return {'Hu_1':Hu_1,'Hv_1':Hv_1,'Hu_2':Hu_2,'Hv_2':Hv_2, 'H_d':H_d, 'H_psi':H_psi}



