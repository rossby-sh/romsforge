import numpy as np



def compute_sigma_3d(h, hc, theta_s, theta_b, N, Vtransform=2, Vstretching=4):


    s = (np.arange(1, N+1) - N - 0.5) / N

    # Song & Haidvogel (1994)
    if Vstretching == 1:  
        C = (1 - theta_b) * (np.sinh(theta_s * s) / np.sinh(theta_s)) + \
            theta_b * (-0.5 + 0.5 * np.tanh(theta_s * (s + 0.5)) / np.tanh(0.5 * theta_s))
    
    # Shchepetkin (2005)
    elif Vstretching == 2:  
        Csur = (1 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1)
        Cbot = -1 + (1 - np.sinh(theta_b * (s + 1))) / np.sinh(theta_b)
        Cweight = (s + 1) ** 2 * (1 + (2 / theta_b) * (1 - (s + 1) ** theta_b))
        C = Cweight * Csur + (1 - Cweight) * Cbot
    
    # Geyer (shallow sediment applications)
    elif Vstretching == 3:  
        Hscale = 3  
        Csur = -np.log(np.cosh(Hscale * np.abs(s) ** theta_s)) / np.log(np.cosh(Hscale))
        Cbot = np.log(np.cosh(Hscale * (s + 1) ** theta_b)) / np.log(np.cosh(Hscale)) - 1
        Cweight = 0.5 * (1 - np.tanh(Hscale * (s + 0.5)))
        C = Cweight * Cbot + (1 - Cweight) * Csur

    # Shchepetkin improved
    elif Vstretching == 4:  
        if theta_s > 0:
            Cs = (1 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1)
        else:
            Cs = -s ** 2
        if theta_b > 0:
            weight = np.exp(theta_b * Cs)
            C = (weight - 1) / (1 - np.exp(-theta_b))
        else:
            C = Cs
    else:
        raise ValueError("Unsupported Vstretching method. Choose between 1-4.")

    
    s_3d = s[:, np.newaxis, np.newaxis]
    C_3d = C[:, np.newaxis, np.newaxis]
    
    h_3d = h[np.newaxis, :, :]

    if Vtransform == 1:
        Z = hc * s_3d + (h_3d - hc) * C_3d
    elif Vtransform == 2:
        Z = (hc * s_3d + h_3d * C_3d) / (hc + h_3d)
    else:
        raise ValueError("Unsupported Vtransform method. Choose between 1 or 2.")

    return Z 