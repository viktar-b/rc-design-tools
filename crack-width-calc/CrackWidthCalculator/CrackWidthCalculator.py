import numpy as np

import xlwings as xw

@xw.func
def get_crack_width(
    cnom_tension,
    cnom_compression, 
    c_dur, 
    moment 
):
    """
    To activate a function in the excel file, press "Import Functions" 
    in the xlwings tab 
    """
    crack_width = 1 
    return crack_width 

class ConcreteSection: 
    """
    Standard units: 
        dimensions - mm
        stresses - N/mm^2 or MPa 
        time - days 
  
    TODO 
    Standard parameters are assumed to be constant at this stage 
    Write setters and getters for the fields below to verify dimesions are valid 
    Calcs are for slab only, refactor to include beams 
    """
    _aggregate_size = 20 
    _width = 1000



    def __init__(self, depth):
        """
        depth - thickness of the section 
        width - to be taken as 1000 mm if slab is 

        """
        self.depth = depth 
        

class CrackWidth: 
    pass 


class SectionProperties: 
    pass


class ConcreteProperties: 

    gamma_c = 1.5
    alpha_cc_flexure = 0.85
    alpha_cc_shear = 1 
    alpha_ct = 1

    #TODO values are assumed as context. write a class for creep calculations
    relative_humidity = 0.75
    exposed_perimter = 1000
    age_of_concrete_at_moment_considered = 2557 
    age_of_concrete_at_loading = 28
    creep = 1.26

    def __init__(self, f_ck):
        self.f_ck = f_ck
        self.f_av = self.get_f_av()
        self.f_cm = self.get_f_cm()
        self.f_ctm = self.get_f_ctm()
        self.E_cm = self.get_E_cm()
        self.f_cd_flexural = self.get_f_cd_flexural()
        self.f_cd_shear = self.get_f_cd_shear()
        self.f_ctd = self.get_f_ctd()
   
    def get_f_av(self):
        return 0.459*self.f_ck

    def get_f_cm(self):
        return self.f_ck + 8

    def get_f_ctm(self):
        if self.f_ck < 50:
            return 0.3*self.f_ck**(2/3)
        else:
            return 2.12*np.log(1 + self.f_cm/10)            
    
    def get_E_cm(self):
        return 22*((self.f_ck+8)/10)**0.3*1000

    def get_f_cd_flexural(self):
        return self.f_ck * self.alpha_cc_flexure / self.gamma_c
    
    def get_f_cd_shear(self):
        return self.f_ck * self.alpha_cc_shear / self.gamma_c

    def get_f_ctd(self):
            return 0.7 * self.f_ctm * self.alpha_ct / self.gamma_c

    def to_string(self):
        return  f"f_ck = {self.f_ck:0.2f}\n" + \
                f"f_av = {self.f_av:0.2f}\n" + \
                f"f_cm = {self.f_cm:0.2f}\n" + \
                f"f_ctm = {self.f_ctm:0.2f}\n" + \
                f"E_cm = {self.E_cm:0.2f}\n" + \
                f"f_cd_flexural = {self.f_cd_flexural:0.2f}\n" + \
                f"f_cd_shear = {self.f_cd_shear:0.2f}\n" + \
                f"f_ctd = {self.f_ctd:0.2f}\n"    

    



    

class SteelProperties: 
    pass







"""
Preparation - bottom-up calcualtions for the crack width 

w_k = s_rmax*(e_sm - e_cm)

s_rmax for long term loading
    if spacing > 5 (c + diam/2)
        s_rmax = 1.3*(h-x)      comm: x - d_c cocnrete depth in compression 
    else 
        s_rmax = k_3*c + k_1*k_2*k_4*diam/rho_eff   comm: 
                                                        c - cover to bar controlling crack for crackwidth calculations 
                                                        rho_eff = A_s/A_c,eff

(e_sm - e_cm) for long term loading  

    max(0,6*(sigma_s/E_s ), (sigma_s - k_t*f_ct_eff/rho_p_eff*(1 + alpha_e*rho_eff))/E_s)     comm:

                                    k_t = 0.4
                                    f_ct_eff = 3.2 MPa for f_ck = 35 
                                    alpha_e = E_s / E_cm = 200,000 / 34.077 for f_ck = 35 
                                    d - depth to centroid F77 = sum(A_s_layer*d_layer)/A_s_total 
                                    
                                    For pure bending:
                                        sigma_s = M_sls_permanent*10^6*(d - d_c_lt)/I_cracked + M_sls_permanent*10^6*(d - d_c_st)/I_cracked
                                    
                                    To be determined: 
                                        rho_eff = A_s/A_ceff
                                            A_ceff = h_ceff*width 
                                                h_ceff = max( min(
                                                                    (h' - d_c)/3,  
                                                                    2.5*(h-d), 
                                                                    h/2)
                                                                    ), 
                                                            bar_group)  =   c_dur + d_shear_bar + d_first_layer + 
                                                                            dist_between_layers (if second is present) +  
                                                                            d_second_layer 
                                                                            + aggregate_size + 5 TODO??
                                        d_c_lt - F192 = 
                                        d_c_st - F175
                                        I_cracked 
                                       
                                        
    """