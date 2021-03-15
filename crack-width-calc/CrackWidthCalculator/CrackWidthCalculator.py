import numpy as np

import xlwings as xw

@xw.func
def get_crack_width(c_nom_tension, c_nom_compression, c_dur, moment,
    rf_diameter_first_layer = 0,
    rf_spacing_first_layer = 0,
    rf_diameter_second_layer = 0,
    rf_spacing_second_layer = 0,
    rf_diameter_third_layer = 0,
    rf_spacing_third_layer = 0):

    """
    To activate a function in the excel file, press "Import Functions" 
    in the xlwings tab 
    """

    section_properrties = SectionProperties(c_nom_tension, c_nom_compression, c_dur,
        rf_diameter_first_layer,
        rf_spacing_first_layer,
        rf_diameter_second_layer,
        rf_spacing_second_layer,
        rf_diameter_third_layer,
        rf_spacing_third_layer
    )
    

    serviceability_checks = ServiceabilityChecks(section_properrties, moment)
    crack_width_calc = CrackWidthCalc(serviceability_checks)

    crack_width = crack_width_calc.get_crack_width()
    return crack_width 


class CrackWidthCalc: 

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

    def __init__(self, section_properrties):
        self.section_properties = section_properrties
        
    def get_crack_width(self, c_nom_tension, c_nom_compression, c_dur, moment):
        self.c_nom_tension = c_nom_tension
        self.c_nom_compression = c_nom_compression 
        self.c_dur = c_dur
        self.moment = moment

        s_rmax = self.get_s_rmax()
        e_sm = self.get_e_sm()
        e_cm = self.get_e_cm()

        crack_width = s_rmax*(e_sm - e_cm)

        return crack_width

    def get_s_rmax():
        pass

    def get_e_sm():
        pass

    def get_e_cm():
        pass

class ServiceabilityChecks: 

    def __init__(self, section_properties, moment):
        self.section_properties = section_properties
        self.sls_moment = moment
        
        self.E_ceff_short_term = self.section_properties.concrete_properties.E_cm
        self.E_ceff_long_term = self.E_ceff_short_term / (1 + self.section_properties.concrete_properties.creep)

        self.d_c_st = self.get_concrete_depth_in_compression(self.E_ceff_short_term)
        self.d_c_lt = self.get_concrete_depth_in_compression(self.E_ceff_long_term)

        self.I_cracked_short_term = self.get_cracked_second_moment_of_area(self.d_c_st, self.E_ceff_short_term)
        self.I_cracked_long_term = self.get_cracked_second_moment_of_area(self.d_c_lt, self.E_ceff_long_term)

        self.longterm_sigma_s = self.get_longterm_sigma_s()

    def get_longterm_sigma_s(self):
        """                       
        For pure bending:
            sigma_s = M_sls_permanent*10^6*(d - d_c_lt)/I_cracked + M_sls_permanent*10^6*(d - d_c_st)/I_cracked
        """
        d = self.section_properties.depth_to_centroid

        sigma_s = self.sls_moment*10**6 * (d - self.d_c_lt) / self.I_cracked_long_term + \
                    self.sls_moment*10*6 * (d - self.d_c_st) / self.I_cracked_short_term

        return sigma_s  

    def get_cracked_second_moment_of_area(self, d_c, E_ceff):
        A_s  = self.section_properties.steel_area_tension 
        d = self.section_properties.depth_to_centroid
        E_s = self.section_properties.steel_properties.E_s 
        b = self.section_properties.width
        
        I_cracked = A_s * (d - d_c)**2 + 1/3 * (E_ceff / E_s) * (b * d_c**3)
  
        return I_cracked  


    def get_concrete_depth_in_compression(self, effective_modulus_of_concrete): 
        """
        Concrete depth in compression in pure bending
        """
        E_ceff = effective_modulus_of_concrete

        A_s  = self.section_properties.steel_area_tension 
        E_s = self.section_properties.steel_properties.E_s 
        d = self.section_properties.depth_to_centroid
        b = self.section_properties.width

        d_c  = (-A_s*E_s + np.sqrt((A_s*E_s)**2 + 2*b*A_s*E_s*E_ceff*d)) / (b * E_ceff)  
        
        return d_c

    def to_string(self):
        return f"E_ceff_short_term = {self.E_ceff_short_term:0.2f}\n" + \
                f"d_c_st = {self.d_c_st:0.2f}\n" + \
                f"I_cracked_short_term = {self.I_cracked_short_term:0.2f}\n" + \
                f"E_ceff_long_term = {self.E_ceff_long_term:0.2f}\n" + \
                f"d_c_lt = {self.d_c_lt:0.2f}\n" + \
                f"I_cracked_long_term = {self.I_cracked_long_term:0.2f}\n" + \
                f"longterm_sigma_s = {self.longterm_sigma_s:0.2f}\n"


class SectionProperties: 
    
    #PROJECT SPECIFIC SECTION PROPERTIES (to be moved to init)
    depth = 750 
    width = 1000 
    f_yd = 500
    f_ck = 35
    spacing_between_rf = 25

    def __init__(self, c_nom_tension, c_nom_compression, c_dur,
                        rf_diameter_first_layer, rf_spacing_first_layer,
                        rf_diameter_second_layer, rf_spacing_second_layer,
                        rf_diameter_third_layer, rf_spacing_third_layer):
        
        self.c_nom_tension = c_nom_tension
        self.c_nom_compression = c_nom_compression
        self.c_dur = c_dur

        self.rf_diameter_first_layer = rf_diameter_first_layer
        self.rf_spacing_first_layer = rf_spacing_first_layer
        self.rf_diameter_second_layer = rf_diameter_second_layer
        self.rf_spacing_second_layer = rf_spacing_second_layer
        self.rf_diameter_third_layer = rf_diameter_third_layer
        self.rf_spacing_third_layer = rf_spacing_third_layer
        
        self._area_1 = self.get_single_layer_steel_area(self.rf_diameter_first_layer, self.rf_spacing_first_layer)
        self._area_2 = self.get_single_layer_steel_area(self.rf_diameter_second_layer, self.rf_spacing_second_layer) 
        self._area_3 = self.get_single_layer_steel_area(self.rf_diameter_third_layer, self.rf_spacing_third_layer)
        
        self.steel_area_tension = self.get_steel_area_tension()
        self.depth_to_centroid = self.get_depth_to_centroid()

        self.steel_properties = SteelProperties(self.f_yd)
        self.concrete_properties = ConcreteProperties(self.f_ck)

        self.depth_to_neutral_axis = self.depth / 2
        self.second_moment_of_area =  self.width * self.depth**3 / 12 

    def get_steel_area_tension(self):
        return self._area_1 + self._area_2 + self._area_3

    def get_depth_to_centroid(self):
        d_1 = self.depth - self.c_nom_tension - self.rf_diameter_first_layer/2   
        d_2 = d_1 - self.rf_diameter_first_layer/2  - self.spacing_between_rf - self.rf_diameter_second_layer/2 
        d_3 = d_2 - self.rf_diameter_second_layer/2 - self.spacing_between_rf - self.rf_diameter_third_layer/2

        return (d_1*self._area_1 + \
                d_2*self._area_2 + \
                d_3*self._area_3) / self.steel_area_tension
    
    def get_single_layer_steel_area(self, diameter, spacing):
        if diameter != 0:
            return (np.pi*diameter**2 / 4) * (self.width / spacing)
        else:
            return 0
    
    def to_string(self):
        return f"steel_area_tension = {self.steel_area_tension:0.2f}\n" + \
                f"depth_to_centroid = {self.depth_to_centroid:0.2f}\n"


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
    creep = 1.256

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
    
    E_s = 200000
    gamma_s = 1.15 

    def __init__(self, f_yk):
        self.f_yk = f_yk 

        self.f_yd = f_yk / self.gamma_s 

    def to_string(self):
        return f"E_s = {self.E_s:0.2f}\n" + \
                f"gamma_s = {self.gamma_s:0.2f}\n" + \
                f"f_yk = {self.f_yk:0.2f}\n" + \
                f"f_yd = {self.f_yd:0.2f}\n"