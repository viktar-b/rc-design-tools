import pandas as pd 
from pandas.api.types import is_string_dtype
from pandas.api.types import is_bool_dtype
import numpy as np


def get_nominal_covers(df: pd.DataFrame)->pd.DataFrame:
    """
    Convert fixed-format table with initial parameters for concrete panels 
    into cover thicknesses defined by projects requirements
    """

    fixed_output_df = pd.DataFrame(columns=[
        "XD", "XC", "XF",
        "critical_exposure",
        "c_dur_min",
        "c_dur_y",
        "c_dur_ci",
        "c_dur_cpf",
        "c_dur_racc",
        "c_dur_drcm",
        "sum_c_dur",
        "bar_diameter",
        "maximum",
        "c_dev",
        "final_cover"
    ])

    i = 0
    for row in df.iterrows():
        
        single_row_data = row[1:][0]     

        XD, XC, XF = ExposureClasses.get_exposure_classes(
            single_row_data.is_buried,
            single_row_data.road_zone,
            single_row_data.road_type,
            single_row_data.humidity_level,
            single_row_data.element_type
        )

        concrete_covers = NominalCover.get_covers(
            XD, XC, XF,
            single_row_data.concrete_class,
            single_row_data.is_prestressed,
            single_row_data.is_consequence_class_three,
            single_row_data.is_racc,
            single_row_data.is_drcm,
            single_row_data.is_corrosion_inhibitor,
            single_row_data.is_cpf,
            single_row_data.is_stainless_steel
        )

        fixed_output_df.loc[i] =    [XD, XC, XF] + \
                                    [NominalCover.critical_exposure] + \
                                    list(concrete_covers) + \
                                    [0, 0, 0, 0]
        i += 1
    return fixed_output_df
       

class FormatChecker:

    input_columns = [
        'element_name',
        'is_buried',
        'road_zone',
        'road_type',
        'element_type',
        'humidity_level',
        'concrete_class',
        'is_consequence_class_three',
        'is_racc',
        'is_drcm',
        'is_prestressed',
        'is_corrosion_inhibitor',
        'is_cpf',
        'is_stainless_steel'
    ]

    @classmethod
    def is_columns(self, dfcolumns):
        return self.input_columns == dfcolumns.tolist()
    
    @staticmethod
    def is_road_zone(dfcolumn):
        check = dfcolumn.isin(['inner', 'outer']) | pd.isna(dfcolumn)
        return all(check)

    @staticmethod
    def is_road_type(dfcolumn):
        check = dfcolumn.isin(['m', 'a']) | pd.isna(dfcolumn)
        return all(check)

    @staticmethod
    def is_element_type(dfcolumn):
        check = dfcolumn.isin(['horizontal', 'vertical']) 
        return all(check)
    
    @staticmethod
    def is_humidity_level(dfcolumn):
        check = dfcolumn.isin(['cyclic', 'moderate']) | pd.isna(dfcolumn)
        return all(check)

    @staticmethod
    def is_concrete_class(dfcolumn):
        check = dfcolumn.isin(['c28/35', 'c35/45', 'c40/50', 'c50/60']) 
        return all(check)

    @classmethod
    def format(cls, df):
        df.columns = df.columns.str.strip()
        df.columns = df.columns.str.lower()
        df = df.apply(lambda x: x.str.lower() if is_string_dtype(x) else x)
        df = df.apply(lambda x: x.str.strip() if is_string_dtype(x) else x)
    
        check_matrix = [
            cls.is_columns(df.columns),
            is_bool_dtype(df.is_buried),
            cls.is_road_zone(df.road_zone),
            cls.is_road_type(df.road_type),
            cls.is_element_type(df.element_type),
            cls.is_humidity_level(df.humidity_level),
            cls.is_concrete_class(df.concrete_class),
            is_bool_dtype(df.is_consequence_class_three),
            is_bool_dtype(df.is_racc),
            is_bool_dtype(df.is_drcm),
            is_bool_dtype(df.is_prestressed),
            is_bool_dtype(df.is_corrosion_inhibitor),
            is_bool_dtype(df.is_cpf),
            is_bool_dtype(df.is_stainless_steel)
        ]

        if all(check_matrix):
            return df
        else: 
            raise ValueError(
                "Invalid format. Please check", [
                    cls.input_columns[i] for i in range(len(check_matrix)) \
                    if check_matrix[i] == False \
                ] \
            )


class ExposureClasses:

    is_buried = None
    road_zone = None
    road_type = None
    humidity_level = None
    element_type = None

    @classmethod 
    def get_exposure_classes(cls, is_buried, 
                                    road_zone, 
                                    road_type,
                                    humidity_level,
                                    element_type ):
        """
        is_buried:      True/False          \n
        road_zone:      inner/outer/np.nan  \n
        road_type:      "m", "a" or np.nan  \n
        humidity_level: "moderate", "cyclic", np.nan \n
        element_type:   "vertical", "horizontal"     \n
        """
        cls.is_buried = is_buried
        cls.road_zone =  road_zone
        cls.road_type = road_type
        cls.humidity_level = humidity_level
        cls.element_type = element_type
        XD = cls.determine_XD(road_type)
        XC = cls.determine_XC(humidity_level)
        XF = cls.determine_XF(element_type)
        return XD, XC, XF

    @classmethod
    def determine_XD(self, road_type):
        if self.is_buried:
            if self.road_zone == "inner":
                return  "XD3" if road_type == ("m" or "a") else "XD1"
            elif self.road_zone == "outer":
                return  "XD2" if road_type == ("m" or "a") else "XD1"
            elif pd.isna(self.road_zone):
                return np.nan
        else:
            if self.road_zone == "inner":
                return  "XD3" if road_type == ("m" or "a") else "XD1"
            elif self.road_zone == "outer":
                return "XD1"
            elif pd.isna(self.road_zone):
                return np.nan

    @classmethod
    def determine_XC(self, humidity_level):
        if self.is_buried: 
            return "XC2"
        else: 
            if humidity_level == "moderate":
                return "XC3"
            elif humidity_level == "cyclic":
                return "XC4"
            elif pd.isna(humidity_level):
                return "XC1"

    @classmethod
    def determine_XF(self, element_type):
        if self.is_buried:
            return np.nan
        else: 
            if pd.isna(self.road_zone):
                return "XF3" if element_type == "horizontal" else "XF1"
            else: 
                return "XF4" if element_type == "horizontal" else "XF2"


class NominalCover: 

    num_category = {
        "XD1": 4,
        "XD2": 5, 
        "XD3": 5, 
        "XC1": 3, 
        "XC2": 2, 
        "XC3": 1, 
        "XC4": 1, 
        "XF": 3
    }
    critical_exposure = None
    critical_category = None
    concrete_class = None 
    c_dur = None

    @classmethod 
    def get_covers(cls, XD, XC, XF, \
                        concrete_class,\
                        is_prestressed, \
                        is_consequence_class_three, \
                        is_racc,\
                        is_drcm,\
                        is_corrosion_inhibitor, \
                        is_cpf,\
                        is_stainless_steel):

        cls.concrete_class = concrete_class

        c_dur_min = cls.get_c_min_dur(XD,XC,XF)
        c_dur_y = cls.get_c_dur_y(is_prestressed, is_consequence_class_three)
        c_dur_ci = cls.get_c_dur_ci(is_corrosion_inhibitor)
        c_dur_cpf = cls.get_c_dur_cpf(is_cpf)
        c_dur_racc = cls.get_c_dur_racc(is_racc)
        c_dur_drcm = cls.get_c_dur_drcm(is_drcm)
        
        cls.c_dur = c_dur_min + \
                    c_dur_y + \
                    c_dur_ci + \
                    c_dur_cpf + \
                    c_dur_racc + \
                    c_dur_drcm
        
        #if stainless is used, c_dur is overwritten 
        cls.c_dur = cls.get_c_ss(is_stainless_steel, XD)

        return  c_dur_min, \
                c_dur_y, \
                c_dur_ci, \
                c_dur_cpf, \
                c_dur_racc, \
                c_dur_drcm, \
                cls.c_dur

    @classmethod
    def get_c_min_dur(cls, XD, XC, XF):
        #TODO create JSON instead 
        c_min_dur_values = {
        "XD1":  {
                    "c28/35": 50,
                    "c35/45": 40,
                    "c40/50": 35,
                    "c50/60": 30
                },
        "XD2":  {
                    "c28/35": 60,
                    "c35/45": 50,
                    "c40/50": 50,
                    "c50/60": 50
                },
        "XD3":  {
                    "c28/35": 60,
                    "c35/45": 50,
                    "c40/50": 50,
                    "c50/60": 50
                },

        "XC1":  {
                    "c28/35": 20,
                    "c35/45": 20,
                    "c40/50": 20,
                    "c50/60": 20
                },
        "XC2":  {
                    "c28/35": 25,
                    "c35/45": 25,
                    "c40/50": 25,
                    "c50/60": 25
                },
        "XC3":  {
                    "c28/35": 45,
                    "c35/45": 35,
                    "c40/50": 30,
                    "c50/60": 30
                },
        "XC4":  {
                    "c28/35": 45,
                    "c35/45": 35,
                    "c40/50": 30,
                    "c50/60": 30
                },
        "XF":  {
                    "c28/35": 20,
                    "c35/45": 20,
                    "c40/50": 20,
                    "c50/60": 20
                },
        }
        #TODO fix. bad architecture
        try:
            xd_cover_category = (c_min_dur_values[XD][cls.concrete_class], XD)
        except:
            xd_cover_category = tuple()
        try:
            xc_cover_category = (c_min_dur_values[XC][cls.concrete_class], XC)
        except:
            xc_cover_category = tuple()
        try:
            xf_cover_category = (c_min_dur_values[XF[:2]][cls.concrete_class], XF[:2])
        except:
            xf_cover_category = tuple()
        c_min_dur, cls.critical_exposure = max(
            xd_cover_category, 
            xc_cover_category,  
            xf_cover_category
        )
        cls.critical_category =  cls.num_category[cls.critical_exposure]
        return c_min_dur

    @classmethod
    def get_c_dur_y(self, is_prestressed, is_consequence_class_three):
        if is_consequence_class_three:
            c_dur_values = {
                1: 0,
                2: 0,
                3: 0,
                4: 10,
                5: 10
            }
            exposure =  10 if self.critical_exposure == "XC3" else \
                        c_dur_values[self.critical_category]
            climate_change = 5 
            prestress = 10 if is_prestressed else 0
            total_dur_y = exposure+climate_change+prestress
            return total_dur_y
        else: 
            return 0

    @classmethod
    def get_c_dur_ci(cls, is_corrosion_inhibitor):
        if is_corrosion_inhibitor:
            c_dur_corrosion_values = {
                1: -5,
                2: 0,
                3: 0,
                4: 0,
                5: 0
            }
            c_dur_ci = -5 if cls.critical_exposure == "XD3" else \
                                c_dur_corrosion_values[cls.critical_category]
            return c_dur_ci
        else:
            return 0

    @staticmethod
    def get_c_dur_cpf(is_cpf):
        c_dur_cpf = -5 if is_cpf else 0  
        return c_dur_cpf

    #TODO add is_racc method
    @classmethod
    def get_c_dur_racc(cls, is_racc):
        if is_racc:
            c_dur_racc_values = {
                1: -5,
                2: -5,
                3: -5,
                4: 0,
                5: 0
            }
            c_dur_racc = c_dur_racc_values[cls.critical_category]
            return c_dur_racc
        else: 
            return 0

    @classmethod
    def get_c_dur_drcm(cls, is_drcm):
        if is_drcm:
            c_dur_drcm_values = {
                1: 0,
                2: 0,
                3: 0,
                4: -5,
                5: -5
            }
            c_dur_drcm = c_dur_drcm_values[cls.critical_category]
            return c_dur_drcm
        else: 
            return 0

    @classmethod
    def get_c_ss(cls, is_stainless_steel, XD):
        if is_stainless_steel:
            if cls.critical_category == 1 or XD == "XD3":
                return 10  
            else:
                return cls.c_dur
        else:
            return cls.c_dur