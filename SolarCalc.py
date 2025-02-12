import sys
import pytz
import numpy as np
import pandas as pd
from tzfpy import get_tz
from dataclasses import dataclass, field

@dataclass
class TzFuncs():
    DST: bool = False
    UTC: bool = False
    TimeZone: str = None
    Timestamp: object = None
    Latitude: float = None
    Longitude: float = None

    def __post_init__(self):
        if self.TimeZone is None:
            self.findTZ()
        self.convert(self.Timestamp.copy(),self.UTC,not self.UTC)  

    def findTZ(self):
        self.DST=True
        self.TimeZone = pytz.timezone(get_tz(self.Longitude,self.Latitude))
        
    def convert(self,Input_Time,from_UTC=False,to_UTC=False):
        if isinstance(Input_Time,pd.DatetimeIndex):
            Input_Time=Input_Time.to_series()
        if from_UTC == False:
            self.Local_Time=Input_Time
            self.to_StandardTime(to_UTC)
        else:
            self.UTC_Time=Input_Time
            self.fromUTC()

    def to_StandardTime(self,to_UTC=False):
        offset = self.Local_Time.apply(lambda x: self.TimeZone.dst(x,is_dst=self.DST))
        if self.DST == True:
            self.Standard_Time = pd.DatetimeIndex(self.Local_Time-offset)
            self.Local_Time = self.Local_Time.apply(lambda x: self.TimeZone.localize(x,is_dst=self.DST))
        else:
            self.Standard_Time = pd.DatetimeIndex(self.Local_Time)
            self.Local_Time = (self.Local_Time+offset).apply(lambda x: self.TimeZone.localize(x,is_dst=self.DST))
        if to_UTC == True:
            self.toUTC()
        self.Local_Time = pd.DatetimeIndex(self.Local_Time)

    def toUTC(self):
        self.UTC_Time = pd.DatetimeIndex(self.Local_Time.apply(lambda x: x.astimezone(pytz.utc)))

    def fromUTC(self):
        self.UTC_Time = self.UTC_Time.apply(lambda x: pytz.utc.localize(x, is_dst=self.DST))
        self.Local_Time = self.UTC_Time.apply(lambda x: x.astimezone(self.TimeZone).replace(tzinfo=None))
        self.to_StandardTime(to_UTC=False)
        self.UTC_Time = pd.DatetimeIndex(self.UTC_Time)

@dataclass
class SunStats(TzFuncs):
    # Equations Adapted From:
    # https://gml.noaa.gov/grad/solcalc/calcdetails.html
    # More rearouses can be found here: https://squarewidget.com/solar-coordinates/
    # Meeus, J. (1991). Astronomical algorithms. Richmond.
    # Timezone determination done using custom script
    Timestamp: object = None
    Latitude: float = None
    Longitude: float = None
    Slope: float = 0
    Aspect: float = 180
    transmissivity: float = 0.85
    Io: float = 1361.0

    def __post_init__(self):
        if type(self.Timestamp) == list and len(self.Timestamp)==3 and len(self.Timestamp[2])<8:
            self.Timestamp = pd.date_range(self.Timestamp[0],self.Timestamp[1],freq=self.Timestamp[2])
        elif type(self.Timestamp) == list:
            self.Timestamp = pd.to_datetime(self.Timestamp)
        elif type(self.Timestamp) == str:
            self.Timestamp = pd.to_datetime([self.Timestamp])
        super().__post_init__()
        JD = self.UTC_Time.to_julian_date().values
        self.JC = (JD-2451545)/36525 # Convert to Julian Century
        self.TIME = (JD-.5)%1 # Local Time - as fraction of day
        self.orb_Pos()
        self.siteSpecific()
        self.fluxDensity()

    def orb_Pos(self):
        # Get orbital position
        GM_lon_sun = (280.46646+self.JC*(36000.76983 + self.JC*0.0003032))%360 #Geometric Mean Longitude of the Sun (deg)
        GM_anom_sun = 357.52911+self.JC*(35999.05029 - 0.0001537*self.JC)#Geometric Mean Anomaly of the Sun (deg)
        EC_orb = 0.016708634-self.JC*(0.000042037+0.0000001267*self.JC)#Eccentricity of Earth's Orbit
        eq_Ctr = np.sin(np.radians(GM_anom_sun))*(1.914602-self.JC*(0.004817+0.000014*self.JC))+np.sin(np.radians(2*GM_anom_sun))*(0.019993-0.000101*self.JC)+np.sin(np.radians(3*GM_anom_sun))*0.000289 #Sun Eq of Ctr
        True_lon_sun = GM_lon_sun+eq_Ctr #Sun True Longitude (deg)
        True_anom_sun = GM_anom_sun+eq_Ctr #Sun True Anomaly (deg)
        self.rad_Vector = (1.000001018*(1-EC_orb*EC_orb))/(1+EC_orb*np.cos(np.radians(True_anom_sun))) #Center to Center Distance from Sun to Earth (AUs)
        AP_lon_sun = True_lon_sun-0.00569-0.00478*np.sin(np.radians(125.04-1934.136*self.JC)) #Apparent Longitude of the Sun (deg)
        Obliquqe = 23+(26+((21.448-self.JC*(46.815+self.JC*(0.00059-self.JC*0.001813))))/60)/60 #Mean Oblique Ecliptic (deg)
        corr_Oblique = Obliquqe+0.00256*np.cos(np.radians(125.04-1934.136*self.JC)) #Corrected Oblique Ecliptic (deg)
        var_y = np.tan(np.radians(corr_Oblique/2))*np.tan(np.radians(corr_Oblique/2))# var y
        self.right_Ascension = np.radians(np.arctan2(np.cos(np.radians(corr_Oblique))*np.sin(np.radians(AP_lon_sun)),np.cos(np.radians(AP_lon_sun)))) #Right Ascension (deg)
        self.Declination = np.degrees(np.arcsin(np.sin(np.radians(corr_Oblique))*np.sin(np.radians(AP_lon_sun)))) #Declination (deg)
        self.eq_Time = 4*np.degrees(var_y*np.sin(2*np.radians(GM_lon_sun))-2*EC_orb*np.sin(np.radians(GM_anom_sun))+
                                    4*EC_orb*var_y*np.sin(np.radians(GM_anom_sun))*np.cos(2*np.radians(GM_lon_sun))-
                                    0.5*var_y**2*np.sin(4*np.radians(GM_lon_sun))-1.25*EC_orb**2*np.sin(2*np.radians(GM_anom_sun))
                                    ) #Equations of Time (minutes)

    def siteSpecific(self):
        with np.errstate(invalid='ignore'):
            self.HA_sunrise = np.degrees(np.arccos(np.cos(np.radians(90.833))/(np.cos(np.radians(self.Latitude))*np.cos(np.radians(self.Declination)))-np.tan(np.radians(self.Latitude))*np.tan(np.radians(self.Declination)))) #HA Sunrise (deg)
        self.Noon_LST = (720-4*self.Longitude-self.eq_Time+60)/1440 #Solar Noon (LST)
        self.Sunrise_LST = (self.Noon_LST*1440-self.HA_sunrise*4)/1440 #Sunrise Time (LST)
        self.Sunset_LST = (self.Noon_LST*1440+self.HA_sunrise*4)/1440 #Sunset Time (LST)
        self.Day_Length=8*self.HA_sunrise #Sunlight Duration (minutes)
        self.Solar_Time=(self.TIME*1440+self.eq_Time+4*self.Longitude-60)%1440 # True Solar Time (min)
        self.hour_Angle = np.copy(self.Solar_Time)
        self.hour_Angle[self.Solar_Time/4<0] = self.Solar_Time[self.Solar_Time/4<0]/4+180
        self.hour_Angle[self.Solar_Time/4>0] = self.Solar_Time[self.Solar_Time/4>0]/4-180
        self.Zenith = np.degrees(np.arccos(np.sin(np.radians(self.Latitude))*np.sin(np.radians(self.Declination))+np.cos(np.radians(self.Latitude))*np.cos(np.radians(self.Declination))*np.cos(np.radians(self.hour_Angle))))# Solar Zenith Angle (deg)
        Solar_elevation = 90-self.Zenith #Solar Elevation Angle (deg)
        Refraction = np.zeros(Solar_elevation.shape)
        Refraction[Solar_elevation<=85] = 58.1/np.tan(np.radians(Solar_elevation[Solar_elevation<85]))-0.07/((np.tan(np.radians(Solar_elevation[Solar_elevation<85])))**3)+0.000086/(np.tan(np.radians(Solar_elevation[Solar_elevation<85]))**5)
        Refraction[Solar_elevation<=5] = 1735+Solar_elevation[Solar_elevation<5]*(-518.2)+(Solar_elevation[Solar_elevation<5]**2)*(103.4)+(Solar_elevation[Solar_elevation<5]**3)*(-12.79)+((Solar_elevation[Solar_elevation<5])**4)*0.711
        Refraction[Solar_elevation<=-0.575] = -20.772/np.tan(np.radians(Solar_elevation[Solar_elevation<=-0.575]))
        Refraction/=3600#pprox Atmospheric Refraction (deg)
        self.Solar_elevation = Solar_elevation+Refraction #Solar Elevation corrected for atm refraction (deg)
        self.Azimuth = self.hour_Angle*0 # Solar Azimuth Angle (deg cw from N)
        self.Azimuth[self.hour_Angle>0] = (np.degrees(
            np.arccos(((np.sin(np.radians(self.Latitude))*np.cos(np.radians(self.Zenith[self.hour_Angle>0])))-
                                  np.sin(np.radians(self.Declination[self.hour_Angle>0])))/(np.cos(np.radians(self.Latitude))*np.sin(np.radians(self.Zenith[self.hour_Angle>0]))))
                                  )+180)%360
        self.Azimuth[self.hour_Angle<=0] = (540-np.degrees(
            np.arccos(((np.sin(np.radians(self.Latitude))*np.cos(np.radians(self.Zenith[self.hour_Angle<=0])))-
                                  np.sin(np.radians(self.Declination[self.hour_Angle<=0])))/(np.cos(np.radians(self.Latitude))*np.sin(np.radians(self.Zenith[self.hour_Angle<=0]))))
                                  ))%360

    def fluxDensity(self):
        # Get incident angle on surface, amounting for slope and bearing
        # Use that to estimate direct irradiance
        # This will ignore diffuse and reflected irradiance, so will give an underestimate of total irradiance

        alpha = np.radians(self.Slope)
        zeta = np.radians(self.Zenith)
        a = np.radians(self.Azimuth)
        b = np.radians(self.Aspect)
        cosi = np.cos(alpha)*np.cos(zeta)+np.sin(alpha)*np.sin(zeta)*np.cos(a-b)
        self.incidence = np.degrees(np.arccos(cosi))
        self.I = self.Io*self.rad_Vector*np.cos(np.arccos(cosi))
        self.I[self.I<0] = 0
        m = (1/(np.cos(np.radians(self.Zenith))))
        m[m<0]=np.nan
        self.SW_in = self.I*self.transmissivity**m
        self.SW_in[np.isnan(self.SW_in)]=0


@dataclass
class SolarArray(SunStats):
    # Solar Panel Specs
    Panel_Watts: float = 160 # Factory Rating
    Panel_area: float = 1.482*0.674
    Panel_Voltage: float = 21
    Charger_Voltage: float = 14.5
    Panel_per_Sequence: int = 1 # Can wire in sequence first to increase voltage - set to 1 = not in sequence
    Panel_Sequences_in_Parallel: int = 3 #	Can wire in parallel to increase amperage
    Degredation_coeff: float = .005 # Assume a continuous % drop / year
    Panel_Age: int = 0
    Slope: float = 0
    Bearing: float = 180

    # Charger Specs
    Charger_Max_Current:float = 60 # Amps
    Charging_Efficiency:float = 0.8 # Assume charging is 80% efficient

    def __post_init__(self):
        super().__post_init__()
        self.N_Panels = self.Panel_per_Sequence*self.Panel_Sequences_in_Parallel
        self.Panel_Array_Voltage = self.Panel_Voltage*self.Panel_per_Sequence
        self.Panel_Watts = self.Panel_Watts*(1-(self.Degredation_coeff*self.Panel_Age))
        self.Panel_Efficiency = self.Panel_Watts/self.Panel_area/1000
        self.total_area = self.Panel_area*self.N_Panels
        self.Panel_Max_Output_coeff = self.Panel_Efficiency*self.Panel_area * self.Panel_per_Sequence * self.Panel_Sequences_in_Parallel


@dataclass
class BatteryBank(SolarArray):
    Batt_Voltage:float = 12
    Batt_Amp_Hours:float = 50	
    Batt_Weight_KG:float = 18.3 # Estimated from here - https://www.powerstream.com/Size_SLA.htm
    Batt_per_Sequence:int = 1 # Can wire in sequence first to increase voltage - set to 1 = not in sequence
    Batt_Sequences_in_Parallel:int = 6 #	Can wire in parallel to increase amperage
    Topping_Charge:float = 0.75
    # Rough estimate of relation between charge pct and voltage for AGM batteries
    # from: https://bluettipower.co.uk/blogs/news/full-battery-voltage-chart
    # AGM table, reduced to a per-volt basis
    Voltage_Curve: np.array = field(default_factory=lambda:np.array([[1.0,1.071],[0.99,1.067],[0.9,1.063],[0.8,1.042],[0.7,1.025],[0.6,1.013],[0.5,1.004],[0.4,0.996],[0.3,0.984],[0.2,0.972],[0.1,0.959],[0.0,0.921]]))

    def __post_init__(self):
        super().__post_init__()
        self.N_Batt = self.Batt_Sequences_in_Parallel*self.Batt_per_Sequence
        self.Total_Batt_Weight = self.Batt_Weight_KG*self.N_Batt
        self.Supply_Voltage = self.Batt_Voltage*self.Batt_per_Sequence
        self.Voltage_Curve[:,1]*= self.Supply_Voltage
        self.Supply_Amp_Hours_Max = self.Batt_Amp_Hours*self.Batt_Sequences_in_Parallel
        self.Supply_Amp_Hours = np.zeros(self.Timestamp.shape)
        self.Supply_Amp_Hours[0] = self.Supply_Amp_Hours_Max
        self.Supply_pct_Available = self.Supply_Amp_Hours/self.Supply_Amp_Hours[0]
        self.Topping_Decay_Rate = -self.Charger_Max_Current/(1-self.Topping_Charge)
        self.Voltage_Curve = np.poly1d(np.polyfit(self.Voltage_Curve[:,0],self.Voltage_Curve[:,1],3))
        self.Voltage = self.Voltage_Curve(self.Supply_pct_Available)


@dataclass
class System(BatteryBank):
    Full_Load_Watts: float = 1.0
    Base_Load_Watts: float = 0.1
    Discharge_Limit: float = 0.5
    Restart_Threshold: float = 0.9
    System_On: int = 1
    System_Reset_Cycles = 24

    def __post_init__(self):
        if type(self.Full_Load_Watts) is list:
            self.Full_Load_Watts = sum(self.Full_Load_Watts)
        elif type(self.Full_Load_Watts) is np.ndarray:
            self.Full_Load_Watts = self.Full_Load_Watts.sum()
        self.Demand = self.Full_Load_Watts#.copy()
        super().__post_init__()

        # Charging        
        # Assume linear decrease in charge efficiency  

        # Output of Panels (Watts pre panel x number of panels x efficiency)
        self.Max_Input_Watts=self.SW_in*self.Panel_Max_Output_coeff
        self.Amps_In = self.Max_Input_Watts/self.Panel_Array_Voltage
        self.Charge_Potential_Watts = self.Amps_In*self.Charger_Voltage
        self.Charge_Balance_Watts = self.Charge_Potential_Watts*0

        if np.nanmax(self.Amps_In) > self.Charger_Max_Current:
            print('Waring - Current Input May Exceed Capacity of Charger!!')

    def TrackState(self):
        for i,t in enumerate(self.Timestamp):
            self.Supply_pct_Available[i] = self.Supply_Amp_Hours[i]/self.Supply_Amp_Hours_Max
            # Scale for topping charge with assumed linear decrease
            Charge_Rate_adj = min(1,(self.Supply_pct_Available[i]*self.Topping_Decay_Rate-self.Topping_Decay_Rate)/self.Charger_Max_Current)*self.Charging_Efficiency
            self.Charge_Balance_Watts[i] = self.Charge_Potential_Watts[i]*Charge_Rate_adj-self.Demand
            if self.Charge_Balance_Watts[0]>0:
                self.Voltage[i] = self.Charger_Voltage
            else:
                self.Voltage[i] = self.Voltage_Curve(self.Supply_pct_Available[i])
            if i < len(self.Timestamp)-1:
                self.Supply_Amp_Hours[i+1] = self.Supply_Amp_Hours[i]+self.Charge_Balance_Watts[i]/self.Voltage[i]
            if self.Supply_pct_Available[i]<self.Discharge_Limit:
                self.Demand = self.Base_Load_Watts
                self.System_On = -1*self.System_Reset_Cycles
            elif self.Supply_pct_Available[i]>self.Restart_Threshold and self.System_On<=0:
                self.System_On += 1
            elif self.System_On == 1:
                self.Demand = self.Full_Load_Watts



            # print(self.Supply_pct_Available[i])
            # print(t)


            
        # # Calculate energy output from charger after voltage conversion
        # self.Watts_Out = self.Amps_In*self.Charger_Voltage*self.Charging_Efficiency
        
        # self.Watts_In_Out = self.Watts_Out-self.Full_Load_Watts
        # self.Amps_Available = np.ones(self.Watts_In_Out.shape)*self.Supply_Amp_Hours
        # self.Amps = (self.Watts_In_Out/self.Voltage).cumsum()
        # for i,a in enumerate(self.Amps):
        #     self.Amps_Available[i] = min(self.Amps_Available[i]+a,self.Supply_Amp_Hours)
        
        # self.PowerSupply = pd.DataFrame(
        #     index = self.Timestamp,
        #     data = {'Watts_In_Out':self.Watts_In_Out,
        #             'Amps_Available':self.Amps_Available}
        # )