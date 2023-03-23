
# Written by June Skeeter
# Equations Adapted From:
# https://gml.noaa.gov/grad/solcalc/calcdetails.html
# More rearouses can be found here: https://squarewidget.com/solar-coordinates/

# Meeus, J. (1991). Astronomical algorithms. Richmond.


import numpy as np

class sunPosition():

    def __init__(self,JD,GMT_offset=0):
        self.tz = GMT_offset
        self.JD = JD-self.tz/24
        self.JC = (self.JD-2451545)/36525 # Convert to Julian Century
        self.TIME = (JD-.5)%1 # Local Time - as fraction of day
        self.orb_Pos()

    def orb_Pos(self):
        GM_lon_sun = (280.46646+self.JC*(36000.76983 + self.JC*0.0003032))%360 #Geometric Mean Longitude of the Sun (deg)
        GM_anom_sun = 357.52911+self.JC*(35999.05029 - 0.0001537*self.JC)#Geometric Mean Anomaly of the Sun (deg)
        EC_orb = 0.016708634-self.JC*(0.000042037+0.0000001267*self.JC)#Eccentricity of Earth's Orbit
        eq_Ctr = np.sin(np.radians(GM_anom_sun))*(1.914602-self.JC*(0.004817+0.000014*self.JC))+np.sin(np.radians(2*GM_anom_sun))*(0.019993-0.000101*self.JC)+np.sin(np.radians(3*GM_anom_sun))*0.000289 #Sun Eq of Ctr
        True_lon_sun = GM_lon_sun+eq_Ctr #Sun True Longitude (deg)
        True_anom_sun = GM_anom_sun+eq_Ctr #Sun True Anomaly (deg)
        rad_Vector = (1.000001018*(1-EC_orb*EC_orb))/(1+EC_orb*np.cos(np.radians(True_anom_sun))) #Center to Center Distance from Sun to Earth (AUs)
        AP_lon_sun = True_lon_sun-0.00569-0.00478*np.sin(np.radians(125.04-1934.136*self.JC)) #Apparent Longitude of the Sun (deg)
        Obliquqe = 23+(26+((21.448-self.JC*(46.815+self.JC*(0.00059-self.JC*0.001813))))/60)/60 #Mean Oblique Ecliptic (deg)
        corr_Oblique = Obliquqe+0.00256*np.cos(np.radians(125.04-1934.136*self.JC)) #Corrected Oblique Ecliptic (deg)
        var_y = np.tan(np.radians(corr_Oblique/2))*np.tan(np.radians(corr_Oblique/2))# var y
        self.right_Ascension = np.radians(np.arctan2(np.cos(np.radians(corr_Oblique))*np.sin(np.radians(AP_lon_sun)),np.cos(np.radians(AP_lon_sun)))) #Right Ascension (deg)
        self.Declination = np.degrees(np.arcsin(np.sin(np.radians(corr_Oblique))*np.sin(np.radians(AP_lon_sun)))) #Declination (deg)
        self.eq_Time = 4*np.degrees(
            var_y*np.sin(2*np.radians(GM_lon_sun))-
                                    2*EC_orb*np.sin(np.radians(GM_anom_sun))+4*EC_orb*var_y*np.sin(np.radians(GM_anom_sun))*np.cos(2*np.radians(GM_lon_sun))-
                                    0.5*var_y**2*np.sin(4*np.radians(GM_lon_sun))-
                                    1.25*EC_orb**2*np.sin(2*np.radians(GM_anom_sun))
                                    ) #Equations of Time (minutes)
                
    def calc_for_Coords(self,LAT,LON):
        self.Lat = LAT
        self.Lon = LON
        self.HA_sunrise = np.degrees(np.arccos(np.cos(np.radians(90.833))/(np.cos(np.radians(self.Lat))*np.cos(np.radians(self.Declination)))-np.tan(np.radians(self.Lat))*np.tan(np.radians(self.Declination)))) #HA Sunrise (deg)
        
        self.Noon_LST = (720-4*self.Lon-self.eq_Time+self.tz*60)/1440 #Solar Noon (LST)
        self.Sunrise_LST = (self.Noon_LST*1440-self.HA_sunrise*4)/1440 #Sunrise Time (LST)
        self.Sunset_LST = (self.Noon_LST*1440+self.HA_sunrise*4)/1440 #Sunset Time (LST)
        
        self.Day_Length=8*self.HA_sunrise #Sunlight Duration (minutes)
        self.Solar_Time=(self.TIME*1440+self.eq_Time+4*self.Lon-60*self.tz)%1440 # True Solar Time (min)
        
        #  AB2=(TIME*1440+V2+4*LON-60*TZ)%1440 # True Solar Time (min)
        
        hour_Angle = np.copy(self.Solar_Time)

        hour_Angle[self.Solar_Time/4<0] = self.Solar_Time[self.Solar_Time/4<0]/4+180
        hour_Angle[self.Solar_Time/4>0] = self.Solar_Time[self.Solar_Time/4>0]/4-180
        
        self.Zenith = np.degrees(np.arccos(np.sin(np.radians(self.Lat))*np.sin(np.radians(self.Declination))+np.cos(np.radians(self.Lat))*np.cos(np.radians(self.Declination))*np.cos(np.radians(hour_Angle))))# Solar Zenith Angle (deg)

        Solar_elevation = 90-self.Zenith #Solar Elevation Angle (deg)
        Refraction = np.zeros(Solar_elevation.shape)
        
        Refraction[Solar_elevation<=85] = 58.1/np.tan(np.radians(Solar_elevation[Solar_elevation<85]))-0.07/((np.tan(np.radians(Solar_elevation[Solar_elevation<85])))**3)+0.000086/(np.tan(np.radians(Solar_elevation[Solar_elevation<85]))**5)
        Refraction[Solar_elevation<=5] = 1735+Solar_elevation[Solar_elevation<5]*(-518.2)+(Solar_elevation[Solar_elevation<5]**2)*(103.4)+(Solar_elevation[Solar_elevation<5]**3)*(-12.79)+((Solar_elevation[Solar_elevation<5])**4)*0.711
        Refraction[Solar_elevation<=-0.575] = -20.772/np.tan(np.radians(Solar_elevation[Solar_elevation<=-0.575]))
        
        Refraction/=3600#pprox Atmospheric Refraction (deg)
        self.Solar_elevation = Solar_elevation+Refraction #Solar Elevation corrected for atm refraction (deg)

        
        self.Azimuth = hour_Angle*0 # Solar Azimuth Angle (deg cw from N)
        self.Azimuth[hour_Angle>0] = (np.degrees(
            np.arccos(((np.sin(np.radians(self.Lat))*np.cos(np.radians(self.Zenith[hour_Angle>0])))-
                                  np.sin(np.radians(self.Declination[hour_Angle>0])))/(np.cos(np.radians(self.Lat))*np.sin(np.radians(self.Zenith[hour_Angle>0]))))
                                  )+180)%360
        self.Azimuth[hour_Angle<=0] = (np.degrees(
            np.arccos(((np.sin(np.radians(self.Lat))*np.cos(np.radians(self.Zenith[hour_Angle<=0])))-
                                  np.sin(np.radians(self.Declination[hour_Angle<=0])))/(np.cos(np.radians(self.Lat))*np.sin(np.radians(self.Zenith[hour_Angle<=0]))))
                                  )+180)%360


