import ugradio
import snap_spec
from ugradio import sdr 
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.coordinates import SkyCoord
import time
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
import glob
import argparse 
parser = argparse.ArgumentParser()
parser.add_argument('--filename', '-n', help='name to give file')
parser.add_argument('--LO', '-l', help='value of LO')
parser.add_argument('--gain', '-g', help='gain of sdrs')

#parser to name files more conviniently
args = parser.parse_args()
file = args.filename
lo = float(args.LO)
g = float(args.gain)
#record_time = float(args.record_time)

#set the LO 

#initialize the sdr object
sdr0 = ugradio.sdr.SDR(device_index=0, direct = False, center_freq = lo, sample_rate = 3.2e6, gain = g)
sdr1 = ugradio.sdr.SDR(device_index=1, direct = False, center_freq = lo, sample_rate = 3.2e6, gain = g)

#general process maybe
#convert galactic lingitudes to topocentric coordinates
#point the leusch 
#take the data 
#fft and square and average to get power spectrum 
#time.sleep(440)   tell the guy to chill for 440 secs idk cos max slew time is 220 secs which is rlly slow
#and then start pointing again ig 

# lat, lon and alt for leushcner
lat = 37.9183
lon = -122.1067
alt = 304. # m

#to get power spec 
def fft(signal):
	return np.fft.fft(signal)
	
def perform_power(signal):
	return np.abs((signal))**2
	
def shift(signal):
	return np.fft.fftshift(signal)
	
# -----------------------------------------
# -----------------------------------------
freqs = np.fft.fftshift(np.fft.fftfreq(2048, 1/3.2))

def power(collection_run): 
    collection_run = collection_run - np.mean(collection_run)
    
    
    signal_agg = collection_run[...,0]+1j*collection_run[...,1]
    
    
    #signal_agg_final = signal_agg - np.mean(signal_agg)
    pwr = shift(np.mean(perform_power(fft(signal_agg)), axis=0))
    
    return pwr


#now pls lets start 
gal_lat = 0 
start_lon = -10 
end_lon = 250  #everything is in galactic coords 

#declare the dish object now
dish = ugradio.leusch.LeuschTelescope()

#this is an array of all of the different longitudes in galactic coords
g_lons = np.arange(-10, 252, 2)
s_galactics = [SkyCoord(l= i, b=0, frame = 'galactic', unit='deg') for i in g_lons]
s_topos = [s.transform_to('icrs') for s in s_galactics]

#now lets try pointing and stuff 

#lets start our list of data values ig
point =0 
flops = {}

try: 
    for i in np.arange(len(s_topos)):
        point += 1 
	long = point*2 -10
        jd = ugradio.timing.julian_date()
        alt,az = ugradio.coord.get_altaz(ra = s_topos[i].ra, dec= s_topos[i].dec, jd = jd, lat = lat, lon = lon, alt = alt)
        if (15 < alt <85) and (5 < az < 350):
            #now we can point the big boi to the given alt az 
            dish.point(alt, az)
            print("I do be pointing bro")
            #lets get the data and then make them power specs 
            data = ugradio.sdr.capture_data([sdr0, sdr1], 32768, 300)

            # check gains maybe 

            first = power(data[0])
            second = power(data[1])

            #now lets try to save the data 

            np.savez(f'{file}point{point}', pol0 = first, pol1 = second, ra = s_topos[i].ra, dec = s_topos[i].dec,long = long, alt=alt, az=az, date = jd, missed = flops)
        else:
            flops.update({point:[alt,az]})
            continue   
except KeyboardInterrupt:
     pass 
# l is longitiude and b is latitude 
