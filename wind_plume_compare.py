###########################################
#
# WindPlume.py
#
###########################################

from numpy import *
from pandas import *
import time

def WindVector(intense_bin, orient_bin, orient, compass):
    # posit wind velocity magnitude (m/sec)
    if intense_bin == 'CALM':
        min_v = 0.
        max_v = 0.5
    elif intense_bin == '>46':
        min_v = 46.
        max_v = 46.
    else:
        delin = intense_bin.index('-')
        min_v = float(intense_bin[:delin]) - 0.5
        max_v = float(intense_bin[delin+1:]) + 0.5
    v = random.uniform(min_v, max_v) * 0.447   # convert MPH to m/sec
    # posit origin direction (radians)
    if orient_bin != 'CALM':    
        delin = orient.index(orient_bin)
        center_angle = compass[delin]
        min_angle = center_angle - 22.5/2.
        max_angle = center_angle + 22.5/2.
    else:
        min_angle = 0.
        max_angle = 360.    
    angle = random.uniform(min_angle, max_angle)
    if angle < 0.: angle = 360. + angle
    return v, angle*pi/180.                 # convert angle to radians

def Stability(v, day_f, stability):
    # assign atmospheric stability
    if v<2:                                 # wind velocity component
        row = 0
    elif v>=2 and v<3:
        row = 1
    elif v>=3 and v<4:
        row = 2        
    elif v>=4 and v<6:
        row = 3
    else:
        row = 4
    if day_f <= 0.5:                        # solar irradiation component
        col = 0                             # sun high in sky
    else:
        col = 1                             # sun low in sky or cloudy
    return stability[row][col]

def Dispersion(x, stable_class):
    # assign transverse and vertical dispersion coefficients (standard terrain)
    if stable_class == 'A':
        sigma_y = 0.22*x / sqrt(1.0 + 0.0001*x)
        sigma_z = 0.2*x
    elif stable_class == 'B':
        sigma_y = 0.16*x / sqrt(1.0 + 0.0001*x)
        sigma_z = 0.12*x        
    elif stable_class == 'C':
        sigma_y = 0.11*x / sqrt(1.0 + 0.0001*x)
        sigma_z = 0.08*x / sqrt(1.0 + 0.0002*x)   
    elif stable_class == 'D':
        sigma_y = 0.08*x / sqrt(1.0 + 0.0001*x)
        sigma_z = 0.06*x / sqrt(1.0 + 0.0015*x) 
    elif stable_class == 'E':
        sigma_y = 0.06*x / (1.0 + 0.0001*x)
        sigma_z = 0.03*x / (1.0 + 0.0003*x) 
    else:
        sigma_y = 0.04*x / (1.0 + 0.0001*x)
        sigma_z = 0.016*x / (1.0 + 0.0003*x)         
    return sigma_y, sigma_z

def Rotate(x, y, theta):
    # rotate coordinate system about theta
    x_prime = x*cos(theta) - y*sin(theta)
    y_prime = x*sin(theta) + y*cos(theta)
    return x_prime, y_prime

def C(x, y, z, u, H, F, stable_class):
    # calculate time-integrated concentration at x, y, z
    if x > 0:
        sigma_y, sigma_z = Dispersion(x, stable_class)
        f = exp(-(y**2./(2.*sigma_y**2)))
        g1 = exp(-((z - H)**2./(2.*sigma_z**2)))
        g2 = exp(-((z + H)**2./(2.*sigma_z**2)))
        return F/u * f/(sigma_y*sqrt(2.*pi)) * (g1 + g2)/(sigma_z*sqrt(2.*pi))
    else:
        return 0.

### main script ###


def WindPlume(num_trials, d_min, d_max, H, z, Q, log_avg_C0, log_stdev_C0):

    # working parameter sets
    compass = array([90., 67.5, 45., 22.5, 0., 337.5, 315., 292.5, 270., 247.5, 225., 202.5, 180., 157.5, 135., 112.5, 0.])
    stability_table = [['A', 'B', 'F'], ['A', 'C', 'E'], ['B', 'C', 'D'], ['C', 'D', 'D'], ['C', 'D', 'D']]

    # read in meteorology data
    wind_df = read_csv('wind.txt', index_col=0, sep='\t')

    # record bin labels for velocity ranges and direction
    orient = list(wind_df.columns)
    intense = list(wind_df.index)

    # set up wind bins lookup table
    wind_bins_df = wind_df.stack()
    wind_bins_df = wind_bins_df.reset_index()
    wind_bins_df.columns = ['intensity', 'direction', 'freq'] 
    wind_bins_df['cumul'] = wind_bins_df['freq'].cumsum()
    wind_bins_df = wind_bins_df[wind_bins_df['freq'] > 0]
    wind_bins_df = wind_bins_df.reset_index(drop=True)

    # run trials

    # placeholders for results arrays
    theta_out = empty(num_trials, dtype=float)
    wind_origin = empty(num_trials, dtype='|S10')
    x_prime_out = empty(num_trials, dtype=float)
    y_prime_out = empty(num_trials, dtype=float)
    v_out = empty(num_trials, dtype=float)
    s_out = empty(num_trials, dtype='|S10')
    c_out = empty(num_trials, dtype=float)

    # posit source term vector
    log_C0 = random.normal(log_avg_C0, log_stdev_C0, num_trials)
    C0 = 10.**log_C0
    F = Q * C0 

    # posit receptor location vector
    d = random.uniform(d_min, d_max, num_trials)
    psi = random.uniform(0., 2.*pi, num_trials)
    x = d * cos(psi)
    y = d * sin(psi)

    for i in xrange(num_trials):

        # select random numbers to choose wind vector and solar irradiation conditions
        r = random.uniform(0.,1.)
        day_f = random.uniform(0.,1.)
        floor = (r > wind_bins_df['cumul'])
        bins_index = floor.sum() - 1
        bins_index = max(0, bins_index)

        # posit meteorology
        v, theta = WindVector(wind_bins_df['intensity'][bins_index], wind_bins_df['direction'][bins_index], orient, compass)
        if theta>pi: theta_point = theta - pi    # point theta 180 degrees away for orientation of local positive x-axis
        else: theta_point = theta + pi
        x_prime, y_prime = Rotate(x[i], y[i], -theta_point)
        s = Stability(v, day_f, stability_table)

        # run Gaussian plume model
        c = C(x_prime, y_prime, z, v, H, F[i], s)

        # record output from this trial
        wind_origin[i] = wind_bins_df['direction'][bins_index]
        theta_out[i] = theta
        x_prime_out[i] = x_prime
        y_prime_out[i] = y_prime
        v_out[i] = v
        s_out[i] = s
        c_out[i] = c

    # summarize results
    trials_summary = {'x':x, 'y':y, 'theta':theta_out, 'wind_origin':wind_origin,
        'x_prime':x_prime_out, 'y_prime':y_prime_out, 'v':v_out, 's':s_out, 'C':c_out, 'C0':C0}
    results_df = DataFrame(trials_summary)
    results_df.to_csv('results.csv')
    
### run script ###

start = time.clock()

num_trials = 10000
d_min = 5.
d_max = 60.
H = 5.
z = 2.
Q = 2.
log_avg_C0 = 1.0 
log_stdev_C0 = 0.3

WindPlume(num_trials, d_min, d_max, H, z, Q, log_avg_C0, log_stdev_C0)

end = time.clock()
print(end - start)
print 'Done.'






