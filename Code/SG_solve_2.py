# Import needed libraries
from astropy.coordinates import SkyCoord
from shutil import copy
import threading
from queue import Empty
from multiprocessing import Queue
import multiprocessing as mp
import os
import time
import pickle
import sqlite3
from watchdog.events import PatternMatchingEventHandler
from watchdog.observers import Observer
from photutils.centroids import centroid_sources, centroid_com
from astropy.coordinates import AltAz
from astropy.coordinates import EarthLocation
from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import NoConvergence
from astropy import units as u
import circle_fit
import re
import numpy as np
import subprocess
import pandas as pd
from datetime import datetime
import os.path
from multiprocessing import Pool
import warnings
from astropy.utils.exceptions import AstropyUserWarning, AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)


def check_initial_conditions(event_path):
    """ Check whether tracking is ON, otherwise skip solving that particular image """
    if ("tracking-0" in os.path.splitext(event_path)[0] or
          "tracking=0" in os.path.splitext(event_path)[0]):
        print("Tracking off, skipping: ", event_path)
        return 0
    else:
        return 1


def check_LEDs_ON(event_path):
    """ Check if the LEDs are ON, otherwise solve only starfield """
    if ("LED-0" in os.path.splitext(event_path)[0]) or ("LED=0" in os.path.splitext(event_path)[0]):
        print("LEDs turned off, only solving starfield: ", event_path)
        return 0
    else:
        return 1
    

def check_if_parked(event_path):
    """ Check if the telescope is parked, if so, skip solving the image """
    if ("parked=1" in os.path.splitext(event_path)[0]) or ("parked-1" in os.path.splitext(event_path)[0]):
        print("Telescope is parked, skip solving!", event_path)
        return 0
    else:
        return 1


def get_earth_location_of_a_site(earth_location):
    """ In case of network connection error, it loads site information from the local storage.
    Not needed anymore, soon to be removed. """
    try:
        ORM = EarthLocation.of_site(earth_location)
    except:
        with open("/home/ccddev/Starguider/Code/ORM_backup", "rb") as f:
            ORM = pickle.load(f)
    return ORM

def get_telescope_coordinates():
    """ Get LST1 coordinates """
    # LST1 location at LaPalma
    telescope_coordinates = {"longitude" : 342.108612,
                             "lattitude" : 28.761389,
                             "elevation" : 2147}
    lst1_location = EarthLocation(
        lon=telescope_coordinates["longitude"] * u.deg,
        lat=telescope_coordinates["lattitude"] * u.deg,
        height=telescope_coordinates["elevation"] * u.m)
    return lst1_location
    

def copy_file(event_path):
    """ Create a safe copy of SG image in Backup folder in case something goes wrong when solving."""
    os.makedirs(os.path.dirname(event_path) + "/Backup", exist_ok=True)
    copy(event_path, os.path.dirname(event_path) + "/Backup")


def check_header_names(event_path):
    """ Get some data from the header of SG image"""
    num_of_tries = 5
    for i in range(num_of_tries):
        try:
            hdul = fits.open(event_path)
            header = hdul[0].header
            break
        except OSError:
            time.sleep(1)
            continue
    try:
        linux_time = hdul[0].header["TIME"]
    except:
        try:
            linux_time = hdul[0].header["UNIXTIME"]
        except:
            pass

    try:
        zenith = header["ZENITH"]
    except:
        try:
            zenith = header["ELEVA"]
        except:
            pass

    try:
        rotation = header["ROTATION"]
    except:
        rotation = "False"

    return hdul, header, linux_time, zenith, rotation


def rotate_image(event_path):
    """ Rotate SG image before solving """
    with fits.open(event_path, mode='update') as hdul:
        hdul[0].data = hdul[0].data[::-1, ::-1]  # Ako okrecem sliku
        hdul[0].header["ROTATION"] = "True"
        hdul.flush()


def crop_image_for_astrometry(event_path):
    """ Discard part of SG image containing PMT camera, so astrometry will solve only part of the image where stars are contained."""
    with fits.open(event_path, mode='update') as hdul:
        height, width = hdul[0].data.shape
        hdul[0].data = hdul[0].data[:, :1268]
        hdul.flush()


def run_astrometry_plate_solver(event_path):
    """ Running astrometry for solving starfield """
    x_cen = 964.5
    y_cen = 726.5

    command = [
        "/usr/local/astrometry/bin/solve-field",
        #         "--sigma", str(100),
        "--overwrite",
        "--scale-units", "degwidth",
        "--scale-low", str(6.5), "--scale-high", str(6.62),
        str(event_path), "--dir", os.path.dirname(event_path),
        "--rdls", "none", 
        "--temp-axy",
        "--match", "none", "--index-xyls", "none",
        "--no-plots",
        "--new-fits", "none",
        "--tweak-order", "4",
        "--pixel-error", "0.2",
        "--crpix-x", str(x_cen), "--crpix-y", str(y_cen),
        "--continue",
        "--skip-solved",
        "--uniformize", str(0),
        "--no-remove-lines",
        # "--config", str(config_file)
    ]
    result = subprocess.run(command, stdout=subprocess.PIPE)
    return result


def run_astrometry_plate_solver_verify(event_path):
    """ Running astrometry for solving starfield with verify option """
    x_cen = 964.5
    y_cen = 726.5

    command = [
        "/usr/local/astrometry/bin/solve-field",
#         "--sigma", str(200),
        "--overwrite",
        "--scale-units","degwidth",
        "--scale-low", str(6.5), "--scale-high", str(6.62),
        str(event_path), "--dir", os.path.dirname(event_path),
        "--rdls", "none", "--temp-axy",
        "--match", "none", "--index-xyls", "none",
        "--no-plots",
        "--new-fits", "none",
        "--tweak-order", "4",
        "--pixel-error", "0.2",
        "--crpix-x", str(x_cen), "--crpix-y", str(y_cen),
        "--continue",
#         "--config", str(config_file),
        "--no-remove-lines",
        "--uniformize", str(0),
        "--verify", str(event_path)
    ]

    result = subprocess.run(command, stdout=subprocess.PIPE)
    return result


def save_astrometry_terminal_output(event_path, result):
    """ Save astrometry output from terminal in Output_astrometry.txt file """
    if (os.path.exists(os.path.dirname(event_path) + "/Output_astrometry.txt") == False):
        with open(os.path.dirname(event_path) + "/Output_astrometry.txt", 'w') as f:
            f.write(result.stdout.decode('utf-8'))
    else:
        with open(os.path.dirname(event_path) + "/Output_astrometry.txt", 'a') as f:
            f.write(result.stdout.decode('utf-8'))


def check_if_solved_successfully(event_path):
    """ Check if the image is solved successfully i.e. whether .corr exissts, if not write empty row in database and .csv file"""
    now = datetime.utcnow()
    if not (os.path.exists(os.path.splitext(event_path)[0] + ".corr")): 
        print("{0} -- event {1} SOLVED NOK".format(
                now.strftime("%Y/%m/%d %H:%M:%S"), event_path))
        return 0
    else:
        print("{0} -- event {1} SOLVED OK".format(
                now.strftime("%Y/%m/%d %H:%M:%S"), event_path))
        return 1


# Mask sector around each of the 6 LEDs
def sector_mask(shape, centre, radius, angle_range):
    """ Mask sector around each of the 6 LEDs 
    Return a boolean mask for a circular sector. The start/stop angles in  
    `angle_range` should be given in clockwise order.
    """

    x, y = np.ogrid[:shape[0], :shape[1]]
    cx, cy = centre
    tmin, tmax = np.deg2rad(angle_range)

    # ensure stop angle > start angle
    if tmax < tmin:
        tmax += 2*np.pi

    # convert cartesian --> polar coordinates
    r2 = (x-cx)*(x-cx) + (y-cy)*(y-cy)
    theta = np.arctan2(x-cx, y-cy) - tmin

    # wrap angles between 0 and 2*pi
    theta %= (2*np.pi)
    # circular mask
    circmask = r2 <= radius*radius
    # angular mask
    anglemask = theta <= (tmax-tmin)

    return circmask*anglemask


def solve_leds(event_path, w, target_source, target_X, target_Y, drive_source, scale):
    """ Solve the other part of the image(the one with the Cherenkov camera)
       1. Find centroid of each LED, 
       2. Fit a circle on the centroids, 
       3. Find a xy coordinates of the center of the fitted circle
       4. Shift the center due to SG shift from the center of the mirror carrier
       4. Convert xy coords of the center to RADec 
       5. Calculate misspointing of the center of LEDs from the nominal position of the telescope """
    # ORM = get_earth_location_of_a_site('Roque de los Muchachos')
    square_size = (45, 45)
    footprint = np.zeros(square_size, dtype=bool)
    mask = sector_mask(footprint.shape, (22, 22), 22.5, (0, 360))
    footprint[~mask] = True
    x_init = (418, 211, 79, 84, 223, 435)
    y_init = (75, 208, 446, 713, 945, 1068)

    hdul,_, linux_time,_,_ = check_header_names(event_path)
    image = hdul[0].data
    image = image[:, 1368:]
    time = datetime.utcfromtimestamp(linux_time)
    lst1_location = get_telescope_coordinates()
    aa = AltAz(obstime=time,location=lst1_location, obswl=0.35*u.micron, relative_humidity=0.5,temperature=10*u.deg_C,pressure=790*u.hPa)

    x, y = centroid_sources(image, x_init, y_init,
                            footprint=footprint, centroid_func=centroid_com)
    
    xc_1, yc_1, R_1 = circle_fit.fit_circle2(x, y)
    # x1, x2 = circle_fit.create_circle(xc_1, yc_1, R_1)
#     ax.plot(x1, x2, label = "Fitted circle")
    ra = np.tan(np.deg2rad(R_1*scale/3600))*28*1000
    # ra_arcsec = R_1*scale
    
    # LED center before fix
    
    LED_centerX = xc_1 + 1368
    LED_centerY = yc_1
    LED_center = w.all_pix2world(LED_centerX, LED_centerY, 1)
    LED_center_RADec = SkyCoord(ra=LED_center[0], dec=LED_center[1], unit="deg", frame='icrs')
    LED_center_aa = LED_center_RADec.transform_to(aa)

    shift_x_mm = 90
    shift_y_mm = 496
    # Fixing LEDs
    shift_x_pix = (3600/scale) * np.rad2deg(np.arctan(shift_x_mm/28000))
    shift_y_pix = (3600/scale) * np.rad2deg(np.arctan(shift_y_mm/28000))

    LED_centerX_fix = xc_1 + 1368 - shift_x_pix
    LED_centerY_fix = yc_1 + shift_y_pix

    LED_center_fix = w.all_pix2world(LED_centerX_fix, LED_centerY_fix, 1)
    LED_center_fix_RADec = SkyCoord(ra=LED_center_fix[0], dec=LED_center_fix[1], unit="deg", frame='icrs')
    LED_center_fix_aa = LED_center_fix_RADec.transform_to(aa)

    # Convert drive RADec coordinates to pixels
    try:
        drive_X, drive_Y = w.all_world2pix(drive_source.ra.value, drive_source.dec.value, 
                                              1, detect_divergence=True, quiet=False, maxiter=30)
    except NoConvergence as e:
        print("Warning: Astropy all_world2pix failed to converge within specified accuracy/number of iterations."+ 
              " Accepting best solution!")
        drive_X = e.best_solution[0][0]
        drive_Y = e.best_solution[0][1]
        
    target_source_RADec = SkyCoord(
        target_source.ra.value, target_source.dec.value, frame="icrs", unit="deg")
    target_source_aa = target_source_RADec.transform_to(aa)

    # Difference between nominal and unfixed
    dx_unfixed = np.tan(np.deg2rad(
        np.abs(LED_centerX - target_X.item())*scale/3600))*28000
    dy_unfixed = np.tan(np.deg2rad(
        np.abs(LED_centerY - target_Y.item())*scale/3600))*28000
    
    # Difference between drive and unfixed
    dx_unfixed_drive = np.tan(np.deg2rad(
        np.abs(LED_centerX - drive_X.item())*scale/3600))*28000
    dy_unfixed_drive = np.tan(np.deg2rad(
        np.abs(LED_centerY - drive_Y.item())*scale/3600))*28000

    #     LST misspointing (drive and nominal) before shift
    dRADec_LST_drive_LED_center = drive_source.separation(LED_center_RADec).arcmin
    dRADec_LST_sim_LED_center = target_source_RADec.separation(
        LED_center_RADec).arcmin
    

    # Difference between nominal and fixed
    dx_fixed = np.tan(np.deg2rad(
        np.abs(LED_centerX_fix - target_X.item())*scale/3600))*28000
    dy_fixed = np.tan(np.deg2rad(
        np.abs(LED_centerY_fix - target_Y.item())*scale/3600))*28000
    
    # Difference between nominal and fixed
    dx_fixed_drive = np.tan(np.deg2rad(
        np.abs(LED_centerX_fix - drive_X.item())*scale/3600))*28000
    dy_fixed_drive = np.tan(np.deg2rad(
        np.abs(LED_centerY_fix - drive_Y.item())*scale/3600))*28000

    # LST misspointing after correction
    dRADec_LST_sim_LED_center_fix = target_source_RADec.separation(
        LED_center_fix_RADec).arcmin
    dRADec_LST_drive_LED_center_fix = drive_source.separation(
        LED_center_fix_RADec).arcmin

    temp_dict = {
        **{"LED_centerX[px]": LED_centerX,
           "LED_centerY[px]": LED_centerY,
           "LED_center_RA[deg]": LED_center_RADec.ra.value,
           "LED_center_Dec[deg]": LED_center_RADec.dec.value,
           "LED_center_Alt[deg]": LED_center_aa.alt.value,
           "LED_center_Az[deg]": LED_center_aa.az.value,
           "Radius[px]": R_1,
           "Radius[mm]": ra,
           "LST_drive_X[px]": drive_X.item(),
           "LST_drive_Y[px]": drive_Y.item(),
           "LED_centerX_fix[px]": LED_centerX_fix,
           "LED_centerY_fix[px]": LED_centerY_fix,
           "LST_sim_RA[deg]": target_source_RADec.ra.value,
           "LST_sim_Dec[deg]": target_source_RADec.dec.value,
           "LST_sim_Alt[deg]": target_source_aa.alt.value,
           "LST_sim_Az[deg]": target_source_aa.az.value,
           "LED_center_fix_RA[deg]": LED_center_fix_RADec.ra.value,
           "LED_center_fix_Dec[deg]": LED_center_fix_RADec.dec.value,
           "LED_center_fix_Alt[deg]": LED_center_fix_aa.alt.value,
           "LED_center_fix_Az[deg]": LED_center_fix_aa.az.value,
           "dx_unfixed[mm]": dx_unfixed,
           "dy_unfixed[mm]": dy_unfixed,
           "dx_fixed[mm]": dx_fixed,
           "dy_fixed[mm]": dy_fixed,
           "dx_unfixed_drive[mm]": dx_unfixed_drive,
           "dy_unfixed_drive[mm]": dy_unfixed_drive,
           "dx_fixed_drive[mm]": dx_fixed_drive,
           "dy_fixed_drive[mm]": dy_fixed_drive,
           "dRADec_LST_sim_LED_center[arcmin]": dRADec_LST_sim_LED_center,
           "dRADec_LST_sim_LED_center_fix[arcmin]": dRADec_LST_sim_LED_center_fix,
          "dRADec_LST_drive_LED_center[arcmin]": dRADec_LST_drive_LED_center,
           "dRADec_LST_drive_LED_center_fix[arcmin]": dRADec_LST_drive_LED_center_fix}
    }
    return temp_dict


def save_data_of_solved_images_csv(list_of_dfs, src_path):
    """ Save Dataframe with aquired and analysed data of the SG image in the csv file, 
        in the directory of the particular day of observing """
    for lst in list_of_dfs:
        df = pd.DataFrame(lst)
        day = df.iloc[0].Datetime.day
        month = df.iloc[0].Datetime.month
        year = df.iloc[0].Datetime.year
        if df.iloc[0].Datetime.time().hour <= 8:
            if day == 1 and month in [5, 7, 10, 12]:
                day = 30
                month = month - 1
            elif day == 1 and month == 3:
                day = 28
                month = 2
            elif day == 1 and month in [2, 4, 6, 8, 9, 11]:
                day = 31
                month = month - 1
            elif day == 1 and month == 1:
                day = 31
                month = 12
                year = year - 1
            else:
                day = day - 1
        if not(
            os.path.exists(os.path.dirname(src_path) + "/" +
            '%d' % year + "/" + 
            '%02d' % month + "/" + 
            '%02d' % day + "/" + "SG_solved.csv")):
            df.to_csv(
                os.path.dirname(src_path) + "/" +
                '%d' % year + "/" +
                '%02d' % month + "/" +
                '%02d' % day + "/" +
                "SG_solved.csv", header=True, index=False)
        else:
            df.to_csv(
                os.path.dirname(src_path) + "/" +
                '%d' % year + "/" +
                '%02d' % month + "/" +
                '%02d' % day + "/"
                "SG_solved.csv", mode='a', header=False, index=False)
    print("Output in .csv saved successfully!")


def save_data_of_solved_images(list_of_dfs, src_path):
    """ Save Dataframe with aquired and analysed data of the SG image in the database """
    # Create your connection.
    #cnx = sqlite3.connect("/local/home/ccddev/SG_plots/SG_solved.db")
    ####cnx = sqlite3.connect("/project/SG_astrometry/Data/SG_solved.db")
    cnx = sqlite3.connect("/home/ccddev/Starguider/Data/SG_solved.db")
    # cur = con.cursor()
    for lst in list_of_dfs:
        df = pd.DataFrame(lst)
        df.to_sql("SG_solved", con=cnx, if_exists='append', index=False)
    # Save (commit) the changes
    cnx.commit()

    # We can also close the connection if we are done with it.
    # Just be sure any changes have been committed or they will be lost.
    cnx.close()
    print("Database updated successfully!")


def starfield_solve(event):
    """ Main function for solving SG images """
    # ORM = get_earth_location_of_a_site('Roque de los Muchachos')
    temp_list = []
    event_path = event.src_path
    # event_path = event
    columns_from_astrometry = [
        "Scale","Sep_Mean[arcsec]", "Sep_STD[arcsec]", "Sep_SE[arcsec]", "#Matched_stars", "Avg_brightness", 
        "LST_drive_RA[deg]", "LST_drive_Dec[deg]", "LST_drive_Alt[deg]", "LST_drive_Az[deg]",
        "SG_astr_RA[deg]", "SG_astr_Dec[deg]", "SG_astr_Alt[deg]", "SG_astr_Az[deg]", "LST_sim_X[px]","LST_sim_Y[px]"
    ]
    columns_from_leds = [
        "LED_centerX[px]", "LED_centerY[px]", "LED_center_RA[deg]", "LED_center_Dec[deg]",
        "LED_center_Alt[deg]", "LED_center_Az[deg]", "Radius[px]", "Radius[mm]","LST_drive_X[px]", "LST_drive_Y[px]",
        "LED_centerX_fix[px]", "LED_centerY_fix[px]",
        "LST_sim_RA[deg]", "LST_sim_Dec[deg]", "LST_sim_Alt[deg]", "LST_sim_Az[deg]",
        "LED_center_fix_RA[deg]", "LED_center_fix_Dec[deg]", "LED_center_fix_Alt[deg]", "LED_center_fix_Az[deg]",
        "dx_unfixed[mm]", "dy_unfixed[mm]", "dx_fixed[mm]", "dy_fixed[mm]",
        "dx_unfixed_drive[mm]", "dy_unfixed_drive[mm]","dx_fixed_drive[mm]","dy_fixed_drive[mm]",
        "dRADec_LST_sim_LED_center[arcmin]", "dRADec_LST_sim_LED_center_fix[arcmin]",
        "dRADec_LST_drive_LED_center[arcmin]","dRADec_LST_drive_LED_center_fix[arcmin]"]

    # Backup fits files
    copy_file(event_path)

    hdul, header, linux_time, _, _ = check_header_names(event_path)
    # if verify == False:
    rotate_image(event_path)
#     if rotation == "False":
#         rotate_image(event_path)
#     hdul, header, linux_time, zenith, rotation = check_header_names(event_path)
    saved_header = hdul[0].header
    # Prije croppanja zbog rekonstrukcije slike poslije
    saved_data = hdul[0].data

    crop_image_for_astrometry(event_path)

    LST_drive_RA = saved_header["RA_LST"]
    LST_drive_Dec = saved_header["DEC_LST"]
    time = datetime.utcfromtimestamp(linux_time)

    avg_b = hdul[0].data.sum() / (saved_header["NAXIS1"]
                                  * saved_header["NAXIS2"])
    
    lst1_location = get_telescope_coordinates()
    aa = AltAz(obstime=time,location=lst1_location, obswl=0.35*u.micron, relative_humidity=0.5,temperature=10*u.deg_C,pressure=790*u.hPa)

    target_source = SkyCoord(ra=saved_header["RA_TRGT"], dec=saved_header["DEC_TRGT"],
                             unit="deg", frame='icrs')
    
    ###################################################
    drive_source = SkyCoord(ra=saved_header["RA_LST"], dec=saved_header["DEC_LST"], unit="deg", frame="icrs")
    altaz_from_radec_from_header = drive_source.transform_to(aa)

    try:
        # Start astrometry.net plate solving
        result = run_astrometry_plate_solver(event_path)
        #   Output in .txt file
        save_astrometry_terminal_output(event_path, result)
    except :
        print("Error with astrometry.net solving, skipping image")

    LEDs = 1 if check_LEDs_ON(event_path) else 0

    # Check if solved successfully i.e. whether .wcs exists, if not write empty row
    if not check_if_solved_successfully(event_path):
        return ({**{"Datetime": [time], "Unixtime": [linux_time], "Successful": 0, "LEDs" : int(LEDs)},
                 **{col: None for col in columns_from_astrometry}, **{col: None for col in columns_from_leds}})

    table = Table(fits.getdata(os.path.splitext(event_path)[0] + ".corr"))
    hdul.close()

    hdul = fits.open(os.path.splitext(event_path)[0] + ".wcs")  # open wcs
    header = hdul[0].header
    header["IMAGEW"] = 1928

    del header["HISTORY"]
    for item in saved_header:
        if item != "NAXIS1" and item != "NAXIS2" and item not in header:
            header.set(item, saved_header[item])
        # fits.writeto(event, saved_data, header, overwrite=True)
    fits.update(event_path, saved_data, header)
    hdul.close()

    w = WCS(header)
    SG_astr = w.all_pix2world(header["CRPIX1"], header["CRPIX2"], 1)
    SG_astr = SkyCoord(ra=SG_astr[0], dec=SG_astr[1], unit="deg", frame='icrs')
    
    try:
        target_X, target_Y = w.all_world2pix(target_source.ra.value, target_source.dec.value, 
                                              1, detect_divergence=True, quiet=False, maxiter=30)
    except NoConvergence as e:
        print("Warning: Astropy all_world2pix failed to converge within specified accuracy/number of iterations."+ 
              " Accepting best solution!")
        target_X = e.best_solution[0][0]
        target_Y = e.best_solution[0][1]
        
        
    skycoord_field = SkyCoord(ra=table['field_ra'],
                              dec=table['field_dec'],
                              unit="deg", frame="icrs")
    skycoord_index = SkyCoord(ra=table["index_ra"],
                              dec=table["index_dec"],
                              unit="deg", frame="icrs")
    separation = skycoord_index.separation(skycoord_field)
    table["sep"] = separation.arcsec

    SG_astr_aa = SG_astr.transform_to(aa)
    string = header["COMMENT"][-8]
    scale = float(re.findall("\d+\.\d+", string)[0])

    temp = {
        "Datetime": time,
        "Unixtime": linux_time,
        "Successful": 1,
        "LEDs": int(LEDs),
        "Scale" : scale,
        "Sep_Mean[arcsec]": table["sep"].mean(),
        "Sep_STD[arcsec]": table["sep"].std(ddof=1),
        "Sep_SE[arcsec]": table['sep'].std(ddof=1) / (np.sqrt(len(table))),
        "#Matched_stars": len(table),
        "Avg_brightness": avg_b,
        "LST_drive_RA[deg]": LST_drive_RA,
        "LST_drive_Dec[deg]": LST_drive_Dec,
        "LST_drive_Alt[deg]": altaz_from_radec_from_header.alt.value,
        "LST_drive_Az[deg]": altaz_from_radec_from_header.az.value,
        "SG_astr_RA[deg]": SG_astr.ra.value,
        "SG_astr_Dec[deg]": SG_astr.dec.value,
        "SG_astr_Alt[deg]": SG_astr_aa.alt.value,
        "SG_astr_Az[deg]": SG_astr_aa.az.value,
        "LST_sim_X[px]": target_X.item(),
        "LST_sim_Y[px]": target_Y.item(),
    }
    # Call function to solve other part of the image - LEDs
    if LEDs == 0:
        temp_list.append({**temp,**{col: None for col in columns_from_leds}})
    else:
        led_return = solve_leds(event_path, w, target_source, target_X, target_Y, drive_source, scale)
        temp_list.append({**temp, **led_return})
    return temp_list


class FileLoaderWatchdog(PatternMatchingEventHandler):
    '''
    Watchdog module(open-source python API library) is used to monitor filesystem for any changes
    (i.e. creation of a new SG image) to the given directory.
    - watchdog.observers.Observer is a class that will watch for a change and then dispatch the event to specified handler.
    - watchdog.events.PatternMatchingEventHandler is the class that will take the event dispatched by the observer and 
    perform specified action - run script for solving image and determining LST misspointing.
    '''

    def __init__(self, queue, patterns):
        PatternMatchingEventHandler.__init__(self, patterns=patterns)
        self.queue = queue

    def process(self, event):
        '''
        event.event_type
            'modified' | 'created' | 'moved' | 'deleted'
        event.is_directory
            True | False
        event.src_path
            path/to/observed/file
        '''
        self.queue.put(event)

    def on_created(self, event):
        """ The backup directory keeps a safe copy of SG images in case something goes wrong when solving.
            Therefore, files in it will not trigger watchdog and will not be analysed.
            Also, images with tracking OFF will be skipped. """
        if ("Backup" not in event.src_path and check_initial_conditions(event.src_path) and check_if_parked(event.src_path)):
            self.process(event)
            now = datetime.utcnow()
            # print("{0} -- event {1} off the queue ...".format(
            #     now.strftime("%Y/%m/%d %H:%M:%S"), event.src_path))


def queue_get_all(q):
    """ Get images in queue """
    items = []
    maxItemsToRetrieve = 10
    for numOfItemsRetrieved in range(0, maxItemsToRetrieve):
        try:
            if numOfItemsRetrieved == maxItemsToRetrieve:
                break
            items.append(q.get_nowait())
        except Empty:
            break
    return items


def process_load_queue(q, path_watch):
    """ This is the worker thread function. It is run as a daemon 
       threads that only exit when the main thread ends.

       Args
       ==========
         q:  Queue() object
    """
    PROCESSES = mp.cpu_count() - 1
    while True:
        if not q.empty():
            # mp.set_start_method('spawn')

            with Pool(processes=PROCESSES) as pool:
                events = queue_get_all(q)
                multiple_results = [pool.apply_async(
                    starfield_solve, (event,)) for event in events]
                output = [res.get() for res in multiple_results if res.get()]
                save_data_of_solved_images(output, path_watch)
                save_data_of_solved_images_csv(output, path_watch)
        else:
            time.sleep(1)



def main():
    """ Main function for starting watchdog """
    # create queue
    watchdog_queue = Queue()

    # Set up a worker thread to process database load

    # setup watchdog to monitor directory for trigger files
    # defining the patterns attribute to watch only for files with .fits extension.
    pattern = ["*.fits"]
    
    # path_watch is /fefs/onsite/data/aux/lst1/cdm/SG_images/Sorted/ 
    # but it is mounted to /home/ccddev/Starguider/Data/
    
    ###path_watch = "/project/SG_astrometry/Data/"
    path_watch = "/home/ccddev/Starguider/Data/"
    # path_watch = "/home/toni/projects/astrometry/Test_2.4/SG_images/"
    event_handler = FileLoaderWatchdog(watchdog_queue, patterns=pattern)
    # create an instance of Observer that will monitor the specified folder for any events.
    observer = Observer()
    observer.schedule(event_handler, path=path_watch, recursive=True)
    observer.start()
    worker = threading.Thread(
        target=process_load_queue, args=(watchdog_queue, path_watch))

    worker.setDaemon(True)
    worker.start()
    print("SG_solve script has started!\n")

    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()


if __name__ == '__main__':
    main()
