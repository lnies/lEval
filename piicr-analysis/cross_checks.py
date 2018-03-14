# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

from cross_checks_functions import *

# ---------------------------------------------------------------------------
# Please change before use to the correct folder -- Instructions in REAM-ME-first
# ---------------------------------------------------------------------------

upper_folder_path = '/Volumes/dfs/DATA/2017/2017-06-Cd-run_prep/Cross-checks/test'
reference_isotope = '87Rb'
measurement_isotope = '87Rb'

#upper_folder_path_reference_ion = '/Users/jonaskarthein/cernbox/Python/2017-05-14_85Rb_88Sr_88Rb_300ms_005_root/88Rb' # Windoof:'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr preperation\\2017-05_piicr_cross-checks\\2017-05-02_50ms_ion-number-tests\\87Rb' # Mac:'/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/2017-05-02_50ms/87Rb'
#upper_folder_path_measurement = '/Users/jonaskarthein/cernbox/Python/2017-05-14_85Rb_88Sr_88Rb_300ms_005_root/88Sr' # Windoof:'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr preperation\\2017-05_piicr_cross-checks\\2017-05-02_50ms_ion-number-tests\\87Rb' # Mac:'/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/2017-05-02_50ms/87Rb'
#upper_folder_path_cross_checks = '/Users/jonaskarthein/cernbox/Python/2017-05-14_85Rb_88Sr_88Rb_300ms_005_root/cross-checks/88Rb-88Sr-88Rb' # Windoof: 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr preperation\\2017-05_piicr_cross-checks\\2017-05-02_50ms_meas_number_test\\cross-checks\\85Rb-87Rb-85Rb' #MAC:'/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/2017-05-02_50ms/cross-checks/85Rb-87Rb-85Rb'


# ---------------------------------------------------------------------------
# Run skript
# ---------------------------------------------------------------------------

cross_checks(upper_folder_path, reference_isotope, measurement_isotope)
