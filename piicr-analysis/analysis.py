# ---------------------------------------------------------------------------
# Written by Jonas Karthein in 2016/2017. Questions to jonas.karthein@cern.ch
# ---------------------------------------------------------------------------

from analysis_functions import *
import platform

# ---------------------------------------------------------------------------
# Please change before use to the correct folder -- Instructions in REAM-ME-first
# ---------------------------------------------------------------------------

upper_folder_path = '/Users/jonaskarthein/cernbox/Python/test/2017-05-15-85Rb-87Rb-85Rb_1000ms_008/' # 'G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\2017-6-test' # '/Users/jonaskarthein/cernbox/Python/2017-05-15-85Rb-87Rb-85Rb_1000ms_007_root/' # Windoof:'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-05_ArKr preperation\\2017-05_piicr_cross-checks\\2017-05-02_50ms_ion-number-tests\\87Rb' # Mac:'/Volumes/Macintosh HD/Users/jonaskarthein/cernbox/Python/2017-05-02_50ms/87Rb'
isotopes = ['85Rb', '87Rb']

# ---------------------------------------------------------------------------
# Run skript
# ---------------------------------------------------------------------------

for i in isotopes:
	if platform.system() == 'Windows':
		analysis_main('%s\\%s' % (upper_folder_path, str(i)))
	elif platform.system() == 'Darwin': # MacOS
		analysis_main('%s%s' % (upper_folder_path, str(i)))
