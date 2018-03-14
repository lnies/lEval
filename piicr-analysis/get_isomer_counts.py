import numpy as np
import pandas as pd
import sys, os
import glob
import json#simplejson as json
# np.sum(df[(df['run'] == 1)]['#ground'])


def get_isomer_counts(upper_z_folder, isotopes):
    # function to merge all count info and calculate isomer/ground ratio
    merge_counts_dict = {}
    merge_counts_list = [['isotope', 'z-class', 'run', 'filename', '#ground', '#isomer', 'contamination%', 'cyc_acc_time']]
    for isotope in isotopes:
        merge_counts_dict[isotope] = {}
        for z in [os.path.basename(z[:-1]) for z in glob.glob(upper_z_folder+'*/') if 'Z0-0' not in z and '/Z' in z]:
            print z
            merge_counts_dict[isotope][z] = {}
            for run in [os.path.basename(run[:-1]) for run in glob.glob(upper_z_folder+z+'/*/')]:
                merge_counts_dict[isotope][z][run] = {}
                filename = '{}{}/{}/{}/p1p2/counts_info_dict.json'.format(upper_z_folder, z, run, isotope)
                times_filename = '{}{}/{}/{}/freq_config/piicr_excitation.csv'.format(upper_z_folder, z, run, isotope)
                try:

                    with open(filename) as json_data:
                        tmp_dict = json.load(json_data)
                    for file in tmp_dict.keys():
                        cyc_acc_time = pd.read_csv(times_filename, delimiter=',')['Cyclotron accumulation time / microseconds'][0]
                        merge_counts_dict[isotope][z][run][file] = tmp_dict[file]
                        merge_counts_list.append([isotope, z, int(run[-3:]), file,
                                                 tmp_dict[file]['#manual-spots-ground'],
                                                 tmp_dict[file]['#manual-spots-isomer'],
                                                 tmp_dict[file]['percent-contamination']*100, cyc_acc_time])
                except:
                    print('File in isotope={} , Z={}, run={} not found!'.format(isotope, z, run))

    save_filename = '{}merge_counts_dict.json'.format(upper_z_folder)
    with open(save_filename, 'w') as f:
        json.dump(merge_counts_dict, f, sort_keys=True, indent=4)
    # calculation of rates
    df = pd.DataFrame(merge_counts_list, columns=merge_counts_list.pop(0))
    df['#ground-unc'] = np.sqrt(df['#ground'])
    df['#isomer-unc'] = np.sqrt(df['#isomer'])
    df['ground/isomer'] = df['#ground']/df['#isomer']
    df['ground/isomer-unc'] = np.sqrt( (df['#ground-unc']/df['#isomer'])**2 + (df['#ground']*df['#isomer-unc']/df['#isomer']**2)**2 )
    # back-extrapolation for 129Cd
    half_life_ground = [154, 2] # time and unc. in ms
    half_life_isomer = [146, 8]
    isoltrap_time_dict = {4:222.6,      # whole time from ISCOOL to UT w/o t_acc
                      5:222.6,
                      6:322.6,
                      7:222.6,
                      8:182.6,
                      9:182.6,
                      10:182.6,
                      11:202.6,
                      12:202.6,
                      13:202.6,
                      14:130.1,
                      15:130.1}

    df['isoltrap_time'] = df.apply(lambda x: isoltrap_time_dict[x['run']], axis=1)
    df['#ground-ISCOOL'] = df['#ground'] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0])
    df['#ground-ISCOOL-unc'] = np.sqrt(
                                       (df['#ground-unc'] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       + (df['#ground'] * 2 / half_life_ground[0] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       + (df['#ground'] * df['isoltrap_time'] * half_life_ground[1] / half_life_ground[0]**2 * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       )
    df['#isomer-ISCOOL'] = df['#isomer'] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_isomer[0])
    df['#isomer-ISCOOL-unc'] = np.sqrt(
                                       (df['#isomer-unc'] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       + (df['#isomer'] * 2 / half_life_ground[0] * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       + (df['#isomer'] * df['isoltrap_time'] * half_life_ground[1] / half_life_ground[0]**2 * np.exp((df['isoltrap_time']+df['cyc_acc_time']/1000) / half_life_ground[0]))**2
                                       )

    print df.head

# summing counts
    total_g = np.sum(df['#ground'])
    print total_g
    total_g_unc = np.sqrt(total_g)
    print total_g_unc

    total_m = np.sum(df['#isomer'])
    total_m_unc = np.sqrt(total_m)
    total_ratio = float(total_g)/float(total_m)
    total_ratio_unc = np.sqrt( (total_g_unc/total_m)**2
                              + (-total_g*total_m_unc/total_m/total_m)**2)
    print 'Ratio ISOLTRAP :', total_ratio, '+/-', total_ratio_unc

# summing extrapolated counts
    total_g_iscool = np.sum(df['#ground-ISCOOL'])
    total_g_iscool_unc = np.sum(np.sqrt(df['#ground-ISCOOL-unc']**2))

    total_m_iscool = np.sum(df['#isomer-ISCOOL'])
    total_m_iscool_unc = np.sum(np.sqrt(df['#isomer-ISCOOL-unc']**2))

    total_ratio_iscool = float(total_g_iscool)/float(total_m_iscool)
    total_ratio_iscool_unc = np.sqrt( (total_g_iscool_unc/total_m_iscool)**2
                              + (-total_g_iscool*total_m_iscool_unc/total_m_iscool/total_m_iscool)**2)
    print 'Ratio ISOLTRAP@ISCOOL :', total_ratio_iscool, '+/-', total_ratio_iscool_unc





if __name__ == '__main__':
    # upper_z_folder = '/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/129Cd/'
    # isotopes = ['129mCd', '129gCd']
    upper_z_folder = '/Volumes/Karthein_Jonas/Doktor/2017-11-02_Cd-analysis/PI-ICR/127Cd/'
    isotopes = ['127mCd', '127gCd']
    get_isomer_counts(upper_z_folder, isotopes)
