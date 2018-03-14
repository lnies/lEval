import glob
from read_write_functions import *
import matplotlib as mpl
import matplotlib.pyplot as plt

SIMION_file =0
axial_phase_shift =0
axial_phase_shift_step= 1 #microseconds step width
length = 0
endx = 0
endy = 0
angle = 0
x_line = 0
y_line = 0
angle_step = 7.2
line_counter = 48


def ion_cloud_3d(spot_positions, axis_scaling_factor, dont_show, file_name, number_of_clouds):
    """
    Nice 3D ioncloud plot.

    The function gives a nice combination between the 2D MCP plot and the
    time of flight information. It therefore reconstructs a picture of the
    3D ioncloud.
    """
    mpl.rc('font', family='serif', serif='Utopia')      # Utopia LaTeX font!!
    mpl.rc('text', usetex=False)

    if number_of_clouds == 1:
                    #Changed for SIMION csv files
        if SIMION_file ==0:
            xss = [float(x[1])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
            yss = [float(x[2])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
            zss = [float(x[3])*1 for x in spot_positions]   # in mikroseconds
        else:
            xss = [float(x[4])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
            yss = [float(x[3])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
            zss = [float(x[1])*1 for x in spot_positions]   # in mikroseconds
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xss, yss, zss, zdir='xss', s=2, c='green')
        ax.set_xlabel('X direction / mm', fontsize=22)
        ax.set_ylabel('Y direction / mm', fontsize=22)
        ax.set_zlabel('TOF / $\mu$s', fontsize=22)
        ax.set_xlim([-20, 20])
        ax.set_ylim([-20, 20])

        ax.view_init(15, 5)
    else:
        os.chdir(folder_path)
        csv_files_2 = []
        for file in glob.glob("2017-07-18_85Rb-MagAmp_2350mV-MagPhase_*_spot_positions.csv"):        # searches for all ...polar.csv files
            csv_files_2.append(file)
        n = 0
        fig = plt.figure(figsize=(12,12))

        ax = fig.add_subplot(111, projection='3d')
        # Drawing a line every xy angle
        angle_line_no = 0
        for i in range(0,line_counter):
            x_line = 5
            y_line = 1
            angle_line_no = angle_line_no/axis_scaling_factor+angle_step/axis_scaling_factor
            endy = 1 +length * math.sin(math.radians(angle))
            endx = length * math.cos(math.radians(angle))
            ax.plot([x_line, endx], [y_line, endy],c='black')



        m=-(axial_phase_shift_step) #correction factor for axial phase shift
        for i in csv_files_2:

            m=m+axial_phase_shift_step
            if n == 1: #this became the color counter 1 means 2 colors, 2 means 3 colors and so on after the counter is reset.
                n=0
            else:
                n+=1
            colors = ['green', 'red', 'blue', 'black', 'magenta', 'gold', 'lightseagreen', 'darkviolet', 'violet']
            #colors = ['green', 'red', 'blue', 'black']
            color = colors[n]
            #n+=1

            import_list_2 = read_csv(folder_path, os.path.splitext(i)[0])
            #Changed for SIMION csv files
            if SIMION_file ==0:
                xss = [float(x[1])*axis_scaling_factor for x in import_list_2]  # in nanoseconds
                yss = [float(x[2])*axis_scaling_factor for x in import_list_2]  # in nanoseconds
                zss = [float(x[3])*1 for x in import_list_2]   # in mikroseconds
            else:
                xss = [float(x[4])*axis_scaling_factor for x in import_list_2]  # # in nanoseconds
                yss = [float(x[3])*axis_scaling_factor for x in import_list_2]  # in nanoseconds
                if axial_phase_shift ==0:
                    zss = [float(x[1])*1 for x in import_list_2]   # in mikroseconds
                else:
                    zss = [(float(x[1])*1)-float(m) for x in import_list_2]   # in mikroseconds
            ax.scatter(xss, yss, zss, zdir='xss', s=2, c=color)
            ax.set_xlabel('X direction / mm', fontsize=22)
            ax.set_ylabel('Y direction / mm', fontsize=22)
            ax.set_zlabel('TOF / $\mu$s', fontsize=22)
            if SIMION_file ==0:
                ax.set_xlim([-10, 10])      #usually it is -20,20
                ax.set_ylim([-10, 10])
                ax.set_zlim([49.0, 51.75])  #for the ToF cut
            else:
                ax.set_xlim([-20, 20])      #usually it is -20,20
                ax.set_ylim([-20, 20])
                ax.set_zlim([300, 301.5])  #for the ToF cut
                #ax.set_xlim([-0.1, 0.1])      #usually it is -20,20
                #ax.set_ylim([0.5, 0.7])
                #ax.set_zlim([-0.0001, 0.0001])  #for the ToF cut


            ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0)) #gives plane colors
            ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
            ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

            #ax.view_init(0, 24)
            #ax.view_init(42, 54)
            ax.view_init(-90, 90)
            #ax.view_init(0, 0)





    #     plt.draw()
    #     plt.pause(.001)
    # for angle in range(0, 360):
    #     ax.view_init(28, angle)
    #     plt.draw()
    #     plt.pause(.001)


    plt.savefig('%s_ioncloud.pdf' % file_name)

    if dont_show == 0:
        plt.show()
    elif dont_show == 1:
        plt.close()


if __name__ == '__main__':
    folder_path = 'E:\\Welker\\Mag_Phase_Scan' # 'E:\\Welker\\Mag_Phase_Scan G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-06-Cd-run\\129Cd_PIICR\\analyzed\\auto-run-004-jonas\\129Cs\\p1p2 G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\harm-new' or G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\Master\\harm-new\\test2
    os.chdir(folder_path)

    number_of_clouds = 2            # Please type the amount of to be collected .csv files here
    if SIMION_file==0:
        axis_scaling_factor = 0.031746032
    else:
        axis_scaling_factor = 1
    dont_show = 0           # if =1 the plots are not displayed, for showing =0! file_name = '133Cs_c_091'  # if only one file should be looked at, uncomment! if peak_finder_plot = 1 this number will automatically set to 0 since they don't work together.
    import_list = []
    csv_files = []
    for file in glob.glob("Harmonic_E_field_r_0_to_1_0mm*.csv"):        # searches for all ...polar.csv files
        csv_files.append(file)
    if number_of_clouds == 1:
        for i in csv_files:
            import_list = read_csv(folder_path, os.path.splitext(i)[0])

            ion_cloud_3d(import_list, axis_scaling_factor, dont_show, i, number_of_clouds)
    else:
        ion_cloud_3d(import_list, axis_scaling_factor, dont_show, 'ion_cloud3d', number_of_clouds)
