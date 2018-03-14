import glob
from read_write_functions import *
import matplotlib as mpl
import matplotlib.pyplot as plt

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
        xss = [float(x[1])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
        yss = [float(x[2])*axis_scaling_factor for x in spot_positions]  # in nanoseconds
        zss = [float(x[3])*0.000025 for x in spot_positions]   # in mikroseconds
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(xss, yss, zss, zdir='xss', s=2, c='green')
        ax.set_xlabel('X direction / mm', fontsize=22)
        ax.set_ylabel('Y direction / mm', fontsize=22)
        ax.set_zlabel('TOF / $\mu$s', fontsize=22)
        ax.set_xlim([-20, 20])
        ax.set_ylim([-20, 20])

        ax.view_init(0, 0)
    else:
        os.chdir(folder_path)
        csv_files_2 = []
        for file in glob.glob("*positions.csv"):        # searches for all ...polar.csv files
            csv_files_2.append(file)
        n = 0
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(111, projection='3d')
        for i in csv_files_2:
            colors = ['green', 'red', 'blue', 'black']
            color = colors[n]
            n+=1
            import_list_2 = read_csv(folder_path, os.path.splitext(i)[0])
            xss = [float(x[1])*axis_scaling_factor for x in import_list_2]  # in nanoseconds
            yss = [float(x[2])*axis_scaling_factor for x in import_list_2]  # in nanoseconds
            zss = [float(x[3])*0.000025 for x in import_list_2]   # in mikroseconds
            ax.scatter(xss, yss, zss, zdir='xss', s=2, c=color)
            ax.set_xlabel('X direction / mm', fontsize=22)
            ax.set_ylabel('Y direction / mm', fontsize=22)
            ax.set_zlabel('TOF / $\mu$s', fontsize=22)
            ax.set_xlim([-20, 20])
            ax.set_ylim([-20, 20])

            ax.view_init(0, 0)





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
    folder_path = 'G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-06-Cd-run\\129Cd_PIICR\\analyzed\\auto-run-004-jonas\\129Cs\\p1p2' # 'G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\harm-new' or 129Cs_005_p2_mcp G:\\Experiments\\ISOLTRAP\\USER\\Karthein_Jonas\\Master\\harm-new\\test2
    os.chdir(folder_path)

    number_of_clouds = 2

    axis_scaling_factor = 0.031746032
    dont_show = 0           # if =1 the plots are not displayed, for showing =0! file_name = '133Cs_c_091'  # if only one file should be looked at, uncomment! if peak_finder_plot = 1 this number will automatically set to 0 since they don't work together.
    import_list = []
    csv_files = []
    for file in glob.glob("129Cs_005_p*_spot_.csv"):        # searches for all ...polar.csv files
        csv_files.append(file)
    if number_of_clouds == 1:
        for i in csv_files:
            import_list = read_csv(folder_path, os.path.splitext(i)[0])

            ion_cloud_3d(import_list, axis_scaling_factor, dont_show, i, number_of_clouds)
    else:
        ion_cloud_3d(import_list, axis_scaling_factor, dont_show, 'ion_cloud3d', number_of_clouds)
