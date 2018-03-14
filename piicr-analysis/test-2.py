from python_plotter_functions import *
from read_write_functions import *

spot_positions = read_csv('G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-06-Cd-run\\129Cd_PIICR\\analyzed\\auto-run-007-jonas\\129gCd\\p1p2', '129gCd_007_p2_spot_positions')
center = read_csv('G:\\Experiments\\ISOLTRAP\\DATA\\2017\\2017-06-Cd-run\\129Cd_PIICR\\analyzed\\auto-run-007-jonas\\129gCd\\p1p2', '133Cs_c_000_spot_positions')

# spot_positions = read_csv('/Users/jonaskarthein/cernbox/Python/2017-07-02-129Cd-auto-run-007/129gCd/p1p2/', 'spot_positions_127Cd-p2-155')

only_spot_positions = [['X / ch.', 'Y / ch.']]

for i in spot_positions:
    print float(i[1]),float(i[2])
    if (float(i[2]) > 70 and float(i[1]) < -100  and float(i[1]) > -270):
        only_spot_positions.append([float(i[1])-100,float(i[2])+100])
        only_spot_positions.append([float(i[1])-100,-float(i[2])-100])

for i in center:
    only_spot_positions.append([float(i[1])+20,float(i[2])+20])

# python_plot(only_spot_positions, '$^{127_{g&m}}$Cd, IS-574, 2016', '%s_mcp_4' % ('mcp_test'), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, '2dhistogram-mcp', 20, 1, 'on', 'off', 'mcp', [37.17, 4.42, 35.06, 4.49], [-177.16, 5.91, -277.29, 5.79])

# python_plot(only_spot_positions, '$^{127_{g&m}}$Cd, IS-574', '%s_mcp_4' % ('mcp_test'), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, '2dhistogram-mcp', 50, 1, 'on', 'off', 'mcp', [37.17, 4.42, 35.06, 4.49], [-177.16, 5.91, -277.29, 5.79])

python_plot(only_spot_positions, '$^{129}$Cd, PI-ICR excitation patterns', '%s_piicr_excitation' % ('mcp_test'), 'no', '', 'no', 'gauss', 'Utopia', 'red', 'black', 'full', 0, 8, 6, '2dhistogram-mcp', 60, 1, 'on', 'off', 'mcp', [37.17, 4.42, 35.06, 4.49], [-177.16, 5.91, -277.29, 5.79])
