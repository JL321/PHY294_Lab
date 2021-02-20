import matplotlib.pyplot as plt
import numpy as np
import os

def main():
    data_folder = '/home/jamesl/Downloads/BlackbodyRadiationData6'
    file_list = os.listdir(data_folder)
    path_dict = {}
    volt_list = ['{}V'.format(i) for i in range(4, 11)]
    spread_list = []
    for volt in volt_list:
        path_dict[volt] = []
    for file in file_list:
        if file[0:2] == '30':
            voltageVal = file[3:6]
            if voltageVal[0:2] != '10':
                voltageVal = voltageVal[0:2]
            path_dict[voltageVal].append(os.path.join(data_folder, file))
        elif file[0:2] == '80':
            spread_list.append(os.path.join(data_folder, file))

    init_list = []
    init_std = []
    for path in spread_list:
        line_counter = 0
        x_list = []
        y_list = []
        with open(path, 'r') as ofile:
            for line in ofile:
                if line_counter > 1:
                    num_list = line.split('\t')
                    x = float(num_list[0])
                    y = float(num_list[1][:-1])
                    x_list.append(x)
                    y_list.append(y)
                line_counter += 1
            ofile.close()
        x_list = np.array(x_list)
        y_list = np.array(y_list)
        scan_idx = 0
        for i, x in enumerate(x_list):
            if x > 60:
                scan_idx = i
                break
        sub_x = x_list[scan_idx:]
        sub_y = y_list[scan_idx:]
        ind = np.argsort(sub_y)
        top_five = sub_x[ind[-5:]]
        avg_init = np.mean(top_five)
        std = np.std(top_five)
        init_list.append(avg_init)
        init_std.append(std)

    init_theta = np.mean(init_list)
    print("Initial theta: {}".format(init_theta))
    print("Initial uncertainty: {}".format(np.sqrt(np.sum(np.square(std)))))
    avg_theta = []
    std_theta = []
    for volt in volt_list:
        for display_path in path_dict[volt]:
            line_counter = 0
            x_list = []
            y_list = []
            with open(display_path, 'r') as ofile:
                for line in ofile:
                    if line_counter > 1:
                        num_list = line.split('\t')
                        x = float(num_list[0])
                        y = float(num_list[1][:-1])
                        x_list.append(x)
                        y_list.append(y)
                    line_counter += 1
            x_list = np.array(x_list)
            y_list = np.array(y_list)
            ind = np.argsort(y_list)
            top_angles = x_list[ind[-10:]]
            avg = np.mean(top_angles)
            std = np.std(top_angles)
            avg_theta.append(avg)
            std_theta.append(std)
            ofile.close()

    theta_list = init_theta-np.array(avg_theta)
    lambda_list = []
    temp_list = []

    #plt.scatter(x_list, y_list, s=2, c='blue')
    #plt.show()

if __name__ == '__main__':
    main()