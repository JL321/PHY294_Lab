import matplotlib.pyplot as plt
import numpy as np
import os
from sklearn.linear_model import LinearRegression

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
            current = float(file.split('-')[2][:-1]) 
            path_dict[voltageVal].append([os.path.join(data_folder, file), current])
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
    T_knot = 20+273 #  Accounting for unit conversion for alpha
    R_knot = 1.1
    alpha_knot = 4.5*(10**-3)
    print("Initial theta: {}".format(init_theta))
    print("Initial uncertainty: {}".format(np.sqrt(np.sum(np.square(std)))))
    avg_theta = []
    std_theta = []
    
    #  Note that temperature will be in K
    temp_list = []
    intensity_list = []
    intensity_uncertainty = []
    abs_min = -0.38716
    for volt in volt_list:
        for tup in path_dict[volt]:
            display_path, current = tup
            volt_val = int(volt[:-1])
            line_counter = 0
            x_list = []
            y_list = []
            temp_list.append(T_knot+(volt_val/(current*R_knot)-1)/alpha_knot)
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
            start_idx = 0
            end_idx = 0
            #  Calculations pertaining to Stefan-Boltzmann law
            unit_uncertainty = 0.000005
            for i, x in enumerate(x_list):
                if x >= 10 and start_idx == 0:
                    start_idx = i
                if x >= 26 and end_idx == 0:
                    end_idx = i
            
            peak_x = x_list[start_idx:end_idx]
            # Subtract min to account for negative intensity
            print(np.min(y_list[start_idx:end_idx]))
            peak_y = y_list[start_idx:end_idx]-abs_min
            width_list = []
            area = 0
            for i in range(len(peak_x)-1):
                width = peak_x[i+1]-peak_x[i]
                width_list.append(width)
                midpoint = (peak_y[i+1]+peak_y[i])/2

                area += width*midpoint
            print("Final area: {}".format(area))
            intensity_list.append(area)
            width_list = np.array(width_list)
            intensity_uncertainty.append(np.sqrt(unit_uncertainty*(np.square(np.sum(width_list)))))
            
            #  Calculations pertaining to Wien's displacement law
            ind = np.argsort(y_list)
            top_angles = x_list[ind[-10:]]
            avg = np.mean(top_angles)
            std = np.std(top_angles)
            avg_theta.append(avg)
            std_theta.append(std)
            ofile.close()

    intensity_list = np.array(intensity_list)
    temp_list = np.array(temp_list)
    theta_list = init_theta-np.array(avg_theta)
    A = 13900
    B = 1.689
        
    def get_lambda(theta):
        denom = np.sqrt(np.square(2*np.sin(theta*np.pi/180)/np.sqrt(3)+0.5)+.75)-B
        A_dup = np.ones(denom.shape)*A
        lamb = np.sqrt(A_dup/denom)
        return lamb 
    
    #  Peak wavelength - in nm
    lambda_list = get_lambda(theta_list)

    expected_const = 2.898*(10**-3)
    empirical_const = lambda_list*(10**-9)*temp_list
    print("Empirical constants: {}".format(empirical_const))
    print("Avg constant: {}".format(np.mean(empirical_const)))
    print("Intensity uncertainties: {}".format(intensity_uncertainty))
    t_4 = np.power(temp_list, 4)
    plt.scatter(t_4, intensity_list, s=2, c='blue')
    #plt.show()
    
    reg = LinearRegression().fit(np.expand_dims(t_4, 1), np.expand_dims(intensity_list, 1))
    m = np.squeeze(reg.coef_)
    b = np.squeeze(reg.intercept_)
    x_space = np.linspace(t_4[0], t_4[-1])
    plt.plot(x_space, x_space*m+b)
    print("Slope: {} Theoretical slope: {}".format(m, 1.714*10**-9))
    plt.show()
    

if __name__ == '__main__':
    main()