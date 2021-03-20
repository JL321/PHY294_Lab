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

    init_std = np.std(init_list)
    init_theta = np.mean(init_list)
    T_knot = 20+273 #  Accounting for unit conversion for alpha
    R_knot = 1.1
    alpha_knot = 4.5*(10**-3)
    print("Initial theta: {}".format(init_theta))
    print("Initial uncertainty: {}".format(np.sqrt(np.sum(np.square(std)))))
    avg_theta = []

    #  Note that temperature will be in K
    temp_list = []
    temp_uncertainty = []
    intensity_list = []
    intensity_uncertainty = []
    abs_min = -0.38716
    v_u = 0.05
    I_u = 0.0005
    for volt in volt_list:
        for tup in path_dict[volt]:
            display_path, current = tup
            volt_val = int(volt[:-1])
            line_counter = 0
            x_list = []
            y_list = []
            temp_list.append(T_knot+(volt_val/(current*R_knot)-1)/alpha_knot)
            temp_uncertainty.append(np.sqrt(np.square(I_u*volt_val*np.power(current, -2)*np.power(R_knot*alpha_knot, -1))+np.square(v_u*np.power(current*R_knot*alpha_knot, -1))))
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
            #print(np.min(y_list[start_idx:end_idx]))
            peak_y = y_list[start_idx:end_idx]-abs_min
            width_list = []
            midpoint_list = []
            area = 0
            for i in range(len(peak_x)-1):
                width = peak_x[i+1]-peak_x[i]
                width_list.append(width)
                midpoint = (peak_y[i+1]+peak_y[i])/2
                midpoint_list.append(midpoint)
                area += width*midpoint
            #print("Final area: {}".format(area))
            intensity_list.append(area)
            width_list = np.array(width_list)
            midpoint_list = np.array(midpoint_list)
            intensity_uncertainty.append(np.sqrt(np.sum(2*np.square(unit_uncertainty*(width_list))+2*np.square(unit_uncertainty*midpoint_list))))
            
            #  Calculations pertaining to Wien's displacement law
            ind = np.argsort(y_list)
            top_angles = x_list[ind[-10:]]
            avg = np.mean(top_angles)
            avg_theta.append(avg)
            ofile.close()


    std_theta = np.std(np.array(avg_theta))
    theta_u = np.sqrt(np.square(std_theta)+np.square(init_std))
    intensity_uncertainty = np.array(intensity_uncertainty)
    temp_uncertainty = np.array(temp_uncertainty)

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
    print("L List: {}".format(lambda_list))
    print("Temp list: {}".format(temp_list))
    lamb_uncertainty = (1/3)*theta_u*A*(1/lambda_list)*(np.power(np.square((2/3)*np.sin(theta_list)+1)+(3/4), -1/2))*(2*np.sin(theta_list*np.pi/180)/3+1)*np.cos(theta_list*np.pi/180)
    print("Lamb uncertainty: {}".format(lamb_uncertainty))
    print("Temp uncertainty: {}".format(temp_uncertainty))
    w_uncertainty = np.mean(np.sqrt(np.square(lamb_uncertainty*temp_list)+np.square(temp_uncertainty*lambda_list)))*(10**-9)
    expected_const = 2.898*(10**-3)
    empirical_const = lambda_list*(10**-9)*temp_list
    print("Empirical constants: {}".format(empirical_const))
    print("Avg constant: {}".format(np.mean(empirical_const)))
    print("Intensity list: {}".format(intensity_list))
    print("Intensity uncertainties: {}".format(intensity_uncertainty))
    print("Wien's uncertainty: {}".format(w_uncertainty))
    t_4 = np.power(temp_list, 4)
    plt.scatter(t_4, intensity_list, s=2, c='black')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Intensities")
    #plt.show()
    plt.errorbar(t_4, intensity_list, yerr=intensity_uncertainty, fmt='.', c='black')

    reg = LinearRegression().fit(np.expand_dims(t_4, 1), np.expand_dims(intensity_list, 1))
    m = np.squeeze(reg.coef_)
    b = np.squeeze(reg.intercept_)
    chi_square = np.sum(np.square(intensity_list-(m*t_4+b))/intensity_uncertainty)
    reduced_c = chi_square/(intensity_list.shape[0]-2)
    x_space = np.linspace(t_4[0], t_4[-1])
    plt.plot(x_space, x_space*m+b)
    print("Chi square value: {}".format(chi_square))
    print("Reduced c: {}".format(reduced_c))
    print("Slope: {} Theoretical slope: {}".format(m, 5.67*10**-8))
    plt.show()

    residuals = intensity_list - (m*t_4+b)
    plt.scatter(t_4, residuals, s=4, c='red')
    plt.xlabel("Temperature (K)")
    plt.ylabel("Intensity Residuals")
    plt.show()
    
    h = 6.062*10**-34
    c = 3*10**8
    k = 1.38064*10**-23
    new_temp = [1000, 2000, 2250, 2500, 2750, 3000, np.round(np.max(temp_list)), 4000]
    color = ['b', 'g', 'r', '--r', '-.r', ':r', 'b', '--b', "c"]
    for t, col in zip(new_temp, color):
        x = np.linspace(0, 2000, 40)
        x_revised = x*(10**-9)
        i_bb = 2*h*c**2/(x_revised**5*(np.exp(h*c/(x_revised*k*t))-1))
        plt.plot(x, i_bb, col, label=str(t)+"k")

    plt.axvline(380)
    plt.axvline(700)
    plt.xlabel("Wavelength (nanometers)")
    plt.ylabel("Intensity")
    plt.legend(loc="upper left")
    plt.show()

if __name__ == '__main__':
    main()