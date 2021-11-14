#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 22:22:02 2021

@author: Yuxin Chen, hfan15
"""

import assignment3_utilities as pa3
import cismath as cis 
from registration import registration

first_body_filename = 'Problem3-BodyA.txt'
N_A, a_list, a_tip =pa3.read_body(first_body_filename) 


second_body_filename = 'Problem3-BodyB.txt'
N_B, b_list, b_tip =pa3.read_body(second_body_filename)

first_reading_file = 'PA3-A-Debug-SampleReadingsTest.txt'
A_list, B_list, D_list, N_D, N_samples = pa3.read_sample_readings(first_reading_file, N_A, N_B)

d_list = []
for i in range(N_samples):
    A_sublist = A_list[i*N_A:(i+1)*N_A]
    B_sublist = B_list[i*N_B:(i+1)*N_B]
    Fa = registration(a_list, A_sublist)
    Fb = registration(b_list, B_sublist)
    
    
    d = Fb.inv()*(Fa * a_tip)
    d_list.append(d)
    


Mesh = pa3.load_mesh_from_file('Problem3Mesh.sur')














"""
def main(dataset):
    
    print("\n\n\n*********************DATASET {}*************************".format(dataset.upper()))
    name = dataset
    
    if dataset in 'abcdefg':
        calbody_filename = 'pa1-debug-'+ name +'-calbody.txt'
        calreading_filename = 'pa1-debug-'+ name +'-calreadings.txt'
        empivot_filename = 'pa1-debug-' + name + '-empivot.txt'
        optpivot_filename = 'pa1-debug-' + name +'-optpivot.txt'
        output_filename = 'pa1-debug-'+ name +'-output1.txt'
    
    elif dataset in 'hijk':
        calbody_filename = 'pa1-unknown-'+ name +'-calbody.txt'
        calreading_filename = 'pa1-unknown-'+ name +'-calreadings.txt'
        empivot_filename = 'pa1-unknown-' + name + '-empivot.txt'
        optpivot_filename = 'pa1-unknown-' + name +'-optpivot.txt'
        output_filename = 'pa1-unknown-'+ name +'-output1.txt'
    
    
    #=========================Registration, Problem 4==============================
    
    print("Problem 4: Compute Expected C")
    C_expected, N_C, N_Frames = pa1.compute_C_expected(calbody_filename, calreading_filename)
    
    # read the C values from given output file
    if dataset in 'abcdefg':
        C_output = pa1.get_C_output(output_filename)
        # Compare the C_expected value with C values from given output file
        average_error_x, average_error_y, average_error_z = pa1.compare_c(C_expected, C_output)
    
        print('Average error of C_Expected = ')
        print( average_error_x, average_error_y, average_error_z)
        print('\n')
    
    
    #====================Pivot Calibration, Problem 5 and Problem 6===============
    
    print("Problem 5: EM Problem pivot Calibration")
    pt_EM, p_pivot_EM = pivot_calibration.EM_Pivot_Calibration(empivot_filename)
    print("tp = \n",pt_EM)
    
    print("Problem 6: Optical Probe pivot Calibration")
    pt_OP, p_pivot_OP = pivot_calibration.OP_Pivot_Calibration(optpivot_filename, calbody_filename)
    print("tp = \n",pt_OP)
    
    #Error analysis
    if dataset in 'abcdefg':
        pt_EM_answer, pt_OP_answer = pa1.get_calibration_output(output_filename)
        print("Error Analysis:{} \n".format(dataset))
        print("Expected EM probe tip:\n{}".format(pt_EM_answer))
        print("Expected Optical probe tip:\n{}".format(pt_OP_answer))
        print("")
        print("Error of EM probe tip:\n{}".format(pt_EM-pt_EM_answer))
        print("Error of Optical probe tip:\n{}".format(pt_OP-pt_OP_answer))
    
    
    #====================Export Output Text File==================================
    
    with open('../OUTPUT/{}'.format(output_filename), 'w') as f:
        f.write(str(N_C) +', ' + str(N_Frames) + ', ' + 'pa1-debug-'+ name +'-output1.txt')
        f.write('\n')
        f.write('  ' + str(round(pt_EM.x, 2)) +',   '+str(round(pt_EM.y, 2)) +',   '+ str(round(pt_EM.z, 2)))
        f.write('\n')
        f.write('  ' + str(round(pt_OP.x,2)) + ',   '+str(round(pt_OP.y, 2)) + ',   ' + str(round(pt_OP.z, 2)))
        f.write('\n')
    
        for row in range(len(C_expected)):
            f.write('  ' + str(round(C_expected[row, 0], 2)) +',   ' + str(round(C_expected[row, 1], 2)) + ',   ' + str(round(C_expected[row, 2], 2)))
            f.write('\n')


if __name__ == '__main__':
    for dataset in 'abcdefghijk':
        main(dataset)
"""