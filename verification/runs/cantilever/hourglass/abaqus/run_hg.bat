@echo off
rem C3D8R hourglass-control variants
call abaqus job=c3d8r_n1_default interactive
call abaqus job=c3d8r_n1_stiff1 interactive
call abaqus job=c3d8r_n1_stiff10 interactive
call abaqus job=c3d8r_n1_enhanced interactive
call abaqus job=c3d8r_n2_default interactive
call abaqus job=c3d8r_n2_stiff1 interactive
call abaqus job=c3d8r_n2_stiff10 interactive
call abaqus job=c3d8r_n2_enhanced interactive
