@echo off
rem Run every Hertz contact case that fits the Learning Edition node cap
call abaqus job=hertz_CAX4_n1 interactive
call abaqus job=hertz_CAX4_n2 interactive
call abaqus job=hertz_CAX4R_n1 interactive
call abaqus job=hertz_CAX4R_n2 interactive
call abaqus job=hertz_CAX4I_n1 interactive
call abaqus job=hertz_CAX4I_n2 interactive
call abaqus job=hertz_CAX8_n1 interactive
call abaqus job=hertz_CAX8R_n1 interactive
