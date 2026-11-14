@echo off
rem Run every plate-with-hole case that fits the Learning Edition node cap
call abaqus job=kirsch_CPS4_n1 interactive
call abaqus job=kirsch_CPS4_n2 interactive
call abaqus job=kirsch_CPS4_n4 interactive
call abaqus job=kirsch_CPS4R_n1 interactive
call abaqus job=kirsch_CPS4R_n2 interactive
call abaqus job=kirsch_CPS4R_n4 interactive
call abaqus job=kirsch_CPS4I_n1 interactive
call abaqus job=kirsch_CPS4I_n2 interactive
call abaqus job=kirsch_CPS4I_n4 interactive
call abaqus job=kirsch_CPS8_n1 interactive
call abaqus job=kirsch_CPS8_n2 interactive
call abaqus job=kirsch_CPS8R_n1 interactive
call abaqus job=kirsch_CPS8R_n2 interactive
call abaqus job=plate_CPS4_n1 interactive
call abaqus job=plate_CPS4_n2 interactive
call abaqus job=plate_CPS4_n4 interactive
call abaqus job=plate_CPS4R_n1 interactive
call abaqus job=plate_CPS4R_n2 interactive
call abaqus job=plate_CPS4R_n4 interactive
call abaqus job=plate_CPS4I_n1 interactive
call abaqus job=plate_CPS4I_n2 interactive
call abaqus job=plate_CPS4I_n4 interactive
call abaqus job=plate_CPS8_n1 interactive
call abaqus job=plate_CPS8_n2 interactive
call abaqus job=plate_CPS8R_n1 interactive
call abaqus job=plate_CPS8R_n2 interactive
