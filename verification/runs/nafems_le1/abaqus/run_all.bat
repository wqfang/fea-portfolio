@echo off
rem Run every NAFEMS LE1 case that fits the Learning Edition node cap
call abaqus job=le1_CPS4_n1 interactive
call abaqus job=le1_CPS4_n2 interactive
call abaqus job=le1_CPS4_n4 interactive
call abaqus job=le1_CPS4_n8 interactive
call abaqus job=le1_CPS4R_n1 interactive
call abaqus job=le1_CPS4R_n2 interactive
call abaqus job=le1_CPS4R_n4 interactive
call abaqus job=le1_CPS4R_n8 interactive
call abaqus job=le1_CPS4I_n1 interactive
call abaqus job=le1_CPS4I_n2 interactive
call abaqus job=le1_CPS4I_n4 interactive
call abaqus job=le1_CPS4I_n8 interactive
call abaqus job=le1_CPS8_n1 interactive
call abaqus job=le1_CPS8_n2 interactive
call abaqus job=le1_CPS8_n4 interactive
call abaqus job=le1_CPS8R_n1 interactive
call abaqus job=le1_CPS8R_n2 interactive
call abaqus job=le1_CPS8R_n4 interactive
