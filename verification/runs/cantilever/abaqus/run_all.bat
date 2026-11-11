@echo off
rem Run every cantilever case that fits the Learning Edition node cap
call abaqus job=cant_C3D8_n1 interactive
call abaqus job=cant_C3D8_n2 interactive
call abaqus job=cant_C3D8R_n1 interactive
call abaqus job=cant_C3D8R_n2 interactive
call abaqus job=cant_C3D8I_n1 interactive
call abaqus job=cant_C3D8I_n2 interactive
call abaqus job=cant_C3D20_n1 interactive
call abaqus job=cant_C3D20R_n1 interactive
