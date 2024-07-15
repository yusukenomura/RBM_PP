#!/bin/sh -x
ifort -debug full -check all -warn all -traceback -g  sub.f90 post_process.f90 
