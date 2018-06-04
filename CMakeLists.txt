cmake_minimum_required(VERSION 3.8)
project(tristan_mp_pitp)
enable_language (Fortran)

set(SOURCE_FILES
        aux1.F90
        communications.F90
        domain.F90
        dynamic_domain.F90
        fieldboundaries.F90
        fields.F90
        filter.F90
        fparser.F90
        globaldata.F90
        initialize.F90
        inputparser.c
        inputparserf.F90
        mpidummy.F90
        optimized_filters.F90
        output.F90
        overload.F90
        par.c
        par.h
        particles.F90
        particles_movedeposit.F90
        restart.F90
        selectprt.F90
        system.c
        systemf.F90
        tristan
        tristan.F90
        tristanmainloop.F90
        user_shock.F90)

add_custom_target(tristan_mp_pitp COMMAND make -C "/home/vadim/tristan-mp-pitp")