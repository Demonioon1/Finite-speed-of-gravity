# Finite-speed-of-gravity
This is a simulation of N body system under a Newtonian Force that propagates at finite velocity. 
To run install FPM (Fortran Pakage Manager) 

You can change the propagation velocity in the (app/main.f90) is the variable (vp) vp=1 means 1 c (c = light speed)
You can change (N number of particles, nt iteartion number, and dt iteration time step) in the archive (src/par.f90) 

the second library (scr/p2.f90) is a Polinomial solver wrote by: 

((((( WRITTEN using F90 intrinsics by
  !        Alan Miller
  !        amiller @ bigpond.net.au
  !     Latest revision - 1 February 1997))))))
  
  
  The program use real(kind=18) pricision, you can change it, in both library modules in the folder (scr)
  
 HOW TO RUN:
 
 1) Install FPM (Fortran Pakage Manager) 
 In the terminal:
2) : cd "folder p2" 
3) : fpm run
4) : in3.dat 

in3.dat are the initial conditions (position and velocity) for Sun, Earth and Moon. 

Units : position: AU (astronomic units)
        time:     y (years)
        
        
 <img width="522" alt="Captura de pantalla 2022-12-10 a la(s) 04 19 01" src="https://user-images.githubusercontent.com/93294992/206842982-c5350db5-922f-4524-afb0-0136f3a17e05.png">
<img width="594" alt="Captura de pantalla 2022-12-10 a la(s) 04 17 56" src="https://user-images.githubusercontent.com/93294992/206843063-5c49782f-1aa9-49ef-9f94-187f5ef0b6f7.png">
<img width="599" alt="Captura de pantalla 2022-12-10 a la(s) 04 18 09" src="https://user-images.githubusercontent.com/93294992/206843068-ee8300dc-ef7e-4fff-8f0c-0a5af64e7c24.png">
     
     
 

