grim
====

General Relativistic Implicit Magnetohydrodynamics

Run using 
time ./grim -ts_monitor -snes_monitor -snes_converged_reason -snes_rtol 1e-50 -snes_atol 1e-4 -snes_lag_jacobian 100 -snes_lag_jacobian_persists TRUE -ts_final_time 1000 -ts_max_snes_failures -1 -ts_max_steps 50000 -snes_stol 1e-50
