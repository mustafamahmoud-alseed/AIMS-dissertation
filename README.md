# AIMS-dissertation

1. LF-simulation:

 Run the script in two steps:

 - First, run without any interventions, such that
`season <- ( 1 + amp*cos(2*pi*times/180))`.

- Second, run with interventions when
`season <- ( 1   amp*cos(2*pi*times/180))*exp(-m)` 
where m is the intervention varying from 0.1 to 0.4.
Run the script when `m = 0.1, 0.2, 0.3` and `0.4`, thus you can observe the impact of the interventions.

The focus should be on the human's plottings.


2. LF-fitting:

- Run the script without interventions for `I_{ha},I_{hc}` and `I_{ha} + I_{hc}`

- Repeat the process with interventions
