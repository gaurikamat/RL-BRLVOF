
# Iteration parameters
nsim = 100
ngibbs=1e3
postburnin = 1e2
endval=1e3

# Error levels
if(er=="P"){prob <- 0}
if(er=="L"){prob <- 2}
if(er=="H"){prob <- 4}

# File sizes
N_A=500
N_B=1000
n_m=300


# Association between XA and XB
cond_assoc=0.1 #links
non_link_assoc = 0.5#non-links

if(a_par==0.5){
  a_m1=0.5
  a_m2=0.5
  a_m3=0.5
  a_m4=0.5
  a_u1=0.5
  a_u2=0.5
  a_u3=0.5
  a_u4=0.5
}


if(a_par==1){
  a_m1=1
  a_m2=1
  a_m3=1
  a_m4=1
  a_u1=0.5
  a_u2=0.5
  a_u3=0.5
  a_u4=0.5
}


if(a_par==2){
  a_m1=2
  a_m2=2
  a_m3=2
  a_m4=2
  a_u1=0.5
  a_u2=0.5
  a_u3=0.5
  a_u4=0.5
}


if(a_par==5){
  a_m1=5
  a_m2=5
  a_m3=5
  a_m4=5
  a_u1=0.5
  a_u2=0.5
  a_u3=0.5
  a_u4=0.5
}

if(a_par==10){
  a_m1=10
  a_m2=10
  a_m3=10
  a_m4=10
  a_u1=0.5
  a_u2=0.5
  a_u3=0.5
  a_u4=0.5
}




## Regression parameters

# covariate mean
x_mn = 1
x_var = 4

# regression error sd
if(er_sd == "H"){error_sd = 0.5}
if(er_sd == "L"){error_sd = 0.1}










