setwd("C:/Users/cani/OneDrive - Australian Institute of Marine Science/Desktop/Buoyancy")

col_len = seq(0.001,0.002,1e-5)  # colony length
col_wid = seq(0.00005,0.00015,1e-6)  # colony width
alpha_vert = 2.0 / ((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 *    #shape factor independent of orientation for a vertically oriented colony
            atanh(((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^2.0) ^ 0.5 / (col_len / 2.0))

beta_vert = 2.0 * (col_len / 2.0) ^ 2.0 / ((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 1.5 * #shape factor dependent on orientation for a vertically oriented colony
            (atanh(((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 / (col_len / 2.0)) - 
            ((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 / (col_len / 2.0))

Phi_vert = 16.0 / (3.0 * 2.0 * ((col_len / 2.0) * (col_wid / 2.0) ^ 2.0) ^ (1.0 / 3.0) * (alpha_vert + beta_vert)) #Vertically oriented form resistance

alpha_hor = 2.0 / (((col_len / 2.0) ^ 2.0 - ((col_wid / 2.0) ^ 2.0) ^ 0.5)) *  #shape factor independent of orientation for a horizontally oriented colony
            atanh(((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 / (col_len / 2.0))
            
beta_hor = (col_wid / 2.0) ^ 2.0 / (((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 1.5) * #shape factor dependent on orientation for a horizontally oriented colony
           ((col_len / 2.0) * ((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 / ((col_wid / 2.0) ^ 2.0) * 
              atanh(((col_len / 2.0) ^ 2.0 - (col_wid / 2.0) ^ 2.0) ^ 0.5 / (col_len / 2.0)))

Phi_hor = 16.0 / (3.0 * 2.0 * (col_len / 2.0 * (col_wid / 2.0) ^ 2.0) ^ (1.0 / 3.0) * (alpha_hor + beta_hor)) #Horizontally oriented form resistance
            
Phi_total = (Phi_vert + Phi_hor) / 2

colrad_eff = (3.0 * (col_wid / 2.0) ^ 2.0 * col_len / 4.0) ^ (1.0 / 3.0)  #effective radius
sink_vel = 2 * 9.81 * colrad_eff ^ 2 *(1030-1020) / (9 * 0.001 * Phi_total) #sinking velocity

plot(colrad_eff * 10 ^ 4, sink_vel * 10 ^ 4,
     #ylab=parse(text="form~resistance~factor"),
     ylab=parse(text="sinking~velocity~(10^-4~m~s^-1)"),
     xlab=parse(text="effective~radius~(10^-4~m)"),
     #xlab = "form resistance factor",
     type="l",
     cex.lab=1.4,  # increase label font size
     lwd=3,
     #xaxt='n',
     #yaxt='n',
     mgp=c(2.4,1,.0), las=1, 
     col="blue") 


