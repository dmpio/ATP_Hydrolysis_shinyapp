#-------------------------------------------------------------------#
#-------------------------CONSTANTS---------------------------------#
#-------------------------------------------------------------------#

# Eq., Acid, and Mg-Binding Constants at T=25°C, I=0
#---------------------------------------------------#
# kref = Kref @ T=25°C, I=0; 
# dH = ΔH° (kJ/mol) @ T=25°C, I=0;
# prod = array of product charges (abs. val.);
# react = array of reactant charges (abs. val.);
# name =  string for associative use

# ATP Acid Diss. Const.
is0_Ka_atp <- list(
        kref = 2.512e-08,
        dH = -6.30,
        prod = c(1, 4),
        react = c(3),
        name = 'Ka_atp'
)

# ATP Mg Bind. Const.
is0_Kb_mg_atp <- list(
        kref = 1.514e+06,
        dH = 22.90,
        prod = c(NA),
        react = c(4),
        name = 'Kb_mg_atp'
)

# Protonated ATP Mg Bind. Const.
is0_Kb_mg_hatp <- list(
        kref = 4.266e+03,
        dH = 16.90,
        prod = c(1),
        react = c(2, 3),
        name = 'Kb_mg_hatp'
)

# ADP Acid Diss. Const.
is0_Ka_adp <- list(
        kref = 6.607e-08,
        dH = -5.60,
        prod = c(1, 3),
        react = c(2),
        name = 'Ka_adp'
)

# ADP Mg Bind. Const.
is0_Kb_mg_adp <- list(
        kref = 4.466e+04,
        dH = 19,
        prod = c(1),
        react = c(2, 3),
        name = 'Kb_mg_adp'
)

# AMP Acid Diss. Const.
is0_Ka_amp <- list(
        kref = 1.862e-07,
        dH = -5.40,
        prod = c(2),
        react = c(NA),
        name = 'Ka_amp'
)

# AMP Mg Bind. Const.
is0_Kb_mg_amp <- list(
        kref = 6.165e+02,
        dH = 11.30,
        prod = c(NA),
        react = c(2, 2),
        name = 'Kb_mg_amp'
)

# Protonated ADP Mg Bind. Const. 
is0_Kb_mg_hadp <- list(
        kref = 3.163e+02,
        dH = 12.50,
        prod = c(NA),
        react = c(2, 2),
        name = 'Kb_mg_hadp'
)

# Phosphocreatine Acid Diss. Const.
is0_Ka_pcr <- list(
        kref = 8.854e-06,
        dH = 2.66,
        prod = c(2),
        react = c(NA),
        name = 'Ka_pcr'
)

# Phosphocreatine Mg Bind. Const.
is0_Kb_mg_pcr <- list(
        kref = 2.320e+02,
        dH = 8.19,
        prod = c(NA),
        react = c(2, 2),
        name = 'Kb_mg_pcr'
)

# Phosphate Acid Diss. Const.
is0_Ka_pho <- list(
        kref = 6.026e-8,
        dH = 12.2,
        prod = c(2),
        react = c(NA),
        name = 'Ka_pho'
)

# Phosphate Mg Bind. Const. 
is0_Kb_mg_pho <- list(
        kref = 5.128e+02,
        dH = 8.19,
        prod = c(NA),
        react = c(2, 2),
        name = 'Kb_mg_pho'
)

# Creatine Kinase
is0_Kref_CK <- list(
        kref = 2.580461e+08,
        dH = -17.55,
        prod = c(4),
        react = c(1, 2, 3),
        name = 'Keq_CK'
)

# ATP hydrolysis
is0_Kref_atpHyd <- list(
        kref = 2.946e-1,
        dH = -20.50,
        prod = c(3, 2, 1),
        react = c(4),
        name = 'Keq_atpHyd'
)

# Adenylate Kinase
is0_Kref_AK <- list(
        kref = 2.248e-1,
        dH = -1.50,
        prod = c(4, 2),
        react = c(3, 3),
        name = 'Keq_AK'
)


# Other Constants & Lists of Constants
#---------------------------------------------------#

#Ideal Gas Const. in J/K•mol
const_R <- 8.3144598

# Temp in Kelvin for all reference constants
temp25 <- 25 + 273.15

# Constants used to calculate the Keq of Creatine Kinase at new Temp and I.S.
is0_CK_constants <- list(is0_Ka_atp, is0_Kb_mg_atp, is0_Kb_mg_hatp, 
                         is0_Ka_adp, is0_Kb_mg_adp, is0_Kb_mg_hadp,
                         is0_Ka_pcr, is0_Kb_mg_pcr)

# Constants used to calculate the Keq of ATP hydrolysis at new Temp and I.S.
is0_ATP_constants <- list(is0_Ka_atp, is0_Kb_mg_atp, is0_Kb_mg_hatp, 
                          is0_Ka_adp, is0_Kb_mg_adp, is0_Kb_mg_hadp,
                          is0_Ka_pho, is0_Kb_mg_pho)

# Constants used to calculate the Keq of Adenylate Kinase at new Temp and I.S.
is0_AK_constants <- list(is0_Ka_atp, is0_Kb_mg_atp, is0_Kb_mg_hatp, 
                         is0_Ka_adp, is0_Kb_mg_adp, is0_Kb_mg_hadp,
                         is0_Ka_amp, is0_Kb_mg_amp)

# UI Constants & Lists of Constants
#---------------------------------------------------#


ui_tempC_const <-  list(min=26, 
                        max=46, 
                        step=0.1, 
                        default_range=c(35, 39), 
                        default_value=37)

ui_bufferIS_const <- list(min=15, 
                          max=300, 
                          step=1, 
                          default_range=c(150, 190), 
                          default_value=170)

ui_pH_const <- list(min=5, 
                    max=9, 
                    step=0.1, 
                    default_range=c(6.8, 7.8), 
                    default_value=7.1)

ui_mg_const <- list(min=0.05, 
                    max=2.0, 
                    step=0.01, 
                    default_range=c(0.2, 0.8),
                    default_value=0.5)

ui_atp_const <- list(min=1, 
                    max=10, 
                    step=0.01, 
                    default_range=c(3, 7),
                    default_value=4.99)

ui_adp_const <- list(min=0.001, 
                    max=0.03, 
                    step=0.001, 
                    default_range=c(0.002, 0.015),
                    default_value=0.01)

ui_phosphate_const <- list(min=0.05, 
                           max=20, 
                           step=0.05, 
                           default_range=c(5, 12),
                           default_value=10)

