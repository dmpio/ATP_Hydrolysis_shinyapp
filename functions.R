source(file = "constants.R")

#-------------------------------------------------------------------#
#-------------------------FUNCTIONS---------------------------------#
#-------------------------------------------------------------------#
solve_gamma <- function(products, reactants, temp, ionic_strength){
  # products = array of all charged product ions
  # reactants = array of all charged reactant ions
  # temp = temperature in Kelvins
  # ionic_strength = ionic strength in Molarity
  
  # Function to solve for Γ
  # Debye-Hückel
  
  Am <- 3*(-16.39023 + (261.3371/temp) + 3.3689633*log(temp) - 1.437167*(temp/100) + 0.111995*((temp/100)^2))
  
  # constant with units of sqrt(kg)/sqrt(mol)
  B <-  1.6
  
  len_r = length(reactants)
  reactant_gammas = vector(length=len_r)
  
  if (is.na(reactants[[1]])){
    reactant_gammas[[1]] = 1
  } else {
    for (i in 1:len_r){
      pwr <- -Am * sqrt(ionic_strength) * reactants[[i]]^2 / (1 + B * sqrt(ionic_strength))
      reactant_gammas[[i]] <- exp(pwr)
    }
  }
  
  len_p = length(products)
  product_gammas = vector(length=len_p)
  
  if (is.na(products[[1]])){
    product_gammas[[1]] = 1
  } else {
    for (i in 1:len_p){
      pwr <- -Am * sqrt(ionic_strength) * products[[i]]^2 / (1 + B * sqrt(ionic_strength))
      product_gammas[[i]] <- exp(pwr)
    }
  }
  
  return(prod(product_gammas)/prod(reactant_gammas))
}


vant_hoff <- function(temp1, temp2, Kref1, dH){
  # temp1 = temperature for constant Kref1 (Kelvins)
  # Kref1 = constant at temp1 and ionic strength 0 (the original constant from which Kref2 will be adjusted from and returned)
  # temp2 = temperature for constant Kref2 (the new temperature to which Kref1 is to be adjusted)(Kelvins)
  # dH = ΔH° associated with Kref1 at temp1 and ionic strength 0 (in kilojoules)
  
  # return: modified constant at ionic strength 0 temp2
  
  # Multiply dH by 1000 to convert from kJ to J
  pwr <- (-1000*dH/const_R)*(1/temp2 - 1/temp1) + log(Kref1)
  
  return(exp(pwr))
}


adjust_Kabs <- function(ionic_strength, tempK, is0_constants){
  # Adjust a vector of Acid-base and Mag-binding dissociation constants for temp and ionic strength
  
  # A temporary dictionary to hold the modified acid and Mg-binding Krefs
  new_Kabs = list()
  
  # iterate through the acid and Mg-binding Krefs in the reaction
  for(dict in is0_constants) {
    
    # new kref adjusted from I=0, T=25°C to I=0, T=tempK
    kref2 <- vant_hoff(temp1 = temp25, 
                       temp2 = tempK, 
                       Kref1 = dict$kref, 
                       dH = dict$dH
    )
    
    # new kref adjusted from I=0, T=25°C to I=i.s. of the step, T=tempK
    kref3 <- kref2/solve_gamma(products = dict$prod, 
                               reactants = dict$react, 
                               temp = tempK, 
                               ionic_strength = ionic_strength
    )
    # save the adjusted kref for later use in eq reaction
    new_Kabs[[dict$name]] <- kref3
  }
  return(new_Kabs)
}


solve_ATP_Hyd_Keq <- function(new_Kabs, Kref, proton, freeMg){
  
  # Solve the Equilibrium reaction for ATP Hydrolysis accounting for pH & Free-Mg
  
  # Assembling the modified constants into the equilibrium reaction
  
  
  eq_num_adp <- 1 + (proton/new_Kabs$Ka_adp) + (new_Kabs$Kb_mg_adp*freeMg) + (new_Kabs$Kb_mg_hadp*proton*freeMg/new_Kabs$Ka_adp)
  eq_num_pho <- 1 + (proton/new_Kabs$Ka_pho) + (new_Kabs$Kb_mg_pho*freeMg)
  eq_den_atp <- 1 + (proton/new_Kabs$Ka_atp) + (new_Kabs$Kb_mg_atp*freeMg) + (new_Kabs$Kb_mg_hatp*proton*freeMg/new_Kabs$Ka_atp)
  
  return((Kref * eq_num_adp * eq_num_pho) / (proton * eq_den_atp))
}


solve_CK_Keq <- function(new_Kabs, Kref, proton, freeMg){
  # Solve the Equilibrium reaction for Creatine Kinase accounting for pH & Free-Mg
  
  # Assembling the modified constants into the equilibrium reaction
  eq_num_atp <- 1 + (proton/new_Kabs$Ka_atp) + (new_Kabs$Kb_mg_atp*freeMg) + (new_Kabs$Kb_mg_hatp*proton*freeMg/new_Kabs$Ka_atp)
  eq_den_adp <- 1 + (proton/new_Kabs$Ka_adp) + (new_Kabs$Kb_mg_adp*freeMg) + (new_Kabs$Kb_mg_hadp*proton*freeMg/new_Kabs$Ka_adp)
  eq_den_pcr <- 1 + (proton/new_Kabs$Ka_pcr) + (new_Kabs$Kb_mg_pcr*freeMg)
  
  return((Kref * proton * eq_num_atp) / (eq_den_adp * eq_den_pcr))
}


solve_AK_Keq <- function(new_Kabs, Kref, proton, freeMg){
  # Solve the Equilibrium reaction for Adenylate Kinase accounting for pH & Free-Mg
  
  # Assembling the modified constants into the equilibrium reaction
  eq_num_atp <- 1 + (proton/new_Kabs$Ka_atp) + (new_Kabs$Kb_mg_atp*freeMg) + (new_Kabs$Kb_mg_hatp*proton*freeMg/new_Kabs$Ka_atp);
  eq_num_amp <- 1 + (proton/new_Kabs$Ka_amp) + (new_Kabs$Kb_mg_amp*freeMg);
  eq_den_adp <- 1 + (proton/new_Kabs$Ka_adp) + (new_Kabs$Kb_mg_adp*freeMg) + (new_Kabs$Kb_mg_hadp*proton*freeMg/new_Kabs$Ka_adp);
  
  return((Kref * eq_num_atp * eq_num_amp) / (eq_den_adp * eq_den_adp))
}


calc_apparent_Keq <- function(reaction, is0_Kref, tempC, buffer_IS, pH, mg, ...){
  
  ### --- Function Variables --- ###
  # reaction  = (string) equilibrium reaction to adjust
  # is0_kref  = (named list) reaction information. 
  #              Valid options include is0_Kref_atpHyd, is0_Kref_CK, is0_Kref_AK
  # tempC     = (number) temperature in Celsius to adjust reaction to
  # buffer_IS = (number) ionic strength of the buffer
  # pH        = (number)
  # mg        = (number) free magnesium in the buffer, ideally empirically determined
  
  # ---
  # return - list containing adjusted equilibrium constants (KJ/mol)
  
  # additional parameters
  optional <- list(...)
  
  
  # Temperature in Kelvin
  # Data entered as Celsius
  tempK <- tempC + 273.15
  
  # Ionic Strength in Molarity
  # Data entered as millimolar
  ionic_strength <- buffer_IS/1000
  
  # Molar Concentration of H+ ions
  proton <-  10 ^ (-pH)
  
  # Free Mg in Molarity
  # Data entered as millimolar
  freeMg <- mg/1000
  
  
  # Adjust the Kref of reaction @ I=0, T=25°C to tempK
  Kref_newT_is0 <- vant_hoff(temp1 = temp25, 
                             temp2 = tempK, 
                             Kref1 = is0_Kref$kref, 
                             dH = is0_Kref$dH
  )
  
  # Adjust the Kref of reaction @ I=0, T=tempK to Ionic Strength of assay
  Kref_adj = Kref_newT_is0/solve_gamma(products = is0_Kref$prod, 
                                       reactants = is0_Kref$react, 
                                       temp = tempK, 
                                       ionic_strength = ionic_strength)
  
  # Adjust the Kb and Ka values of the Eq. reaction
  # and multiply them by the adjusted Krefs 
  # to calculate the apparent Keq
  
  # Adenylate Kinase 
  if (reaction=='AK'){
    
    new_Kabs <- adjust_Kabs(ionic_strength = ionic_strength, 
                            tempK = tempK, 
                            is0_constants = is0_AK_constants
    )
    
    Keq <- solve_AK_Keq(new_Kabs = new_Kabs, 
                        Kref = Kref_adj, 
                        proton = proton, 
                        freeMg = freeMg
    )
    results <- list("Keq"=Keq)
  }
  
  # Creatine Kinase
  if (reaction=='CK'){
    
    new_Kabs <- adjust_Kabs(ionic_strength = ionic_strength, 
                            tempK = tempK, 
                            is0_constants = is0_CK_constants
    )
    
    Keq <- solve_CK_Keq(new_Kabs = new_Kabs, 
                        Kref = Kref_adj, 
                        proton = proton, 
                        freeMg = freeMg
    )
    
    results <- list("Keq"=Keq)
  }
  
  # ATP Hydrolysis
  if (reaction=='ATP'){
    
    # ADP in Molarity
    # Data entered as millimolar
    if(exists("adp", optional)){
      adp <- optional$adp/1000
    } else {
      return(c("Error - No ADP Specified"))
    }
    
    # ATP in Molarity
    # Data entered as millimolar
    if(exists("atp", optional)){
      atp <- optional$atp/1000
    } else {
      return(c("Error - No ATP Specified"))
    }
    
    # Phosphate in Molarity
    # Data entered as millimolar
    if(exists("phosphate", optional)){
      phosphate <- optional$phosphate/1000
    } else {
      return(c("Error - No Phosphate Specified"))
    }
    
    new_Kabs <- adjust_Kabs(ionic_strength = ionic_strength, 
                            tempK = tempK, 
                            is0_constants = is0_ATP_constants
    )
    
    # The Modified ATP Hydrolysis Eq. Constant
    Keq <- solve_ATP_Hyd_Keq(new_Kabs = new_Kabs, 
                             Kref = Kref_adj, 
                             proton = proton, 
                             freeMg = freeMg
    )
    
    
    # The ΔG°' of ATP hydrolysis (Joules/mol)
    dG_ATP_keq  <-  -1*const_R*tempK*log(Keq)
    
    # The ΔG' of ATP hydrolysis (Joules/mol)
    dG_ATP <-  dG_ATP_keq + const_R*tempK*log((adp * phosphate)/atp);
    
    results <- list(
      "dG_ATP" = dG_ATP,
      "dG_ATP_keq"= dG_ATP_keq,
      "Keq" = Keq
    )
  }
  
  return(results)
}
