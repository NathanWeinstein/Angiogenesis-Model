library(BoolNet)
load("wild_type.Rdata")

environmentsByNumberOfAttractors <- function(classed){
	types <- classed$types
	ltypes <- lengths(types)
	nt <- names(types)
	max_attractors <- max(ltypes)
	envsByNofAttractors <- vector(mode = "list", length = max_attractors)
	for(i in 1:length(ltypes)){
		l <- ltypes[i]
		envsByNofAttractors[[l]] <- append(envsByNofAttractors[[l]], nt[i])
	}
	return(envsByNofAttractors)
}

analyzeECBehavior <- function(net, Attractors, sortedEnvironments, classed, microenvironment_variables, kind){
	# First get the right environments
	if(kind == "Phalanx"){
		environments <- sortedEnvironments$Phalanx
	}else if(kind == "Stalk"){
		environments <- sortedEnvironments$Stalk
	}else if(kind == "Tip"){
		environments <- sortedEnvironments$Tip
	}else if(kind == "Atypical"){
		environments <- sortedEnvironments$Atypical
	}else{
		print("No such kind")
	}
	#print(environments)
	# Build the environment matrix
	np <- length(environments)
	mel <- length(microenvironment_variables)
	environment_matrix <- matrix(NA, nrow = np, ncol = mel)
	for(i in 1:np){
		pbin <- decimal_to_binary(environments[i], microenvironment_variables)
		environment_matrix[i,] <- pbin
	}
	colnames(environment_matrix) <- microenvironment_variables
	rownames(environment_matrix) <- environments
	# Chack which molecules have fixed values
	for(i in 1:length(microenvironment_variables)){
		#print(microenvironment_variables[i])
		f <- factor(environment_matrix[,i])
		#print(summary(f))
	}
	return(environment_matrix)
}

# Analize a list of environments
environmentsToMatrix <- function(microenvironment_variables, environments){
	# Build the environment matrix
	np <- length(environments)
	mel <- length(microenvironment_variables)
	environment_matrix <- matrix(NA, nrow = np, ncol = mel)
	for(i in 1:np){
		pbin <- decimal_to_binary(environments[i], microenvironment_variables)
		environment_matrix[i,] <- pbin
	}
	colnames(environment_matrix) <- microenvironment_variables
	rownames(environment_matrix) <- environments
	# Chack which molecules have fixed values
	for(i in 1:length(microenvironment_variables)){
		#print(microenvironment_variables[i])
		f <- factor(environment_matrix[,i])
		#print(summary(f))
	}
	return(environment_matrix)
}

# Dynamic data structures in R are slow but here I do not know how many attractor indexes will be in the end
getAttractorsFromeEnvironments <- function(envlist, classed){
	classed_indexes <- classed$indexes
	attractor_indexes <- c()
	for(env in envlist){
		env_attr_indexes <- unlist(classed_indexes[env])
		attractor_indexes <- append(attractor_indexes, env_attr_indexes)
	}
	return(attractor_indexes)
}

# Check which attractors divide and call mural cells
getDividingMuralFixed <- function(envm, classed, Attractors, net){
	all_divide <- c()
	some_divide <- c()
	none_divide <- c()
	all_pdgfb <- c()
	some_pdgfb <- c()
	none_pdgfb <- c()
	attractor_lengths <- c()
	atris <- c()
	natr <- 0
	classed_indexes <- classed$indexes
	fixed <- rep(1, length(net$genes))
	names(fixed) <- net$genes
	envlist <- rownames(envm)
	env0_attr_indexes <- unlist(classed_indexes[envlist[1]])
	atr0_frame <- getAttractorSequence(Attractors, env0_attr_indexes[1])
	st0 <- atr0_frame[1,]
	for(env in envlist){
		divides <- 1
		no_division <- 1
		mural <- 1
		no_mural <- 1
		# First get the attractor indexes for the environment
		env_attr_indexes <- unlist(classed_indexes[env])
		for(atri in env_attr_indexes){
			natr <- natr + 1
			atr_frame <- getAttractorSequence(Attractors, atri)
			l <- nrow(atr_frame)
			atris <- append(atris, atri)
			attractor_lengths <- append(attractor_lengths, l)
			pdgfb <- 1
			cyclind <- 1
			for(s in 1:l){
				st <- atr_frame[s, ]
				if(st["FOXO1"] == 0 && st["SMAD2"] == 0){
					pdgfb <- 0
				}
				if(st["Bcatenin"] == 0 || st["LEF1"] == 0){
					cyclind <- 0
				}
				for(j in 1:length(st0)){
					if(st[j] != st0[j]){
						fixed[j] <- 0
					}
				}
			}
			if(pdgfb == 1){
				no_mural <- 0
			}else{
				mural <- 0
			}
			if(cyclind == 1){
				no_division <- 0
			}else{
				divides <- 0
			}
		}
		if(divides == 1){
			all_divide <- append(all_divide, env)
		}else if(no_division == 1){
			none_divide <- append(none_divide, env)
		}else{
			some_divide <- append(some_divide, env)
		}
		if(mural == 1){
			all_pdgfb <- append(all_pdgfb, env)
		}else if(no_mural == 1){
			none_pdgfb <- append(none_pdgfb, env)
		}else{
			some_pdgfb <- append(some_pdgfb, env)
		}
	}
	names(attractor_lengths) <- atris
	results <- list(all_divide = all_divide,
			some_divide = some_divide,
			none_divide = none_divide,
			all_pdgfb = all_pdgfb,
			some_pdgfb = some_pdgfb,
			none_pdgfb = none_pdgfb,
			number_of_attractors = natr,
			attractor_lengths = attractor_lengths,
			fixed = st0[fixed == 1])
	return(results)
}


## First analyze Phalanx cell behavior
sortedEnvironments <- behavior$Sorted
net <- behavior$Network
classed <- behavior$Classed
Attractors <- behavior$Attractors
microenvironment_variables <- c("VEGFC_Dp", "VEGFAxxxP", "ANG1", "Oxygen", "ShearStress", "JAGp", "DLL4p", "WNT5a", "WNT7a", "FGF","IGF", "BMP9", "BMP10", "TGFB1", "VEGFC_D", "AMPATP")
phalanx_env_matrix <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Phalanx")
phalanx_env_analisis <- getDividingMuralFixed(phalanx_env_matrix, classed, Attractors, net)
Phalanx <- list("Environment Matrix" = phalanx_env_matrix,"Analisis" = phalanx_env_analisis)
## Second analyze Tip cell behavior

tip_env_matrix <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Tip")
# Get the list of VEGF active environments that cause Tip cell behavior
vegf_tip_matrix <- tip_env_matrix[which(tip_env_matrix[, "VEGFAxxxP"] == 1 | tip_env_matrix[, "VEGFC_Dp"] == 1),]
vegf_tip_analisis <- getDividingMuralFixed(vegf_tip_matrix, classed, Attractors, net)
TipI <- list("Environment Matrix" = vegf_tip_matrix,"Analisis" = vegf_tip_analisis)
no_vegf_tip_matrix <- tip_env_matrix[which(tip_env_matrix[, "VEGFAxxxP"] == 0 & tip_env_matrix[, "VEGFC_Dp"] == 0),]
no_vegfs1 <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 1 | no_vegf_tip_matrix[, "WNT5a"] == 1 | no_vegf_tip_matrix[, "WNT7a"] == 1),]
no_vegfs1_analisis <- getDividingMuralFixed(no_vegfs1, classed, Attractors, net)
TipII <- list("Environment Matrix" = no_vegfs1,"Analisis" = no_vegfs1_analisis)
o12_tip_matrix <- no_vegf_tip_matrix [which(no_vegf_tip_matrix[, "ShearStress"] == 0 & no_vegf_tip_matrix[, "WNT5a"] == 0 & no_vegf_tip_matrix[, "WNT7a"] == 0),]
no_vegfs2_analisis <- getDividingMuralFixed(o12_tip_matrix, classed, Attractors, net)
TipIII <- list("Environment Matrix" = o12_tip_matrix,"Analisis" = no_vegfs2_analisis)
Tip <- list("TipI" = TipI,"TipII" = TipII, "TipIII" = TipIII)
## Now lets analyze stalk cell behavior

stalk_env_matrix <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Stalk")

# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = *, "ShearStress" = *, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = *, "WNT7a" = *, "FGF" = *,"IGF" = *, "BMP9" = *, 
#   "BMP10" = *, "TGFB1" = *, "VEGFC_D" = *, "AMPATP" = *) 2^14 = 16384 dim(stalk_env_matrix) -> 12096 16

st_env_summary <- summary(stalk_env_matrix)

# How many 0f the 16384 have WNT active?
# 2^14 - 2^12 = 3 * 2^12 = 12288 ... not all are stalk cells
wnt_stalk_m <- stalk_env_matrix[which(stalk_env_matrix[, "WNT5a"] == 1 | stalk_env_matrix[, "WNT7a"] == 1),]
# Now how many of those 12288 would also lack sufficient Oxygen or IGF? 12288*(3/4) = 9216 All of those are Stalk cells!
wnt_nigf_noxy <- wnt_stalk_m[which(wnt_stalk_m[, "Oxygen"] == 0 | wnt_stalk_m[, "IGF"] == 0),]
wnt_nigf_noxy_analisis <- getDividingMuralFixed(wnt_nigf_noxy, classed, Attractors, net)
StalkI <- list("Environment Matrix" = wnt_nigf_noxy,"Analisis" = wnt_nigf_noxy_analisis)
wnt_igf_oxy <- wnt_stalk_m[which(wnt_stalk_m[, "Oxygen"] == 1 & wnt_stalk_m[, "IGF"] == 1),] # In this group AMPATP = 1
wnt_igf_oxy_analisis <- getDividingMuralFixed(wnt_igf_oxy, classed, Attractors, net)
StalkII <- list("Environment Matrix" = wnt_igf_oxy,"Analisis" = wnt_igf_oxy_analisis)
# Now the other stalk cells
nwnt_stalk <- stalk_env_matrix[which(stalk_env_matrix[, "WNT5a"] == 0 & stalk_env_matrix[, "WNT7a"] == 0),]
nwnt_stalk_analisis <- getDividingMuralFixed(nwnt_stalk, classed, Attractors, net)
StalkIII <- list("Environment Matrix" = nwnt_stalk,"Analisis" = nwnt_stalk_analisis)
Stalk <- list("StalkI" = StalkI,"StalkII" = StalkII, "StalkIII" = StalkIII)
# Are they all TGF active
nwnt_tgfs <- nwnt_stalk[which(nwnt_stalk[,"BMP9"] == 1 | nwnt_stalk[,"BMP10"] == 1| nwnt_stalk[,"TGFB1"] == 1),] # All of them
# Are they all either JAGp == 1 or DLL4p == 0, NOTCH inactive
nwnt_notchs <- nwnt_stalk[which(nwnt_stalk[,"JAGp"] == 1 | nwnt_stalk[,"DLL4p"] == 0),] # All of them
# c("VEGFC_Dp" = 0, "VEGFAxxxP" = 0, "ANG1" = *, "Oxygen" = *, "ShearStress" = *, 
#   "JAGp" = 0, "DLL4p" = 1, "WNT5a" = *, "WNT7a" = *, "FGF" = *,"IGF" = *, "BMP9" = *, 
#   "BMP10" = *, "TGFB1" = *, "VEGFC_D" = *, "AMPATP" = *) 
# Now lets analyze the other behaviors
otherm <- analyzeECBehavior(net, Attractors, sortedEnvironments, classed, microenvironment_variables, "Atypical")
# c("VEGFC_Dp" = 1, "VEGFAxxxP" = 1, "ANG1" = *, "Oxygen" = *, "ShearStress" = 0, 
#   "JAGp" = *, "DLL4p" = *, "WNT5a" = 0, "WNT7a" = 0, "FGF" = 0,"IGF" = *, "BMP9" = 0, 
#   "BMP10" = 0, "TGFB1" = 0, "VEGFC_D" = *, "AMPATP" = *) 2^7, 2^7 = 128, 128*3 = 384
other_vegf <- otherm[which(otherm[, "VEGFAxxxP"] == 1 | otherm[, "VEGFC_Dp"] == 1),]
analisis_other_vegf <- getDividingMuralFixed(other_vegf, classed, Attractors, net)
AtypicalI <- list("Environment Matrix" = other_vegf,"Analisis" = analisis_other_vegf)
other_nvegf <- otherm[which(otherm[, "VEGFAxxxP"] == 0 & otherm[, "VEGFC_Dp"] == 0),]

# Using AMPATP and Oxygen first
o_nv_amp_no <- other_nvegf[which(other_nvegf[, "AMPATP"] == 1 | other_nvegf[, "Oxygen"] == 0),]
o_nv_amp_no_nss <- other_nvegf[which(other_nvegf[, "AMPATP"] == 1 | other_nvegf[, "Oxygen"] == 0 | other_nvegf[, "ShearStress"] == 0),]

# Using IGF first
other_nv_igf <- other_nvegf[which(other_nvegf[, "IGF"] == 1),]

# Using TGF
other_nv_tgf <- other_nvegf[which(other_nvegf[, "BMP9"] == 1 | other_nvegf[, "BMP10"] == 1 | other_nvegf[, "TGFB1"] == 1 ),]
o_nv_tgf_amp_no_nss <- other_nv_tgf[which(other_nv_tgf[, "AMPATP"] == 1 | other_nv_tgf[, "Oxygen"] == 0 | other_nv_tgf[, "ShearStress"] == 0 ),]
o_nv_tgf_amp_no_nss_igf <- o_nv_tgf_amp_no_nss[which(o_nv_tgf_amp_no_nss[, "IGF"] == 1),]
analisis_o2 <- getDividingMuralFixed(o_nv_tgf_amp_no_nss_igf, classed, Attractors, net)
AtypicalII <- list("Environment Matrix" = o_nv_tgf_amp_no_nss_igf,"Analisis" = analisis_o2)
o_nv_tgf_amp_no_nss_nigf <- o_nv_tgf_amp_no_nss[which(o_nv_tgf_amp_no_nss[, "IGF"] == 0),]
analisis_o3 <- getDividingMuralFixed(o_nv_tgf_amp_no_nss_nigf, classed, Attractors, net)
AtypicalIII <- list("Environment Matrix" = o_nv_tgf_amp_no_nss_nigf,"Analisis" = analisis_o3)
o_nv_tgf_namp_o_ss <- other_nv_tgf[which(other_nv_tgf[, "AMPATP"] == 0 & other_nv_tgf[, "Oxygen"] == 1 & other_nv_tgf[, "ShearStress"] == 1 ),]
analisis_o4 <- getDividingMuralFixed(o_nv_tgf_namp_o_ss, classed, Attractors, net)
AtypicalIV <- list("Environment Matrix" = o_nv_tgf_namp_o_ss,"Analisis" = analisis_o4)
o_nv_ntgf <- other_nvegf[which(other_nvegf[, "BMP9"] == 0 & other_nvegf[, "BMP10"] == 0 & other_nvegf[, "TGFB1"] == 0 ),]
o_nv_ntgf_ss <- o_nv_ntgf[which(o_nv_ntgf[, "ShearStress"] == 1),]
o_nv_ntgf_nigf_ss <- o_nv_ntgf[which(o_nv_ntgf[, "IGF"] == 0 & o_nv_ntgf[, "ShearStress"] == 1),]
analisis_o6 <- getDividingMuralFixed(o_nv_ntgf_nigf_ss, classed, Attractors, net)
AtypicalVI <- list("Environment Matrix" = o_nv_ntgf_nigf_ss,"Analisis" = analisis_o6)
o_nv_ntgf_igf_ss <- o_nv_ntgf[which(o_nv_ntgf[, "IGF"] == 1 & o_nv_ntgf[, "ShearStress"] == 1),]
o_nv_ntgf_igf_ss_no_adp <- o_nv_ntgf_igf_ss[which(o_nv_ntgf_igf_ss[, "AMPATP"] == 1 | o_nv_ntgf_igf_ss[, "Oxygen"] == 0),]
analisis_o7 <- getDividingMuralFixed(o_nv_ntgf_igf_ss_no_adp, classed, Attractors, net)
AtypicalVII <- list("Environment Matrix" = o_nv_ntgf_igf_ss_no_adp,"Analisis" = analisis_o7)
o_nv_ntgf_nss <- o_nv_ntgf[which(o_nv_ntgf[, "ShearStress"] == 0),]
o_nv_ntgf_nigf_nss <- o_nv_ntgf[which(o_nv_ntgf[, "IGF"] == 0 & o_nv_ntgf[, "ShearStress"] == 0),]
analisis_o5 <- getDividingMuralFixed(o_nv_ntgf_nigf_nss, classed, Attractors, net)
AtypicalV <- list("Environment Matrix" = o_nv_ntgf_nigf_nss,"Analisis" = analisis_o5)
o_nv_nt_i_nss <- o_nv_ntgf[which(o_nv_ntgf[, "IGF"] == 1 & o_nv_ntgf[, "ShearStress"] == 0),]
o_nv_nt_i_nss_ofa <- o_nv_nt_i_nss[which(o_nv_nt_i_nss[,"Oxygen"] == 0 | o_nv_nt_i_nss[,"FGF"] == 0 | o_nv_nt_i_nss[,"AMPATP"] == 1),]
analisis_o8 <- getDividingMuralFixed(o_nv_nt_i_nss_ofa, classed, Attractors, net)
AtypicalVIII <- list("Environment Matrix" = o_nv_nt_i_nss_ofa,"Analisis" = analisis_o8)
o_nv_nt_i_nss_nofa <- o_nv_nt_i_nss[which(o_nv_nt_i_nss[,"Oxygen"] == 1 & o_nv_nt_i_nss[,"FGF"] == 1 & o_nv_nt_i_nss[,"AMPATP"] == 0),]
analisis_o9 <- getDividingMuralFixed(o_nv_nt_i_nss_nofa, classed, Attractors, net)
AtypicalIX <- list("Environment Matrix" = o_nv_nt_i_nss_nofa,"Analisis" = analisis_o9)
Atypical <- list("AtypicalI" = AtypicalI, "AtypicalII" = AtypicalII, "AtypicalIII" = AtypicalIII, "AtypicalIV" = AtypicalIV, "AtypicalV" = AtypicalV, "AtypicalVI" = AtypicalVI, "AtypicalVII" = AtypicalVII, "AtypicalVIII" = AtypicalVIII, "AtypicalIX" = AtypicalIX)
microenvironments_by_behavior <- list("Phalanx" = Phalanx, "Stalk" = Stalk, "Tip" = Tip, "Atypical" = Atypical)
# Now lets analize the environments that cause cell division
pdivide <- phalanx_env_analisis$all_divide
tdivide <- c(vegf_tip_analisis$all_divide, no_vegfs1_analisis$all_divide, no_vegfs2_analisis$all_divide)
sdivide <- c(wnt_nigf_noxy_analisis$all_divide, wnt_igf_oxy_analisis$all_divide, nwnt_stalk_analisis$all_divide)
odivide <- c(analisis_other_vegf$all_divide,
		analisis_o2$all_divide,
		analisis_o3$all_divide,
		analisis_o4$all_divide,
		analisis_o5$all_divide,
		analisis_o6$all_divide,
		analisis_o7$all_divide,
		analisis_o8$all_divide,
		analisis_o9$all_divide)
all_divide <- c(pdivide, tdivide, sdivide, odivide)
all_divide_m <- environmentsToMatrix(microenvironment_variables, all_divide)
analisis_d_m <- getDividingMuralFixed(all_divide_m, classed, Attractors, net)

# Now lets analize the environments that inhibit cell division
pndivide <- phalanx_env_analisis$none_divide
tndivide <- c(vegf_tip_analisis$none_divide, no_vegfs1_analisis$none_divide, no_vegfs2_analisis$none_divide)
sndivide <- c(wnt_nigf_noxy_analisis$none_divide, wnt_igf_oxy_analisis$none_divide, nwnt_stalk_analisis$none_divide)
ondivide <- c(analisis_other_vegf$none_divide,
		analisis_o2$none_divide,
		analisis_o3$none_divide,
		analisis_o4$none_divide,
		analisis_o5$none_divide,
		analisis_o6$none_divide,
		analisis_o7$none_divide,
		analisis_o8$none_divide,
		analisis_o9$none_divide)
none_divide <- c(pndivide, tndivide, sndivide, ondivide)
none_divide_m <- environmentsToMatrix(microenvironment_variables, none_divide)
analisis_nd_m <- getDividingMuralFixed(none_divide_m, classed, Attractors, net)

# Now lets analize the environments that neither cause, nor inhibit cell division
psdivide <- phalanx_env_analisis$some_divide
tsdivide <- c(vegf_tip_analisis$some_divide, no_vegfs1_analisis$some_divide, no_vegfs2_analisis$some_divide)
ssdivide <- c(wnt_nigf_noxy_analisis$some_divide, wnt_igf_oxy_analisis$some_divide, nwnt_stalk_analisis$some_divide)
osdivide <- c(analisis_other_vegf$some_divide,
		analisis_o2$some_divide,
		analisis_o3$some_divide,
		analisis_o4$some_divide,
		analisis_o5$some_divide,
		analisis_o6$some_divide,
		analisis_o7$some_divide,
		analisis_o8$some_divide,
		analisis_o9$some_divide)
some_divide <- c(psdivide, tsdivide, ssdivide, osdivide)
some_divide_m <- environmentsToMatrix(microenvironment_variables, some_divide)
analisis_sd_m <- getDividingMuralFixed(some_divide_m, classed, Attractors, net)
environments_by_division <- list("All_divide" = all_divide_m,"None_divide"=none_divide_m,"Some_divide"=some_divide_m)
grouped_environments <- list("grouped_by_behavior" = microenvironments_by_behavior,"grouped_by_division" = environments_by_division)
save(grouped_environments, file = "grouped_environments.Rdata")

# Now lets analize the environments that cause mural cell recruitment
ppdgfb <- phalanx_env_analisis$all_pdgfb
tpdgfb <- c(vegf_tip_analisis$all_pdgfb, no_vegfs1_analisis$all_pdgfb, no_vegfs2_analisis$all_pdgfb)
spdgfb <- c(wnt_nigf_noxy_analisis$all_pdgfb, wnt_igf_oxy_analisis$all_pdgfb, nwnt_stalk_analisis$all_pdgfb)
opdgfb <- c(analisis_other_vegf$all_pdgfb,
		analisis_o2$all_pdgfb,
		analisis_o3$all_pdgfb,
		analisis_o4$all_pdgfb,
		analisis_o5$all_pdgfb,
		analisis_o6$all_pdgfb,
		analisis_o7$all_pdgfb,
		analisis_o8$all_pdgfb,
		analisis_o9$all_pdgfb)
all_pdgfb <- c(ppdgfb, tpdgfb, spdgfb, opdgfb)
all_pdgfb_m <- environmentsToMatrix(microenvironment_variables, all_pdgfb)
all_pdgfb_m_nigf <- all_pdgfb_m[which(all_pdgfb_m[, "IGF"] == 0),] # All (3/4)*2^12 = 3072 are in
all_pdgfb_m_igf <- all_pdgfb_m[which(all_pdgfb_m[, "IGF"] == 1),]
all_pdgfb_m_w_igf <- all_pdgfb_m_igf[which(all_pdgfb_m_igf[, "WNT5a"] == 1 | all_pdgfb_m_igf[, "WNT7a"] == 1),]
# All (3/4)^3*2^12 = 1728 are in
all_pdgfb_m_w_igf_no_amp <- all_pdgfb_m_w_igf[which(all_pdgfb_m_w_igf[, "Oxygen"] == 0 | all_pdgfb_m_w_igf[, "AMPATP"] == 1),]
analisis_p_m <- getDividingMuralFixed(all_pdgfb_m, classed, Attractors, net)

# Now lets analize the environments that inhibit mural cell recruitment
pnpdgfb <- phalanx_env_analisis$none_pdgfb
tnpdgfb <- c(vegf_tip_analisis$none_pdgfb, no_vegfs1_analisis$none_pdgfb, no_vegfs2_analisis$none_pdgfb)
snpdgfb <- c(wnt_nigf_noxy_analisis$none_pdgfb, wnt_igf_oxy_analisis$none_pdgfb, nwnt_stalk_analisis$none_pdgfb)
onpdgfb <- c(analisis_other_vegf$none_pdgfb,
		analisis_o2$none_pdgfb,
		analisis_o3$none_pdgfb,
		analisis_o4$none_pdgfb,
		analisis_o5$none_pdgfb,
		analisis_o6$none_pdgfb,
		analisis_o7$none_pdgfb,
		analisis_o8$none_pdgfb,
		analisis_o9$none_pdgfb)
none_pdgfb <- c(pnpdgfb, tnpdgfb, snpdgfb, onpdgfb)
none_pdgfb_m <- environmentsToMatrix(microenvironment_variables, none_pdgfb)
n_p_m_o <- none_pdgfb_m[which(none_pdgfb_m[,"Oxygen"] == 1),]
n_p_m_o_na <- none_pdgfb_m[which(none_pdgfb_m[,"Oxygen"] == 1 & none_pdgfb_m[,"AMPATP"] == 0),]
n_p_m_no_a <- none_pdgfb_m[which(none_pdgfb_m[,"Oxygen"] == 0 | none_pdgfb_m[,"AMPATP"] == 1),]
n_p_m_no_a1 <- n_p_m_no_a[which(n_p_m_no_a[,"ShearStress"] == 1),]
n_p_m_no_a2 <- n_p_m_no_a1[which(n_p_m_no_a1[,"IGF"] == 1 | n_p_m_no_a1[,"BMP9"] == 0),]
n_p_m_no_a3 <- n_p_m_no_a1[which(n_p_m_no_a1[,"IGF"] == 0 & n_p_m_no_a1[,"BMP9"] == 1),]
n_p_m_no_a4 <- n_p_m_no_a[which(n_p_m_no_a[,"ShearStress"] == 0),]
n_p_m_no_a5 <- n_p_m_no_a4[which(n_p_m_no_a4[,"BMP9"] == 0),]
n_p_m_no_a6 <- n_p_m_no_a5[which(n_p_m_no_a5[,"FGF"] == 1 |n_p_m_no_a5[,"BMP10"] == 1|n_p_m_no_a5[,"TGFB1"] == 1),] # ALL 448
n_p_m_no_a7 <- n_p_m_no_a4[which(n_p_m_no_a4[,"BMP9"] == 1),]
n_p_m_o_na_v <- n_p_m_o_na[which(n_p_m_o_na[,"VEGFC_Dp"] == 1 | n_p_m_o_na[,"VEGFAxxxP"] == 1),] # All (3/4)*2^14 = 12288 are here!
n_p_m_o_na_nv <- n_p_m_o_na[which(n_p_m_o_na[,"VEGFC_Dp"] == 0 & n_p_m_o_na[,"VEGFAxxxP"] == 0),]
npmonanv_w <- n_p_m_o_na_nv[which(n_p_m_o_na_nv[,"WNT5a"] == 1 | n_p_m_o_na_nv[,"WNT7a"] == 1),]
npmonanv_w_inb <- npmonanv_w[which(npmonanv_w[,"IGF"] == 1 | npmonanv_w[,"BMP9"] == 0),] #(9/16)*2^12 = 2304
npmonanv_w_nib <- npmonanv_w[which(npmonanv_w[,"IGF"] == 0 & npmonanv_w[,"BMP9"] == 1),] #(3/4)*2^8 = 192
npmonanv_nw <- n_p_m_o_na_nv[which(n_p_m_o_na_nv[,"WNT5a"] == 0 & n_p_m_o_na_nv[,"WNT7a"] == 0),]
npmonanv_nw_nb9 <- npmonanv_nw[which(npmonanv_nw[,"BMP9"] == 0),] # All 512!
npmonanv_nw_b9 <- npmonanv_nw[which(npmonanv_nw[,"BMP9"] == 1),]
npmonanv_nw_b91 <- npmonanv_nw_b9[which(npmonanv_nw_b9[,"DLL4p"] == 1 | npmonanv_nw_b9[,"JAGp"] == 0 | npmonanv_nw_b9[,"FGF"] == 1),]
npmonanv_nw_b92 <- npmonanv_nw_b9[which(npmonanv_nw_b9[,"IGF"] == 0 | npmonanv_nw_b9[,"ShearStress"] == 0),]
npmonanv_nw_b93 <- npmonanv_nw_b9[which(npmonanv_nw_b9[,"IGF"] == 1 & npmonanv_nw_b9[,"ShearStress"] == 1),] #128 all in
analisis_np_m <- getDividingMuralFixed(none_pdgfb_m, classed, Attractors, net)

# Now lets analize the environments that neither cause, nor inhibit mural cell recruitment
pspdgfb <- phalanx_env_analisis$some_pdgfb
tspdgfb <- c(vegf_tip_analisis$some_pdgfb, no_vegfs1_analisis$some_pdgfb, no_vegfs2_analisis$some_pdgfb)
sspdgfb <- c(wnt_nigf_noxy_analisis$some_pdgfb, wnt_igf_oxy_analisis$some_pdgfb, nwnt_stalk_analisis$some_pdgfb)
ospdgfb <- c(analisis_other_vegf$some_pdgfb,
		analisis_o2$some_pdgfb,
		analisis_o3$some_pdgfb,
		analisis_o4$some_pdgfb,
		analisis_o5$some_pdgfb,
		analisis_o6$some_pdgfb,
		analisis_o7$some_pdgfb,
		analisis_o8$some_pdgfb,
		analisis_o9$some_pdgfb)
some_pdgfb <- c(pspdgfb, tspdgfb, sspdgfb, ospdgfb)
some_pdgfb_m <- environmentsToMatrix(microenvironment_variables, some_pdgfb)
analisis_sp_m <- getDividingMuralFixed(some_pdgfb_m, classed, Attractors, net)

# Using Notch
other_nv_notch <- other_nvegf[which(other_nvegf[, "DLL4p"] == 1 & other_nvegf[, "JAGp"] == 0),]
o_nv_nnotch <- other_nvegf[which(other_nvegf[, "DLL4p"] == 0 | other_nvegf[, "JAGp"] == 1),]
o_nv_nn_igf <- o_nv_nnotch[which(o_nv_nnotch[, "IGF"] == 1),]
o_nv_nn_tgf <- o_nv_nnotch[which(o_nv_nnotch[, "BMP9"] == 1 | o_nv_nnotch[, "BMP10"] == 1 | o_nv_nnotch[, "TGFB1"] == 1 ),]
o_nv_nn_tgf_amp_no_nss = o_nv_nn_tgf[which(o_nv_nn_tgf[, "AMPATP"] == 1 | o_nv_nn_tgf[, "Oxygen"] == 0 | o_nv_nn_tgf[, "ShearStress"] == 0 ),]

# Now lets sort the environments by number of attractors
envsByNofAttr <- environmentsByNumberOfAttractors(classed)




