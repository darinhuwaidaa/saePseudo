#' @title Small Area Estimation using Averaging Pseudo Area Level Model
#'
#' @description Provides function for small area estimation at area level
#'     using averaging pseudo area level model for variables of interest. A
#'     dataset produced by data generation are also provided. This package
#'     estimates small areas at the village level and then aggregates them to
#'     the sub-district, region, and provincial levels.
#'
#' @param prov Vector containing information of province
#' @param reg Vector containing information of region
#' @param sub Vector containing information of subdistrict
#' @param vill Vector containing information of village
#' @param y Direct estimation for each area
#' @param x Auxiliary variable for each area
#' @param var Sampling variances of direct estimators for each domain
#' @param N Number of population in each area
#' @param method Method used to fit the Fay-Herriot model, which can be either "ML", "REML" or "FH" methods. Default is "REML" method
#'
#' @export avgPseudo
#'
#' @import dplyr
#' @import sae
#'
#' @return This function returns a list of the following objects:
#'  \item{Est_Area3}{A dataframe with the values of Small Area Estimation with averaging pseudo area level model for sub-district level}
#'  \item{Est_Area2}{A dataframe with the values of Small Area Estimation with averaging pseudo area level model for region level}
#'  \item{Est_Area1}{A dataframe with the values of Small Area Estimation with averaging pseudo area level model for provincial level}
#'
#' @examples
#' # Load Dataset
#' data(dataVill)
#' saeAVG.Pseudo <- avgPseudo(prov = dataVill$Area1, reg = dataVill$Area2, sub = dataVill$Area3,
#'                  vill = dataVill$Area4, y = dataVill$ydir_area4, x = dataVill$X1,
#'                  var = dataVill$vardir_area4, N = dataVill$N, method="REML")
#'
#' # Result
#' saeAVG.Pseudo$Est_Area3
#' saeAVG.Pseudo$Est_Area2
#' saeAVG.Pseudo$Est_Area1
#'

avgPseudo <- function(prov, reg, sub, vill, y, x, var, N, method="REML") {

  result <- list(Est_Area3=NA, Est_Area2=NA, Est_Area1=NA)
  Province <- NA
  Region <- NA
  Subdistrict <- NA
  N_vill <- NA
  Y_eblup <- NA
  MSE_eblup <- NA
  Weight_agr_villsub <- NA
  y_agr_villsub <- NA
  rse_agr_villsub <- NA
  Weight_agr_subreg <- NA
  y_agr_subreg <- NA
  rse_agr_subreg <- NA
  Weight_agr_regprov <- NA
  y_agr_regprov <- NA
  rse_agr_regprov <- NA

  # Check the NA value in the var and N variables
  if (any(is.na(var)))
    stop("Argument var contains NA values.")
  if (any(is.na(N)))
    stop("Argument N contains NA values.")

  # Select method for SAE
  if (method=="REML") {
    sae_result <- mseFH(y ~ x, var, method = "REML")
  } else if (method=="ML") {
    sae_result <- mseFH(y ~ x, var, method = "ML")
  } else if (method=="FH") {
    sae_result <- mseFH(y ~ x, var, method = "FH")
  }

  # Extract eblup est, mse, and rse
  y_eblup <- sae_result$est$eblup
  mse_eblup <- sae_result$mse
  rse_eblup <- (sqrt(mse_eblup) / y_eblup)*100

  # Create the eblup_data dataset
  eblup_data <- data.frame(Village = vill, Y_eblup = y_eblup, MSE_eblup = mse_eblup, RSE_eblup = rse_eblup)

  # Create the sae_vill dataset
  sae_vill <- data.frame(Province = prov, Region = reg, Subdistrict = sub, N_vill = N)

  # Join the sae_vill dataset with eblup_data
  sae_vill <- bind_cols(sae_vill, eblup_data)

  # Calculate N_sub (amount of data (population) in each Subdistrict)
  N_sub <- sae_vill %>%
    group_by(Province, Region, Subdistrict) %>%
    summarize(N_sub = sum(N_vill))

  # Join N_sub into the sae_vill dataset
  sae_vill <- left_join(sae_vill, N_sub, by = c("Province", "Region", "Subdistrict"))

  # Calculate Weight_agr_villsub (Nij/Ni) value
  sae_vill <- sae_vill %>%
    mutate(Weight_agr_villsub = N_vill / N_sub)

  # Calculate y_agr_villsub per Subdistrict
  vill_agr_sub <- sae_vill %>%
    group_by(Province, Region, Subdistrict) %>%
    summarize(y_agr_villsub = sum(Weight_agr_villsub * Y_eblup))

  # Calculate variance per Subdistrict
  var_agr_villsub <- sae_vill %>%
    group_by(Province, Region, Subdistrict) %>%
    summarise(var_agr_villsub = sum(Weight_agr_villsub^2 * MSE_eblup))

  # Left join for the var_agr_villsub column in the vill_agr_sub dataset (to make it easier to calculate rse)
  vill_agr_sub <- left_join(vill_agr_sub, dplyr::select(var_agr_villsub, Province, Region, Subdistrict, var_agr_villsub), by = c("Province", "Region", "Subdistrict"))

  # Cakculate rse_agr_villsub
  vill_agr_sub <- vill_agr_sub%>%
    mutate(rse_agr_villsub = (sqrt(var_agr_villsub) / y_agr_villsub)*100)

  # Join N_sub into vill_agr_sub dataset
  vill_agr_sub <- left_join(vill_agr_sub, N_sub, by = c("Province", "Region", "Subdistrict"))

  # Calculate N_reg (amount of data (population) in each Region)
  N_reg <- vill_agr_sub %>%
    group_by(Province, Region) %>%
    summarize(N_reg = sum(N_sub))

  # Join N_reg into the vill_agr_sub dataset
  vill_agr_sub <- left_join(vill_agr_sub, N_reg, by = c("Province", "Region"))

  # Calculate Weight_agr_subreg (Nij/Ni) value
  vill_agr_sub <- vill_agr_sub%>%
    mutate(Weight_agr_subreg = N_sub / N_reg)

  # Show vill_agr_sub results for selected columns (user output)
  result_villsub <- vill_agr_sub %>%
    dplyr::select(Province, Region, Subdistrict, y_agr_villsub, rse_agr_villsub)

  # Calculate y_agr_subreg per Region
  sub_agr_reg <- vill_agr_sub %>%
    group_by(Province, Region) %>%
    summarize(y_agr_subreg = sum(Weight_agr_subreg * y_agr_villsub))

  # Calculate var_agr_subreg
  var_agr_subreg <- vill_agr_sub %>%
    group_by(Province, Region) %>%
    summarise(var_agr_subreg = sum(Weight_agr_subreg^2 * var_agr_villsub))

  # Left join for the var_agr_subreg column in the sub_agr_reg dataset (to make it easy to calculate rse)
  sub_agr_reg <- left_join(sub_agr_reg, dplyr::select(var_agr_subreg, Province, Region, var_agr_subreg), by = c("Province", "Region"))

  # Calculate rse_agr_subreg
  sub_agr_reg <- sub_agr_reg%>%
    mutate(rse_agr_subreg = (sqrt(var_agr_subreg) / y_agr_subreg)*100)

  # Join N_reg into sub_agr_reg dataset
  sub_agr_reg <- left_join(sub_agr_reg, N_reg, by = c("Province", "Region"))

  # Calculate N_prov (amount of data (population) in each Province)
  N_prov <- sub_agr_reg %>%
    group_by(Province) %>%
    summarize(N_prov = sum(N_reg))

  # Join N_prov into sub_agr_reg dataset
  sub_agr_reg <- left_join(sub_agr_reg, N_prov, by = "Province")

  # Calculate Weight_agr_regprov (Nij/Ni) value
  sub_agr_reg <- sub_agr_reg%>%
    mutate(Weight_agr_regprov = N_reg / N_prov)

  # Show sub_agr_reg results for selected column (user output)
  result_subreg <- sub_agr_reg %>%
    dplyr::select(Province, Region, y_agr_subreg, rse_agr_subreg)

  # Calculate y_agr_regprov per Province
  reg_agr_prov <- sub_agr_reg %>%
    group_by(Province) %>%
    summarize(y_agr_regprov = sum(Weight_agr_regprov * y_agr_subreg))

  # Calculate var_agr_regprov
  var_agr_regprov <- sub_agr_reg %>%
    group_by(Province) %>%
    summarise(var_agr_regprov = sum(Weight_agr_regprov^2 * var_agr_subreg))

  # Left join for the var_agr_regprov column in the reg_agr_prov dataset (to make it easy to calculate rse)
  reg_agr_prov <- left_join(reg_agr_prov, dplyr::select(var_agr_regprov, Province, var_agr_regprov), by = "Province")

  # Calculate rse_agr_regprov
  reg_agr_prov <- reg_agr_prov%>%
    mutate(rse_agr_regprov = (sqrt(var_agr_regprov) / y_agr_regprov)*100)

  # Show reg_agr_prov results for selected column (user output)
  result_regprov <- reg_agr_prov %>%
    dplyr::select(Province, y_agr_regprov, rse_agr_regprov)

  result$Est_Area3 = result_villsub
  result$Est_Area2 = result_subreg
  result$Est_Area1 = result_regprov

  return(result)
}


