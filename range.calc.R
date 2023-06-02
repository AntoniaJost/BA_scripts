# function for calculating ylim for plots
range_calc = function(eval_value) {
  ylim = c()
  is_naN_inf = is.na(eval_value) | is.nan(eval_value) | is.infinite(eval_value)
  ylim = range(eval_value[!is_naN_inf])
  return(ylim)
}