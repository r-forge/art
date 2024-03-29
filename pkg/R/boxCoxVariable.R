#-------------------------------------------------------------------------------
# Revision history:
# checked in 2009-09-29 by J. Fox (moved/renamed from car)
#-------------------------------------------------------------------------------

# constructed variable for Box-Cox transformation (J. Fox)

boxCoxVariable <- function(y) {
	geo.mean <- exp(mean(log(y), na.rm=TRUE))
	y*(log(y/geo.mean) - 1)
}
