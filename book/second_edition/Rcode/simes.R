#-------------------------------------------------------------

simes = function(n,p)
{
# This function returns the Simes adjusted p-value
# based on n p-values in the vector p
#
 p.sort=sort(p)
 psimes=min(p.sort*n/c(1:n))

 psimes
}

#-------------------------------------------------------------