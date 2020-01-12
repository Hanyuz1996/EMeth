dlaplace <-
function(y,m=0,s=1){
    return(exp(-abs(y-m)/s)/(2*s))
}
