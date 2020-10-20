using gcvspl, Test

# Dummy data
x = collect(0.:1:100)
y = 2*x

spline_coeff = 1.0
SplOrder = 3

c, WK, IER = gcvspl.gcv_spl(x,y,ese,spline_coeff;SplineOrder = Int32(SplOrder-1)) # with cubic spline as the default
y_calc = gcvspl.spl_der(x,x,c,SplineOrder= Int32(SplOrder-1))

@test isapprox(y,y_calc,atol=1e-1)
