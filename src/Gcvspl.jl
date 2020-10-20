#############################################################################
#Copyright (c) 2020 Charles Le Losq
#
#The MIT License (MIT)
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the #Software without restriction, including without limitation the rights to use, copy, #modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, #and to permit persons to whom the Software is furnished to do so, subject to the #following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, #INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#############################################################################
module Gcvspl

unixpath = "../deps/src/gcvspline/libgcvspl"
const gcvlib = joinpath(dirname(@__FILE__), @static Sys.isunix() ? unixpath : winpath)

#function __init__()
#	if (Libdl.dlopen_e(gcvlib) == C_NULL)
#        error("GCVSPL not properly compiled. Run Pkg.build(\"gcvspl\"). Windows auto-build is not setup, you might want to build the library manually.")
#	end
#end

include("gcvspl_wrapper.jl")

export gcv_spl, spl_der

end # module
