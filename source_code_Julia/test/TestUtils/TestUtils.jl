"""
# TestUtils file

- Julia version: 
- Author: Mars Semenova
- Date: 2024-06-19

# Examples

```jldoctest
julia>
```
"""

using Test

coefTypes = [Base.BitInteger_types..., BigInt, Float32, Float64]
floatCoefTypes = [Float32, Float64, BigFloat] # TODO: add BigRationals

"""
Functions which enable logging on fail.
"""
_onFail(body, x) = error("I might have overlooked something: $x")
_onFail(body, _::Test.Pass) = nothing
_onFail(body, _::Test.Fail)= body()

"""
Helper function which logs test details on fail.
"""
function _logInfo(args...)
    for x = 1:length(args)
        @info args[x]
    end
end