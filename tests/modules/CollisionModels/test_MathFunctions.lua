function erf(z)
    local z2 = math.abs(z)
    local t = 1 / (1 + 0.32759109962 * z2)
    local res = (    - 1.061405429 ) * t
    res = (res + 1.453152027 ) * t
    res = (res - 1.421413741 ) * t
    res = (res + 0.2844966736) * t
    res =((res - 0.254829592 ) * t) * math.exp(-z2*z2)
    res = res + 1
    if z < 0 then res = -res end
    return res
end

print(erf(2.0))
print(erf(1.0))
print(erf(0.5))