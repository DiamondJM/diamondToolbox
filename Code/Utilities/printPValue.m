function pValString = printPValue(pVal)

if pVal < 0.0001
    pValString = '$p < .0001$';
elseif pVal < 0.001
    pValString = sprintf('$p = %.4f$',pVal);
elseif pVal < 0.05
    pValString = sprintf('$p = %.3f$',pVal);
else
    pValString = sprintf('$p = %.2f$',pVal);
end

pValString = strrep(pValString, '0.', '.');