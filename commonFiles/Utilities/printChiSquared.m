function pVal = printChiSquared(n1,n2,N1,N2)

[pVal,chi2stat]= chi2_fun(n1,n2,N1,N2);
pValString = printPValue(pVal);
fprintf('\\chi^2 (1, N = %d) = %.2f, P %s \n',N1 + N2,chi2stat,pValString);


end