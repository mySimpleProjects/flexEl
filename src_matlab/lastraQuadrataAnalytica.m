function wa = lastraQuadrataAnalytica(data)

DD = (data.E*data.thickness^3)/(12*(1-data.nu^2));
alfa = 0.004062;
wa = alfa *(data.q*data.L^4)/DD;
